"""
Spatial binning of multi-frequency cubes before spectral fitting.

Used by `specidx_map(..., binning=...)` (radio_sed.py) to group pixels into
bins that reach a target signal-to-noise in every band, so that one spectrum
per bin (instead of per pixel) is fitted.

Standalone on purpose: imports only numpy and scipy, never `mlibs`.

Methods are registered in `BINNING_METHODS`; to add a new one (e.g. Voronoi),
write a function with the same call signature as `contour_binning` that
returns an integer label map, and add it to the dictionary.

Currently available:
- 'contour': contour binning, following Sanders (2006, MNRAS 371, 829).
  Bins grow along contours of a (smoothed) reference brightness map, so they
  follow the morphology of the emission instead of being compact/round.
"""
import heapq
import numpy as np
from scipy import ndimage


def bin_snr_per_band(sums, n_pix, rms, beam_area_px=1.0):
    """
    Signal-to-noise of a bin, per band.

    The bin signal is the mean pixel value; its noise is rms/sqrt(N_eff),
    where N_eff = max(1, n_pix / beam_area_px) is the number of independent
    beams in the bin (pixels within one beam have correlated noise).

    Parameters
    ----------
    sums : array (nfreq,)
        Sum of the pixel values of the bin, per band.
    n_pix : int
        Number of pixels in the bin.
    rms : array (nfreq,)
        Noise rms per pixel, per band.
    beam_area_px : float
        Beam area in pixels.
    """
    n_eff = max(1.0, n_pix / beam_area_px)
    return (sums / n_pix) / (rms / np.sqrt(n_eff))


def contour_binning(ref_map, cube, mask, rms, target_snr=1.0,
                    beam_area_px=1.0, smooth_sigma_px=1.0, constraint=3.0,
                    verbose=0):
    """
    Contour binning (Sanders 2006) of the pixels inside `mask`.

    Algorithm
    ---------
    1. The reference map is (optionally) Gaussian-smoothed; its values define
       the "contours" bins grow along.
    2. A new bin starts at the brightest pixel not yet binned.
    3. The bin grows by adding the neighbouring (4-connected) unbinned pixel
       whose reference value is closest to the value of the starting pixel,
       until the bin reaches `target_snr` in *every* band (the weakest band
       sets the bin size).
    4. Geometric constraint: a pixel is only accepted if its distance to the
       bin centroid is <= constraint * sqrt(n_pix / pi), which avoids long
       thin bins. Rejected pixels stay free for other bins.
    5. Bins that cannot reach the target (they run out of free neighbours)
       are merged into the touching bin whose mean reference value is the
       closest. Isolated ones (no touching bin) are kept as they are.

    Parameters
    ----------
    ref_map : array (ny, nx)
        Map whose brightness contours the bins follow (e.g. the mean over
        bands of the masked cube).
    cube : array (ny, nx, nfreq)
        Image cube used for the signal-to-noise of each bin.
    mask : bool array (ny, nx)
        Pixels to bin.
    rms : array (nfreq,)
        Noise rms per pixel, per band (same units as cube).
    target_snr : float
        Signal-to-noise each bin must reach in every band.
    beam_area_px : float
        Beam area in pixels (noise correlation, see `bin_snr_per_band`).
    smooth_sigma_px : float or None
        Sigma (pixels) of the Gaussian smoothing of `ref_map`. 0/None: none.
    constraint : float or None
        Geometric constraint factor (Sanders uses ~2-3). None: no constraint.
    verbose : int
        If > 0, print a summary.

    Returns
    -------
    labels : int array (ny, nx)
        Bin index (0..nbins-1) of each pixel; -1 outside `mask`.
    """
    mask = np.asarray(mask, dtype=bool)
    ny, nx = mask.shape
    rms = np.maximum(np.asarray(rms, dtype=float), np.finfo(float).tiny)
    cube = np.asarray(cube)
    valid = mask & np.all(np.isfinite(cube), axis=2)

    ref = np.where(valid, np.asarray(ref_map, dtype=float), 0.0)
    if smooth_sigma_px:
        # normalised convolution, so pixels outside the mask do not dilute it
        w = ndimage.gaussian_filter(valid.astype(float), smooth_sigma_px)
        ref = ndimage.gaussian_filter(ref, smooth_sigma_px)
        with np.errstate(invalid='ignore', divide='ignore'):
            ref = np.where(valid, ref / w, 0.0)

    labels = np.full((ny, nx), -1, dtype=np.int64)
    # pixels in descending order of reference value (bin seeds)
    order = np.argsort(-ref[valid], kind='stable')
    seeds = np.column_stack(np.nonzero(valid))[order]

    reached = []
    nbins = 0
    for si, sj in seeds:
        if labels[si, sj] != -1:
            continue
        b = nbins
        nbins += 1
        v0 = ref[si, sj]
        labels[si, sj] = b
        n = 1
        sums = cube[si, sj, :].astype(float).copy()
        ci, cj = float(si), float(sj)          # running sums for centroid
        heap = []
        seen = {(si, sj)}

        def push_neighbours(i, j):
            for di, dj in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                a, c = i + di, j + dj
                if (0 <= a < ny and 0 <= c < nx and valid[a, c]
                        and labels[a, c] == -1 and (a, c) not in seen):
                    seen.add((a, c))
                    heapq.heappush(heap, (abs(ref[a, c] - v0), a, c))

        push_neighbours(si, sj)
        ok = np.all(bin_snr_per_band(sums, n, rms, beam_area_px) >= target_snr)
        while not ok and heap:
            _, a, c = heapq.heappop(heap)
            if labels[a, c] != -1:
                continue
            if constraint is not None:
                r_max = constraint * np.sqrt(n / np.pi)
                if np.hypot(a - ci / n, c - cj / n) > max(r_max, 1.0):
                    continue
            labels[a, c] = b
            n += 1
            sums += cube[a, c, :]
            ci += a
            cj += c
            push_neighbours(a, c)
            ok = np.all(bin_snr_per_band(sums, n, rms, beam_area_px) >= target_snr)
        reached.append(ok)

    reached = np.array(reached, dtype=bool)
    n_failed = int((~reached).sum())

    # merge bins that did not reach the target into the closest-valued
    # touching bin (smallest first)
    if n_failed:
        sizes = np.bincount(labels[labels >= 0], minlength=nbins)
        mean_ref = (np.bincount(labels[labels >= 0], weights=ref[labels >= 0],
                                minlength=nbins) / np.maximum(sizes, 1))
        objects = ndimage.find_objects(labels + 1)
        for b in sorted(np.nonzero(~reached)[0], key=lambda k: sizes[k]):
            sl = objects[b]
            if sl is None:
                continue
            sl = tuple(slice(max(s.start - 1, 0), s.stop + 1) for s in sl)
            sub = labels[sl]
            this = sub == b
            if not this.any():
                continue                        # already merged away
            ring = ndimage.binary_dilation(this) & ~this & (sub >= 0)
            neighbours = np.unique(sub[ring])
            if neighbours.size == 0:
                continue                        # isolated: keep as is
            target = neighbours[np.argmin(np.abs(mean_ref[neighbours] - mean_ref[b]))]
            sub[this] = target                  # view into labels
            # keep the merged bin's statistics up to date
            sizes[target] += sizes[b]
            mean_ref[target] = (mean_ref[target] * (sizes[target] - sizes[b])
                                + mean_ref[b] * sizes[b]) / sizes[target]
            sizes[b] = 0
            objects[target] = _union_slices(objects[target], objects[b])

    labels = _relabel_consecutive(labels)
    if verbose > 0:
        n_final = labels.max() + 1
        npx = np.bincount(labels[labels >= 0])
        print(f" > contour binning: {int(valid.sum())} pixels -> {n_final} bins "
              f"(target SNR={target_snr} in every band, "
              f"beam area={beam_area_px:.1f} px); pixels per bin: "
              f"median {np.median(npx):.0f}, max {npx.max()}; "
              f"{n_failed} bins below target were merged or kept.")
    return labels


def _union_slices(s1, s2):
    if s1 is None:
        return s2
    if s2 is None:
        return s1
    return tuple(slice(min(a.start, b.start), max(a.stop, b.stop))
                 for a, b in zip(s1, s2))


def _relabel_consecutive(labels):
    out = np.full_like(labels, -1)
    used = labels >= 0
    _, inv = np.unique(labels[used], return_inverse=True)
    out[used] = inv
    return out


BINNING_METHODS = {
    'contour': contour_binning,
}


def make_bins(method, ref_map, cube, mask, rms, **kwargs):
    """
    Dispatch to a binning method in `BINNING_METHODS`.

    Returns an int label map (ny, nx): bin index per pixel, -1 outside.
    """
    if method not in BINNING_METHODS:
        raise ValueError(f"Unknown binning method '{method}'. "
                         f"Available: {sorted(BINNING_METHODS)}")
    return BINNING_METHODS[method](ref_map, cube, mask, rms, **kwargs)


def bin_spectra(labels, cube, cube_res):
    """
    Mean image and residual spectrum of each bin.

    The residual is averaged the same way as the image, so the bin's
    residual is the noise of the bin's mean spectrum (the same convention
    specidx_map uses for single pixels).

    Returns
    -------
    spectra, spectra_res : arrays (nbins, nfreq)
    """
    nb = labels.max() + 1
    nf = cube.shape[2]
    used = labels >= 0
    lab = labels[used]
    counts = np.bincount(lab, minlength=nb).astype(float)
    spectra = np.empty((nb, nf))
    spectra_res = np.empty((nb, nf))
    for k in range(nf):
        spectra[:, k] = np.bincount(lab, weights=cube[:, :, k][used], minlength=nb) / counts
        spectra_res[:, k] = np.bincount(lab, weights=cube_res[:, :, k][used], minlength=nb) / counts
    return spectra, spectra_res
