"""
Astrometric alignment of radio/optical images across frequency and instrument.

Unlike the other files loaded by ``mlibs.py``, this module is **self-contained**:
it carries its own imports and depends only on numpy/scipy/astropy. It can
therefore be imported directly (``import image_alignment``) for testing, without
pulling in CASA, JAX, bilby and the rest of the ``mlibs`` preamble. It is still
listed in ``mlibs.py`` so every function below is also reachable as
``mlibs.<name>``.

Why this exists
---------------
The previous estimators (``peak_image_alignment``,
``structural_image_alignment`` in ``data_io.py``) were both built on *phase*
correlation -- the cross-power spectrum divided by its own modulus. For
beam-convolved radio images that is the wrong statistic: the restoring beam is a
low-pass filter, so whitening the spectrum boosts exactly those frequencies that
contain nothing but noise. Measured on NGC 7469 sub-band data, the same skimage
routine gave (2.48, -0.54) with ``normalization='phase'`` and (2.380, -1.620)
with ``normalization=None``, against an injected truth of (2.37, -1.62).

``structural_image_alignment`` had a second, fatal problem: it multiplied *both*
images by the *same* mask before correlating. The mask's autocorrelation then
dominates the correlation surface and pins the peak at zero shift, which is why
it returned exactly (0, 0) for every input ever given to it. Nothing here
multiplies both images by a shared mask before an FFT; masks select pixels for
the chi2/MI sums, and the cross-correlation path uses a soft taper.

Sign convention
---------------
Throughout this module ``(dy, dx)`` is **the correction to add to the target's
pixel coordinates to bring it onto the reference** -- i.e. the value you hand
straight to ``scipy.ndimage.shift`` or to :func:`apply_pixel_shift`. This
matches what ``peak_image_alignment`` returned, so call sites do not flip.
"""

import os
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage, optimize, signal
from astropy.io import fits
from astropy.stats import mad_std

__all__ = [
    'ShiftResult',
    'estimate_image_shift',
    'apply_pixel_shift',
    'split_shift',
    'align_image_list',
]

# Methods that can be requested via `method=`.
_METHODS = ('chi2', 'mi', 'xcorr', 'chi2_shift', 'auto', 'ensemble')

# Legacy `shift_correction_mode` values -> new method names.
_LEGACY_METHODS = {
    'peak': 'chi2',
    'structural': 'chi2',
    'image_diff': 'chi2',
}


# ---------------------------------------------------------------------------
# result container
# ---------------------------------------------------------------------------
@dataclass
class ShiftResult:
    """Outcome of a single alignment estimate.

    Attributes
    ----------
    dy, dx : float
        The correction to apply to the target, in pixels (see the module
        docstring for the sign convention). These are *after* the QA gate, so
        they are what should actually be applied.
    dy_err, dx_err : float
        1-sigma uncertainties, or NaN when the backend cannot provide them.
    dy_raw, dx_raw : float
        The estimate before the QA gate zeroed or clamped it.
    method : str
        Backend that produced the final answer.
    flags : list of str
        Diagnostics worth acting on: ``'consistent_with_zero'``,
        ``'clamped'``, ``'component_flip_suspected'``, ``'no_uncertainty'``,
        ``'beam_matched'``, ``'not_converged'``, ``'empty_mask'``.
    diagnostics : dict
        Per-backend estimates, mask size, rms values, beams -- everything needed
        to understand why the answer came out the way it did.
    """

    dy: float = 0.0
    dx: float = 0.0
    dy_err: float = np.nan
    dx_err: float = np.nan
    dy_raw: float = 0.0
    dx_raw: float = 0.0
    method: str = ''
    flags: list = field(default_factory=list)
    diagnostics: dict = field(default_factory=dict)

    def __iter__(self):
        """So ``dy, dx = estimate_image_shift(...)`` keeps working."""
        return iter((self.dy, self.dx))

    @property
    def magnitude(self):
        return float(np.hypot(self.dy, self.dx))

    def __repr__(self):
        err = ''
        if np.isfinite(self.dy_err):
            err = f' +/- ({self.dy_err:.3f}, {self.dx_err:.3f})'
        flg = f' flags={self.flags}' if self.flags else ''
        return (f'ShiftResult(dy={self.dy:.3f}, dx={self.dx:.3f}){err} '
                f'[{self.method}]{flg}')


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _as_plane(image):
    """Return a float 2D array from an array or a FITS path.

    WSClean writes ``(1, 1, ny, nx)``; CASA exports sometimes ``(1, ny, nx)``.
    Both collapse to the celestial plane.
    """
    if isinstance(image, str):
        data = fits.getdata(image)
    else:
        data = np.asarray(image)
    data = np.squeeze(np.asarray(data, dtype=float))
    if data.ndim != 2:
        raise ValueError(f'expected a 2D image, got shape {data.shape}')
    return data


def _header_of(image):
    """Header for a FITS path, else None."""
    if isinstance(image, str):
        return fits.getheader(image)
    return None


def _beam_of(header):
    """(bmaj, bmin, bpa) in degrees/degrees/degrees, or None if absent."""
    if header is None:
        return None
    bmaj, bmin = header.get('BMAJ'), header.get('BMIN')
    if bmaj is None or bmin is None:
        return None
    if not (np.isfinite(bmaj) and np.isfinite(bmin)) or bmaj <= 0:
        return None
    return float(bmaj), float(bmin), float(header.get('BPA', 0.0) or 0.0)


def _cell_of(header):
    """Pixel scale in degrees, or None."""
    if header is None:
        return None
    for key in ('CDELT1', 'CD1_1'):
        if key in header and header[key]:
            return abs(float(header[key]))
    return None


def _noise_of(image, data):
    """Noise estimate: the companion WSClean residual when one exists, else
    ``mad_std`` of the image itself.

    A residual is the honest noise estimate for a CLEANed image; ``mad_std`` of
    the restored image is biased high wherever the source fills a good fraction
    of the frame.
    """
    if isinstance(image, str):
        for a, b in (('-image', '-residual'), ('.image', '.residual')):
            if a in image:
                cand = image.replace(a, b)
                if os.path.exists(cand):
                    try:
                        return float(mad_std(_as_plane(cand), ignore_nan=True))
                    except Exception:
                        break
    return float(mad_std(data, ignore_nan=True))


def _clean(data):
    """NaN/inf -> 0, so the FFTs below stay finite."""
    return np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)


# ---------------------------------------------------------------------------
# beam matching
# ---------------------------------------------------------------------------
def _deconvolve_beams(target, source):
    """Gaussian `k` such that ``source (*) k == target``.

    Standard elliptical-beam difference (the same algebra CASA's
    ``deconvolvefrombeam`` uses). Beams are ``(maj, min, pa)``, pa in degrees.
    Returns ``(maj, min, pa)`` or None when `target` is not broader than
    `source` in every direction, i.e. when the difference is not a real
    Gaussian.
    """
    maj2, min2, pa2 = target
    maj1, min1, pa1 = source
    p2, p1 = np.deg2rad(pa2), np.deg2rad(pa1)

    alpha = ((maj2 * np.cos(p2)) ** 2 + (min2 * np.sin(p2)) ** 2
             - (maj1 * np.cos(p1)) ** 2 - (min1 * np.sin(p1)) ** 2)
    beta = ((maj2 * np.sin(p2)) ** 2 + (min2 * np.cos(p2)) ** 2
            - (maj1 * np.sin(p1)) ** 2 - (min1 * np.cos(p1)) ** 2)
    gamma = 2.0 * ((min2 ** 2 - maj2 ** 2) * np.sin(p2) * np.cos(p2)
                   - (min1 ** 2 - maj1 ** 2) * np.sin(p1) * np.cos(p1))

    s = alpha + beta
    t = np.hypot(alpha - beta, gamma)
    if s < t or s <= 0:
        return None
    maj = np.sqrt(0.5 * (s + t))
    bmin = np.sqrt(0.5 * (s - t))
    if abs(gamma) + abs(alpha - beta) == 0:
        pa = 0.0
    else:
        pa = 0.5 * np.rad2deg(np.arctan2(-gamma, alpha - beta))
    return float(maj), float(bmin), float(pa)


def _gaussian_kernel(bmaj_px, bmin_px, bpa_deg, truncate=4.0):
    """Unit-sum elliptical Gaussian kernel; FWHM inputs in pixels.

    `bpa` follows the FITS convention (degrees east of north), so the major
    axis runs along +y when bpa = 0.
    """
    smaj = bmaj_px / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    smin = bmin_px / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    half = int(np.ceil(truncate * max(smaj, smin)))
    half = max(half, 1)
    y, x = np.mgrid[-half:half + 1, -half:half + 1]
    th = np.deg2rad(bpa_deg)
    # rotate into the beam frame: major axis along +y at bpa = 0
    xr = x * np.cos(th) + y * np.sin(th)
    yr = -x * np.sin(th) + y * np.cos(th)
    k = np.exp(-0.5 * ((xr / smin) ** 2 + (yr / smaj) ** 2))
    return k / k.sum()


def _match_beams(ref, tgt, ref_beam, tgt_beam, cell, tol=0.05):
    """Bring both images to a common resolution before estimating a shift.

    No intensity-based estimator is unbiased when the two images have different
    PSFs: a resolved source convolved with a broader beam has its light
    redistributed asymmetrically wherever the source is asymmetric, which
    displaces the correlation peak. Measured degradation on NGC 7469 was 0.05 px
    at +2 px extra FWHM but 0.95 px at +8 px, so this matters well before the
    beams look obviously different.

    Both images are convolved to a common beam broad enough to be reachable
    from either, which is more robust than deconvolving one towards the other.
    Returns ``(ref, tgt, applied)``.
    """
    if ref_beam is None or tgt_beam is None or cell is None:
        return ref, tgt, False
    # already matched?
    if (abs(ref_beam[0] - tgt_beam[0]) <= tol * max(ref_beam[0], tgt_beam[0])
            and abs(ref_beam[1] - tgt_beam[1]) <= tol * max(ref_beam[1], tgt_beam[1])):
        return ref, tgt, False

    # a circular common beam, 2% broader than the broadest axis present
    common_ax = 1.02 * max(ref_beam[0], ref_beam[1], tgt_beam[0], tgt_beam[1])
    common = (common_ax, common_ax, 0.0)

    out = []
    for img, beam in ((ref, ref_beam), (tgt, tgt_beam)):
        k = _deconvolve_beams(common, beam)
        if k is None or max(k[0], k[1]) / cell < 0.5:
            out.append(img)          # already at (or above) the common beam
            continue
        kern = _gaussian_kernel(k[0] / cell, k[1] / cell, k[2])
        out.append(signal.fftconvolve(_clean(img), kern, mode='same'))
    return out[0], out[1], True


# ---------------------------------------------------------------------------
# masking
# ---------------------------------------------------------------------------
def _build_mask(ref, tgt, rms_ref, rms_tgt, nsigma=5.0, grow_px=6,
                mask=None):
    """Pixels that carry signal in *either* image, dilated by `grow_px`.

    The union rather than the intersection: a component that is bright at one
    frequency and faint at the other still constrains the shift, and dropping it
    would bias the fit towards whichever component survives the cut.

    A caller-supplied `mask` (e.g. from ``mlibs.mask_dilation``) is used as-is;
    this local fallback exists so the module stays importable without
    ``image_morphometry``.
    """
    if mask is not None:
        m = np.asarray(mask).astype(bool)
        if m.shape != ref.shape:
            raise ValueError(f'mask shape {m.shape} != image shape {ref.shape}')
        return m
    m = (ref > nsigma * rms_ref) | (tgt > nsigma * rms_tgt)
    if not m.any():                       # nothing above the cut -- back off
        m = (ref > 3.0 * rms_ref) | (tgt > 3.0 * rms_tgt)
    if m.any() and grow_px > 0:
        m = ndimage.binary_dilation(m, iterations=int(grow_px))
    return m


def _taper(mask, sigma=3.0):
    """Soft-edged mask, so the FFT path does not ring off a hard boundary."""
    return ndimage.gaussian_filter(mask.astype(float), sigma)


# ---------------------------------------------------------------------------
# shifting
# ---------------------------------------------------------------------------
def apply_pixel_shift(data, dy, dx):
    """Shift `data` by ``(dy, dx)`` pixels with a Fourier phase ramp.

    Flux conserving and free of the interpolation-kernel bias that
    ``scipy.ndimage.shift`` carries at sub-pixel offsets. NaNs are zero-filled
    before the transform and restored (shifted) afterwards.
    """
    data = np.asarray(data, dtype=float)
    if dy == 0 and dx == 0:
        return data.copy()
    bad = ~np.isfinite(data)
    filled = np.where(bad, 0.0, data)
    spec = ndimage.fourier_shift(np.fft.fft2(filled), (dy, dx))
    # For an even-length axis the Nyquist term picks up a phase factor that has
    # no conjugate partner, so `np.real` silently discards part of it and the
    # shift stops being invertible. Radio images are band-limited far below
    # Nyquist, so that term carries no signal -- dropping it costs nothing and
    # makes the operation exactly reversible.
    ny, nx = spec.shape
    if ny % 2 == 0:
        spec[ny // 2, :] = 0.0
    if nx % 2 == 0:
        spec[:, nx // 2] = 0.0
    out = np.real(np.fft.ifft2(spec))
    if bad.any():
        moved = ndimage.shift(bad.astype(float), (dy, dx), order=1,
                              mode='constant', cval=0.0) > 0.5
        out[moved] = np.nan
    return out


def split_shift(dy, dx):
    """Split a shift into an exact integer part and a sub-pixel remainder.

    The integer part is meant to be absorbed by moving a cutout's centre --
    which pulls in real neighbouring pixels instead of padding with zeros, and
    involves no interpolation at all -- leaving only the remainder, always
    ``<= 0.5`` px, to be applied as a phase ramp.

    Returns ``(dy_int, dx_int, dy_frac, dx_frac)``.
    """
    dy_i = int(np.round(dy))
    dx_i = int(np.round(dx))
    return dy_i, dx_i, float(dy - dy_i), float(dx - dx_i)


# ---------------------------------------------------------------------------
# backends -- each returns the *displacement* D of the target w.r.t. the
# reference; `estimate_image_shift` negates once, at the end, to get the
# correction.
# ---------------------------------------------------------------------------
def _coarse_xcorr(ref, tgt, mask, max_shift):
    """Integer-pixel un-whitened cross-correlation peak.

    Deliberately *not* phase-normalised (see the module docstring). Used only
    to seed the refiners, so that a source with several components cannot drop
    them into the wrong local optimum.
    """
    w = _taper(mask)
    r = _clean(ref) * w
    t = _clean(tgt) * w
    r = r - r.mean()
    t = t - t.mean()
    cc = np.fft.fftshift(np.real(np.fft.ifft2(
        np.fft.fft2(t) * np.conj(np.fft.fft2(r)))))
    ny, nx = cc.shape
    cy, cx = ny // 2, nx // 2
    b = int(min(max_shift, cy - 1, cx - 1))
    sub = cc[cy - b:cy + b + 1, cx - b:cx + b + 1]
    j, i = np.unravel_index(np.argmax(sub), sub.shape)
    return float(j - b), float(i - b)


def _starts(x0, max_shift=None):
    """Seed points for the refiners: the coarse peak, plus the origin.

    Deduplicated, and clipped to the search box when one is given.
    """
    cand = [np.asarray(x0, dtype=float), np.zeros(2)]
    if max_shift is not None:
        cand = [np.clip(c, -max_shift, max_shift) for c in cand]
    out = []
    for c in cand:
        if not any(np.allclose(c, o, atol=1e-6) for o in out):
            out.append(c)
    return out


def _chi2_cost_factory(ref, tgt, mask, rms_ref, rms_tgt):
    """chi2 of ``tgt ~= a * shift(ref, dy, dx) + b`` over `mask`.

    `a` and `b` are solved by linear least squares at every step, so only
    ``(dy, dx)`` are non-linear. `a` absorbs the global spectral index -- which
    is what lets a 5 GHz image be matched against a 33 GHz one at all -- and `b`
    absorbs any residual background offset.

    The weight is ``sqrt(rms_t^2 + a^2 rms_r^2)``: the model carries the
    reference's noise too, scaled by the fitted flux ratio. Without this a
    high-SNR reference paired with a low-SNR target reports an over-optimistic
    uncertainty.
    """
    w = mask
    y = tgt[w]
    ones = np.ones_like(y)

    def cost(p):
        dy, dx = float(p[0]), float(p[1])
        model = apply_pixel_shift(ref, dy, dx)[w]
        A = np.vstack([model, ones]).T
        coef, *_ = np.linalg.lstsq(A, y, rcond=None)
        a = coef[0]
        sig = np.sqrt(rms_tgt ** 2 + (a * rms_ref) ** 2)
        return float(np.sum(((y - A @ coef) / sig) ** 2))

    return cost


def _hessian_errors(cost, p, step=0.05):
    """1-sigma errors from the numerical Hessian of chi2 at ``delta chi2 = 1``."""
    H = np.zeros((2, 2))
    for i in range(2):
        for j in range(2):
            acc = 0.0
            for si, sj, sign in ((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)):
                q = np.array(p, dtype=float)
                q[i] += si * step
                q[j] += sj * step
                acc += sign * cost(q)
            H[i, j] = acc / (4.0 * step * step)
    try:
        cov = 2.0 * np.linalg.inv(H)
        d = np.diag(cov)
        if np.any(d <= 0):
            return np.nan, np.nan
        return float(np.sqrt(d[0])), float(np.sqrt(d[1]))
    except np.linalg.LinAlgError:
        return np.nan, np.nan


def _refine_chi2(ref, tgt, mask, rms_ref, rms_tgt, x0, max_shift):
    cost = _chi2_cost_factory(ref, tgt, mask, rms_ref, rms_tgt)
    # Two starts: the coarse cross-correlation peak, and the origin. The coarse
    # peak is usually the better seed, but it is exactly what goes wrong when
    # one component dominates in only one of the two bands -- in that case it
    # can land 10+ px away and a local optimiser will happily stay there.
    best = None
    for start in _starts(x0):
        res = optimize.minimize(cost, start, method='Nelder-Mead',
                                options=dict(xatol=1e-3, fatol=1e-3,
                                             maxiter=2000))
        if best is None or res.fun < best.fun:
            best = res
    dy, dx = float(best.x[0]), float(best.x[1])
    ey, ex = _hessian_errors(cost, best.x)
    return dy, dx, ey, ex, bool(best.success)


def _mi_cost_factory(ref, tgt, mask, rms_ref, rms_tgt, bins=48, smooth=1.0):
    """Negative mutual information of the two images' joint histogram.

    MI makes no assumption that one image is a scaled copy of the other -- only
    that they share structure. That is the one thing that survives a component
    whose brightness ratio *inverts* between bands, where every least-squares or
    correlation estimator is pulled towards whichever component dominates.
    Intensities are arcsinh-scaled in units of the rms so the histogram spans
    the dynamic range instead of piling into the lowest bin, and the joint
    histogram is Gaussian-smoothed to give the optimiser a differentiable
    surface.
    """
    R0 = np.arcsinh(_clean(ref) / rms_ref)
    T0 = np.arcsinh(_clean(tgt) / rms_tgt)
    tvals = T0[mask]
    tedges = np.linspace(np.nanpercentile(tvals, 0.5),
                         np.nanpercentile(tvals, 99.9), bins + 1)

    def neg_mi(p):
        Rs = apply_pixel_shift(R0, float(p[0]), float(p[1]))[mask]
        redges = np.linspace(np.nanpercentile(Rs, 0.5),
                             np.nanpercentile(Rs, 99.9), bins + 1)
        H, _, _ = np.histogram2d(Rs, tvals, bins=[redges, tedges])
        H = ndimage.gaussian_filter(H, smooth)
        tot = H.sum()
        if tot <= 0:
            return 0.0
        P = H / tot
        px = P.sum(1)[:, None]
        py = P.sum(0)[None, :]
        nz = P > 0
        return -float(np.sum(P[nz] * np.log(P[nz] / (px * py)[nz])))

    return neg_mi


def _refine_mi(ref, tgt, mask, rms_ref, rms_tgt, x0, max_shift,
               n_bootstrap=6, rng=None):
    neg_mi = _mi_cost_factory(ref, tgt, mask, rms_ref, rms_tgt)
    bounds = [(-max_shift, max_shift), (-max_shift, max_shift)]
    # Multi-start matters more here than for chi2: MI is the fallback precisely
    # for the case where the coarse seed is untrustworthy, so it must not be
    # anchored to it.
    res = None
    for start in _starts(x0, max_shift):
        r = optimize.minimize(neg_mi, start, method='Powell', bounds=bounds,
                              options=dict(xtol=1e-3, ftol=1e-5))
        if res is None or r.fun < res.fun:
            res = r
    dy, dx = float(res.x[0]), float(res.x[1])

    ey = ex = np.nan
    if n_bootstrap and n_bootstrap > 1:
        # MI has no chi2-like curvature scale, so the uncertainty comes from
        # re-fitting against independent noise realisations of the target.
        rng = np.random.default_rng(0) if rng is None else rng
        samples = []
        for _ in range(int(n_bootstrap)):
            noisy = tgt + rng.normal(0.0, rms_tgt, tgt.shape)
            f = _mi_cost_factory(ref, noisy, mask, rms_ref, rms_tgt)
            r = optimize.minimize(f, res.x, method='Powell', bounds=bounds,
                                  options=dict(xtol=1e-2, ftol=1e-4))
            samples.append(r.x)
        s = np.asarray(samples)
        ey, ex = float(np.std(s[:, 0], ddof=1)), float(np.std(s[:, 1], ddof=1))
    return dy, dx, ey, ex, bool(res.success)


def _refine_xcorr(ref, tgt, mask, rms_ref, rms_tgt, x0, max_shift,
                  upsample=100):
    """Upsampled *plain* cross-correlation via skimage.

    ``normalization=None`` is not a detail -- it is the whole point. See the
    module docstring.
    """
    try:
        from skimage.registration import phase_cross_correlation
    except ImportError as exc:                                # pragma: no cover
        raise ImportError(
            "method='xcorr' needs scikit-image (conda install scikit-image)."
        ) from exc
    w = _taper(mask)
    out = phase_cross_correlation(_clean(ref) * w, _clean(tgt) * w,
                                  upsample_factor=upsample,
                                  normalization=None)
    sh = out[0] if isinstance(out, tuple) else out
    # skimage returns the shift that registers `moving` onto `reference`,
    # i.e. the correction; we want the displacement, hence the sign flip.
    return -float(sh[0]), -float(sh[1]), np.nan, np.nan, True


def _refine_chi2_shift(ref, tgt, mask, rms_ref, rms_tgt, x0, max_shift):
    """Adam Ginsburg's ``image_registration.chi2_shift``.

    An optional external cross-check: an independently written chi2 registration
    with its own error model. Kept optional because a core cutout routine should
    not acquire a hard dependency.
    """
    try:
        from image_registration import chi2_shift
    except ImportError as exc:
        raise ImportError(
            "method='chi2_shift' needs the `image_registration` package:\n"
            "    pip install image_registration\n"
            "It is optional -- the built-in method='chi2' needs nothing extra."
        ) from exc
    w = _taper(mask)
    dx, dy, edx, edy = chi2_shift(_clean(ref) * w, _clean(tgt) * w,
                                  err=float(rms_tgt), return_error=True,
                                  upsample_factor='auto')
    # chi2_shift returns x first; its (dx, dy) is already the displacement of
    # image2 relative to image1, which is this module's internal convention, so
    # only the axis order changes here. Both are pinned by the injection test --
    # the sign is genuinely not what the upstream docstring implies.
    return float(dy), float(dx), float(edy), float(edx), True


_BACKENDS = {
    'chi2': _refine_chi2,
    'mi': _refine_mi,
    'xcorr': _refine_xcorr,
    'chi2_shift': _refine_chi2_shift,
}


# ---------------------------------------------------------------------------
# the public estimator
# ---------------------------------------------------------------------------
def estimate_image_shift(reference, target,
                         method='chi2',
                         mask=None,
                         rms_ref=None, rms_target=None,
                         match_beam=True,
                         ref_beam=None, target_beam=None, cell_size=None,
                         nsigma=5.0, grow_px=6,
                         max_shift=None,
                         n_sigma_reject=2.0,
                         mi_bootstrap=6,
                         flip_tolerance=1.0,
                         verbose=0):
    """Estimate the shift that brings `target` onto `reference`.

    Parameters
    ----------
    reference, target : str or 2D ndarray
        FITS paths or image planes. Paths are preferred: they let the beam, the
        pixel scale and the companion residual's noise be read automatically.
    method : {'chi2', 'mi', 'xcorr', 'chi2_shift', 'auto', 'ensemble'}
        ``'chi2'``   forward model with a free flux scale (default; best
                     general-purpose accuracy, good to SNR ~ 3);
        ``'mi'``     mutual information (the only one robust to a component
                     whose brightness ratio inverts between bands, but noisier);
        ``'xcorr'``  upsampled plain cross-correlation via scikit-image;
        ``'chi2_shift'`` external ``image_registration`` backend;
        ``'auto'``   run chi2 and mi, and fall back to mi when they disagree --
                     which is the signature of a component flip;
        ``'ensemble'`` run everything available, take the median, and use the
                     inter-method scatter as the uncertainty.
    mask : 2D bool ndarray, optional
        Pixels to fit over. Defaults to the union of the >`nsigma` regions of
        both images, dilated by `grow_px`. Pass ``mlibs.mask_dilation(...)[0]``
        here to reuse morphen's beam-aware masking.
    rms_ref, rms_target : float, optional
        Noise levels. Default: the companion ``-residual`` image if one sits
        next to the input, else ``mad_std`` of the image.
    match_beam : bool
        Convolve both images to a common resolution before estimating when
        their beams differ by more than 5%. Only affects the estimate; the
        returned shift applies to the native images.
    max_shift : float, optional
        Largest shift to search for and to accept, in pixels. Defaults to one
        reference beam major axis (or 15 px when no beam is known).
    n_sigma_reject : float
        A shift smaller than this many sigma is reported as zero, with a
        ``'consistent_with_zero'`` flag. This is what prevents a genuinely
        sub-pixel offset from being rounded up into a whole-pixel "correction"
        that leaves the images worse aligned than doing nothing. Set to 0 to
        always apply the raw estimate.
    flip_tolerance : float
        In ``'auto'``, the pixel disagreement between chi2 and mi beyond which
        the component-flip fallback triggers.

    Returns
    -------
    ShiftResult
        Unpacks as ``dy, dx`` for backwards compatibility.
    """
    method = _LEGACY_METHODS.get(method, method)
    if method not in _METHODS:
        raise ValueError(f'unknown method {method!r}; choose from {_METHODS}')

    ref_hdr, tgt_hdr = _header_of(reference), _header_of(target)
    ref = _as_plane(reference)
    tgt = _as_plane(target)
    if ref.shape != tgt.shape:
        raise ValueError(f'reference {ref.shape} and target {tgt.shape} must '
                         'be on the same pixel grid -- cut both to the same '
                         'centre and size first.')

    if rms_ref is None:
        rms_ref = _noise_of(reference, ref)
    if rms_target is None:
        rms_target = _noise_of(target, tgt)
    rms_ref = max(float(rms_ref), 1e-30)
    rms_target = max(float(rms_target), 1e-30)

    if ref_beam is None:
        ref_beam = _beam_of(ref_hdr)
    if target_beam is None:
        target_beam = _beam_of(tgt_hdr)
    if cell_size is None:
        cell_size = _cell_of(ref_hdr) or _cell_of(tgt_hdr)

    flags = []
    if max_shift is None:
        if ref_beam is not None and cell_size:
            max_shift = max(3.0, ref_beam[0] / cell_size)
        else:
            max_shift = 15.0
    max_shift = float(max_shift)

    # --- masking happens on the *native* images, before any beam matching, so
    #     the mask reflects where the real signal is.
    m = _build_mask(ref, tgt, rms_ref, rms_target, nsigma=nsigma,
                    grow_px=grow_px, mask=mask)
    if m.sum() < 16:
        flags.append('empty_mask')
        return ShiftResult(method=method, flags=flags,
                           diagnostics=dict(mask_px=int(m.sum())))

    ref_e, tgt_e, matched = ref, tgt, False
    if match_beam:
        ref_e, tgt_e, matched = _match_beams(ref, tgt, ref_beam, target_beam,
                                             cell_size)
        if matched:
            flags.append('beam_matched')

    x0 = _coarse_xcorr(ref_e, tgt_e, m, max_shift)

    diagnostics = dict(mask_px=int(m.sum()), rms_ref=rms_ref,
                       rms_target=rms_target, coarse=x0,
                       ref_beam=ref_beam, target_beam=target_beam,
                       cell_size=cell_size, max_shift=max_shift,
                       beam_matched=matched, per_method={})

    def _run(name):
        kw = {}
        if name == 'mi':
            kw['n_bootstrap'] = mi_bootstrap
        dy, dx, ey, ex, ok = _BACKENDS[name](ref_e, tgt_e, m, rms_ref,
                                             rms_target, x0, max_shift, **kw)
        diagnostics['per_method'][name] = dict(dy=dy, dx=dx, dy_err=ey,
                                               dx_err=ex, converged=ok)
        return dy, dx, ey, ex, ok

    used = method
    if method in _BACKENDS:
        dy, dx, ey, ex, ok = _run(method)
        if not ok:
            flags.append('not_converged')

    elif method == 'auto':
        c_dy, c_dx, c_ey, c_ex, c_ok = _run('chi2')
        m_dy, m_dx, m_ey, m_ex, m_ok = _run('mi')
        # The trigger must clear *both* estimators' noise, not just chi2's.
        # MI is the coarser statistic by a wide margin, so keying the threshold
        # off chi2's (very small) formal error alone made `auto` discard good
        # chi2 answers in favour of noisier MI ones whenever the two happened to
        # differ by half a pixel.
        def _n(v):
            return float(v) if np.isfinite(v) else 0.0
        scale = max(flip_tolerance,
                    3.0 * np.hypot(_n(c_ey), _n(c_ex)),
                    3.0 * np.hypot(_n(m_ey), _n(m_ex)))
        if np.hypot(c_dy - m_dy, c_dx - m_dx) > scale:
            # chi2 and MI only diverge like this when the intensity relation
            # between the bands is not a single global scaling -- a component
            # flip. MI is the one that survives that.
            flags.append('component_flip_suspected')
            dy, dx, ey, ex, ok, used = m_dy, m_dx, m_ey, m_ex, m_ok, 'mi'
        else:
            dy, dx, ey, ex, ok, used = c_dy, c_dx, c_ey, c_ex, c_ok, 'chi2'
        if not ok:
            flags.append('not_converged')

    else:   # ensemble
        got = []
        for name in ('chi2', 'mi', 'xcorr', 'chi2_shift'):
            try:
                r = _run(name)
                got.append((r[0], r[1]))
            except ImportError:
                continue
            except Exception as exc:                          # pragma: no cover
                warnings.warn(f'alignment backend {name!r} failed: {exc}')
        if not got:                                           # pragma: no cover
            raise RuntimeError('every alignment backend failed')
        arr = np.asarray(got)
        dy, dx = float(np.median(arr[:, 0])), float(np.median(arr[:, 1]))
        if len(arr) > 1:
            ey, ex = float(np.std(arr[:, 0], ddof=1)), float(np.std(arr[:, 1], ddof=1))
        else:
            ey = ex = np.nan
        ok, used = True, 'ensemble'

    # displacement -> correction
    dy, dx = -dy, -dx
    dy_raw, dx_raw = dy, dx

    # --- QA gate -----------------------------------------------------------
    if not (np.isfinite(ey) and np.isfinite(ex)):
        flags.append('no_uncertainty')
    elif n_sigma_reject > 0:
        if np.hypot(dy, dx) < n_sigma_reject * np.hypot(ey, ex):
            flags.append('consistent_with_zero')
            dy = dx = 0.0

    if np.hypot(dy, dx) > max_shift:
        flags.append('clamped')
        s = max_shift / np.hypot(dy, dx)
        dy, dx = dy * s, dx * s

    result = ShiftResult(dy=dy, dx=dx, dy_err=ey, dx_err=ex,
                         dy_raw=dy_raw, dx_raw=dx_raw, method=used,
                         flags=flags, diagnostics=diagnostics)
    if verbose >= 1:
        print(f'        > {result}')
    return result


# ---------------------------------------------------------------------------
# convenience: a whole band at once
# ---------------------------------------------------------------------------
def align_image_list(reference, images, **kwargs):
    """Estimate a shift for every image in `images` against one `reference`.

    Returns a list of :class:`ShiftResult`. Useful for the consistency check
    that matters most in practice: sub-bands from a single imaging run should
    all return the *same* shift, because the offset is a property of the
    dataset (array, self-cal solution), not of the sub-band. Scatter between
    them is a direct measure of the estimator's noise on your data.
    """
    return [estimate_image_shift(reference, im, **kwargs) for im in images]
