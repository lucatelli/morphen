"""
Visual inspection of astrometric alignment, in the notebook.

Companion to :mod:`image_alignment`: that module *measures and applies* the
shift, this one lets you look at the result. It replaces
``mlibs.imview(imagelist_c)`` -- i.e. ``casaviewer.imview``, which needs a GUI
and is therefore unusable on a remote notebook server.

The entry point is :func:`visualise_alignment`. It blinks through a list of
cutouts with the reference's contours frozen on top, and adds three static
diagnostic panels, then returns an :class:`AlignmentReport` so the check is
scriptable and not only visual.

Why every comparison here is amplitude-free
-------------------------------------------
The obvious panel to draw is an image difference, ``t - (a * r + b)`` with a
fitted flux scale ``a``. It is deliberately *not* here. A single scalar ``a``
only flattens the residual if the whole source shares one spectral index. In
multi-component sources -- exactly the case this pipeline targets, where one
component can inverts its spectrum relative to another -- no scalar works, and
the residual would be dominated by spectral structure masquerading as
misalignment.

Everything below therefore compares *positions*, never flux ratios:

* reference contours frozen across frames (contours mark where emission is,
  not how bright it is);
* each band contoured at its **own** sigma levels, so a band that is ten times
  fainter still contours correctly;
* intensity-weighted centroids per component, which do not move when a
  component changes brightness;
* 1-D cuts individually normalised to their own peak.

A component can flip brightness between bands without any of these moving.

Independence
------------
Like :mod:`image_alignment`, this module is **self-contained**: numpy, scipy,
astropy and matplotlib only. It can be imported directly for testing without
the ``mlibs`` preamble (no CASA, no JAX). It is listed in ``mlibs.py`` so
:func:`visualise_alignment` is also reachable as ``mlibs.visualise_alignment``.

The raster is drawn with the same convention as ``plotting.eimshow`` -- a grey
linear underlay at alpha 0.5 beneath an ``asinh`` layer with ``asinh_a=0.075``
-- so frames look like the rest of the package's figures, without inheriting
``eimshow``'s CASA-backed beam lookup or its forced ``projection='offset'``
when an axis is supplied.
"""

import os
import io
import gc
import base64
import uuid
import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy import ndimage
from astropy.io import fits
from astropy.stats import mad_std
from astropy.visualization import simple_norm

import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from matplotlib import gridspec

try:  # sibling module; present both standalone and inside the mlibs namespace
    from image_alignment import (estimate_image_shift, ShiftResult,
                                 _as_plane, _beam_of, _cell_of, _noise_of)
    _HAVE_ALIGN = True
except Exception:  # pragma: no cover - only when the sibling is missing
    _HAVE_ALIGN = False

__all__ = [
    'AlignmentReport',
    'visualise_alignment',
]

# Frequency ramp used for every "one colour per band" panel, so P2/P3/P4 agree.
_FREQ_CMAP = 'turbo'
# Colour of the frozen reference contours -- deliberately outside magma_r's range.
_REF_COLOUR = '#00d5ff'


# ---------------------------------------------------------------------------
# small header / data readers
#
# These duplicate a few lines of data_io/radio_utils on purpose: those files
# are only importable through the mlibs exec preamble, and this module has to
# stay standalone-importable for the test suite. Where the real helper *is*
# reachable (i.e. we are running inside mlibs) it is preferred, via _mlibs().
# ---------------------------------------------------------------------------
def _mlibs(name):
    """The repo's own helper `name` if we are running inside mlibs, else None."""
    fn = globals().get(name)
    return fn if callable(fn) else None


def _plane(image):
    """2D float array from a path or an array, WSClean 4D shapes included."""
    if _HAVE_ALIGN:
        return _as_plane(image)
    if isinstance(image, str):
        data = fits.getdata(image)
    else:
        data = np.asarray(image)
    data = np.squeeze(np.asarray(data, dtype=float))
    if data.ndim != 2:
        raise ValueError(f'expected a 2D image, got shape {data.shape}')
    return data


def _header(image):
    return fits.getheader(image) if isinstance(image, str) else None


def _beam_arcsec(header):
    """(bmaj, bmin, bpa) in arcsec/arcsec/deg, or None."""
    if header is None:
        return None
    beam = _beam_of(header) if _HAVE_ALIGN else None
    if beam is None:
        bmaj, bmin = header.get('BMAJ'), header.get('BMIN')
        if bmaj is None or bmin is None or not bmaj:
            return None
        beam = (float(bmaj), float(bmin), float(header.get('BPA', 0.0) or 0.0))
    return beam[0] * 3600.0, beam[1] * 3600.0, beam[2]


def _cell_arcsec(image, header):
    """Pixel scale in arcsec, or None."""
    get_cell_size = _mlibs('get_cell_size')
    if get_cell_size is not None and isinstance(image, str):
        try:
            cell = float(get_cell_size(image))
            if np.isfinite(cell) and cell > 0:
                return cell
        except Exception:
            pass
    if header is None:
        return None
    cell = _cell_of(header) if _HAVE_ALIGN else header.get('CDELT1')
    return abs(float(cell)) * 3600.0 if cell else None


def _freq_hz(image, header):
    """Observing frequency in Hz, or NaN."""
    getfreqs = _mlibs('getfreqs')
    if getfreqs is not None and isinstance(image, str):
        try:
            freq = float(np.atleast_1d(getfreqs([image]))[0])
            if np.isfinite(freq) and freq > 0:
                return freq
        except Exception:
            pass
    if header is None:
        return np.nan
    for i in (3, 4):
        if str(header.get(f'CTYPE{i}', '')).upper().startswith('FREQ'):
            val = header.get(f'CRVAL{i}')
            if val:
                unit = str(header.get(f'CUNIT{i}', 'Hz')).strip().lower()
                scale = {'hz': 1.0, 'khz': 1e3, 'mhz': 1e6, 'ghz': 1e9}.get(unit, 1.0)
                return float(val) * scale
    for key in ('RESTFRQ', 'RESTFREQ', 'FREQ'):
        if header.get(key):
            return float(header[key])
    return np.nan


def _rms_of(image, data, residual=None):
    """Noise from an explicit residual, else the WSClean companion, else mad_std."""
    if residual is not None:
        try:
            return float(mad_std(_plane(residual), ignore_nan=True))
        except Exception:
            pass
    if _HAVE_ALIGN:
        return _noise_of(image, data)
    return float(mad_std(data, ignore_nan=True))


def _applied_shift(header):
    """What t_cutout_2D_radec recorded, as a dict. Empty when never aligned."""
    if header is None or 'ASTRSHFY' not in header:
        return {}
    out = {'applied_dy': float(header.get('ASTRSHFY', np.nan)),
           'applied_dx': float(header.get('ASTRSHFX', np.nan)),
           # ASTRMETH is the *requested* mode, not the backend that won when
           # 'auto' fell through to 'mi' -- reported separately from the live
           # ShiftResult.method for exactly that reason.
           'applied_method': str(header.get('ASTRMETH', '')),
           'applied_ref': str(header.get('ASTRREF', '')),
           # Set by t_cutout_2D_radec when the estimate hit max_shift. The
           # applied number is then a bound, not a measurement, so the frame
           # is expected to still be off -- say so rather than let the reader
           # blame the estimator.
           'applied_clamped': bool(header.get('ASTRCLMP', False))}
    for key in ('CRVAL1', 'CRVAL2'):
        if key in header and ('O' + key) in header:
            out['d' + key] = (float(header[key]) - float(header['O' + key])) * 3600.0
    return out


# ---------------------------------------------------------------------------
# containers
# ---------------------------------------------------------------------------
@dataclass
class _Frame:
    """Everything read or measured for one image in the list."""
    image: object                 # path or array, as handed in
    name: str = ''
    data: np.ndarray = None
    header: object = None
    freq: float = np.nan          # Hz
    rms: float = np.nan
    peak: float = np.nan
    peak_yx: tuple = (0, 0)
    beam: tuple = None            # (bmaj, bmin, bpa) arcsec/arcsec/deg
    cell: float = None            # arcsec/px
    applied: dict = field(default_factory=dict)
    result: object = None         # ShiftResult, or None
    flags: list = field(default_factory=list)
    centroids: list = field(default_factory=list)   # per component (y, x, ey, ex)

    @property
    def freq_ghz(self):
        return self.freq / 1e9 if np.isfinite(self.freq) else np.nan

    @property
    def beam_px(self):
        """Beam major axis in pixels -- the yardstick a shift is judged against."""
        if self.beam and self.cell:
            return self.beam[0] / self.cell
        return np.nan


@dataclass
class AlignmentReport:
    """Result of :func:`visualise_alignment`.

    Attributes
    ----------
    table : pandas.DataFrame
        One row per image: frequency, beam, rms, the shift that was applied
        (from the header) and the residual shift still present (measured now).
    frames : list of bytes
        The rendered PNG of each blink frame, in display order.
    results : list
        The :class:`~image_alignment.ShiftResult` behind each residual row, or
        None where the measurement was skipped.
    reference : str
        The image everything was compared against.
    components : list
        Boolean masks of the components tracked in the centroid panel.
    summary_png : bytes or None
        The rendered static diagnostics figure (contour overlay, centroids,
        cuts). Bytes rather than a Figure because the figure is closed as soon
        as it is rasterised, to keep memory bounded.
    verdict : str
        ``'aligned'``, ``'marginal'``, ``'misaligned'`` or ``'not measured'``.
    """

    table: object = None
    frames: list = field(default_factory=list)
    results: list = field(default_factory=list)
    reference: str = ''
    components: list = field(default_factory=list)
    summary_png: object = None
    verdict: str = 'not measured'
    html: str = ''

    def __repr__(self):
        n = len(self.frames)
        return (f'AlignmentReport({n} frames, ref={os.path.basename(str(self.reference))!r}, '
                f'verdict={self.verdict!r})')


# ---------------------------------------------------------------------------
# components and centroids
# ---------------------------------------------------------------------------
def _components_from(data, rms, nsigma=5.0, min_area=1.0, n_components=None,
                     min_separation=None, verbose=0):
    """Deblend the reference's emission into components.

    Connected-component labelling alone is not enough here. Two components
    joined by a bridge of diffuse emission above the threshold come back as a
    *single* blob, and the centroid of a blob that contains two components
    moves when their brightness ratio changes -- which is exactly the failure
    mode the centroid panel is supposed to be immune to. So each connected
    region is watershed-split on its own local maxima, with `min_separation`
    (one beam by default) as the smallest resolvable separation.

    Returns boolean masks, brightest first.
    """
    d = np.nan_to_num(data, nan=0.0)
    detection = d > nsigma * rms
    if not detection.any():
        return []
    labels, nlab = ndimage.label(detection, structure=np.ones((3, 3)))

    split = None
    try:
        from skimage.feature import peak_local_max
        from skimage.segmentation import watershed
        sep = int(max(2, round(min_separation or 3)))
        peaks = peak_local_max(d, min_distance=sep,
                               threshold_abs=max(3.0, nsigma) * rms,
                               labels=detection)
        if len(peaks) > 1:
            markers = np.zeros(d.shape, dtype=int)
            markers[tuple(peaks.T)] = np.arange(1, len(peaks) + 1)
            split = watershed(-d, markers, mask=detection)
    except ImportError:  # pragma: no cover - scikit-image is a declared dep
        if verbose:
            warnings.warn('scikit-image is unavailable, so touching components '
                          'cannot be deblended; centroids of a merged blob do '
                          'respond to spectral-index changes')
    if split is not None:
        labels, nlab = split, int(split.max())

    comps = []
    for lab in range(1, nlab + 1):
        m = labels == lab
        if m.sum() < max(1.0, min_area):
            continue
        comps.append((float(np.nanmax(d[m])), m))
    comps.sort(key=lambda t: -t[0])
    if n_components:
        comps = comps[:int(n_components)]
    return [m for _, m in comps]


def _centroid(data, mask, rms, nsigma=3.0):
    """Intensity-weighted centroid inside `mask`, with first-order noise errors.

    For ``c = sum(I_j r_j) / sum(I_j)``, ``dc/dI_j = (r_j - c) / S``, so
    ``sigma_c^2 = rms^2 * sum((r_j - c)^2) / S^2``. Only pixels above
    ``nsigma * rms`` contribute, which keeps the faint tail from dragging the
    centroid around between bands of different depth.
    """
    d = np.nan_to_num(data, nan=0.0)
    sel = mask & (d > nsigma * rms)
    if sel.sum() < 3:
        return np.nan, np.nan, np.nan, np.nan
    yy, xx = np.nonzero(sel)
    w = d[sel].astype(float)
    total = w.sum()
    if not np.isfinite(total) or total <= 0:
        return np.nan, np.nan, np.nan, np.nan
    yc = float((w * yy).sum() / total)
    xc = float((w * xx).sum() / total)
    ey = float(rms * np.sqrt(np.sum((yy - yc) ** 2)) / total)
    ex = float(rms * np.sqrt(np.sum((xx - xc) ** 2)) / total)
    return yc, xc, ey, ex


# ---------------------------------------------------------------------------
# rendering primitives
# ---------------------------------------------------------------------------
def _view_slice(ref, box_size=None, center=None, nsigma=5.0, margin_beams=3.0):
    """The (row, col) slice every panel is drawn through.

    Cutouts are often much larger than the source, and a source occupying a
    tenth of the frame is useless for judging a sub-pixel offset. By default
    the view is the emission's bounding box grown by `margin_beams`; pass
    ``box_size='full'`` to disable, or an integer for a fixed box.
    """
    ny, nx = ref.data.shape
    if isinstance(box_size, str) and box_size == 'full':
        return (slice(0, ny), slice(0, nx))
    cy, cx = center if center is not None else ref.peak_yx
    if box_size is None:
        det = np.nan_to_num(ref.data, nan=0.0) > nsigma * ref.rms
        beam_px = ref.beam_px if np.isfinite(ref.beam_px) else 4.0
        if det.any():
            ys, xs = np.nonzero(det)
            pad = margin_beams * beam_px
            half = max(ys.max() - ys.min(), xs.max() - xs.min()) / 2.0 + pad
            cy, cx = (ys.min() + ys.max()) / 2.0, (xs.min() + xs.max()) / 2.0
        else:
            half = min(ny, nx) / 2.0
        box_size = 2 * half
    half = max(8.0, float(box_size) / 2.0)
    y0 = int(max(0, round(cy - half)))
    y1 = int(min(ny, round(cy + half)))
    x0 = int(max(0, round(cx - half)))
    x1 = int(min(nx, round(cx + half)))
    if y1 - y0 < 4 or x1 - x0 < 4:
        return (slice(0, ny), slice(0, nx))
    return (slice(y0, y1), slice(x0, x1))


def _cut(data, view):
    return data if view is None else data[view]


def _extent_of(data, view=None):
    """Pixel-offset extent, matching eimshow's array-centred convention.

    When a `view` is given the extent still refers to the *full* frame's
    centre, so a crosshair or centroid marker computed on the uncropped array
    lands in the right place.
    """
    ny, nx = data.shape
    if view is None:
        return [-nx / 2.0, nx / 2.0, -ny / 2.0, ny / 2.0]
    ys, xs = view
    return [xs.start - nx / 2.0, xs.stop - nx / 2.0,
            ys.start - ny / 2.0, ys.stop - ny / 2.0]


def _to_offset(y, x, shape):
    """(row, col) -> the offset coordinates used by _extent_of."""
    return x - shape[1] / 2.0, y - shape[0] / 2.0


def _limits(frame, scaling, ref, vmin_factor=3.0, vmax_factor=0.5):
    """(vmin, vmax) for one frame under the requested scaling policy."""
    if scaling == 'absolute':
        # One physical scale for everything. Honest about real flux
        # differences, but a steep spectrum leaves the high-frequency frames
        # nearly blank -- useful mostly within a single band.
        return vmin_factor * ref.rms, vmax_factor * ref.peak
    if scaling == 'peak':
        # eimshow's default: every frame at its own best contrast. Easiest to read
        # per frame, but the apparent size of the source changes between frames
        # purely from rescaling, which the eye can read as motion.
        return vmin_factor * frame.rms, vmax_factor * frame.peak
    # 'sigma' (default): identical dynamic range *in units of each frame's own
    # noise*, so the noise floor looks the same everywhere and any change in
    # the emission's apparent extent is real.
    k = vmax_factor * ref.peak / ref.rms if ref.rms > 0 else 10.0
    return vmin_factor * frame.rms, k * frame.rms


def _raster(ax, data, vmin, vmax, cmap='magma_r', extent=None, asinh_a=0.075):
    """The package's standard two-layer image render (cf. plotting.eimshow)."""
    g = np.nan_to_num(data, nan=0.0)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        norm_lin = simple_norm(g, stretch='linear', max_percent=99.0)
        norm = simple_norm(g, stretch='asinh', asinh_a=asinh_a,
                           vmin=vmin, vmax=max(vmax, vmin * 1.001))
    ax.imshow(g, origin='lower', cmap='gray', norm=norm_lin, alpha=0.5,
              extent=extent, aspect='equal')
    return ax.imshow(g, origin='lower', cmap=cmap, norm=norm, alpha=1.0,
                     extent=extent, aspect='equal')


def _draw_beam(ax, beam, cell, extent, colour='white'):
    """Beam ellipse in the lower-left corner. Pure astropy -- no CASA imhead."""
    if not beam or not cell:
        return
    bmaj, bmin, bpa = beam
    x0, x1, y0, y1 = extent
    pos = (x0 + 0.12 * (x1 - x0), y0 + 0.12 * (y1 - y0))
    el = Ellipse(pos, bmin / cell, bmaj / cell, angle=bpa,
                 facecolor=colour, edgecolor='0.3', lw=0.8, alpha=0.9)
    ax.add_artist(el)
    el.set_clip_box(ax.bbox)


def _ref_levels(ref, floor_sigma=5.0):
    """Geometric contour levels on the reference, as eimshow builds them.

    `floor_sigma` is the same threshold the masks and components use, so
    raising `sigma_mask` to keep sidelobes out of the analysis also keeps them
    out of the frozen contours -- otherwise the panel would still draw a ring
    of sidelobe contours the report deliberately ignores.
    """
    top = 0.9 * float(np.nanmax(ref.data))
    bottom = float(floor_sigma) * ref.rms
    if not np.isfinite(top) or top <= bottom:
        return np.asarray([bottom])
    return np.geomspace(bottom, top, 6)


def _fig_to_png(fig, dpi=110):
    """Rasterise and drop the figure -- the memory-safe idiom used repo-wide.

    eimshow's own docstring warns that keeping many axes alive on large images
    fills RAM, so no figure survives this call.
    """
    buf = io.BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight',
                facecolor=fig.get_facecolor())
    buf.seek(0)
    png = buf.read()
    buf.close()
    fig.clf()
    plt.close(fig)
    gc.collect()
    return png


# ---------------------------------------------------------------------------
# P1: the blinking panel -- image + frozen reference contours + readout
# ---------------------------------------------------------------------------
def _verdict_of(frame):
    """('aligned'|'marginal'|'misaligned'|'not measured', colour) for one frame.

    Judged against the beam, not against a pixel count: 0.3 px means something
    quite different at 6 px/beam than at 2 px/beam.
    """
    res = frame.result
    if res is None:
        return 'not measured', '0.45'
    mag = float(np.hypot(res.dy, res.dx))
    err = float(np.hypot(np.nan_to_num(res.dy_err), np.nan_to_num(res.dx_err)))
    beam_px = frame.beam_px
    tight = max(0.1, 2.0 * err)
    loose = 0.3 if not np.isfinite(beam_px) else max(0.3, 0.05 * beam_px)
    if mag < tight:
        return 'aligned', '#1a9850'
    if mag < loose:
        return 'marginal', '#e08214'
    return 'misaligned', '#d73027'


def _readout_lines(frame, ref):
    """The text block shown beside each frame.

    Deliberately no filename: WSClean names run to ~90 characters, and this
    block is monospace and unwrapped, so the name alone set the width of the
    whole figure. It is on the axes title (truncated) and in `report.table`.
    """
    lines = []
    if np.isfinite(frame.freq_ghz):
        lines.append(f'nu       = {frame.freq_ghz:.3f} GHz')
    if frame.beam:
        lines.append(f'beam     = {frame.beam[0]:.3f}" x {frame.beam[1]:.3f}"'
                     f' @ {frame.beam[2]:.1f} deg')
    if frame.cell:
        lines.append(f'cell     = {frame.cell:.4f}"/px'
                     + (f'  ({frame.beam_px:.2f} px/beam)'
                        if np.isfinite(frame.beam_px) else ''))
    lines.append(f'peak     = {frame.peak:.4g}')
    lines.append(f'rms      = {frame.rms:.4g}   (S/N {frame.peak / frame.rms:.0f})'
                 if frame.rms > 0 else f'rms      = {frame.rms:.4g}')
    lines.append('')

    ap = frame.applied
    if ap:
        lines.append('-- applied by t_cutout_2D_radec --')
        lines.append(f'shift    = ({ap["applied_dy"]:+.3f}, {ap["applied_dx"]:+.3f}) px'
                     + ('  CLAMPED' if ap.get('applied_clamped') else ''))
        if frame.cell:
            mas = np.hypot(ap['applied_dy'], ap['applied_dx']) * frame.cell * 1e3
            lines.append(f'         = {mas:.2f} mas')
        lines.append(f'mode     = {ap["applied_method"]}  (requested)')
        if ap.get('applied_ref'):
            lines.append(f'ref      = {ap["applied_ref"][:34]}')
    else:
        lines.append('-- no ASTRSHF* header: not aligned --')
    lines.append('')

    res = frame.result
    lines.append('-- residual, measured now --')
    if res is None:
        lines.append('not measured')
        if frame.flags:
            lines.append('flags    = ' + ', '.join(frame.flags))
    else:
        lines.append(f'shift    = ({res.dy:+.3f}, {res.dx:+.3f})'
                     f' +/- ({res.dy_err:.3f}, {res.dx_err:.3f}) px')
        mag = float(np.hypot(res.dy, res.dx))
        if frame.cell:
            lines.append(f'         = {mag * frame.cell * 1e3:.2f} mas')
        if np.isfinite(frame.beam_px) and frame.beam_px > 0:
            lines.append(f'         = {mag / frame.beam_px:.3f} beam')
        lines.append(f'backend  = {res.method}  (actually used)')
        if res.flags:
            lines.append('flags    = ' + ', '.join(res.flags))
    return lines


def _render_frame(frame, ref, scaling='sigma', cmap='magma_r',
                  show_ref_contours=True, show_crosshair=True,
                  show_components=True, view=None, figsize=(11, 5),
                  title=None, floor_sigma=5.0):
    """One blink frame: raster + frozen reference overlay + text readout."""
    fig = plt.figure(figsize=figsize)
    # 0.62 rather than 0.85: with the filename gone the widest readout line is
    # ~50 characters, and the spare column was padding the whole figure.
    gs = gridspec.GridSpec(1, 2, width_ratios=[1.0, 0.62], wspace=0.02,
                           figure=fig)
    ax = fig.add_subplot(gs[0, 0])
    ax_txt = fig.add_subplot(gs[0, 1])
    ax_txt.axis('off')

    fview = view if (view is None or frame.data.shape == ref.data.shape) else None
    extent = _extent_of(frame.data, fview)
    vmin, vmax = _limits(frame, scaling, ref)
    _raster(ax, _cut(frame.data, fview), vmin, vmax, cmap=cmap, extent=extent)

    if show_ref_contours:
        # The *reference's* contours, identical in every frame. Anything that
        # moves relative to them is the thing we are looking for.
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ax.contour(_cut(ref.data, fview),
                       levels=_ref_levels(ref, floor_sigma),
                       colors=_REF_COLOUR, linewidths=0.8, alpha=0.9,
                       extent=_extent_of(ref.data, fview))
    if show_crosshair:
        cx, cy = _to_offset(*ref.peak_yx, shape=ref.data.shape)
        ax.axhline(cy, color=_REF_COLOUR, lw=0.5, ls=':', alpha=0.7)
        ax.axvline(cx, color=_REF_COLOUR, lw=0.5, ls=':', alpha=0.7)
    if show_components and ref.centroids:
        for k, (yc, xc, _, _) in enumerate(ref.centroids):
            if not np.isfinite(yc):
                continue
            px, py = _to_offset(yc, xc, ref.data.shape)
            ax.plot(px, py, marker='+', ms=9, mew=1.4, color=_REF_COLOUR)
            ax.annotate(f'C{k + 1}', (px, py), textcoords='offset points',
                        xytext=(6, 4), color=_REF_COLOUR, fontsize=9)

    _draw_beam(ax, frame.beam, frame.cell, extent,
               colour='black' if '_r' in cmap else 'white')

    if frame.cell:
        ticks = ax.get_xticks()
        ax.set_xticks(ticks)
        ax.set_xticklabels([f'{t * frame.cell:.2f}' for t in ticks], fontsize=8)
        ticks = ax.get_yticks()
        ax.set_yticks(ticks)
        ax.set_yticklabels([f'{t * frame.cell:.2f}' for t in ticks], fontsize=8)
        ax.set_xlabel('offset [arcsec]', fontsize=9)
        ax.set_ylabel('offset [arcsec]', fontsize=9)
        ax.set_xlim(extent[0], extent[1])
        ax.set_ylim(extent[2], extent[3])
    else:
        ax.set_xlabel('offset [px]', fontsize=9)
        ax.set_ylabel('offset [px]', fontsize=9)
    ax.grid(False)

    verdict, colour = _verdict_of(frame)
    ax.set_title(title or os.path.basename(str(frame.name))[:52], fontsize=9)

    ax_txt.text(0.0, 1.0, '\n'.join(_readout_lines(frame, ref)),
                transform=ax_txt.transAxes, va='top', ha='left',
                fontsize=8.5, family='monospace')
    ax_txt.text(0.0, 0.0, f'  {verdict.upper()}  ', transform=ax_txt.transAxes,
                va='bottom', ha='left', fontsize=11, family='monospace',
                color='white', weight='bold',
                bbox=dict(boxstyle='round,pad=0.4', facecolor=colour, lw=0))
    return fig


# ---------------------------------------------------------------------------
# P2: per-sigma contour overlay, every band on one panel
# ---------------------------------------------------------------------------
def _panel_overlay(ax, frames, sigmas=(5, 20), view=None,
                   cmap=_FREQ_CMAP):
    """Each band contoured at *its own* sigma levels, coloured by frequency.

    Scaling every image by its own noise is what makes this amplitude-free: a
    band ten times fainter, or a component that inverts its spectrum, still
    contours at the same significance.
    """
    colours = _freq_colours(frames, cmap)
    ref_shape = frames[0].data.shape
    extent = _extent_of(frames[0].data, view)
    for frame, colour in zip(frames, colours):
        if frame.data.shape != ref_shape:
            continue
        levels = np.sort(np.asarray(sigmas, dtype=float)) * frame.rms
        levels = levels[levels < np.nanmax(frame.data)]
        if levels.size == 0:
            continue
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ax.contour(_cut(frame.data, view), levels=levels, colors=[colour],
                       linewidths=1.0, extent=extent)
    sm = plt.cm.ScalarMappable(cmap=cmap)
    freqs = [f.freq_ghz for f in frames]
    if np.all(np.isfinite(freqs)) and len(set(freqs)) > 1:
        sm.set_clim(min(freqs), max(freqs))
        ax.figure.colorbar(sm, ax=ax, fraction=0.046, pad=0.03,
                           label=r'$\nu$ [GHz]')
    ax.set_aspect('equal')
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_title(f'contours at {", ".join(str(s) for s in sigmas)}'
                 r'$\sigma$ of each band', fontsize=9)
    ax.set_xlabel('offset [px]', fontsize=8)
    ax.set_ylabel('offset [px]', fontsize=8)
    ax.grid(False)


def _freq_colours(frames, cmap=_FREQ_CMAP):
    """One colour per frame along a frequency ramp (falls back to index order)."""
    freqs = np.asarray([f.freq_ghz for f in frames], dtype=float)
    usable = freqs.size and np.all(np.isfinite(freqs))
    span = (freqs.max() - freqs.min()) if usable else 0.0
    if usable and span > 0:
        t = (freqs - freqs.min()) / span
    else:
        t = np.linspace(0, 1, len(frames))
    return [plt.get_cmap(cmap)(v) for v in t]


# ---------------------------------------------------------------------------
# P3: component centroid vs frequency
# ---------------------------------------------------------------------------
def _panel_centroids(ax_x, ax_y, frames, ref):
    """Per-component centroid drift with frequency, relative to the reference.

    Immune to a brightness flip: a component can invert its spectrum and its
    centroid does not move. Flat tracks at zero mean the bands agree.

    Caveat, printed rather than hidden: a genuine frequency-dependent core
    shift (synchrotron opacity) is *indistinguishable* from bad alignment on
    this panel. It surfaces the motion; deciding which it is stays with you.
    """
    n_comp = len(ref.centroids)
    if n_comp == 0:
        for ax in (ax_x, ax_y):
            ax.text(0.5, 0.5, 'no components found', ha='center', va='center',
                    transform=ax.transAxes, fontsize=9, color='0.4')
        return
    markers = ['o', 's', '^', 'D', 'v', 'P', 'X', '*']
    cell_mas = (ref.cell * 1e3) if ref.cell else None
    unit = 'mas' if cell_mas else 'px'
    scale = cell_mas if cell_mas else 1.0

    freqs = np.asarray([f.freq_ghz for f in frames], dtype=float)
    xvals = freqs if np.all(np.isfinite(freqs)) else np.arange(len(frames), dtype=float)
    xlabel = r'$\nu$ [GHz]' if np.all(np.isfinite(freqs)) else 'frame'

    for k in range(n_comp):
        y0, x0, _, _ = ref.centroids[k]
        dxs, dys, exs, eys = [], [], [], []
        for frame in frames:
            if k < len(frame.centroids):
                yc, xc, ey, ex = frame.centroids[k]
            else:
                yc = xc = ey = ex = np.nan
            dxs.append((xc - x0) * scale)
            dys.append((yc - y0) * scale)
            exs.append(ex * scale)
            eys.append(ey * scale)
        colour = plt.get_cmap('tab10')(k % 10)
        mk = markers[k % len(markers)]
        ax_x.errorbar(xvals, dxs, yerr=exs, marker=mk, ms=4, lw=1.0,
                      color=colour, label=f'C{k + 1}', capsize=2)
        ax_y.errorbar(xvals, dys, yerr=eys, marker=mk, ms=4, lw=1.0,
                      color=colour, label=f'C{k + 1}', capsize=2)

    for ax, lab in ((ax_x, r'$\Delta x$'), (ax_y, r'$\Delta y$')):
        ax.axhline(0.0, color='0.5', lw=0.8, ls='--')
        ax.set_ylabel(f'{lab} [{unit}]', fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(True, ls=':', alpha=0.4)
    if ref.beam and ref.cell:
        # A tenth of a beam is the scale at which misalignment starts to matter
        # for a spectral-index map.
        tenth = 0.1 * ref.beam[0] * (1e3 if cell_mas else 1.0 / ref.cell)
        for ax in (ax_x, ax_y):
            ax.axhspan(-tenth, tenth, color='0.8', alpha=0.35, zorder=0)
    ax_x.set_title('component centroid vs frequency', fontsize=9)
    ax_x.legend(fontsize=7, ncol=min(n_comp, 4), frameon=False)
    ax_y.set_xlabel(xlabel, fontsize=8)


# ---------------------------------------------------------------------------
# P4: normalised 1-D cuts
# ---------------------------------------------------------------------------
def _panel_cuts(ax_h, ax_v, frames, ref, cut_width=None, view=None,
                cmap=_FREQ_CMAP):
    """Cuts through the reference peak, each normalised to its own maximum.

    Normalising per profile is what keeps this amplitude-free; a residual shift
    then reads directly as a horizontal displacement of the peaks. Averaging
    over ``cut_width`` rows (one beam by default) keeps the faint bands usable.
    """
    colours = _freq_colours(frames, cmap)
    yc, xc = ref.peak_yx
    ny, nx = ref.data.shape
    xlim = None
    if view is not None:
        ys, xs = view
        cell_v = ref.cell or 1.0
        xlim = ((xs.start - nx / 2.0) * cell_v, (xs.stop - nx / 2.0) * cell_v)
    if cut_width is None:
        beam_px = ref.beam_px
        cut_width = int(max(1, round(beam_px))) if np.isfinite(beam_px) else 1
    half = max(0, int(cut_width) // 2)

    cell = ref.cell or 1.0
    unit = 'arcsec' if ref.cell else 'px'
    xs = (np.arange(nx) - nx / 2.0) * cell
    ys = (np.arange(ny) - ny / 2.0) * cell

    for frame, colour in zip(frames, colours):
        d = np.nan_to_num(frame.data, nan=0.0)
        row = d[max(0, yc - half): min(ny, yc + half + 1), :].mean(axis=0)
        col = d[:, max(0, xc - half): min(nx, xc + half + 1)].mean(axis=1)
        label = (f'{frame.freq_ghz:.1f} GHz' if np.isfinite(frame.freq_ghz)
                 else os.path.basename(str(frame.name))[:14])
        for ax, prof, axis in ((ax_h, row, xs), (ax_v, col, ys)):
            top = np.nanmax(prof)
            if not np.isfinite(top) or top <= 0:
                continue
            ax.plot(axis, prof / top, lw=1.0, color=colour, label=label)

    for ax, lab in ((ax_h, 'horizontal'), (ax_v, 'vertical')):
        if xlim:
            ax.set_xlim(*xlim)
        ax.axvline(0.0, color='0.5', lw=0.8, ls='--')
        ax.set_ylabel(f'{lab} cut\n(norm.)', fontsize=8)
        ax.tick_params(labelsize=8)
        ax.grid(True, ls=':', alpha=0.4)
    ax_h.set_title(f'cuts through the reference peak ({cut_width} px wide)',
                   fontsize=9)
    ax_h.legend(fontsize=6.5, ncol=2, frameon=False)
    ax_v.set_xlabel(f'offset [{unit}]', fontsize=8)


def _summary_figure(frames, ref, panels, overlay_sigmas, cut_width, view=None,
                    figsize=None):
    """The static diagnostics figure: P2 | P3 | P4, whichever were requested."""
    cols = [p for p in ('overlay', 'centroids', 'cuts') if p in panels]
    if not cols:
        return None
    widths = {'overlay': 1.15, 'centroids': 1.0, 'cuts': 1.0}
    if figsize is None:
        figsize = (5.2 * len(cols), 6.4)
    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, len(cols), figure=fig,
                           width_ratios=[widths[c] for c in cols],
                           hspace=0.12, wspace=0.32)
    for j, col in enumerate(cols):
        if col == 'overlay':
            ax = fig.add_subplot(gs[:, j])
            _panel_overlay(ax, frames, sigmas=overlay_sigmas, view=view)
        elif col == 'centroids':
            ax_x = fig.add_subplot(gs[0, j])
            ax_y = fig.add_subplot(gs[1, j], sharex=ax_x)
            plt.setp(ax_x.get_xticklabels(), visible=False)
            _panel_centroids(ax_x, ax_y, frames, ref)
        else:
            ax_h = fig.add_subplot(gs[0, j])
            ax_v = fig.add_subplot(gs[1, j])
            _panel_cuts(ax_h, ax_v, frames, ref, cut_width=cut_width, view=view)
    fig.suptitle(f'alignment diagnostics -- reference: '
                 f'{os.path.basename(str(ref.name))[:60]}', fontsize=10)
    return fig


# ---------------------------------------------------------------------------
# grid (contact sheet) mode
# ---------------------------------------------------------------------------
def _grid_figure(frames, ref, scaling='sigma', cmap='magma_r', view=None,
                 ncols=None, panel_size=2.9, floor_sigma=5.0):
    """Static contact sheet -- no JavaScript, no widgets, works anywhere."""
    n = len(frames)
    ncols = ncols or int(min(4, max(1, np.ceil(np.sqrt(n)))))
    nrows = int(np.ceil(n / ncols))
    fig, axes = plt.subplots(nrows, ncols, squeeze=False,
                             figsize=(panel_size * ncols, panel_size * nrows))
    extent = _extent_of(ref.data, view)
    levels = _ref_levels(ref, floor_sigma)
    for k, ax in enumerate(axes.ravel()):
        if k >= n:
            ax.axis('off')
            continue
        frame = frames[k]
        fview = view if frame.data.shape == ref.data.shape else None
        vmin, vmax = _limits(frame, scaling, ref)
        _raster(ax, _cut(frame.data, fview), vmin, vmax, cmap=cmap,
                extent=_extent_of(frame.data, fview))
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ax.contour(_cut(ref.data, view), levels=levels, colors=_REF_COLOUR,
                       linewidths=0.6, alpha=0.9, extent=extent)
        verdict, colour = _verdict_of(frame)
        label = (f'{frame.freq_ghz:.2f} GHz' if np.isfinite(frame.freq_ghz)
                 else os.path.basename(str(frame.name))[:18])
        if frame.result is not None:
            label += f'\n({frame.result.dy:+.2f}, {frame.result.dx:+.2f}) px'
        ax.set_title(label, fontsize=8, color=colour)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.grid(False)
    fig.suptitle('blink frames -- contours are the reference, identical in every panel',
                 fontsize=9)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# display backends
# ---------------------------------------------------------------------------
def _display_backend(requested='auto', verbose=1):
    """Resolve `mode`. Only Python-side availability can be checked here."""
    if requested != 'auto':
        return requested
    try:
        from IPython import get_ipython
        shell = get_ipython()
        shell_name = shell.__class__.__name__ if shell is not None else ''
    except Exception:
        shell_name = ''
    if shell_name != 'ZMQInteractiveShell':
        return 'grid'
    try:
        import ipywidgets  # noqa: F401
        return 'slider'
    except ImportError:
        return 'animate'


def _html_player(pngs, labels, fps=2):
    """A self-contained base64 player over already-rendered frames.

    Preferred over ``FuncAnimation.to_jshtml()`` because the frames exist
    already: no second rendering path, so the slider and the animation are
    guaranteed to show identical pixels. Pure JS, no ffmpeg, and the output is
    baked into the .ipynb so it survives save, restart and nbconvert. The first
    frame is a plain <img>, so something is visible even if scripts are blocked.
    """
    uid = 'mviz' + uuid.uuid4().hex[:10]
    srcs = ['data:image/png;base64,' + base64.b64encode(p).decode('ascii')
            for p in pngs]
    import json
    return f"""
<div id="{uid}" style="font-family:monospace;font-size:12px;max-width:100%">
  <img id="{uid}-img" src="{srcs[0]}" style="max-width:100%;display:block"/>
  <div style="margin-top:4px;display:flex;align-items:center;gap:6px;flex-wrap:wrap">
    <button id="{uid}-prev">&#9664;</button>
    <button id="{uid}-play">&#9654; play</button>
    <button id="{uid}-next">&#9654;</button>
    <input id="{uid}-range" type="range" min="0" max="{len(srcs) - 1}" value="0"
           step="1" style="flex:1;min-width:120px"/>
    <span id="{uid}-label"></span>
    <label>fps <input id="{uid}-fps" type="number" min="1" max="20"
           value="{int(fps)}" style="width:3.5em"/></label>
  </div>
</div>
<script>
(function() {{
  var srcs = {json.dumps(srcs)};
  var labels = {json.dumps(labels)};
  var root = document.getElementById("{uid}");
  var img = document.getElementById("{uid}-img");
  var rng = document.getElementById("{uid}-range");
  var lab = document.getElementById("{uid}-label");
  var play = document.getElementById("{uid}-play");
  var fpsEl = document.getElementById("{uid}-fps");
  var i = 0, timer = null;
  function show(k) {{
    i = ((k % srcs.length) + srcs.length) % srcs.length;
    img.src = srcs[i]; rng.value = i;
    lab.textContent = (i + 1) + "/" + srcs.length + "  " + labels[i];
  }}
  function stop() {{ if (timer) {{ clearInterval(timer); timer = null; }}
                     play.innerHTML = "&#9654; play"; }}
  function start() {{
    stop();
    var f = Math.max(1, Math.min(20, parseInt(fpsEl.value) || 2));
    timer = setInterval(function() {{ show(i + 1); }}, 1000 / f);
    play.innerHTML = "&#10074;&#10074; pause";
  }}
  play.addEventListener("click", function() {{ timer ? stop() : start(); }});
  document.getElementById("{uid}-prev").addEventListener(
      "click", function() {{ stop(); show(i - 1); }});
  document.getElementById("{uid}-next").addEventListener(
      "click", function() {{ stop(); show(i + 1); }});
  rng.addEventListener("input", function() {{ stop(); show(parseInt(rng.value)); }});
  fpsEl.addEventListener("change", function() {{ if (timer) start(); }});
  show(0);
}})();
</script>
"""


def _show_slider(state, labels, rerender, scaling, fps=2):
    """ipywidgets Play/IntSlider over cached PNGs, plus a live scaling toggle."""
    import ipywidgets as widgets
    from IPython.display import display as _display

    n = len(state['frames'])
    img = widgets.Image(value=state['frames'][0], format='png',
                        layout=widgets.Layout(max_width='100%'))
    slider = widgets.IntSlider(value=0, min=0, max=max(0, n - 1), step=1,
                               description='frame', continuous_update=True)
    play = widgets.Play(value=0, min=0, max=max(0, n - 1), step=1,
                        interval=int(1000 / max(1, fps)), description='play')
    widgets.jslink((play, 'value'), (slider, 'value'))
    scale_dd = widgets.Dropdown(options=['sigma', 'peak', 'absolute'],
                                value=scaling, description='scale',
                                layout=widgets.Layout(width='190px'))
    contours = widgets.Checkbox(value=True, description='ref contours',
                                indent=False)
    crosshair = widgets.Checkbox(value=True, description='crosshair',
                                 indent=False)
    caption = widgets.HTML()

    def _refresh(idx=None):
        idx = slider.value if idx is None else idx
        img.value = state['frames'][idx]
        caption.value = (f'<code>{idx + 1}/{n} &nbsp; {labels[idx]}</code>')

    def _on_frame(change):
        _refresh(change['new'])

    def _on_option(_change):
        # Re-render and re-cache; scrubbing afterwards is a byte swap again.
        state['frames'] = rerender(scale_dd.value, contours.value, crosshair.value)
        _refresh()

    slider.observe(_on_frame, names='value')
    for w in (scale_dd, contours, crosshair):
        w.observe(_on_option, names='value')

    _display(widgets.VBox([
        widgets.HBox([play, slider]),
        widgets.HBox([scale_dd, contours, crosshair]),
        caption, img,
    ]))
    _refresh(0)


# ---------------------------------------------------------------------------
# assembly
# ---------------------------------------------------------------------------
def _build_frame(image, residual=None):
    """Read one image into a :class:`_Frame`."""
    data = _plane(image)
    header = _header(image)
    rms = _rms_of(image, data, residual)
    peak_yx = np.unravel_index(int(np.nanargmax(np.nan_to_num(data, nan=-np.inf))),
                               data.shape)
    return _Frame(image=image,
                  name=image if isinstance(image, str) else 'array',
                  data=data, header=header,
                  freq=_freq_hz(image, header),
                  rms=rms if rms > 0 else float(mad_std(data, ignore_nan=True)),
                  peak=float(np.nanmax(data)),
                  peak_yx=(int(peak_yx[0]), int(peak_yx[1])),
                  beam=_beam_arcsec(header),
                  cell=_cell_arcsec(image, header),
                  applied=_applied_shift(header))


def _same_file(a, b):
    """Whether two path-ish objects name the same file.

    A plain ``==`` is not enough: the reference handed to `t_cutout_2D_radec`
    as ``ref_cutout_image=`` and the same file inside `imagelist` routinely
    differ by a ``./`` or a symlinked parent, and treating them as two files
    would put the reference in the blink set twice.
    """
    if not (isinstance(a, str) and isinstance(b, str)):
        return False
    try:
        return os.path.realpath(a) == os.path.realpath(b)
    except (TypeError, ValueError, OSError):
        return a == b


def _shared_mask(ref, nsigma=5.0, grow_px=6):
    """One mask, built from the reference, used for every frame.

    Measuring each frame over an identical pixel set is what makes the residual
    shifts comparable between bands of very different depth -- otherwise a deep
    band and a shallow one are answering slightly different questions.
    """
    m = np.nan_to_num(ref.data, nan=0.0) > nsigma * ref.rms
    if grow_px and m.any():
        m = ndimage.binary_dilation(m, iterations=int(grow_px))
    return m


def visualise_alignment(imagelist,
                        reference=None,
                        ref_image=None,
                        residuallist=None,
                        mode='auto',
                        panels=('image', 'overlay', 'centroids', 'cuts'),
                        check_shift=True,
                        method='auto',
                        scaling='sigma',
                        components=None,
                        n_components=None,
                        nsigma=5.0,
                        sigma_mask=None,
                        overlay_sigmas=None,
                        cut_width=None,
                        box_size=None,
                        center=None,
                        sort_by_freq=True,
                        cmap='magma_r',
                        figsize=(11, 5),
                        dpi=110,
                        fps=2,
                        save_name=None,
                        verbose=1):
    """Inspect the astrometric alignment of a list of cutouts, in the notebook.

    The in-notebook replacement for ``mlibs.imview(imagelist_c)``: blink
    through the images with the reference's contours frozen on top, backed by
    three static, amplitude-free diagnostic panels and a measured residual
    shift per image.

    Parameters
    ----------
    imagelist : list of str or list of ndarray
        The cutouts to inspect, exactly as built in the notebooks (the third
        return value of ``t_cutout_2D_radec``). Arrays work too, but then the
        frequency, beam and provenance annotations are unavailable.
    ref_image : str, optional
        The image everything is compared against -- normally the very file
        passed to ``t_cutout_2D_radec`` as ``ref_cutout_image=``, i.e. the
        grid the cutouts were resampled onto. It does **not** have to be a
        member of `imagelist`: when it is not, it is read in and joined to the
        blink set, since the contours frozen on every frame come from it and
        being able to flip to it is the point. Without this the first image is
        used, which is only the right baseline by accident.
    reference : str, int or None
        Older spelling of the same thing, also accepting an index into
        `imagelist` or ``'astrref'`` to resolve the file named by the
        ``ASTRREF`` header. `ref_image` wins if both are given.
    residuallist : list of str, optional
        Parallel list of residual images, used for the noise. When omitted the
        WSClean ``-image`` -> ``-residual`` companion is tried automatically,
        falling back to ``mad_std`` of the image.
    mode : {'auto', 'slider', 'animate', 'grid'}
        How the blink frames are shown. ``'auto'`` picks ``'slider'`` when
        ipywidgets is importable inside a notebook kernel, else ``'animate'``,
        else ``'grid'``. Python-side availability cannot prove the widget
        *JavaScript* is rendering: **if the slider cell comes up blank, pass
        ``mode='animate'``**, which embeds the frames in the notebook output
        and always works.
    panels : sequence of str
        Any of ``'image'`` (the blink frames), ``'overlay'``, ``'centroids'``,
        ``'cuts'``.
    check_shift : bool
        Re-measure each image against the reference with
        :func:`~image_alignment.estimate_image_shift`. This is the "did it
        actually work" number, as opposed to the ``ASTRSHF*`` headers, which
        only say what was applied. Costs roughly 0.2-1.3 s per image.
    method : str
        Backend passed through to ``estimate_image_shift``.
    scaling : {'sigma', 'peak', 'absolute'}
        Intensity scaling across frames. ``'sigma'`` (default) gives every
        frame the same dynamic range *in units of its own noise*, so the noise
        floor looks identical everywhere and apparent changes in the source's
        extent are real rather than an artefact of rescaling -- the right
        default when the flux changes a lot with frequency. In ``'slider'``
        mode this is a live dropdown.
    components : list of ndarray, optional
        Explicit boolean masks for the centroid panel. By default the
        reference is deblended automatically.
    sigma_mask : float, optional
        Detection threshold, in units of the reference rms, for everything the
        tool derives from the emission: the components tracked on the centroid
        panel, the mask the residual shift is measured over, the frozen
        contours and the auto-zoom. Set it to the same value the pipeline uses
        for its masks (e.g. ``sigma_mask=sigma_mask_ref``). Raising it is the
        cure for low-level sidelobes being picked up as components and
        reported alongside the real ones. Defaults to `nsigma`.
    n_components, nsigma : int, float
        How many components to keep (brightest first), and the default
        threshold when `sigma_mask` is not given.
    overlay_sigmas : sequence of float, optional
        Contour levels for the overlay panel, in units of *each band's own*
        rms. Defaults to ``(sigma_mask, 4 * sigma_mask)``, so the overlay
        honours the same floor as the rest of the analysis.
    cut_width : int, optional
        Rows/columns averaged for the 1-D cuts. Defaults to one beam.
    box_size : int or 'full', optional
        Side of the displayed box, in pixels. By default the view is zoomed to
        the emission's bounding box plus three beams -- a source occupying a
        tenth of the cutout is useless for judging a sub-pixel offset. Pass
        ``'full'`` to show the whole array. The measurement always uses the
        full image regardless.
    center : (int, int), optional
        ``(y, x)`` centre of that box. Defaults to the reference peak.
    sort_by_freq : bool
        Order the frames by frequency.
    save_name : str, optional
        Write ``<save_name>_diagnostics.png`` and, when Pillow is available,
        ``<save_name>_blink.gif``.

    Returns
    -------
    AlignmentReport
        Carries ``.table`` (a DataFrame, one row per image), ``.results``,
        ``.frames`` and ``.verdict``, so the check can be scripted rather than
        only eyeballed.

    Notes
    -----
    Every comparison drawn here is amplitude-free -- see the module docstring
    for why an image-difference panel is deliberately absent. One consequence
    is worth stating plainly: on the centroid panel a genuine
    frequency-dependent core shift looks exactly like bad alignment. The panel
    shows you the motion; deciding which it is remains yours.

    Examples
    --------
    >>> report = mlibs.visualise_alignment(imagelist_c,
    ...                                        ref_image=new_filename,
    ...                                        sigma_mask=sigma_mask_ref,
    ...                                        residuallist=residuallist_c)
    >>> report.table[['freq_ghz', 'residual_dy', 'residual_dx', 'verdict']]
    """
    if isinstance(imagelist, (str, np.ndarray)):
        imagelist = [imagelist]
    imagelist = list(imagelist)
    if not imagelist:
        raise ValueError('imagelist is empty')
    if check_shift and not _HAVE_ALIGN:
        warnings.warn('image_alignment is not importable; skipping the '
                      'residual-shift measurement')
        check_shift = False

    residuals = list(residuallist) if residuallist is not None else [None] * len(imagelist)
    if len(residuals) != len(imagelist):
        raise ValueError(f'residuallist has {len(residuals)} entries for '
                         f'{len(imagelist)} images')

    # One threshold drives the components, the shift mask, the frozen
    # contours, the overlay and the auto-zoom, so raising it to the pipeline's
    # own mask level takes sidelobes out of every one of them at once.
    thresh = float(sigma_mask) if sigma_mask is not None else float(nsigma)
    if overlay_sigmas is None:
        overlay_sigmas = (thresh, 4.0 * thresh)

    frames = [_build_frame(im, res) for im, res in zip(imagelist, residuals)]

    def _sorted(fr):
        if not sort_by_freq:
            return fr
        freqs = np.asarray([f.freq for f in fr], dtype=float)
        if not np.all(np.isfinite(freqs)):
            return fr
        return [fr[i] for i in np.argsort(freqs)]

    frames = _sorted(frames)

    # -- reference ---------------------------------------------------------
    if ref_image is not None:
        if reference is not None:
            warnings.warn('both ref_image and reference were given; '
                          'using ref_image')
        reference = ref_image

    if reference is None:
        ref = frames[0]
    elif isinstance(reference, (int, np.integer)):
        ref = frames[int(reference)]
    elif isinstance(reference, str) and reference == 'astrref':
        named = {f.applied.get('applied_ref') for f in frames} - {None, '', 'precomputed'}
        ref = None
        if len(named) == 1:
            want = named.pop()
            for f in frames:
                if os.path.basename(str(f.name)).startswith(want[:40]):
                    ref = f
                    break
        if ref is None:
            warnings.warn("reference='astrref' could not be resolved from the "
                          'ASTRREF headers; using the first image')
            ref = frames[0]
    else:
        match = [f for f in frames
                 if f.name == reference or _same_file(f.name, reference)]
        if match:
            ref = match[0]
        else:
            # The cutouts' reference grid is usually a file of its own -- the
            # `ref_cutout_image=` handed to t_cutout_2D_radec -- and not one of
            # the cutouts. Join it to the blink set rather than keeping it
            # off-screen: it is the image whose contours are frozen on every
            # other frame.
            ref = _build_frame(reference)
            frames = _sorted(frames + [ref])
            shapes = [f.data.shape for f in frames if f is not ref]
            if shapes and ref.data.shape != max(set(shapes), key=shapes.count):
                warnings.warn(
                    f'ref_image {os.path.basename(str(ref.name))} has shape '
                    f'{ref.data.shape}, which does not match the cutouts; '
                    'the residual shifts cannot be measured against it. Pass '
                    'the reference *cutout*, on the same grid as imagelist.')

    # -- components and centroids -----------------------------------------
    beam_area_px = 1.0
    if ref.beam and ref.cell:
        beam_area_px = np.pi * (ref.beam[0] / ref.cell) * (ref.beam[1] / ref.cell) / 4.0
    if components is None:
        comps = _components_from(ref.data, ref.rms, nsigma=thresh,
                                 min_area=max(1.0, 0.5 * beam_area_px),
                                 n_components=n_components,
                                 min_separation=(ref.beam_px
                                                 if np.isfinite(ref.beam_px)
                                                 else 3.0),
                                 verbose=verbose)
    else:
        comps = [np.asarray(m, dtype=bool) for m in components]

    ref.centroids = [_centroid(ref.data, m, ref.rms) for m in comps]
    for frame in frames:
        frame.centroids = ([_centroid(frame.data, m, frame.rms) for m in comps]
                           if frame.data.shape == ref.data.shape else [])

    # -- residual shift ----------------------------------------------------
    mask = _shared_mask(ref, nsigma=thresh) if check_shift else None
    for frame in frames:
        if not check_shift:
            continue
        if frame.data.shape != ref.data.shape:
            # estimate_image_shift requires a common grid, and t_cutout_2D_radec
            # can return a smaller array when mode='trim' clips at a border.
            frame.flags.append('grid_mismatch')
            continue
        if frame is ref:
            frame.result = ShiftResult(method='reference', flags=['reference'])
            continue
        try:
            frame.result = estimate_image_shift(
                ref.image if isinstance(ref.image, str) else ref.data,
                frame.image if isinstance(frame.image, str) else frame.data,
                method=method, mask=mask, verbose=0)
        except Exception as exc:
            frame.flags.append(f'failed: {exc}')
            if verbose:
                print(f'  ! shift measurement failed for '
                      f'{os.path.basename(str(frame.name))}: {exc}')

    # -- table -------------------------------------------------------------
    import pandas as pd
    rows = []
    for frame in frames:
        res = frame.result
        verdict, _ = _verdict_of(frame)
        row = {
            'image': os.path.basename(str(frame.name)),
            'freq_ghz': frame.freq_ghz,
            'beam_arcsec': frame.beam[0] if frame.beam else np.nan,
            'px_per_beam': frame.beam_px,
            'rms': frame.rms,
            'peak': frame.peak,
            'snr': frame.peak / frame.rms if frame.rms > 0 else np.nan,
            'applied_dy': frame.applied.get('applied_dy', np.nan),
            'applied_dx': frame.applied.get('applied_dx', np.nan),
            'applied_mode': frame.applied.get('applied_method', ''),
            'applied_clamped': frame.applied.get('applied_clamped', False),
            'residual_dy': res.dy if res is not None else np.nan,
            'residual_dx': res.dx if res is not None else np.nan,
            'residual_dy_err': res.dy_err if res is not None else np.nan,
            'residual_dx_err': res.dx_err if res is not None else np.nan,
            'backend': res.method if res is not None else '',
            'verdict': verdict,
            'flags': ','.join(list(frame.flags)
                              + (list(res.flags) if res is not None else [])),
        }
        mag = np.hypot(row['residual_dy'], row['residual_dx'])
        row['residual_mas'] = mag * frame.cell * 1e3 if frame.cell else np.nan
        row['residual_beam'] = mag / frame.beam_px if np.isfinite(frame.beam_px) else np.nan
        rows.append(row)
    table = pd.DataFrame(rows)

    order = {'aligned': 0, 'not measured': 1, 'marginal': 2, 'misaligned': 3}
    worst = max((_verdict_of(f)[0] for f in frames if f is not ref),
                key=lambda v: order.get(v, 0), default='not measured')

    # -- render ------------------------------------------------------------
    view = _view_slice(ref, box_size=box_size, center=center, nsigma=thresh)

    labels = [(f'{f.freq_ghz:.3f} GHz  ' if np.isfinite(f.freq_ghz) else '')
              + os.path.basename(str(f.name))[:46] for f in frames]

    def _rerender(scale=scaling, contours=True, crosshair=True):
        return [_fig_to_png(_render_frame(f, ref, scaling=scale, cmap=cmap,
                                          show_ref_contours=contours,
                                          show_crosshair=crosshair, view=view,
                                          figsize=figsize,
                                          floor_sigma=thresh), dpi=dpi)
                for f in frames]

    pngs = []
    resolved = _display_backend(mode, verbose=verbose)
    if 'image' in panels:
        if resolved == 'grid':
            fig = _grid_figure(frames, ref, scaling=scaling, cmap=cmap,
                               floor_sigma=thresh,
                               view=view)
            pngs = [_fig_to_png(fig, dpi=dpi)]
        else:
            pngs = _rerender(scaling)

    summary_png = None
    if any(p in panels for p in ('overlay', 'centroids', 'cuts')):
        summary_fig = _summary_figure(frames, ref, panels, overlay_sigmas,
                                      cut_width, view=view)
        if summary_fig is not None:
            summary_png = _fig_to_png(summary_fig, dpi=dpi)

    report = AlignmentReport(table=table, frames=pngs,
                             results=[f.result for f in frames],
                             reference=str(ref.name), components=comps,
                             summary_png=summary_png, verdict=worst)

    # -- display -----------------------------------------------------------
    try:
        from IPython.display import display as _display, HTML, Image as _IPImage
        _has_ipython = True
    except Exception:
        _has_ipython = False

    if 'image' in panels and pngs:
        if verbose:
            print(f'display mode: {resolved}'
                  + ('' if mode != 'auto' else "  (mode='auto')"))
        if resolved == 'slider':
            try:
                _show_slider({'frames': pngs}, labels, _rerender, scaling, fps=fps)
            except Exception as exc:
                warnings.warn(f'ipywidgets slider failed ({exc}); '
                              'falling back to the embedded animation')
                resolved = 'animate'
        if resolved == 'animate':
            report.html = _html_player(pngs, labels, fps=fps)
            if _has_ipython:
                _display(HTML(report.html))
        elif resolved == 'grid' and _has_ipython:
            _display(_IPImage(data=pngs[0]))

    if summary_png is not None and _has_ipython:
        _display(_IPImage(data=summary_png))

    # -- report ------------------------------------------------------------
    if save_name:
        base = os.path.splitext(str(save_name))[0]
        if summary_png is not None:
            with open(base + '_diagnostics.png', 'wb') as fh:
                fh.write(summary_png)
        if len(pngs) > 1:
            try:
                from PIL import Image as _PILImage
                imgs = [_PILImage.open(io.BytesIO(p)).convert('RGB') for p in pngs]
                imgs[0].save(base + '_blink.gif', save_all=True,
                             append_images=imgs[1:],
                             duration=int(1000 / max(1, fps)), loop=0)
                if verbose:
                    print(f'wrote {base}_blink.gif')
            except ImportError:
                warnings.warn('Pillow not available; skipping the GIF')

    if verbose:
        cols = ['image', 'freq_ghz', 'snr', 'applied_dy', 'applied_dx',
                'residual_dy', 'residual_dx', 'residual_beam', 'backend',
                'verdict']
        with pd.option_context('display.width', 200,
                               'display.max_colwidth', 34):
            print(table[cols].to_string(index=False,
                                        float_format=lambda v: f'{v:.3f}'))
        print(f'\nreference : {os.path.basename(str(ref.name))}')
        snrs = [float(np.nanmax(ref.data[m])) / ref.rms for m in comps]
        print(f'components: {len(comps)} at >{thresh:g} sigma'
              + ('  (peak SNR: '
                 + ', '.join(f'{v:.0f}' for v in snrs) + ')' if snrs else ''))
        if len(snrs) > 1 and snrs[-1] < 0.02 * snrs[0]:
            print('            the faintest component is <2% of the brightest '
                  '-- if that is a\n'
                  '            sidelobe, raise sigma_mask or set n_components.')
        print(f'verdict   : {worst}')
        if comps and 'centroids' in panels:
            print('note      : a real frequency-dependent core shift is '
                  'indistinguishable from\n'
                  '            misalignment on the centroid panel.')
    return report
