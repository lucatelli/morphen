"""
 ____________________________________________________________________________
|                                                                            |
|                        Field / Mosaic Source Extraction                    |
|____________________________________________________________________________|

Detection and cataloguing of *many* sources across a mosaic or wide-field image,
plus the bridge that turns that catalogue into per-source cutouts.

Everything else in morphen assumes a postage stamp: one source (simple or
complex) roughly centred in a small image. `sort_by='distance'` means distance
from the image centre, `mad_std` over the whole array is used as *the* noise
level, and the morphometry routines build flux growth curves over the full
frame. None of that survives contact with a mosaic, where the noise varies
across the field with primary-beam response and uv coverage, there is no
privileged central source, and the edges are NaN.

The pipeline here is therefore deliberately split in two:

    1. detection -> catalogue      (`field_source_ext` + `build_field_catalogue`)
    2. catalogue -> cutouts        (`cutouts_from_catalogue`)

and stage 3 (the actual analysis) is the *existing* morphen machinery, run
unchanged on the cutouts, because a cutout is exactly the postage stamp those
routines were written for.

Notes
-----
This file is loaded by `mlibs.py` via `exec()` into a shared namespace, like its
siblings. It therefore relies on names already imported there (`np`, `pd`,
`plt`, `WCS`, `SkyCoord`, `u`, `fits`, `Background2D`, `MedianBackground`,
`MADStdBackgroundRMS`, `SigmaClip`, `convolve`, `detect_sources`, `simple_norm`,
`tqdm`, `os`, ...) as well as on morphen functions defined in `data_io.py`
(`load_fits_data`, `get_cell_size`, `t_cutout_2D_radec`), `radio_utils.py`
(`beam_area2`, `get_beam_size_px`) and `plotting.py` (`eimshow`). It is not
independently importable.
"""


# ---------------------------------------------------------------------------
# background / noise
# ---------------------------------------------------------------------------

def field_beam_size(imagename, default=5.0, verbose=0):
    """
    Beam size in pixels for a field image, with a graceful fallback.

    Returns
    -------
    (beam_fwhm_px, beam_area_px) : tuple of float
        Geometric-mean beam FWHM in pixels and the beam area in pixels. For
        images with no restoring beam in the header (optical/IR), falls back to
        `default` pixels and a beam area of one resolution element.
    """
    try:
        beam_fwhm_px, _, _ = get_beam_size_px(imagename)
    except Exception:
        beam_fwhm_px = default
        if verbose >= 1:
            print(f' -->> No restoring beam found; assuming '
                  f'{default} px resolution element.')
    try:
        beam_area_px = beam_area2(imagename)
    except Exception:
        beam_area_px = np.pi * (beam_fwhm_px / 2.0) ** 2.0

    if not np.isfinite(beam_fwhm_px) or beam_fwhm_px <= 0:
        beam_fwhm_px = default
    if not np.isfinite(beam_area_px) or beam_area_px <= 0:
        beam_area_px = np.pi * (beam_fwhm_px / 2.0) ** 2.0

    return float(beam_fwhm_px), float(beam_area_px)


def field_background(data_2D, coverage_mask=None, box_size=None,
                     filter_size=(3, 3), beam_fwhm_px=5.0,
                     bkg_estimator='median', sigma_clip_sigma=3.0,
                     exclude_percentile=20.0, rms_map=None, verbose=0):
    """
    Estimate a 2D background and a 2D RMS map for a field/mosaic image.

    Unlike `multiscale_segmentation_background` (used by the postage-stamp
    paths), this is `photutils.Background2D` with a beam-scaled mesh: cheap
    enough for a 10k x 10k mosaic and, crucially, it returns a *spatially
    varying* RMS so that the detection threshold follows the primary-beam /
    coverage noise instead of a single global `mad_std`.

    Parameters
    ----------
    data_2D : np.ndarray
        2D image data.
    coverage_mask : np.ndarray, optional
        Boolean array, True where pixels must be ignored (photutils
        convention). Defaults to the non-finite pixels of `data_2D`.
    box_size : int, optional
        Mesh box size in pixels. Defaults to ~10 beams, clipped to [32, 256].
    filter_size : tuple, optional
        Median-filter size applied to the low-resolution background mesh.
    beam_fwhm_px : float, optional
        Beam FWHM in pixels, used to pick the default `box_size`.
    bkg_estimator : str, optional
        'median' (default), 'sextractor' or 'mmm'.
    exclude_percentile : float, optional
        Boxes with more than this percentage of masked pixels are excluded.
        Raised above the photutils default of 10 because mosaic edges are
        typically NaN-heavy.
    rms_map : np.ndarray or str, optional
        Externally supplied noise map (array or FITS path). If given, the
        background is still estimated but this map is used as the RMS, which is
        the honest option when the imager already produced a noise map.

    Returns
    -------
    (bkg_image, bkg_rms) : tuple of np.ndarray
    """
    if coverage_mask is None:
        coverage_mask = ~np.isfinite(data_2D)

    if box_size is None:
        box_size = int(np.clip(10.0 * beam_fwhm_px, 32, 256))
    box_size = int(min(box_size,
                       max(8, min(data_2D.shape[0], data_2D.shape[1]) // 4)))

    if bkg_estimator == 'sextractor':
        _bkg_estimator = SExtractorBackground()
    elif bkg_estimator == 'mmm':
        _bkg_estimator = MMMBackground()
    else:
        _bkg_estimator = MedianBackground()

    if verbose >= 1:
        print(f' ++==>> Background2D: box_size={box_size}, '
              f'filter_size={filter_size}, estimator={bkg_estimator}')

    bkg = Background2D(data_2D, box_size=(box_size, box_size),
                       filter_size=filter_size,
                       coverage_mask=coverage_mask,
                       fill_value=np.nan,
                       exclude_percentile=exclude_percentile,
                       sigma_clip=SigmaClip(sigma=sigma_clip_sigma, maxiters=10),
                       bkg_estimator=_bkg_estimator,
                       bkgrms_estimator=MADStdBackgroundRMS())

    bkg_image = bkg.background
    bkg_rms = bkg.background_rms

    if rms_map is not None:
        if isinstance(rms_map, str):
            rms_map = load_fits_data(rms_map)
        rms_map = np.squeeze(np.asarray(rms_map))
        if rms_map.shape != data_2D.shape:
            raise ValueError(f'rms_map shape {rms_map.shape} does not match '
                             f'image shape {data_2D.shape}.')
        bkg_rms = rms_map
        if verbose >= 1:
            print(' ++==>> Using externally supplied RMS map.')

    return bkg_image, bkg_rms


# ---------------------------------------------------------------------------
# catalogue construction
# ---------------------------------------------------------------------------

def _sexagesimal_strings(sky_coords):
    """
    Sexagesimal RA/Dec strings that round-trip through `conver_str_coords`.

    `conver_str_coords` parses with `unit=(u.hourangle, u.deg)`, so RA must be
    in hours and the separator must be whitespace. Keeping this exact format is
    what allows catalogues produced here to be dropped into the older
    RA_LIST/DEC_LIST notebook loops unchanged.
    """
    ra_str = sky_coords.ra.to_string(unit=u.hourangle, sep=' ',
                                     precision=4, pad=True)
    dec_str = sky_coords.dec.to_string(unit=u.deg, sep=' ', precision=3,
                                       pad=True, alwayssign=True)
    return np.atleast_1d(ra_str), np.atleast_1d(dec_str)


def _iau_names(sky_coords, prefix='J'):
    """
    IAU-style designations (e.g. ``J131534.94+620728.6``), safe as filenames.
    """
    ra_c = sky_coords.ra.to_string(unit=u.hourangle, sep='', precision=2,
                                   pad=True)
    dec_c = sky_coords.dec.to_string(unit=u.deg, sep='', precision=1, pad=True,
                                     alwayssign=True)
    ra_c = np.atleast_1d(ra_c)
    dec_c = np.atleast_1d(dec_c)
    return np.asarray([f'{prefix}{r}{d}' for r, d in zip(ra_c, dec_c)])


def _safe_values(quantity):
    """Return a plain float array from a photutils property (Quantity or not)."""
    values = np.atleast_1d(quantity)
    if hasattr(values, 'value'):
        values = values.value
    return np.asarray(values, dtype=float)


def build_field_catalogue(cat, cat_parent, parent_labels, wcs,
                          cell_size=1.0, beam_area_px=1.0,
                          bkg_rms=None, image_shape=None,
                          edge_margin=0, cutout_size_factor=4.0,
                          name_prefix='J', verbose=0):
    """
    Turn a pair of photutils `SourceCatalog` objects into a tidy DataFrame.

    One row per *deblended component*; the parent (blended island) it belongs
    to is carried alongside in the `parent_*` columns, so nothing is lost and
    cutouts can be keyed on either.

    Parameters
    ----------
    cat : photutils.segmentation.SourceCatalog
        Catalogue of the deblended components.
    cat_parent : photutils.segmentation.SourceCatalog
        Catalogue of the parent (pre-deblend) islands.
    parent_labels : np.ndarray
        For each component in `cat`, the label of the parent island it sits in.
    wcs : astropy.wcs.WCS
        Celestial WCS (naxis=2) of the detection image.
    cell_size : float
        Pixel scale in arcsec, from `get_cell_size`.
    beam_area_px : float
        Beam area in pixels, from `beam_area2`.
    bkg_rms : np.ndarray, optional
        2D RMS map, sampled at each centroid to give `local_rms` / `snr_peak`.
    image_shape : tuple, optional
        (ny, nx), needed for the `edge_flag` column.
    edge_margin : int
        A source whose parent bbox comes within this many pixels of the border
        is flagged (not removed).
    cutout_size_factor : float
        `suggested_cutout_size` = factor x parent bbox size, rounded up to even.

    Returns
    -------
    pandas.DataFrame
    """
    n_src = len(cat)

    xc = _safe_values(cat.xcentroid)
    yc = _safe_values(cat.ycentroid)
    sky = cat.sky_centroid
    ra = np.atleast_1d(sky.ra.degree)
    dec = np.atleast_1d(sky.dec.degree)
    ra_str, dec_str = _sexagesimal_strings(sky)

    # --- parent quantities, broadcast onto the component rows ---------------
    parent_lookup = {}
    p_labels = np.atleast_1d(cat_parent.labels)
    p_xc = _safe_values(cat_parent.xcentroid)
    p_yc = _safe_values(cat_parent.ycentroid)
    p_sky = cat_parent.sky_centroid
    p_ra = np.atleast_1d(p_sky.ra.degree)
    p_dec = np.atleast_1d(p_sky.dec.degree)
    p_xmin = _safe_values(cat_parent.bbox_xmin)
    p_xmax = _safe_values(cat_parent.bbox_xmax)
    p_ymin = _safe_values(cat_parent.bbox_ymin)
    p_ymax = _safe_values(cat_parent.bbox_ymax)
    p_flux = _safe_values(cat_parent.segment_flux)
    p_area = _safe_values(cat_parent.area)

    for k, lab in enumerate(p_labels):
        parent_lookup[int(lab)] = k

    idx_parent = np.asarray([parent_lookup.get(int(lab), -1)
                             for lab in parent_labels])
    valid = idx_parent >= 0
    if not np.all(valid) and verbose >= 1:
        print(f' -->> WARNING: {np.sum(~valid)} component(s) could not be '
              f'matched to a parent island; parent columns set to the '
              f'component values.')

    def _pick(parent_array, fallback):
        out = np.where(valid, parent_array[np.clip(idx_parent, 0, None)],
                       fallback)
        return np.asarray(out, dtype=float)

    xc_parent = _pick(p_xc, xc)
    yc_parent = _pick(p_yc, yc)
    ra_parent = _pick(p_ra, ra)
    dec_parent = _pick(p_dec, dec)
    bbox_xmin = _pick(p_xmin, xc)
    bbox_xmax = _pick(p_xmax, xc)
    bbox_ymin = _pick(p_ymin, yc)
    bbox_ymax = _pick(p_ymax, yc)
    parent_flux = _pick(p_flux, _safe_values(cat.segment_flux))
    parent_area = _pick(p_area, _safe_values(cat.area))

    parent_size_px = np.maximum(bbox_xmax - bbox_xmin + 1.0,
                                bbox_ymax - bbox_ymin + 1.0)
    suggested = np.ceil(cutout_size_factor * parent_size_px / 2.0) * 2.0

    # --- local noise --------------------------------------------------------
    peak_flux = _safe_values(cat.max_value)
    if bkg_rms is not None:
        yi = np.clip(np.round(yc).astype(int), 0, bkg_rms.shape[0] - 1)
        xi = np.clip(np.round(xc).astype(int), 0, bkg_rms.shape[1] - 1)
        local_rms = np.asarray(bkg_rms[yi, xi], dtype=float)
    else:
        local_rms = np.full(n_src, np.nan)
    with np.errstate(divide='ignore', invalid='ignore'):
        snr_peak = peak_flux / local_rms

    # --- flags --------------------------------------------------------------
    if image_shape is not None:
        ny, nx = image_shape
        edge_flag = ((bbox_xmin <= edge_margin) |
                     (bbox_ymin <= edge_margin) |
                     (bbox_xmax >= nx - 1 - edge_margin) |
                     (bbox_ymax >= ny - 1 - edge_margin))
    else:
        edge_flag = np.zeros(n_src, dtype=bool)

    area_px = _safe_values(cat.area)
    equivalent_radius = _safe_values(cat.equivalent_radius)

    df = pd.DataFrame({
        'label': np.atleast_1d(cat.labels).astype(int),
        'parent_label': np.asarray(parent_labels).astype(int),
        'xc': xc,
        'yc': yc,
        'ra': ra,
        'dec': dec,
        'ra_str': ra_str,
        'dec_str': dec_str,
        'iau_name': _iau_names(sky, prefix=name_prefix),
        'xc_parent': xc_parent,
        'yc_parent': yc_parent,
        'ra_parent': ra_parent,
        'dec_parent': dec_parent,
        'peak_flux': peak_flux,
        'segment_flux': _safe_values(cat.segment_flux),
        'segment_fluxerr': _safe_values(cat.segment_fluxerr),
        'kron_flux': _safe_values(cat.kron_flux),
        'kron_fluxerr': _safe_values(cat.kron_fluxerr),
        'local_rms': local_rms,
        'snr_peak': snr_peak,
        'area_px': area_px,
        'area_beams': area_px / beam_area_px,
        'equivalent_radius_px': equivalent_radius,
        'equivalent_radius_arcsec': equivalent_radius * cell_size,
        'semimajor_px': _safe_values(cat.semimajor_sigma),
        'semiminor_px': _safe_values(cat.semiminor_sigma),
        'ellipticity': _safe_values(cat.ellipticity),
        'orientation_deg': _safe_values(cat.orientation),
        'fwhm_px': _safe_values(cat.fwhm),
        'parent_flux': parent_flux,
        'parent_area_px': parent_area,
        'bbox_xmin': bbox_xmin.astype(int),
        'bbox_xmax': bbox_xmax.astype(int),
        'bbox_ymin': bbox_ymin.astype(int),
        'bbox_ymax': bbox_ymax.astype(int),
        'parent_size_px': parent_size_px,
        'suggested_cutout_size': suggested.astype(int),
        'edge_flag': edge_flag,
    })

    # --- per-parent bookkeeping --------------------------------------------
    df['n_components'] = df.groupby('parent_label')['label'].transform('count')
    brightest = df.groupby('parent_label')['segment_flux'].transform('max')
    df['is_parent_primary'] = df['segment_flux'] >= brightest
    # a tie (identical fluxes) would mark two rows; keep only the first.
    dup = df['is_parent_primary'] & df.duplicated(
        subset=['parent_label', 'is_parent_primary'], keep='first')
    df.loc[dup, 'is_parent_primary'] = False

    parent_ids = {lab: i + 1 for i, lab in
                  enumerate(sorted(df['parent_label'].unique()))}
    df['parent_id'] = df['parent_label'].map(parent_ids)

    return df


def sort_field_catalogue(df, sort_by='flux', image_shape=None, verbose=0):
    """
    Sort a field catalogue and (re)assign the 1-based `id` column.

    `sort_by='distance'` is the default of `mp.source_extraction` and means
    "distance from the image centre" - meaningful for a postage stamp, not for
    a field. It is accepted (so the class default does not blow up) but falls
    back to 'flux' with a warning.
    """
    if sort_by == 'distance':
        print(" -->> WARNING: sort_by='distance' is meaningless for a field "
              "image (there is no central source). Using sort_by='flux'.")
        sort_by = 'flux'

    if sort_by == 'flux':
        order = np.argsort(df['segment_flux'].values)[::-1]
    elif sort_by == 'peak':
        order = np.argsort(df['peak_flux'].values)[::-1]
    elif sort_by == 'snr':
        order = np.argsort(df['snr_peak'].values)[::-1]
    elif sort_by == 'area':
        order = np.argsort(df['area_px'].values)[::-1]
    elif sort_by == 'ra':
        order = np.argsort(df['ra'].values)
    elif sort_by == 'dec':
        order = np.argsort(df['dec'].values)
    elif sort_by == 'label':
        order = np.argsort(df['label'].values)
    else:
        raise ValueError(f"Unknown sort_by='{sort_by}'. Use one of "
                         f"'flux', 'peak', 'snr', 'area', 'ra', 'dec', 'label'.")

    df_sorted = df.iloc[order].reset_index(drop=True).copy()
    df_sorted.insert(0, 'id', np.arange(1, len(df_sorted) + 1))
    if verbose >= 1:
        print(f' ++==>> Catalogue sorted by {sort_by}.')
    return df_sorted, np.asarray(order)


def filter_field_catalogue(df, snr_min=None, flux_min=None, peak_min=None,
                           area_min=None, area_max=None, edge_margin=None,
                           image_shape=None, exclude_edge=False,
                           primary_only=False, verbose=1):
    """
    Filter a field catalogue. Returns a *copy*; the input is untouched.

    Kept as a free function (rather than only a method) so it can be applied to
    a CSV loaded back later, without re-running detection.

    Parameters
    ----------
    snr_min, flux_min, peak_min : float, optional
        Lower cuts on `snr_peak`, `segment_flux`, `peak_flux`.
    area_min, area_max : float, optional
        Cuts on `area_px`.
    edge_margin : int, optional
        Recompute `edge_flag` with this margin before applying `exclude_edge`.
        Requires `image_shape`.
    exclude_edge : bool
        Drop rows whose parent bbox touches the border margin.
    primary_only : bool
        Keep only the brightest component of each parent island - i.e. one row
        per source rather than one row per component.
    """
    out = df.copy()
    n0 = len(out)

    if edge_margin is not None and image_shape is not None:
        ny, nx = image_shape
        out['edge_flag'] = ((out['bbox_xmin'] <= edge_margin) |
                            (out['bbox_ymin'] <= edge_margin) |
                            (out['bbox_xmax'] >= nx - 1 - edge_margin) |
                            (out['bbox_ymax'] >= ny - 1 - edge_margin))

    if snr_min is not None:
        out = out[out['snr_peak'] >= snr_min]
    if flux_min is not None:
        out = out[out['segment_flux'] >= flux_min]
    if peak_min is not None:
        out = out[out['peak_flux'] >= peak_min]
    if area_min is not None:
        out = out[out['area_px'] >= area_min]
    if area_max is not None:
        out = out[out['area_px'] <= area_max]
    if exclude_edge:
        out = out[~out['edge_flag'].astype(bool)]
    if primary_only:
        out = out[out['is_parent_primary'].astype(bool)]

    out = out.reset_index(drop=True)
    if 'id' in out.columns:
        out['id'] = np.arange(1, len(out) + 1)

    if verbose >= 1:
        print(f' ++==>> Catalogue filtered: {n0} -> {len(out)} rows.')
    return out


# ---------------------------------------------------------------------------
# detection
# ---------------------------------------------------------------------------

def field_source_ext(imagename, residualname=None, residual=None,
                     sigma=5.0, minarea=100, minarea_factor=1.0, npixels=None,
                     deblend_nthresh=32, deblend_cont=1e-3, deblend=True,
                     connectivity=8, nproc=1,
                     filter_type='conv', psf_data=None,
                     bkg_box_size=None, bkg_filter_size=(3, 3),
                     bkg_estimator='median', rms_map=None,
                     mask=None, apply_mask=False, sigma_mask=6.0,
                     dilation_size=None, iterations=2,
                     edge_margin=0, reject_edge_sources=False,
                     snr_min=None,
                     sort_by='flux', cutout_size_factor=4.0,
                     name_prefix='J',
                     return_masks=False, max_masks=200,
                     show_detection=False, show_bkg_map=False,
                     max_labels_to_plot=500, ell_size_factor=3.0,
                     save_plot=True,
                     save_products=False, products_path=None,
                     products=('data_sub', 'bkg', 'rms'),
                     verbose=1):
    """
    Source extraction over a mosaic / wide-field image.

    Detection is done on a beam-convolved, background-subtracted image against a
    *2D* threshold (`sigma * bkg_rms`), so the effective sensitivity follows the
    noise structure of the mosaic. Sources are detected as parent islands and
    then deblended; the returned catalogue has one row per deblended component,
    with the parent island carried alongside.

    Parameters
    ----------
    imagename : str
        Path to the field/mosaic FITS image.
    residualname, residual : str or np.ndarray, optional
        Residual image, used only as a sanity cross-check on the noise level.
        The detection threshold always comes from the 2D RMS map.
    sigma : float
        Detection threshold in units of the local RMS.
    minarea, minarea_factor, npixels : float / int
        Minimum connected area for a detection, in pixels. If `npixels` is None
        it is `int(minarea * minarea_factor)`, the same convention the
        postage-stamp paths use (one restoring beam by default).
    deblend_nthresh, deblend_cont : int, float
        `nlevels` and `contrast` for `photutils.deblend_sources`.
    nproc : int
        Processes used for deblending. Worth raising for a crowded field.
    filter_type : str or None
        'conv' convolves with a beam-matched Gaussian before detection (better
        SNR for marginal sources); None detects on the unconvolved image.
    bkg_box_size, bkg_filter_size, bkg_estimator, rms_map
        Passed to `field_background`.
    apply_mask, mask, sigma_mask, dilation_size, iterations
        Optional global mask restricting where detection happens. Off by
        default: `mask_dilation` on a full mosaic is expensive and rarely what
        you want.
    edge_margin : int
        Border width, in pixels, used for the `edge_flag` column.
    reject_edge_sources : bool
        If True, labels touching `edge_margin` are removed from the
        segmentation entirely rather than merely flagged.
    snr_min : float, optional
        Drop catalogue rows below this peak SNR.
    sort_by : str
        'flux' (default), 'peak', 'snr', 'area', 'ra', 'dec', 'label'.
    return_masks : bool
        Materialise the per-component boolean masks. Off by default - for a
        field with thousands of detections this is thousands of full-size
        boolean arrays. When True, at most `max_masks` are built.
    show_detection, show_bkg_map, max_labels_to_plot
        Diagnostics. Labelling every source in a crowded field is unusably
        slow, hence the cap.
    save_products : bool
        Write the full-mosaic intermediates (background-subtracted data,
        background, RMS map, and optionally the segmentation images) to FITS.
        The RMS map can be fed straight back in as `rms_map=` on a later run.
    products_path : str, optional
        Where to write them. Defaults to the directory of `imagename`.
    products : tuple of str
        Which to write: any of 'data_sub', 'bkg', 'rms', 'segm', 'segm_parent'.

    Returns
    -------
    tuple
        (masks_deblended, indices, bkg_image, bkg_rms, seg_maps,
         seg_maps_parent, objects_sorted, cat, cat_parent, catalogue, wcs)
    """
    from photutils.segmentation import deblend_sources, SourceCatalog
    from photutils.segmentation import make_2dgaussian_kernel

    # --- load ---------------------------------------------------------------
    data_2D = load_fits_data(imagename)
    data_2D = np.squeeze(np.asarray(data_2D, dtype=float))
    if data_2D.ndim != 2:
        raise ValueError(f'Expected a 2D image after squeezing, got shape '
                         f'{data_2D.shape} for {imagename}.')

    header = fits.getheader(imagename)
    wcs = WCS(header, naxis=2)
    cell_size = get_cell_size(imagename)
    beam_fwhm_px, beam_area_px = field_beam_size(imagename, verbose=verbose)

    if verbose >= 1:
        print(f' ++==>> Field extraction on '
              f'{os.path.basename(imagename)}  {data_2D.shape}')
        print(f'        beam = {beam_fwhm_px:.2f} px '
              f'({beam_area_px:.1f} px/beam), cell = {cell_size:.4f} arcsec')

    coverage_mask = ~np.isfinite(data_2D)
    n_bad = int(np.sum(coverage_mask))
    if n_bad and verbose >= 1:
        print(f'        {n_bad} non-finite pixels '
              f'({100.0 * n_bad / data_2D.size:.1f}%) excluded.')

    # --- background / noise -------------------------------------------------
    bkg_image, bkg_rms = field_background(data_2D,
                                          coverage_mask=coverage_mask,
                                          box_size=bkg_box_size,
                                          filter_size=bkg_filter_size,
                                          beam_fwhm_px=beam_fwhm_px,
                                          bkg_estimator=bkg_estimator,
                                          rms_map=rms_map,
                                          verbose=verbose)

    data_sub = data_2D - bkg_image

    if residual is None and residualname is not None:
        try:
            residual = np.squeeze(load_fits_data(residualname))
        except Exception:
            residual = None
    if residual is not None and verbose >= 1:
        rms_res = mad_std(residual, ignore_nan=True)
        rms_map_med = np.nanmedian(bkg_rms)
        print(f'        RMS: residual = {rms_res:.3e}, '
              f'median(bkg_rms) = {rms_map_med:.3e}')

    if show_bkg_map:
        plot_field_background(data_2D, bkg_image, bkg_rms)

    # --- optional global mask ----------------------------------------------
    if apply_mask and mask is None:
        _, mask = mask_dilation(data_2D, sigma=sigma_mask,
                                iterations=iterations,
                                rms=np.nanmedian(bkg_rms),
                                PLOT=False, show_figure=False,
                                dilation_size=dilation_size)
    if mask is not None:
        coverage_mask = coverage_mask | (~np.asarray(mask, dtype=bool))

    # --- convolution --------------------------------------------------------
    data_filled = np.where(coverage_mask, 0.0, data_sub)
    if filter_type is not None:
        kernel_size = int(2 * np.ceil(2.0 * beam_fwhm_px) + 1)
        kernel = make_2dgaussian_kernel(beam_fwhm_px, size=kernel_size)
        convolved = convolve(data_filled, kernel, normalize_kernel=True)
        convolved = np.where(coverage_mask, 0.0, convolved)
    else:
        convolved = data_filled

    # --- threshold ----------------------------------------------------------
    # NaNs in the RMS map (fully-masked mesh boxes) become +inf so that nothing
    # is ever detected there, rather than raising invalid-comparison warnings.
    threshold = sigma * np.nan_to_num(bkg_rms, nan=np.inf, posinf=np.inf)

    if npixels is None:
        npixels = int(max(1, minarea * minarea_factor))
    if verbose >= 1:
        print(f'        detection: sigma={sigma}, npixels={npixels}, '
              f'deblend nlevels={deblend_nthresh}, contrast={deblend_cont}')

    # --- parent detection ---------------------------------------------------
    seg_parent = detect_sources(convolved, threshold, npixels=npixels,
                                mask=coverage_mask, connectivity=connectivity)

    empty = ([], np.asarray([], dtype=int), bkg_image, bkg_rms, None, None,
             {'xc': np.asarray([]), 'yc': np.asarray([])}, None, None,
             pd.DataFrame(), wcs)
    if seg_parent is None or seg_parent.nlabels == 0:
        print(' -->> WARNING: no sources detected. Try lowering `sigma` or '
              '`minarea_factor`.')
        return empty

    if reject_edge_sources and edge_margin > 0:
        seg_parent.remove_border_labels(border_width=edge_margin,
                                        partial_overlap=True, relabel=True)
        if seg_parent.nlabels == 0:
            print(' -->> WARNING: all sources removed by `edge_margin`.')
            return empty

    if verbose >= 1:
        print(f'        {seg_parent.nlabels} parent island(s) detected.')

    # --- deblending ---------------------------------------------------------
    if deblend and deblend_nthresh > 0:
        seg_maps = deblend_sources(convolved, seg_parent, npixels=npixels,
                                   nlevels=deblend_nthresh,
                                   contrast=deblend_cont,
                                   mode='exponential',
                                   connectivity=connectivity,
                                   relabel=True, nproc=nproc,
                                   progress_bar=bool(verbose >= 1))
        if verbose >= 1:
            print(f'        {seg_maps.nlabels} component(s) after deblending.')
    else:
        seg_maps = seg_parent

    # --- catalogues ---------------------------------------------------------
    cat = SourceCatalog(data_sub, seg_maps, convolved_data=convolved,
                        error=np.nan_to_num(bkg_rms, nan=np.inf),
                        background=bkg_image, mask=coverage_mask, wcs=wcs)
    cat_parent = SourceCatalog(data_sub, seg_parent, convolved_data=convolved,
                               background=bkg_image, mask=coverage_mask,
                               wcs=wcs)

    # component -> parent, sampled at each component's peak pixel. The peak is
    # used rather than the centroid because the centroid of a crescent-shaped
    # component can fall outside its own footprint.
    yidx = np.clip(np.atleast_1d(cat.maxval_yindex).astype(int),
                   0, seg_parent.shape[0] - 1)
    xidx = np.clip(np.atleast_1d(cat.maxval_xindex).astype(int),
                   0, seg_parent.shape[1] - 1)
    parent_labels = seg_parent.data[yidx, xidx]

    # --- dataframe ----------------------------------------------------------
    catalogue = build_field_catalogue(cat, cat_parent, parent_labels, wcs,
                                      cell_size=cell_size,
                                      beam_area_px=beam_area_px,
                                      bkg_rms=bkg_rms,
                                      image_shape=data_2D.shape,
                                      edge_margin=edge_margin,
                                      cutout_size_factor=cutout_size_factor,
                                      name_prefix=name_prefix,
                                      verbose=verbose)
    catalogue['#imagename'] = os.path.basename(imagename)
    catalogue['cell_size'] = cell_size

    catalogue, order = sort_field_catalogue(catalogue, sort_by=sort_by,
                                            image_shape=data_2D.shape,
                                            verbose=verbose)
    indices = np.asarray(order, dtype=int)

    if snr_min is not None:
        catalogue = filter_field_catalogue(catalogue, snr_min=snr_min,
                                           verbose=verbose)

    objects_sorted = {'xc': catalogue['xc'].values,
                      'yc': catalogue['yc'].values,
                      'ra': catalogue['ra'].values,
                      'dec': catalogue['dec'].values}

    # --- masks (opt-in; see docstring) --------------------------------------
    masks_deblended = []
    if return_masks:
        n_masks = min(len(indices), max_masks)
        if n_masks < len(indices):
            print(f' -->> Only the first {n_masks} of {len(indices)} '
                  f'component masks were built (max_masks).')
        labels = np.atleast_1d(seg_maps.labels)
        for k in range(n_masks):
            masks_deblended.append(seg_maps.data == labels[indices[k]])

    if verbose >= 1:
        print(f' ++==>> Field catalogue: {len(catalogue)} component(s) in '
              f'{catalogue["parent_label"].nunique()} island(s).')

    if save_products:
        save_field_products(imagename, bkg_image=bkg_image, bkg_rms=bkg_rms,
                            data_sub=data_sub, seg_maps=seg_maps,
                            seg_maps_parent=seg_parent,
                            products=products, save_path=products_path,
                            verbose=verbose)

    if show_detection:
        plot_field_detections(imagename, catalogue, segm=seg_maps, cat=cat,
                              data_sub=data_sub, bkg_rms=bkg_rms,
                              max_labels_to_plot=max_labels_to_plot,
                              save_name=(imagename + '_field_detection.jpg'
                                         if save_plot else None))

    return (masks_deblended, indices, bkg_image, bkg_rms, seg_maps, seg_parent,
            objects_sorted, cat, cat_parent, catalogue, wcs)


# ---------------------------------------------------------------------------
# plotting / export
# ---------------------------------------------------------------------------

def plot_field_background(data_2D, bkg_image, bkg_rms, figsize=(18, 6)):
    """Three-panel data / background / RMS diagnostic for a field image."""
    finite = np.isfinite(data_2D)
    vmax = np.nanpercentile(data_2D[finite], 99.9) if finite.any() else 1.0
    vmin = np.nanpercentile(data_2D[finite], 1.0) if finite.any() else 0.0
    norm = simple_norm(data_2D[finite] if finite.any() else data_2D,
                       vmax=vmax, vmin=vmin, stretch='asinh', asinh_a=0.05)

    fig = plt.figure(figsize=figsize)
    ax0 = fig.add_subplot(1, 3, 1)
    ax0.imshow(data_2D, origin='lower', norm=norm, cmap='Greys_r')
    ax0.set_title('data')

    ax1 = fig.add_subplot(1, 3, 2)
    ax1.imshow(bkg_image, origin='lower', cmap='Greys_r')
    ax1.set_title('background')

    ax2 = fig.add_subplot(1, 3, 3)
    im2 = ax2.imshow(bkg_rms, origin='lower', cmap='inferno')
    ax2.set_title('background RMS')
    fig.colorbar(im2, ax=ax2, fraction=0.046)

    plt.show()
    plt.close(fig)


def plot_field_detections(imagename, catalogue, segm=None, cat=None,
                          data_sub=None, data_2D=None, bkg_rms=None,
                          panels='segm', show_labels=False,
                          max_labels_to_plot=500,
                          figsize=None, show_cutout_boxes=False,
                          cutout_size=None, vmin_factor=3.0, vmax_factor=0.2,
                          stretch='asinh', asinh_a=0.02,
                          show_kron_apertures=True, cmap_seed=123,
                          save_name=None, dpi=200):
    """
    Field overview in the standard photutils segmentation style.

    Panels follow the photutils segmentation user guide: the deblended
    segmentation image rendered with `SegmentationImage.make_cmap`, and/or the
    background-subtracted data with the sources' Kron apertures overlaid
    (`SourceCatalog.plot_kron_apertures`). This replaces the hand-rolled
    per-source `Ellipse` overlay the postage-stamp paths use - photutils'
    apertures reflect the actual measured extent, and a segmentation image
    shows the real footprints and deblending boundaries rather than an
    idealised ellipse.

    `segm` and `cat` are subset to the labels present in `catalogue`, so the
    display always matches the (possibly filtered) catalogue rather than the
    raw detection.

    Parameters
    ----------
    panels : {'segm', 'data', 'both'}
        Which panels to draw. Default 'segm': one large segmentation map, which
        is the useful view for a field - it shows every footprint and every
        deblending boundary at once. 'data' gives the background-subtracted
        image with Kron apertures; 'both' puts them side by side.
    show_labels : bool
        Annotate each source with its catalogue ID. Off by default - the plot
        is meant to read as segments plus Kron apertures, and on a real field
        the labels are both illegible and slow to draw.
    max_labels_to_plot : int
        Cap applied when `show_labels=True`, counted down the already-sorted
        catalogue.
    figsize : tuple, optional
        Defaults to (12, 12) for a single panel and (20, 10) for 'both'.
    vmin_factor, vmax_factor : float
        Display limits for the data panel: `vmin_factor` in units of the local
        RMS, `vmax_factor` as a fraction of the image peak, following the
        convention of the rest of morphen's plotting. A field is mostly empty
        sky, so percentile-based limits put the display floor inside the noise
        and wash the image out.
    stretch : str
        Defaults to morphen's usual 'asinh' rather than the 'sqrt' of the
        photutils docs, which leaves faint radio sources close to invisible.
    """
    if panels not in ('segm', 'data', 'both'):
        raise ValueError("panels must be one of 'segm', 'data', 'both'.")
    if data_sub is None:
        if data_2D is None:
            data_2D = np.squeeze(load_fits_data(imagename))
        data_sub = data_2D

    if bkg_rms is not None:
        rms = float(np.nanmedian(bkg_rms))
    else:
        rms = float(mad_std(data_sub, ignore_nan=True))
    vmin = vmin_factor * rms
    vmax = max(vmax_factor * float(np.nanmax(data_sub)), 10.0 * vmin)

    # keep the plotted segments/apertures in sync with the catalogue
    labels = None
    if 'label' in catalogue.columns:
        labels = np.asarray(catalogue['label'].values, dtype=int)

    segm_show = None
    if segm is not None:
        segm_show = segm.copy()
        if labels is not None and len(labels) < segm_show.nlabels:
            segm_show.keep_labels(labels)

    cat_show = cat
    if cat is not None and labels is not None and len(labels) < len(cat):
        try:
            cat_show = cat.get_labels(labels)
        except Exception:
            cat_show = cat

    # Fall back to the data panel if there is no segmentation to show.
    if segm_show is None and panels in ('segm', 'both'):
        panels = 'data'

    which = ['segm', 'data'] if panels == 'both' else [panels]
    if figsize is None:
        figsize = (20, 10) if len(which) == 2 else (12, 12)

    fig, axes = plt.subplots(1, len(which), figsize=figsize)
    axes = np.atleast_1d(axes)
    panel_axes = dict(zip(which, axes))

    if 'data' in panel_axes:
        ax_data = panel_axes['data']
        norm_kwargs = {'asinh_a': asinh_a} if stretch == 'asinh' else {}
        norm = simple_norm(np.nan_to_num(data_sub), stretch=stretch,
                           vmin=vmin, vmax=vmax, **norm_kwargs)
        ax_data.imshow(np.nan_to_num(data_sub), origin='lower', cmap='Greys_r',
                       norm=norm, interpolation='nearest')
        ax_data.set_title('Background-subtracted data')

    if 'segm' in panel_axes:
        ax_segm = panel_axes['segm']
        ax_segm.imshow(segm_show, origin='lower',
                       cmap=segm_show.make_cmap(seed=cmap_seed),
                       interpolation='nearest')
        ax_segm.set_title(f'Segmentation image ({segm_show.nlabels} labels)')

    # Kron apertures belong on the data; with no data panel, draw them over the
    # segmentation so the measured extent is still visible against the footprints.
    ax_kron = panel_axes.get('data', panel_axes.get('segm'))
    if show_kron_apertures and cat_show is not None and ax_kron is not None:
        try:
            cat_show.plot_kron_apertures(ax=ax_kron, color='white', lw=0.7)
        except Exception as error:
            print(f' -->> Kron apertures could not be drawn: {error}')

    # Optional annotations, on the first panel. Both are off by default: the
    # figure is meant to read as segments plus Kron apertures.
    ax = axes[0]
    if show_cutout_boxes or show_labels:
        for i in range(len(catalogue)):
            row = catalogue.iloc[i]
            if show_cutout_boxes:
                size = cutout_size if cutout_size is not None \
                    else row['suggested_cutout_size']
                ax.add_patch(plt.Rectangle((row['xc_parent'] - size / 2.0,
                                            row['yc_parent'] - size / 2.0),
                                           size, size, fill=False,
                                           edgecolor='cyan', linewidth=0.5))
            # Labelling every source in a crowded field is unusably slow, so
            # only the first `max_labels_to_plot` rows are annotated.
            if show_labels and i < max_labels_to_plot:
                label = f"ID{int(row['id'])}" if 'id' in catalogue.columns \
                    else f"ID{i + 1}"
                ax.annotate(label, xy=(row['xc'], row['yc']),
                            xytext=(row['xc'] + 20, row['yc'] + 20),
                            color='red', fontsize=6,
                            arrowprops=dict(arrowstyle='-', color='red',
                                            alpha=0.4, lw=0.4))

    suffix = (f' (first {max_labels_to_plot} labelled)'
              if show_labels and len(catalogue) > max_labels_to_plot else '')
    count = f'{len(catalogue)} detections{suffix}'
    if len(which) == 1:
        # One panel: fold the count into the axis title rather than leaving a
        # suptitle floating above it.
        ax.set_title(f'{ax.get_title()} - {count}')
    else:
        fig.suptitle(count)
    for axis in axes:
        axis.axis('off')

    if save_name is not None:
        plt.savefig(save_name, dpi=dpi, bbox_inches='tight')
        print(f' ++==>> Detection map saved to {save_name}')
    plt.show()
    plt.close(fig)


def save_field_products(imagename, bkg_image=None, bkg_rms=None,
                        data_sub=None, seg_maps=None, seg_maps_parent=None,
                        products=('data_sub', 'bkg', 'rms'),
                        save_path=None, prefix=None, overwrite=True,
                        verbose=1):
    """
    Write the full-mosaic intermediate products to FITS.

    Useful both for QA (blink the RMS map against the mosaic to see where the
    primary-beam response is eating your sensitivity) and for reuse: the RMS
    map written here can be fed straight back in as `rms_map=` on a later run
    to skip re-estimating it.

    Parameters
    ----------
    products : tuple of str
        Any of 'data_sub', 'bkg', 'rms', 'segm', 'segm_parent'.
    save_path : str, optional
        Output directory. Defaults to the directory of `imagename`.
    prefix : str, optional
        Output basename stem. Defaults to the image basename without '.fits'.

    Returns
    -------
    dict
        {product: path} for everything written.
    """
    header_in = fits.getheader(imagename)
    wcs_2d = WCS(header_in, naxis=2)

    if save_path is None:
        save_path = os.path.dirname(os.path.abspath(imagename))
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    if prefix is None:
        prefix = os.path.basename(imagename).replace('.fits', '')

    def _header_for(data, bunit=None):
        """2D celestial header, carrying over beam/units from the input."""
        header = wcs_2d.to_header()
        header['NAXIS'] = 2
        header['NAXIS1'] = data.shape[1]
        header['NAXIS2'] = data.shape[0]
        for key in ('BUNIT', 'BMAJ', 'BMIN', 'BPA', 'TELESCOP', 'OBJECT',
                    'RESTFRQ'):
            if key in header_in:
                header[key] = header_in[key]
        if bunit is not None:
            header['BUNIT'] = bunit
        return header

    available = {'data_sub': (data_sub, None, 'data minus background'),
                 'bkg': (bkg_image, None, 'estimated background'),
                 'rms': (bkg_rms, None, 'background RMS map'),
                 'segm': (None if seg_maps is None else seg_maps.data,
                          '', 'deblended segmentation labels'),
                 'segm_parent': (None if seg_maps_parent is None
                                 else seg_maps_parent.data,
                                 '', 'parent segmentation labels')}

    written = {}
    for product in products:
        if product not in available:
            print(f' -->> Unknown product "{product}"; skipping. Choose from '
                  f'{list(available)}.')
            continue
        data, bunit, comment = available[product]
        if data is None:
            print(f' -->> Product "{product}" was not computed; skipping.')
            continue
        outname = os.path.join(save_path, f'{prefix}_{product}.fits')
        hdu = fits.PrimaryHDU(data=np.asarray(data),
                              header=_header_for(np.asarray(data), bunit))
        hdu.header['COMMENT'] = f'morphen field_extraction: {comment}'
        hdu.writeto(outname, overwrite=overwrite)
        written[product] = outname
        if verbose >= 1:
            print(f' ++==>> {product:12s} -> {outname}')

    return written


def catalogue_to_regions(catalogue, filename, coord_frame='fk5',
                         shape='ellipse', cell_size=None, ell_size_factor=3.0,
                         color='green', label_column='id', verbose=1):
    """
    Write a ds9/CARTA region file for a field catalogue.

    The cheapest possible sanity check on a field run: load the region file over
    the mosaic and confirm the ellipses sit on real emission before committing
    CPU to hundreds of cutouts and fits.

    Note the ellipse position angle is photutils' `orientation` (counterclockwise
    from the +x axis), written straight through; it is close enough for visual
    QA but is not a rigorously converted sky position angle.
    """
    if cell_size is None:
        cell_size = float(catalogue['cell_size'].iloc[0]) \
            if 'cell_size' in catalogue.columns else 1.0

    lines = ['# Region file format: DS9 version 4.1',
             f'global color={color} width=1 font="helvetica 10 normal roman" '
             f'select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 '
             f'include=1 source=1',
             coord_frame]

    for _, row in catalogue.iterrows():
        text = ''
        if label_column is not None and label_column in catalogue.columns:
            text = f' # text={{{row[label_column]}}}'
        if shape == 'circle':
            radius = ell_size_factor * row['equivalent_radius_px'] * cell_size
            lines.append(f'circle({row["ra"]:.7f},{row["dec"]:.7f},'
                         f'{radius:.3f}")' + text)
        else:
            a = ell_size_factor * row['semimajor_px'] * cell_size
            b = ell_size_factor * row['semiminor_px'] * cell_size
            lines.append(f'ellipse({row["ra"]:.7f},{row["dec"]:.7f},'
                         f'{a:.3f}",{b:.3f}",{row["orientation_deg"]:.2f})'
                         + text)

    with open(filename, 'w') as handle:
        handle.write('\n'.join(lines) + '\n')

    if verbose >= 1:
        print(f' ++==>> {len(catalogue)} region(s) written to {filename}')
    return filename


# ---------------------------------------------------------------------------
# stage 2: catalogue -> cutouts
# ---------------------------------------------------------------------------

def _even_int(value):
    """Round up to the nearest even integer (Cutout2D likes symmetric boxes)."""
    return int(np.ceil(float(value) / 2.0) * 2)


def _image_shape_from_header(imagename):
    """(ny, nx) of a FITS image without loading the data."""
    header = fits.getheader(imagename)
    return int(header['NAXIS2']), int(header['NAXIS1'])


def _resolve_cutout_size(requested, row, imagename, size_units,
                         cutout_size_factor, cutout_size_min, cutout_size_max,
                         detection_cell_size):
    """
    Cutout box size in pixels *of this image*.

    'auto' sizes the box from the parent island's bounding box in **arcsec**, so
    a heterogeneous imagelist (different frequencies, different cell sizes) gets
    boxes covering the same piece of sky rather than the same pixel count.
    """
    if isinstance(requested, str) and requested == 'auto':
        parent_size = float(row.get('parent_size_px', np.nan))
        if not np.isfinite(parent_size) or parent_size <= 0:
            # An externally supplied source list carries no island size, so
            # there is nothing to scale from; fall back to the minimum box.
            return _even_int(cutout_size_min)
        size_arcsec = cutout_size_factor * parent_size * detection_cell_size
        cell = get_cell_size(imagename)
        size_px = size_arcsec / cell
    elif size_units == 'arcsec':
        cell = get_cell_size(imagename)
        size_px = float(requested) / cell
    else:
        size_px = float(requested)

    size_px = np.clip(size_px, cutout_size_min, cutout_size_max)
    return _even_int(size_px)


def cutouts_from_catalogue(catalogue, imagelist, residuallist=None,
                           modellist=None,
                           cutout_size='auto', size_units='pixel',
                           cutout_size_factor=4.0,
                           cutout_size_min=64, cutout_size_max=2048,
                           cutout_on='parent', edge_policy='shrink',
                           custom_save_path=None, name_column='iau_name',
                           prefix='', suffix='',
                           do_plot=False, rms=None, plot_kwargs=None,
                           skip_existing=False, verbose=1):
    """
    Build per-source cutouts from a field catalogue, for one or many images.

    This generalises the "loop over MFS images x loop over RA/Dec list" pattern
    from the older notebooks. It accepts either a catalogue produced by
    `field_source_ext` or a plain user-supplied source list (see below), so an
    external catalogue goes through exactly the same code path. Cutting is
    delegated to `t_cutout_2D_radec` (data_io.py), the merged/definitive cutout
    routine.

    Parameters
    ----------
    catalogue : pandas.DataFrame or dict
        Either an `SE.catalogue`, or a dict/DataFrame with at least `ra`/`dec`
        in **degrees** (or `ra_str`/`dec_str` in sexagesimal, which are parsed
        with `conver_str_coords`) and optionally a name column.
    imagelist : str or list of str
        Image(s) to cut. A single path is accepted.
    residuallist, modellist : str or list of str, optional
        Matching residual / model image(s), cut on the identical sky position
        and box alongside each image.
    cutout_size : int or 'auto'
        'auto' sizes each box from the parent island's extent (see
        `_resolve_cutout_size`); an int is taken in `size_units`.
    size_units : {'pixel', 'arcsec'}
        Interpretation of an integer `cutout_size`. 'arcsec' is the right choice
        for a multi-resolution imagelist.
    cutout_on : {'parent', 'component'}
        'parent' (default) makes one cutout per deblended *island*, centred on
        the island centroid and large enough to hold all of its components -
        i.e. one postage stamp per astronomical source, which is what the
        downstream morphometry expects. 'component' makes one per component.
    edge_policy : {'shrink', 'skip', 'pad'}
        What to do when the requested box does not fit inside the image.
        'shrink' reduces the box to the largest even one that fits; 'skip'
        drops the source; 'pad' keeps the requested box and NaN-pads it
        (`t_cutout_2D_radec(mode='partial')`).
    custom_save_path : str, optional
        Output directory. Defaults to `<dirname(image)>/cutouts_<image stem>/`.
    name_column : str
        Catalogue column used to name the cutouts (default IAU designation).
    skip_existing : bool
        Do not re-cut a source whose output file already exists.

    Returns
    -------
    dict
        {'imagelist', 'residuallist', 'per_image', 'catalogue', 'skipped',
         'failed'}. `imagelist`/`residuallist` correspond to the *first* image
         of `imagelist` and are ready to be handed to `structural_morphology`.
    """
    if isinstance(imagelist, str):
        imagelist = [imagelist]
    if isinstance(residuallist, str):
        residuallist = [residuallist]
    if isinstance(modellist, str):
        modellist = [modellist]
    if residuallist is None:
        residuallist = [None] * len(imagelist)
    if len(residuallist) != len(imagelist):
        raise ValueError(f'residuallist has {len(residuallist)} entries but '
                         f'imagelist has {len(imagelist)}.')
    if modellist is not None and len(modellist) != len(imagelist):
        raise ValueError(f'modellist has {len(modellist)} entries but '
                         f'imagelist has {len(imagelist)}.')

    table = _normalise_source_table(catalogue, cutout_on=cutout_on,
                                    name_column=name_column, verbose=verbose)

    detection_cell_size = float(table['cell_size'].iloc[0]) \
        if 'cell_size' in table.columns else 1.0

    plot_kwargs = {} if plot_kwargs is None else dict(plot_kwargs)

    per_image = {}
    per_residual = {}
    skipped = []
    failed = []
    out_catalogue = table.copy()

    for j, imagename in enumerate(imagelist):
        residualname = residuallist[j]
        if custom_save_path is None:
            stem = os.path.basename(imagename).replace('.fits', '')
            save_path = os.path.join(os.path.dirname(os.path.abspath(imagename)),
                                     f'cutouts_{stem}')
        else:
            save_path = custom_save_path
        if not os.path.exists(save_path):
            os.makedirs(save_path)

        ny, nx = _image_shape_from_header(imagename)
        image_wcs = WCS(fits.getheader(imagename), naxis=2)

        if verbose >= 1:
            print(f' ++==>> Cutting {len(table)} source(s) from '
                  f'{os.path.basename(imagename)} into {save_path}')

        paths = []
        residual_paths = []
        sizes = []
        for k in tqdm(range(len(table)), disable=(verbose < 1)):
            row = table.iloc[k]
            ra_deg, dec_deg = float(row['ra_cut']), float(row['dec_cut'])
            name = f"{prefix}{row['cutout_name']}{suffix}"

            size_px = _resolve_cutout_size(cutout_size, row, imagename,
                                           size_units, cutout_size_factor,
                                           cutout_size_min, cutout_size_max,
                                           detection_cell_size)

            # --- edge handling, before calling into t_cutout_2D_radec -------
            try:
                x_px, y_px = image_wcs.world_to_pixel(
                    SkyCoord(ra=ra_deg * u.degree, dec=dec_deg * u.degree,
                             frame='icrs'))
                x_px, y_px = float(x_px), float(y_px)
            except Exception as error:
                failed.append((imagename, name, f'WCS: {error}'))
                paths.append(None)
                residual_paths.append(None)
                sizes.append(size_px)
                continue

            fit_half = min(x_px, y_px, nx - 1 - x_px, ny - 1 - y_px)
            if fit_half < size_px / 2.0:
                if edge_policy == 'skip':
                    skipped.append((imagename, name, 'edge'))
                    paths.append(None)
                    residual_paths.append(None)
                    sizes.append(size_px)
                    continue
                elif edge_policy == 'shrink':
                    new_size = _even_int(max(0.0, 2.0 * np.floor(fit_half)))
                    if new_size < cutout_size_min or new_size <= 0:
                        skipped.append((imagename, name,
                                        f'edge (would shrink to {new_size})'))
                        paths.append(None)
                        residual_paths.append(None)
                        sizes.append(size_px)
                        continue
                    if verbose >= 1:
                        print(f'        {name}: box shrunk '
                              f'{size_px} -> {new_size} px (image edge).')
                    size_px = new_size
                # 'pad' keeps the full box via t_cutout_2D_radec(mode='partial')

            if skip_existing:
                existing = glob.glob(os.path.join(save_path, f'{name}*image.fits'))
                if existing:
                    paths.append(existing[0])
                    candidate = existing[0].replace('-image.fits',
                                                    '-residual.fits')
                    residual_paths.append(candidate
                                          if os.path.exists(candidate) else None)
                    sizes.append(size_px)
                    continue

            try:
                result = t_cutout_2D_radec(
                    imagename=imagename,
                    residualname=residualname,
                    modelname=modellist[j] if modellist is not None else None,
                    cutout_size=size_px,
                    ra_f=ra_deg, dec_f=dec_deg,
                    mode=('partial' if edge_policy == 'pad' else 'trim'),
                    custom_save_path=save_path,
                    custom_save_name=name,
                    special_name='',
                    return_paths=True, verbose=0)
                cutout_filename = result['image']
            except Exception as error:
                failed.append((imagename, name, str(error)))
                paths.append(None)
                residual_paths.append(None)
                sizes.append(size_px)
                continue

            paths.append(cutout_filename)
            residual_paths.append(result['residual'])
            sizes.append(size_px)

            if do_plot:
                try:
                    plot_rms = rms
                    if plot_rms is None and 'local_rms' in table.columns:
                        plot_rms = float(row['local_rms'])
                    eimshow(cutout_filename, add_beam=True, vmax_factor=0.9,
                            plot_colorbar=True, rms=plot_rms,
                            plot_title=os.path.basename(
                                cutout_filename).replace('.fits', ''),
                            save_name=cutout_filename.replace('.fits', '_map'),
                            **plot_kwargs)
                    plt.show()
                    plt.close()
                except Exception as error:
                    print(f' -->> Plot failed for {name}: {error}')

        per_image[imagename] = paths
        per_residual[imagename] = residual_paths
        out_catalogue[f'cutout_{j}'] = paths
        out_catalogue[f'cutout_residual_{j}'] = residual_paths
        out_catalogue[f'cutout_size_{j}'] = sizes

    # The first image defines the imagelist handed downstream.
    primary = per_image[imagelist[0]]
    primary_residual = per_residual[imagelist[0]]

    out_catalogue['cutout_image'] = primary
    out_catalogue['cutout_residual'] = primary_residual

    imagelist_out = [p for p in primary if p is not None]
    residuallist_out = [r for p, r in zip(primary, primary_residual)
                        if p is not None]
    if all(r is None for r in residuallist_out):
        residuallist_out = None

    if verbose >= 1:
        print(f' ++==>> {len(imagelist_out)} cutout(s) written per image; '
              f'{len(skipped)} skipped, {len(failed)} failed.')
        for entry in failed[:10]:
            print(f'        FAILED {entry[1]}: {entry[2]}')

    return {'imagelist': imagelist_out,
            'residuallist': residuallist_out,
            'per_image': per_image,
            'per_residual': per_residual,
            'catalogue': out_catalogue,
            'skipped': skipped,
            'failed': failed}


# ---------------------------------------------------------------------------
# stage 2b: multi-band cutouts
# ---------------------------------------------------------------------------

def _wsclean_prefix(path):
    """
    The WSClean run prefix of an image, i.e. everything before its
    ``-MFS-``/``-0000-`` product token.

    ``<prefix>-MFS-image-pb.fits`` and ``<prefix>-0003-residual.fits`` both
    reduce to ``<prefix>``, which is what identifies the band a sub-band image
    belongs to. `prepare_data` returns flat lists with no grouping, so this is
    how the sub-bands are matched back to their parent MFS.
    """
    import re
    return re.sub(r'-(MFS|\d{4})-(image|residual|model|dirty|psf).*$', '',
                  os.path.basename(path))


def _wsclean_token(path):
    """The product token of an image: 'MFS', '0000', ... or None."""
    import re
    match = re.search(r'-(MFS|\d{4})-(image|residual|model|dirty|psf)',
                      os.path.basename(path))
    return match.group(1) if match else None


def build_band_spec(MFS_images, MFS_residuals=None, imagelist=None,
                    residuallist=None, freqlist_MFS=None, freqlist=None,
                    band_names=None, verbose=1):
    """
    Group the flat lists `prepare_data` returns into one entry per band.

    `mlibs.prepare_data` concatenates every band's sub-bands into a single
    `imagelist`, with no record of which sub-band came from which run. They are
    regrouped here by WSClean run prefix (see `_wsclean_prefix`), which is
    exact rather than positional - it does not care about ordering or about how
    many sub-bands each band happens to have.

    Parameters
    ----------
    MFS_images : list of str
        One MFS image per band. Defines the bands.
    MFS_residuals : list of str, optional
        Matching MFS residuals.
    imagelist, residuallist : list of str, optional
        All sub-band images/residuals across all bands, in any order.
    freqlist_MFS, freqlist : array, optional
        Frequencies in Hz for the MFS images and the sub-bands.
    band_names : list of str, optional
        Names used in the cutout filenames. Defaults to the MFS frequency
        rounded to 0.1 GHz (e.g. '6.0GHz'). Override this when two bands would
        collide - two Briggs weightings of the same band, or overlapping
        coverage - since the name has to be unique per band.

    Returns
    -------
    list of dict
        One per band: name, freq, mfs_image, mfs_residual, sub_images,
        sub_residuals, sub_freqs, sub_tokens, cell_size.
    """
    MFS_images = list(np.atleast_1d(MFS_images))
    MFS_residuals = (list(np.atleast_1d(MFS_residuals))
                     if MFS_residuals is not None else [None] * len(MFS_images))
    imagelist = list(np.atleast_1d(imagelist)) if imagelist is not None else []
    residuallist = (list(np.atleast_1d(residuallist))
                    if residuallist is not None else [None] * len(imagelist))
    freqlist_MFS = (np.atleast_1d(freqlist_MFS)
                    if freqlist_MFS is not None else np.full(len(MFS_images), np.nan))
    freqlist = (np.atleast_1d(freqlist)
                if freqlist is not None else np.full(len(imagelist), np.nan))

    # sub-band lookup by run prefix
    sub_by_prefix = {}
    for k, path in enumerate(imagelist):
        sub_by_prefix.setdefault(_wsclean_prefix(path), []).append(k)

    # --- band names ---------------------------------------------------------
    if band_names is None:
        band_names = []
        for j in range(len(MFS_images)):
            freq = freqlist_MFS[j] if j < len(freqlist_MFS) else np.nan
            band_names.append(f'{freq / 1e9:.1f}GHz' if np.isfinite(freq)
                              else f'band{j}')
        if len(set(band_names)) != len(band_names):
            print(' -->> WARNING: frequency-derived band names are not unique '
                  f'({band_names}); disambiguating with an index. Pass '
                  '`band_names=` to control this.')
            band_names = [f'{n}_{j}' for j, n in enumerate(band_names)]
    band_names = list(band_names)

    bands = []
    for j, mfs in enumerate(MFS_images):
        prefix = _wsclean_prefix(mfs)
        idx = sorted(sub_by_prefix.get(prefix, []),
                     key=lambda k: _wsclean_token(imagelist[k]) or '')
        try:
            cell_size = get_cell_size(mfs)
        except Exception:
            cell_size = np.nan
        band = {'name': band_names[j],
                'freq': float(freqlist_MFS[j]) if j < len(freqlist_MFS) else np.nan,
                'prefix': prefix,
                'mfs_image': mfs,
                'mfs_residual': MFS_residuals[j] if j < len(MFS_residuals) else None,
                'sub_images': [imagelist[k] for k in idx],
                'sub_residuals': [residuallist[k] for k in idx],
                'sub_freqs': [float(freqlist[k]) if k < len(freqlist) else np.nan
                              for k in idx],
                'sub_tokens': [_wsclean_token(imagelist[k]) for k in idx],
                'cell_size': cell_size}
        bands.append(band)
        if verbose >= 1:
            print(f" ++==>> band '{band['name']}': "
                  f"{len(band['sub_images'])} sub-band(s), "
                  f"cell = {cell_size:.4f} arcsec")

    matched = sum(len(b['sub_images']) for b in bands)
    if matched != len(imagelist):
        print(f' -->> WARNING: {len(imagelist) - matched} sub-band image(s) '
              f'did not match any MFS run prefix and will be ignored.')

    return bands


def multiband_cutouts_from_catalogue(
        catalogue, bands, workdir=None,
        cutout_size='auto', size_units='pixel', cutout_size_factor=4.0,
        cutout_size_min=64, cutout_size_max=2048,
        cutout_on='parent', edge_policy='shrink',
        correct_shift=False, align=None, ref_band=0,
        shift_correction_mode='chi2_shift', shift_granularity='band',
        include_mfs=True, include_subbands=True,
        name_column='iau_name', skip_existing=False,
        manifest_name='manifest.csv', write_manifest=True, verbose=1):
    """
    Cut every catalogued source out of every band, MFS and sub-band alike.

    One directory per source, so everything needed to analyse a source sits
    together and is deleted in one call when you are done:

        <workdir>/<iau_name>/<iau_name>-<band>-<token>-{image,residual}.fits
        <workdir>/manifest.csv

    The ``<token>`` is WSClean's own ('MFS', '0000', ...), carried through by
    `t_cutout_2D_radec`, so a cutout is traceable back to the exact input.

    Parameters
    ----------
    catalogue : pandas.DataFrame
        `SE.catalogue`, or anything `_normalise_source_table` accepts.
    bands : list of dict
        From `build_band_spec`.
    workdir : str
        Root output directory. Defaults to `<dirname(reference MFS)>/cutouts/`.
    cutout_size : int or 'auto'
        Resolved **once per source, in arcsec**, then converted to pixels per
        band. Bands on different cell sizes therefore cover the same sky rather
        than the same pixel count.
    correct_shift : bool
        **False by default.** The master on/off switch for astrometric
        alignment, matching the argument of the same name on
        `t_cutout_2D_radec` / `cutout_2D_radec`. While it is False nothing is
        ever shifted and every cutout is placed by WCS alone, whatever `align`
        says. Shifting resamples your data, so it has to be asked for.
    align : {None, 'global', 'band'}
        *Which* alignment scheme to use once `correct_shift=True`. Ignored
        while `correct_shift=False`.

        - ``None`` (default) - resolves to ``'global'`` when `correct_shift`
          is enabled.
        - ``'global'`` - the reference band's MFS cutout is the anchor for that
          source; every other band's MFS and every sub-band is shifted onto it.
          This is what you want before comparing bands pixel-by-pixel (spectral
          index, say), but it is only meaningful when the bands share a pixel
          grid. The cell sizes are checked and the mode is refused if they
          differ.
        - ``'band'`` - each band's sub-bands align to that band's own MFS
          cutout; bands are placed by WCS alone. More conservative: the
          cross-correlation is only ever between images of equal resolution.
    shift_granularity : {'band', 'image'}
        Under ``align='global'``, whether the shift is measured once per band or
        once per image.

        - ``'band'`` (default) - measure it on the band's MFS cutout, which is
          the highest-SNR image the band has, and reuse that one value for every
          sub-band of the band. Sub-bands of a single imaging run share a pixel
          grid and a self-cal solution, so there is no real sub-band-to-sub-band
          offset to find; measuring one anyway just adds noise.
        - ``'image'`` - the historical behaviour: estimate independently for
          every sub-band. Useful for diagnosing a band whose sub-bands disagree,
          since the scatter between them measures the estimator's noise on your
          data.
    shift_correction_mode : str
        Estimator passed to `t_cutout_2D_radec` / `estimate_image_shift`:
        'chi2', 'mi', 'xcorr', 'chi2_shift', 'auto' (default) or 'ensemble'.
    ref_band : int
        Index into `bands` of the alignment reference. Default 0.
    include_mfs, include_subbands : bool
        Cut the MFS and/or the sub-band images.
    skip_existing : bool
        Leave already-written cutouts alone (resume a long run).

    Returns
    -------
    dict
        - ``manifest`` : DataFrame, one row per (source, band, token), with
          image/residual paths, frequency, box size, applied shift and status.
          This is the thing to iterate in the analysis loop - not the
          filesystem.
        - ``sources`` : dict keyed by source name - **the unit a multi-band
          analysis works on**. Each entry holds that one source's products as
          flat, index-aligned lists, split MFS vs sub-band exactly the way
          `prepare_data` splits them:

          ===================== =================================================
          ``imagelist``         sub-band images, all bands
          ``residuallist``      matching sub-band residuals
          ``freqlist``          matching frequencies [Hz]
          ``bandlist``          matching band names
          ``tokenlist``         matching WSClean tokens ('0000', ...)
          ``MFS_imagelist``     MFS images, one per band
          ``MFS_residuallist``  matching MFS residuals
          ``freqlist_MFS``      matching MFS frequencies [Hz]
          ``MFS_bandlist``      matching band names
          ``dir``               the source's directory
          ``bands``             the same, grouped per band
          ===================== =================================================

        - ``imagelist`` / ``residuallist`` / ``MFS_imagelist`` /
          ``MFS_residuallist`` : the same, concatenated over *all* sources.
          Rarely what you want - indexing these walks across source boundaries.
        - ``workdir``, ``failed``, ``skipped``.

    Examples
    --------
    One source::

        src = out['sources']['J095551.44+694044.3']   # or list(out['sources'])[0]

        imagelist_c, residuallist_c = src['imagelist'], src['residuallist']
        MFS_images_c = src['MFS_imagelist']
        MFS_residuals_c = src['MFS_residuallist']

        idx = 5
        d = mp.read_data(filename=imagelist_c[idx],
                         residualname=residuallist_c[idx])
        rms = d.rms_res
        print(src['bandlist'][idx], src['tokenlist'][idx], src['freqlist'][idx])

    Both lists run in band order::

        imagelist     : 6.0GHz-0000 ... 6.0GHz-0003, 22.0GHz-0000, ...
        MFS_imagelist : 6.0GHz-MFS, 22.0GHz-MFS, 33.0GHz-MFS
    """
    if align not in (None, 'global', 'band'):
        raise ValueError("align must be None, 'global' or 'band'.")

    # `correct_shift` is the switch; `align` only picks the scheme. Nothing is
    # ever shifted unless the caller explicitly asks for it.
    if not correct_shift:
        if align is not None:
            print(f" -->> WARNING: align={align!r} was given but "
                  "correct_shift=False, so no alignment will be applied. "
                  "Pass correct_shift=True to enable it.")
        align = None
    elif align is None:
        align = 'global'

    # align='band' anchors each band's sub-bands to that band's own MFS cutout,
    # so with no MFS there is nothing to anchor to and every sub-band would be
    # left unshifted without saying so.
    if align == 'band' and not include_mfs:
        print(" -->> WARNING: align='band' needs the MFS cutouts as its "
              "anchors, but include_mfs=False. No alignment will be applied. "
              "Use align='global' or set include_mfs=True.")
        align = None
    if not bands:
        raise ValueError('`bands` is empty; build it with build_band_spec().')
    if not (0 <= ref_band < len(bands)):
        raise ValueError(f'ref_band={ref_band} is out of range for '
                         f'{len(bands)} band(s).')

    # --- grid check: cross-band alignment is only meaningful on a common grid
    if align == 'global' and len(bands) > 1:
        cells = np.asarray([b['cell_size'] for b in bands], dtype=float)
        if not np.allclose(cells, cells[ref_band], rtol=1e-4, equal_nan=False):
            print(' -->> WARNING: bands do not share a pixel scale '
                  f'({np.round(cells, 5)} arcsec). Cross-band peak alignment '
                  "between different grids is not meaningful; falling back to "
                  "align='band'.")
            align = 'band'

    table = _normalise_source_table(catalogue, cutout_on=cutout_on,
                                    name_column=name_column, verbose=verbose)
    detection_cell_size = float(table['cell_size'].iloc[0]) \
        if 'cell_size' in table.columns else float(bands[ref_band]['cell_size'])

    if workdir is None:
        workdir = os.path.join(
            os.path.dirname(os.path.abspath(bands[ref_band]['mfs_image'])),
            'cutouts')
    if not os.path.exists(workdir):
        os.makedirs(workdir)

    # Pre-compute each band's image shape and WCS once, not per source.
    for band in bands:
        band['_shape'] = _image_shape_from_header(band['mfs_image'])
        band['_wcs'] = WCS(fits.getheader(band['mfs_image']), naxis=2)

    rows = []
    sources = {}
    failed = []
    skipped = []

    if verbose >= 1:
        n_prod = sum((1 if include_mfs else 0) +
                     (len(b['sub_images']) if include_subbands else 0)
                     for b in bands)
        print(f' ++==>> {len(table)} source(s) x {n_prod} product(s) '
              f'-> {workdir}')
        if align:
            print(f"        alignment: {align} "
                  f"(reference band '{bands[ref_band]['name']}')")

    for k in tqdm(range(len(table)), disable=(verbose < 1)):
        row = table.iloc[k]
        name = str(row['cutout_name'])
        ra_deg, dec_deg = float(row['ra_cut']), float(row['dec_cut'])
        source_dir = os.path.join(workdir, name)
        if not os.path.exists(source_dir):
            os.makedirs(source_dir)

        # box size, resolved once per source in arcsec
        size_ref = _resolve_cutout_size(cutout_size, row,
                                        bands[ref_band]['mfs_image'],
                                        size_units, cutout_size_factor,
                                        cutout_size_min, cutout_size_max,
                                        detection_cell_size)
        size_arcsec = size_ref * float(bands[ref_band]['cell_size'])

        anchor = None          # global alignment reference for this source
        band_anchor = {}       # per-band MFS cutout, for align='band'
        band_shift = {}        # per-band shift, measured once on the MFS
        rows_by_band = {}      # collected per band, emitted in band order
        source_entry = {'dir': source_dir, 'ra': ra_deg, 'dec': dec_deg,
                        'bands': {}}

        def _cut(band, image, residual, token, freq, ref_cutout,
                 precomputed=None):
            """Cut one product; returns the manifest row."""
            size_px = _even_int(np.clip(size_arcsec / band['cell_size'],
                                        cutout_size_min, cutout_size_max))
            entry = {'source': name, 'ra': ra_deg, 'dec': dec_deg,
                     'band': band['name'], 'band_freq': band['freq'],
                     'token': token, 'freq': freq,
                     'input_image': image, 'cell_size': band['cell_size'],
                     'image': None, 'residual': None,
                     'cutout_size': size_px, 'dx': 0, 'dy': 0,
                     'shift_dy': 0.0, 'shift_dx': 0.0,
                     'shift_err_dy': np.nan, 'shift_err_dx': np.nan,
                     'shift_method': '', 'shift_flags': '',
                     'is_reference': ref_cutout is None and precomputed is None,
                     'status': 'ok'}

            ny, nx = band['_shape']
            try:
                x_px, y_px = band['_wcs'].world_to_pixel(
                    SkyCoord(ra=ra_deg * u.degree, dec=dec_deg * u.degree,
                             frame='icrs'))
            except Exception as error:
                entry['status'] = f'wcs: {error}'
                failed.append((name, band['name'], token, entry['status']))
                return entry

            fit_half = min(float(x_px), float(y_px),
                           nx - 1 - float(x_px), ny - 1 - float(y_px))
            if fit_half < size_px / 2.0:
                if edge_policy == 'skip':
                    entry['status'] = 'skipped: edge'
                    skipped.append((name, band['name'], token, 'edge'))
                    return entry
                if edge_policy == 'shrink':
                    new_size = _even_int(max(0.0, 2.0 * np.floor(fit_half)))
                    if new_size < cutout_size_min or new_size <= 0:
                        entry['status'] = 'skipped: edge'
                        skipped.append((name, band['name'], token, 'edge'))
                        return entry
                    entry['cutout_size'] = size_px = new_size

            save_name = f'{name}-{band["name"]}'
            expected = os.path.join(source_dir, f'{save_name}-{token}-image.fits')
            if skip_existing and os.path.exists(expected):
                entry['image'] = expected
                candidate = expected.replace('-image.fits', '-residual.fits')
                entry['residual'] = candidate if os.path.exists(candidate) else None
                entry['status'] = 'existing'
                return entry

            try:
                result = t_cutout_2D_radec(
                    imagename=image, residualname=residual,
                    ra_f=ra_deg, dec_f=dec_deg, cutout_size=size_px,
                    mode=('partial' if edge_policy == 'pad' else 'trim'),
                    correct_shift=(ref_cutout is not None
                                   or precomputed is not None),
                    ref_cutout_image=ref_cutout,
                    precomputed_shift=precomputed,
                    shift_correction_mode=shift_correction_mode,
                    custom_save_path=source_dir, custom_save_name=save_name,
                    return_paths=True, verbose=0)
            except Exception as error:
                entry['status'] = f'failed: {error}'
                failed.append((name, band['name'], token, str(error)))
                return entry

            entry['image'] = result['image']
            entry['residual'] = result['residual']
            entry['dx'], entry['dy'] = result['offset']
            entry['shift_dy'], entry['shift_dx'] = result['shift']
            sr = result.get('shift_result')
            if sr is not None:
                entry['shift_err_dy'] = sr.dy_err
                entry['shift_err_dx'] = sr.dx_err
                entry['shift_method'] = sr.method
                entry['shift_flags'] = ';'.join(sr.flags)
            return entry

        # --- the reference band's MFS goes first: it is the anchor ----------
        order = [ref_band] + [j for j in range(len(bands)) if j != ref_band]

        for j in order:
            band = bands[j]
            band_rows = []

            if include_mfs and band['mfs_image'] is not None:
                ref = anchor if (align == 'global' and j != ref_band) else None
                entry = _cut(band, band['mfs_image'], band['mfs_residual'],
                             'MFS', band['freq'], ref)
                band_rows.append(entry)
                if entry['image'] is not None:
                    band_anchor[j] = entry['image']
                    # Sub-bands from one imaging run sit on the same grid and
                    # share one self-cal solution, so the offset is a property
                    # of the band, not of the sub-band: measured on NGC 7469,
                    # sub-band-vs-own-MFS is 0.00-0.06 px for every estimator.
                    # Measure it once here, on the highest-SNR image the band
                    # has, and reuse it -- re-estimating per sub-band only adds
                    # scatter (0.09 px in dx on that same data).
                    band_shift[j] = (entry['shift_dy'], entry['shift_dx'])
                    if j == ref_band and anchor is None:
                        anchor = entry['image']

            if include_subbands:
                for i, sub in enumerate(band['sub_images']):
                    token = band['sub_tokens'][i] or f'{i:04d}'
                    pre = None
                    if align == 'global':
                        ref = anchor
                        if shift_granularity == 'band' and j in band_shift:
                            pre, ref = band_shift[j], None
                    elif align == 'band':
                        ref = band_anchor.get(j)
                    else:
                        ref = None
                    entry = _cut(band, sub, band['sub_residuals'][i], token,
                                 band['sub_freqs'][i], ref, precomputed=pre)
                    band_rows.append(entry)

            rows_by_band[j] = band_rows
            ok = [e for e in band_rows if e['image'] is not None]
            source_entry['bands'][band['name']] = {
                'freq': band['freq'],
                'mfs_image': next((e['image'] for e in ok
                                   if e['token'] == 'MFS'), None),
                'mfs_residual': next((e['residual'] for e in ok
                                      if e['token'] == 'MFS'), None),
                'sub_images': [e['image'] for e in ok if e['token'] != 'MFS'],
                'sub_residuals': [e['residual'] for e in ok
                                  if e['token'] != 'MFS'],
                'sub_freqs': [e['freq'] for e in ok if e['token'] != 'MFS']}

        # Emit in band order, not in the order the cuts happened -- the
        # reference band is processed first because it is the alignment anchor,
        # which is an implementation detail and should not leak into the lists.
        ordered = [e for j in range(len(bands)) for e in rows_by_band.get(j, [])]
        rows.extend(ordered)

        # Flat, index-aligned lists for THIS source, split MFS vs sub-band to
        # mirror `prepare_data`'s own convention (MFS_images vs imagelist).
        ok_all = [e for e in ordered if e['image'] is not None]
        ok_mfs = [e for e in ok_all if e['token'] == 'MFS']
        ok_sub = [e for e in ok_all if e['token'] != 'MFS']
        source_entry['imagelist'] = [e['image'] for e in ok_sub]
        source_entry['residuallist'] = [e['residual'] for e in ok_sub]
        source_entry['freqlist'] = [e['freq'] for e in ok_sub]
        source_entry['bandlist'] = [e['band'] for e in ok_sub]
        source_entry['tokenlist'] = [e['token'] for e in ok_sub]
        source_entry['MFS_imagelist'] = [e['image'] for e in ok_mfs]
        source_entry['MFS_residuallist'] = [e['residual'] for e in ok_mfs]
        source_entry['freqlist_MFS'] = [e['freq'] for e in ok_mfs]
        source_entry['MFS_bandlist'] = [e['band'] for e in ok_mfs]

        sources[name] = source_entry

    manifest = pd.DataFrame(rows)
    manifest_path = None
    if write_manifest and len(manifest):
        manifest_path = os.path.join(workdir, manifest_name)
        manifest.to_csv(manifest_path, index=False)

    # The same lists concatenated over all sources, split MFS vs sub-band as
    # above. Rarely what you want -- indexing these walks across source
    # boundaries; the per-source entries in `sources` are the useful unit.
    ok = manifest[manifest['image'].notna()] if len(manifest) else manifest
    sub = ok[ok['token'] != 'MFS'] if len(ok) else ok
    mfs = ok[ok['token'] == 'MFS'] if len(ok) else ok
    imagelist_out = list(sub['image']) if len(sub) else []
    residuallist_out = list(sub['residual']) if len(sub) else []
    mfs_imagelist_out = list(mfs['image']) if len(mfs) else []
    mfs_residuallist_out = list(mfs['residual']) if len(mfs) else []

    if verbose >= 1:
        n_ok = int((manifest['image'].notna()).sum()) if len(manifest) else 0
        print(f' ++==>> {n_ok}/{len(manifest)} cutout(s) written across '
              f'{len(sources)} source(s); {len(skipped)} skipped, '
              f'{len(failed)} failed.')
        if manifest_path:
            print(f' ++==>> manifest -> {manifest_path}')
        for entry in failed[:10]:
            print(f'        FAILED {entry[0]} {entry[1]} {entry[2]}: {entry[3]}')

    return {'manifest': manifest, 'sources': sources, 'workdir': workdir,
            'manifest_path': manifest_path, 'bands': bands,
            'imagelist': imagelist_out, 'residuallist': residuallist_out,
            'MFS_imagelist': mfs_imagelist_out,
            'MFS_residuallist': mfs_residuallist_out,
            'failed': failed, 'skipped': skipped}


def cleanup_multiband_cutouts(result, sources=None, keep_tokens=None,
                              dry_run=True, verbose=1):
    """
    Delete cutouts once a source has been analysed.

    Cutouts are reproducible from the catalogue plus the parent images, so
    there is no reason to hoard them: a few hundred sources across several
    bands runs to tens of GB. The manifest records what existed either way.

    Parameters
    ----------
    result : dict
        The return value of `multiband_cutouts_from_catalogue`.
    sources : str or list of str, optional
        Which sources to clear. Default: all of them.
    keep_tokens : tuple of str, optional
        Product tokens to spare, e.g. ``('MFS',)`` to drop the sub-bands and
        keep the broadband stamps. Default: keep nothing.
    dry_run : bool
        **True by default** - reports what would be removed without touching
        anything. Pass ``dry_run=False`` to actually delete.

    Returns
    -------
    list of str
        The paths removed (or that would be removed).
    """
    import shutil

    manifest = result['manifest']
    if sources is None:
        sources = list(result['sources'])
    elif isinstance(sources, str):
        sources = [sources]

    targets = []
    for name in sources:
        rows = manifest[manifest['source'] == name]
        if keep_tokens is not None:
            rows = rows[~rows['token'].isin(list(keep_tokens))]
        for _, row in rows.iterrows():
            for key in ('image', 'residual'):
                path = row[key]
                if isinstance(path, str) and os.path.exists(path):
                    targets.append(path)
        # drop the directory too when nothing is being kept
        if keep_tokens is None:
            source_dir = result['sources'][name]['dir']
            if os.path.isdir(source_dir):
                targets.append(source_dir)

    if dry_run:
        if verbose >= 1:
            print(f' -->> DRY RUN: {len(targets)} path(s) would be removed. '
                  f'Pass dry_run=False to delete.')
            for path in targets[:10]:
                print(f'        {path}')
        return targets

    removed = []
    for path in targets:
        try:
            if os.path.isdir(path):
                shutil.rmtree(path)
            else:
                os.remove(path)
            removed.append(path)
        except Exception as error:
            print(f' -->> could not remove {path}: {error}')
    if verbose >= 1:
        print(f' ++==>> removed {len(removed)} path(s).')
    return removed


def _normalise_source_table(catalogue, cutout_on='parent',
                            name_column='iau_name', verbose=1):
    """
    Coerce whatever the user passed into a table with `ra_cut`, `dec_cut`,
    `cutout_name` and `parent_size_px`.

    Accepts a full `field_source_ext` catalogue, or a bare dict/DataFrame with
    `ra`/`dec` in degrees or `ra_str`/`dec_str` in sexagesimal - the latter is
    the shape of the hand-maintained RA_LIST/DEC_LIST/SOURCES_LIST used in the
    older notebooks.
    """
    if isinstance(catalogue, dict):
        catalogue = pd.DataFrame(catalogue)
    table = catalogue.copy().reset_index(drop=True)

    if 'ra' not in table.columns or 'dec' not in table.columns:
        if 'ra_str' in table.columns and 'dec_str' in table.columns:
            converted = [conver_str_coords(r, d) for r, d in
                         zip(table['ra_str'], table['dec_str'])]
            table['ra'] = [c[0] for c in converted]
            table['dec'] = [c[1] for c in converted]
        else:
            raise ValueError("Catalogue needs 'ra'/'dec' (degrees) or "
                             "'ra_str'/'dec_str' (sexagesimal) columns.")

    # one row per island, keyed on the brightest component
    if cutout_on == 'parent' and 'is_parent_primary' in table.columns:
        n_before = len(table)
        table = table[table['is_parent_primary'].astype(bool)]
        table = table.reset_index(drop=True)
        if verbose >= 1 and n_before != len(table):
            print(f' ++==>> cutout_on="parent": {n_before} component(s) '
                  f'collapsed to {len(table)} island(s).')

    if cutout_on == 'parent' and 'ra_parent' in table.columns:
        table['ra_cut'] = table['ra_parent']
        table['dec_cut'] = table['dec_parent']
    else:
        table['ra_cut'] = table['ra']
        table['dec_cut'] = table['dec']

    if 'parent_size_px' not in table.columns:
        table['parent_size_px'] = np.nan

    if name_column in table.columns:
        names = table[name_column].astype(str).values
    else:
        sky = SkyCoord(ra=table['ra_cut'].values * u.degree,
                       dec=table['dec_cut'].values * u.degree, frame='icrs')
        names = _iau_names(sky)
    # keep filenames safe
    table['cutout_name'] = [str(n).replace(' ', '_').replace('/', '_')
                            for n in names]

    return table
