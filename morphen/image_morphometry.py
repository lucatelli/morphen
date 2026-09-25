def deprecated(old_name, new_name):
    def decorator(func):
        def wrapper(*args, **kwargs):
            warnings.warn(f"'{old_name}' is deprecated and "
                          f"will be removed in a future version. "
                          f"Use '{new_name}' instead.",
                          category=DeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)
        return wrapper
    return decorator

def circle_area(r):
    """
    Return the pixel area of a circle given the radius.

    Parameters
    ----------
    r : float, ndarray
        Radius of the circle.

    Returns
    -------
    float, ndarray
        Area of the circle.
    """
    return (np.pi * r * r)


def area_to_radii(A):
    return (np.sqrt(A / np.pi))

def radii_to_area(radii,radii_err=None):
    area_from_radii = np.pi*radii*radii
    if radii_err is not None:
        area_from_radii_err = 2.0*np.pi*radii*radii_err
        # return(area_from_radii,area_from_radii_err)
    else:
        area_from_radii_err = 0.0
    return(area_from_radii,area_from_radii_err)

def pix_area_to_kpc_area(area_pix,pix_to_pc):
    kpc_area = (area_pix * (pix_to_pc**2.0)/(1000**2.0))
    return(kpc_area)
    

def cvi(imname):
    try:
        os.system('casaviewer ' + imname)
    except:
        try:
            os.system('~/casaviewer ' + imname)
        except:
            pass


def normalise_in_log(profile):
    def normalise(x):
        return (x- x.min())/(x.max() - x.min())
    y = profile.copy()
    nu_0 = normalise(np.log(y))
    return nu_0


def shuffle_2D(image):
    height, width = image.shape

    # Reshape the image to a 1D array
    image_flat = image.copy().reshape(-1)

    # Shuffle the pixels randomly
    np.random.shuffle(image_flat)

    # Reshape the shuffled 1D array back to a 2D image
    shuffled_image = image_flat.reshape(height, width)

    # Print the shuffled image
    return(shuffled_image)


def estimate_circular_aperture(image, cellsize, std=3):
    if isinstance(image, str) == True:
        g = load_fits_data(image)
    else:
        g = image
    gmask = g[g > std * mad_std(g)]
    npix = len(gmask)
    barea = beam_area2(image, cellsize)
    nbeams = npix / barea
    circ_radii = np.sqrt(npix / np.pi)
    return (circ_radii)





def sort_masks(image, mask_array, sort_by='flux'):
    unique_labels = np.unique(mask_array)
    region_labels = unique_labels[1:]  # 0 is the zero values
    mask_areas = []
    masks = []
    if sort_by == 'area':
        for i in range(len(region_labels)):
            mask_comp = (mask_array == region_labels[i])
            masks.append(mask_comp)
            area_mask = np.sum(mask_comp)
            mask_areas.append(area_mask)
    if sort_by == 'flux':
        for i in range(len(region_labels)):
            mask_comp = (mask_array == region_labels[i])
            masks.append(mask_comp)
            area_mask = np.sum(mask_comp * image)
            mask_areas.append(area_mask)

    mask_areas = np.asarray(mask_areas)
    sorted_indices_desc = np.argsort(mask_areas)[::-1]
    sorted_arr_desc = mask_areas[sorted_indices_desc]
    return (masks, sorted_indices_desc)


def rot180(imagem, x0, y0):
    R = np.matrix([[-1, 0], [0, -1]])
    img180 = nd.affine_transform(imagem, R,
                                 offset=(2. * (y0 - 0.5), 2. * (x0 - 0.5)))
    return img180


def estimate_area(data_mask, cellsize, Omaj, Omin):
    npix = np.sum(data_mask)
    barea = beam_area(Omaj, Omin, cellsize)
    nbeams = npix / barea
    circ_radii = np.sqrt(npix / np.pi)
    # print(npix)
    return (nbeams, circ_radii, npix)


def estimate_area_nbeam(image, mask, cellsize):
    npix = np.sum(mask)
    barea = beam_area2(image, cellsize)
    nbeams = npix / barea
    circ_radii = np.sqrt(npix / np.pi)
    # print(npix)
    return (circ_radii)

def estimate_area2(image,data_mask,cellsize):
    npix = np.sum(data_mask)
    barea= beam_area2(image,cellsize)
    nbeams = npix/barea
    circ_radii = np.sqrt(npix/np.pi)
    # print(npix)
    return(nbeams,circ_radii,npix)

def get_dilation_size(image):
    omaj, omin, _, _, _ = beam_shape(image)
    dilation_size = int(
        np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
    return(dilation_size)

def mask_dilation(image, cell_size=None, sigma=6, rms=None,
                  dilation_size=None, iterations=2, dilation_type='disk',
                  do_filtering=False,
                  PLOT=False, show_figure=False, logger=None,
                  fig_size_factor=1,
                  special_name='', verbose=0,
                  save_fig_name=None,
                  # --- new optional parameters (default off = identical to old behaviour) ---
                  residual=None,
                  sigma_residual=3.0,
                  min_size_beams=0.0,
                  ):
    """
    Mask dilation function.

    Apply a binary dilation to a mask originated from the emission of a source.
    This was originally designed for radio images, where the dilation factor is
    proportional to the beam size. The dilation factor is computed as the geometric
    mean of the major and minor axis of the beam. The dilation factor is then
    converted to pixels using the cell size of the image. Note that the expansion
    occurs in two opposite directions, so the dilation size should be half of the
    size of the beam for one iteration, so that the final dilation will be a unit of
    the beam size. However, the exact size of the dilation will also be determinet by
    the number of iterations.

    Parameters
    ----------
    image : str or np.ndarray
        Image name (FITS path) or 2-D numpy array.
    cell_size : float
        Cell size of the image in arcsec.
    sigma : float
        Sigma level for the mask.
    rms : float
        RMS level for the mask.
    dilation_size : int
        Size of the dilation structuring element in pixels.
    iterations : int
        Number of dilation iterations.
    dilation_type : str
        Type of dilation structuring element: 'disk' or 'square'.
    PLOT : bool
        Plot the mask dilation diagnostics.
    show_figure : bool
        Display the figure interactively.
    logger : logger
        Logger object for verbose messages.
    fig_size_factor : float
        Multiplicative factor for figure size.
    special_name : str
        Optional label appended to the dilated-mask panel title.
    verbose : int
        Verbosity level (0 = silent).
    save_fig_name : str or None
        File path to save the diagnostic figure.
    residual : str or np.ndarray or None
        Residual image corresponding to `image` (FITS path or 2-D array).
        When provided, pixels where ``|residual| > sigma_residual * rms_res``
        are excluded from the mask before dilation.  This rejects sidelobes
        and artefacts that were not properly subtracted by the deconvolver.
        Recommended for self-calibration pipelines.  Default None (disabled).
    sigma_residual : float
        Residual rejection threshold in units of the residual rms
        (``mad_std`` of the residual image).  Default 3.0.  Increase to be
        more inclusive of emission with high residuals (e.g. early self-cal
        steps where cleaning is shallow).
    min_size_beams : float
        Minimum connected-component size in units of beam areas to retain in
        the mask.  Default 0.0 (DISABLED).

        **WARNING - use with extreme caution.**  Setting this > 0 will reject
        components smaller than ``min_size_beams * beam_area_px``.  This WILL
        destroy real compact sources, isolated point sources, secondary nuclei
        (e.g. both nuclei in Arp 220), and faint companions that happen to be
        small in the mask.  Only enable if you are certain the field contains
        a single extended source with no compact companions, and use a very
        small value (< 0.1).  The primary improvement should come from
        ``residual`` gating, not size filtering.
    """

    if isinstance(image, str):
        data = load_fits_data(image)
    else:
        data = image
    if rms is None:
        std = mad_std(data[data > 0])
    else:
        std = rms

    if (dilation_size is None) or (dilation_size == '2x'):
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            if dilation_size == '2x':
                dilation_size = int(
                    2 * np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            else:
                dilation_size = int(
                    np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            if verbose >= 1:
                if logger is not None:
                    logger.debug(f" ==>  Mask dilation size is "
                                 f"{dilation_size} [px]")
                else:
                    print(f" ==>  Mask dilation size is "
                          f"{dilation_size} [px]")
        except:
            if dilation_size is None:
                dilation_size = 5

    if do_filtering:
        sigma_smooth = 3.0
        data_smooth = gaussian_filter(data, sigma=sigma_smooth)
        mask = (data_smooth >= sigma * std)
        mask3 = (data_smooth >= 3 * std)
    else:
        mask = (data >= sigma * std)
        mask3 = (data >= 3 * std)

    data_mask = mask * data

    # ------------------------------------------------------------------
    # Residual gating (new, only active when residual is provided)
    # ------------------------------------------------------------------
    res_data = None
    std_res = None
    mask_for_dilation = mask   # default: use raw sigma mask unchanged

    if residual is not None:
        if isinstance(residual, str):
            res_data = load_fits_data(residual)
        else:
            res_data = np.asarray(residual)
        std_res = mad_std(res_data)
        res_gate = (np.abs(res_data) <= sigma_residual * std_res)
        mask_for_dilation = mask & res_gate
        if verbose >= 1:
            n_rej = int(mask.sum()) - int(mask_for_dilation.sum())
            msg = (f" ==>  Residual gating (sigma_res={sigma_residual}): "
                   f"{mask.sum()} -> {mask_for_dilation.sum()} px "
                   f"({n_rej} rejected)")
            if logger is not None:
                logger.debug(msg)
            else:
                print(msg)

    # ------------------------------------------------------------------
    # Connected-component size filter (new, disabled by default)
    # ------------------------------------------------------------------
    if min_size_beams > 0.0:
        beam_area_px = None
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            cell = get_cell_size(image)
            bmaj_px = omaj / cell
            bmin_px = omin / cell
            beam_area_px = np.pi * bmaj_px * bmin_px / (4 * np.log(2))
        except Exception:
            pass
        if beam_area_px is not None and beam_area_px > 0:
            min_px = max(1, int(min_size_beams * beam_area_px))
            labeled, n_features = ndimage.label(mask_for_dilation)
            filtered = np.zeros_like(mask_for_dilation)
            for comp_id in range(1, n_features + 1):
                comp = (labeled == comp_id)
                if comp.sum() >= min_px:
                    filtered |= comp
            mask_for_dilation = filtered
            if verbose >= 1:
                msg = (f" ==>  Size filter (min={min_size_beams:.3f} beams, "
                       f"{min_px} px): -> {mask_for_dilation.sum()} px")
                if logger is not None:
                    logger.debug(msg)
                else:
                    print(msg)

    # ------------------------------------------------------------------
    # Dilation
    # ------------------------------------------------------------------
    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(
            mask_for_dilation,
            structure=disk(dilation_size),
            iterations=iterations).astype(mask.dtype)
    else:
        data_mask_d = ndimage.binary_dilation(
            mask_for_dilation,
            structure=square(dilation_size),
            iterations=iterations).astype(mask.dtype)

    # ------------------------------------------------------------------
    # Diagnostic plot
    # ------------------------------------------------------------------
    if PLOT:
        _use_residual_panels = (residual is not None)
        n_panels = 7 if _use_residual_panels else 4
        fig = plt.figure(figsize=(int(4 * n_panels * fig_size_factor),
                                  int(4 * fig_size_factor)))

        ax0 = fig.add_subplot(1, n_panels, 1)
        ax0.imshow(mask3, origin='lower', cmap='magma')
        ax0.set_title(r'mask $>$ ' + str(3) + r'$\times\sigma_{\mathrm{rms}}$')
        ax0.axis('off')

        ax1 = fig.add_subplot(1, n_panels, 2)
        ax1.imshow(mask, origin='lower', cmap='magma')
        ax1.set_title(r'mask $>$ ' + str(sigma) + r'$\times\sigma_{\mathrm{rms}}$')
        ax1.axis('off')

        if _use_residual_panels:
            from matplotlib.colors import SymLogNorm
            _linthresh = max(float(std_res), 1e-10)
            _norm_res = SymLogNorm(linthresh=_linthresh, linscale=0.5,
                                   vmin=-5 * _linthresh * 100,
                                   vmax=5 * _linthresh * 100, base=10)
            ax_res = fig.add_subplot(1, n_panels, 3)
            ax_res.imshow(res_data, origin='lower', cmap='RdBu_r',
                          norm=_norm_res)
            ax_res.set_title(f'residual\n(sigma_res={sigma_residual})')
            ax_res.axis('off')

            ax_gated = fig.add_subplot(1, n_panels, 4)
            ax_gated.imshow(mask_for_dilation, origin='lower', cmap='magma')
            ax_gated.set_title(f'gated mask\n({mask_for_dilation.sum()} px)')
            ax_gated.axis('off')

            ax_diff = fig.add_subplot(1, n_panels, 5)
            diff = mask.astype(int) - mask_for_dilation.astype(int)
            ax_diff.imshow(diff, origin='lower', cmap='bwr', vmin=-1, vmax=1)
            ax_diff.set_title(f'rejected pixels\n({(diff==1).sum()} px)')
            ax_diff.axis('off')

            ax2 = fig.add_subplot(1, n_panels, 6)
            ax2.imshow(data_mask_d, origin='lower', cmap='magma')
            ax2.set_title(r'D(mask)' f'[{dilation_size}$|${iterations}]'
                          f'{special_name}')
            ax2.axis('off')

            ax3 = fig.add_subplot(1, n_panels, 7)
        else:
            ax2 = fig.add_subplot(1, n_panels, 3)
            ax2.imshow(data_mask_d, origin='lower', cmap='magma')
            ax2.set_title(r'D(mask)' f'[{dilation_size}$|${iterations}]'
                          f'{special_name}')
            ax2.axis('off')

            ax3 = fig.add_subplot(1, n_panels, 4)

        ax3 = eimshow(data * data_mask_d, ax=ax3, CM='magma',
                      rms=std, add_contours=True,
                      vmin_factor=0.01, vmax_factor=0.1)
        ax3.set_title(r'D(mask) $\times$ data')
        ax3.axis('off')

        fig.subplots_adjust(wspace=0.01, hspace=0.01)
        if show_figure:
            if save_fig_name is not None:
                fig.savefig(save_fig_name, dpi=300, bbox_inches='tight')
            plt.show()
        plt.close(fig)
        plt.close('all')

    return (mask, data_mask_d)


def _region_grow(seed_mask, candidate_mask):
    """Return candidate pixels reachable from seeds via 8-connected paths."""
    labeled, _ = ndimage.label(candidate_mask,
                               structure=np.ones((3, 3), dtype=np.int32))
    seed_labels = set(int(x) for x in labeled[seed_mask & (labeled > 0)].ravel())
    seed_labels.discard(0)
    if not seed_labels:
        return np.zeros_like(candidate_mask, dtype=np.int32)
    return np.isin(labeled, list(seed_labels)).astype(np.int32)


def mask_dilation_snr(image, residual,
                      cell_size=None, sigma=5.0, sigma_seed=10.0,
                      sigma_residual=3.0, rms=None,
                      dilation_size=None, iterations=2, dilation_type='disk',
                      PLOT=False, show_figure=False, logger=None,
                      fig_size_factor=1, special_name='', verbose=0,
                      save_fig_name=None):
    """
    SNR-guided seed-and-grow mask dilation for self-calibration.

    Builds a clean mask by growing outward from high-confidence seeds and
    accepting only emission that is (a) connected to those seeds via a
    continuous path of above-threshold pixels, and (b) passes a residual
    quality gate.  Isolated noise peaks and partially-cleaned sidelobes
    cannot seed the growth and are therefore excluded.

    Algorithm
    ---------
    1. Seed mask: pixels with ``image >= sigma_seed * rms`` - these are
       the high-confidence anchors.
    2. Candidate mask: pixels with ``image >= sigma * rms`` AND
       ``|residual| <= sigma_residual * rms_res`` - lower threshold but
       residual-gated.
    3. Region grow: label all connected components in the candidate mask;
       retain only those components that contain at least one seed pixel.
    4. Dilate the grown mask using a disk or square structuring element.

    Parameters
    ----------
    image : str or np.ndarray
        Cleaned image (FITS path or 2-D numpy array).
    residual : str or np.ndarray
        Residual image corresponding to `image` (FITS path or 2-D array).
        The residual rms is estimated as ``mad_std(residual_data)``.
    cell_size : float, optional
        Cell size in arcsec (used for dilation_size auto-detection from header).
    sigma : float
        Lower threshold for the grow mask in units of image rms.  Default 5.0.
        Pixels at or above this level that are connected to seeds are included.
    sigma_seed : float
        High-confidence seed threshold in units of image rms.  Default 10.0.
        Pixels below this level can still be included if they are connected to
        a seed by continuous paths above `sigma`.  Raise this value if spurious
        seeds appear (e.g., strong artefacts above sigma_seed - use the
        residual gate to suppress them first).
    sigma_residual : float
        Residual rejection threshold applied to the *grow* mask in units of
        ``mad_std(residual)``.  Default 3.0.  Pixels where
        ``|residual| > sigma_residual * rms_res`` are excluded from the
        candidate mask.
    rms : float, optional
        Image rms.  If None, computed as ``mad_std(residual_data)``, which is
        the recommended estimator for self-calibration images.
    dilation_size : int or None
        Structuring element half-size in pixels.  Auto-computed from beam if None.
    iterations : int
        Number of dilation iterations.  Default 2.
    dilation_type : str
        'disk' (default) or 'square'.
    PLOT : bool
        Produce a 7-panel diagnostic figure.
    show_figure : bool
        Display the figure interactively (requires PLOT=True).
    logger : logging.Logger, optional
        Logger for verbose messages.
    fig_size_factor : float
        Multiplicative scale for the diagnostic figure size.
    special_name : str
        Optional suffix for the dilated-mask panel title.
    verbose : int
        Verbosity level (0 = silent, 1 = progress messages).
    save_fig_name : str, optional
        File path to save the diagnostic figure.

    Returns
    -------
    seed_mask : np.ndarray (bool/int)
        The raw high-confidence seed mask (image >= sigma_seed * rms).
        Analogous to the first return of ``mask_dilation``.
    data_mask_d : np.ndarray (int)
        The dilated grown mask.  Pass to WSClean as the cleaning mask.
        Analogous to the second return of ``mask_dilation``.

    Notes
    -----
    Recommended parameters for VLA/e-MERLIN self-calibration:
      - ap1 step (amplitude+phase, most dangerous): sigma=5, sigma_seed=10,
        sigma_residual=3.0
      - For extended sources with shallow cleaning: raise sigma_residual to 4-5
      - For compact sources at high SNR: defaults work well

    This function is safe to call when the residual image does not yet exist
    (it will raise FileNotFoundError for str paths; callers should check
    ``os.path.exists`` before calling).
    """
    if isinstance(image, str):
        data = load_fits_data(image)
    else:
        data = np.asarray(image)

    if isinstance(residual, str):
        res_data = load_fits_data(residual)
    else:
        res_data = np.asarray(residual)

    rms_res = float(mad_std(res_data))
    std = rms if rms is not None else rms_res

    if (dilation_size is None) or (dilation_size == '2x'):
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            if dilation_size == '2x':
                dilation_size = int(
                    2 * np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            else:
                dilation_size = int(
                    np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            if verbose >= 1:
                msg = f" ==>  Mask dilation size is {dilation_size} [px]"
                if logger is not None:
                    logger.debug(msg)
                else:
                    print(msg)
        except Exception:
            if dilation_size is None:
                dilation_size = 5

    # ---- seed-and-grow ------------------------------------------------
    seed_mask = (data >= sigma_seed * std)
    candidate_mask = ((data >= sigma * std) &
                      (np.abs(res_data) <= sigma_residual * rms_res))
    grown = _region_grow(seed_mask, candidate_mask)

    if verbose >= 1:
        labeled_g, n_grown = ndimage.label(grown, structure=np.ones((3, 3)))
        msg = (f" ==>  Seed-and-grow: seeds={seed_mask.sum()} px, "
               f"candidates={candidate_mask.sum()} px, "
               f"grown={grown.sum()} px ({n_grown} components)")
        if logger is not None:
            logger.debug(msg)
        else:
            print(msg)

    # ---- dilation -----------------------------------------------------
    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(
            grown, structure=disk(dilation_size),
            iterations=iterations).astype(grown.dtype)
    else:
        data_mask_d = ndimage.binary_dilation(
            grown, structure=square(dilation_size),
            iterations=iterations).astype(grown.dtype)

    # ---- diagnostic plot ----------------------------------------------
    if PLOT:
        from matplotlib.colors import SymLogNorm as _SymLogNorm
        _linthresh_img = max(float(std), 1e-12)
        _norm_img = _SymLogNorm(linthresh=_linthresh_img, linscale=0.5,
                                vmin=-10 * _linthresh_img,
                                vmax=float(np.nanmax(np.abs(data))), base=10)
        _linthresh_res = max(rms_res, 1e-12)
        _norm_res = _SymLogNorm(linthresh=_linthresh_res, linscale=0.5,
                                vmin=-5 * _linthresh_res * 50,
                                vmax=5 * _linthresh_res * 50, base=10)

        fig = plt.figure(figsize=(int(28 * fig_size_factor),
                                  int(4 * fig_size_factor)))
        axes = [fig.add_subplot(1, 7, i + 1) for i in range(7)]

        axes[0].imshow(data, origin='lower', norm=_norm_img, cmap='inferno')
        axes[0].set_title('image (SymLog)')
        axes[0].axis('off')

        axes[1].imshow(res_data, origin='lower', norm=_norm_res, cmap='RdBu_r')
        axes[1].set_title(f'residual (sigma_res={sigma_residual})')
        axes[1].axis('off')

        axes[2].imshow(seed_mask, origin='lower', cmap='magma')
        axes[2].set_title(f'seed mask\n(sigma_seed={sigma_seed})')
        axes[2].axis('off')

        axes[3].imshow(candidate_mask, origin='lower', cmap='magma')
        axes[3].set_title(f'candidate mask\n(sigma={sigma}, sigma_res={sigma_residual})')
        axes[3].axis('off')

        axes[4].imshow(grown, origin='lower', cmap='magma')
        axes[4].set_title(f'grown\n({grown.sum()} px)')
        axes[4].axis('off')

        axes[5].imshow(data_mask_d, origin='lower', cmap='magma')
        axes[5].set_title(r'D(grown)' f'[{dilation_size}|{iterations}]'
                          f'{special_name}')
        axes[5].axis('off')

        axes[6] = eimshow(data * data_mask_d, ax=axes[6], CM='magma',
                          rms=std, add_contours=True,
                          vmin_factor=0.01, vmax_factor=0.1)
        axes[6].set_title(r'D(grown) $\times$ data')
        axes[6].axis('off')

        fig.subplots_adjust(wspace=0.01, hspace=0.01)
        if show_figure:
            if save_fig_name is not None:
                fig.savefig(save_fig_name, dpi=300, bbox_inches='tight')
            plt.show()
        plt.close(fig)
        plt.close('all')

    return (seed_mask, data_mask_d)


def mask_dilation_v2(image, cell_size=None, sigma=6,rms=None,
                  dilation_size=None,iterations=2, dilation_type='disk',
                  PLOT=False,show_figure=False,logger=None,
                  fig_size_factor = 1,
                  special_name='',verbose=0):
    """
    Mask dilation function.

    Apply a binary dilation to a mask originated from the emission of a source.
    This was originally designed for radio images, where the dilation factor is
    proportional to the beam size. The dilation factor is computed as the geometric
    mean of the major and minor axis of the beam. The dilation factor is then
    converted to pixels using the cell size of the image. Note that the expansion
    occurs in two opposite directions, so the dilation size should be half of the
    size of the beam for one iteration, so that the final dilation will be a unit of
    the beam size. However, the exact size of the dilation will also be determinet by
    the number of iterations.


    Parameters
    ----------
    image : str
        Image name.
    cell_size : float
        Cell size of the image in arcsec.
    sigma : float
        Sigma level for the mask.
    rms : float
        RMS level for the mask.
    dilation_size : int
        Size of the dilation.
    iterations : int
        Number of iterations for the dilation.
    dilation_type : str
        Type of dilation. Options are 'disk' or 'square'.
    PLOT : bool
        Plot the mask dilation.
    show_figure : bool
        Show the figure.
    logger : logger
        Logger object.
    fig_size_factor : float
        Figure size factor.
    special_name : str
        Special name for the mask dilation.
    """
    from photutils.segmentation import detect_sources
    from astropy.stats import sigma_clipped_stats
    
    if isinstance(image, str) == True:
        data = load_fits_data(image)
    else:
        data = image
    if rms is None:
        std = mad_std(data)
    else:
        std = rms

    if (dilation_size is None) or (dilation_size == '2x'):
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            if dilation_size == '2x':
                dilation_size = int(
                    2*np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            else:
                dilation_size = int(
                    np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            if verbose >= 1:
                if logger is not None:
                    logger.debug(f" ==>  Mask dilation size is "
                                f"{dilation_size} [px]")
                else:
                    print(f" ==>  Mask dilation size is "
                        f"{dilation_size} [px]")
        except:
            if dilation_size is None:
                dilation_size = 5
    try:
        beam_size_px = get_beam_size_px(image)[0]
    except:
        beam_size_px = dilation_size
    

    mask_o = (data >= sigma * std)
    mask3 = (data >= 3 * std)
    data_mask = mask_o * data
    
    mean, median, std = sigma_clipped_stats(data, sigma=sigma)
    # threshold = median - (sigma * std)
    threshold = sigma * mad_std(data,ignore_nan=True) 
    segm = detect_sources(data, threshold, npixels=int(2*beam_size_px))
    if segm is None:
        segm_mask = mask_o
        # return (mask,np.zeros(data.shape, dtype=bool))
        # return (mask_o,mask_o)
    else:
        segm_mask = segm.data
    

    
    mask = mask_o * segm_mask

    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=disk(dilation_size),
                                            iterations=iterations).astype(mask.dtype)

    if dilation_type == 'square':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=square(dilation_size),
                                            iterations=iterations).astype(mask.dtype)

    # # Plot the original image, segmentation map, and final mask
    # fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # # Original image with mask overlay
    # axes[0].imshow(data, origin='lower', cmap='gray')
    # axes[0].imshow(mask_o, origin='lower', cmap='Reds', alpha=0.5)
    # axes[0].set_title('Original Image with Mask')
    # axes[0].axis('off')

    # # Segmentation map
    # axes[1].imshow(segm_mask, origin='lower', cmap='viridis')
    # axes[1].set_title('Segmentation Map')
    # axes[1].axis('off')

    # # Final dilated mask
    # axes[2].imshow(data_mask_d, origin='lower', cmap='gray',alpha=0.5)
    # axes[2].imshow(mask, origin='lower', cmap='Reds', alpha=0.5)
    # axes[2].imshow(data*data_mask_d, origin='lower', cmap='gray',alpha=0.9)
    # axes[2].set_title('Final Dilated Mask')
    # axes[2].axis('off')

    # plt.tight_layout()
    # plt.show()

    if PLOT == True:
        fig = plt.figure(figsize=(int(16*fig_size_factor), int(4*fig_size_factor)))
        ax0 = fig.add_subplot(1, 4, 1)
        ax0.imshow((mask3), origin='lower',cmap='magma')
        ax0.set_title(r'mask $>$ ' + str(3) + r'$\times\sigma_{\mathrm{rms}}$')
        ax0.axis('off')
        ax1 = fig.add_subplot(1, 4, 2)
        #         ax1.legend(loc='lower left')
        ax1.imshow((mask_o), origin='lower',cmap='magma')
        ax1.set_title(r'mask $>$ ' + str(sigma) + r'$\times\sigma_{\mathrm{rms}}$')
        ax1.axis('off')
        ax2 = fig.add_subplot(1, 4, 3)
        ax2.imshow(data_mask_d, origin='lower',cmap='magma')
        ax2.set_title(r'D(mask)'f'[{dilation_size}|{iterations}]'f'{special_name}')
        ax2.axis('off')
        ax3 = fig.add_subplot(1, 4, 4)
        ax3 = eimshow(data * data_mask_d, ax=ax3, CM='magma',
                      rms=std,add_contours=True,
                      vmin_factor=0.01,vmax_factor=0.1)
        ax3.set_title(r'D(mask) $\times$ data')
        #         ax3.imshow(np.log(data*data_mask_d))
        fig.subplots_adjust(wspace=0.01, hspace=0.01)
        #         fig.tight_layout()
        ax3.axis('off')
        if show_figure == True:
            plt.show()
        else:
            plt.close()
    #         plt.savefig(image.replace('.fits','_masks.jpg'),dpi=300, bbox_inches='tight')

    # if cell_size is not None:
    #     if isinstance(image, str) == True:
    #         try:
    #             print((data * data_mask_d).sum() / beam_area2(image, cell_size))
    #             print((data * data_mask).sum() / beam_area2(image, cell_size))
    #             print((data).sum() / beam_area2(image, cell_size))
    #         except:
    #             print('Provide a cell size of the image.')
    return (mask.astype('int32'), data_mask_d.astype('int32'))

def t_mask_dilation(image, cell_size=None, sigma=6, rms=None,
                  dilation_size=None, iterations=2, dilation_type='disk',
                  do_filtering=False,
                  PLOT=False, show_figure=False, logger=None,
                  fig_size_factor=1,
                  special_name='', verbose=0,
                  max_distance_fraction=0.25, 
                  use_distance_filter=True,
                  min_pix=None):
    """
    testing-Mask dilation function with connected component filtering.
    Geferson Lucatelli, 2022-2024 (https://doi.org/10.1093/mnras/stae744)
    Improved 2025.
    
    This function creates a mask from astronomical image data, filters connected 
    components based on their proximity to a reference location and their size, 
    and applies binary dilation. The reference location is determined by either 
    the largest component or the component closest to the image center.
    
    Parameters
    ----------
    image : str or ndarray
        Image filename or data array.
    cell_size : float, optional
        Cell size of the image in arcsec (currently unused).
    sigma : float, optional
        Sigma threshold for initial mask creation (default: 6).
    rms : float, optional
        RMS noise level. If None, computed using MAD statistics.
    dilation_size : int or str, optional
        Size of dilation structuring element. If None, computed from beam size.
        If '2x', uses twice the beam-derived size.
    iterations : int, optional
        Number of dilation iterations (default: 2).
    dilation_type : str, optional
        Type of structuring element: 'disk' or 'square' (default: 'disk').
    PLOT : bool, optional
        Whether to create diagnostic plots (default: False).
    show_figure : bool, optional
        Whether to display plots or close them (default: False).
    logger : logger object, optional
        Logger for diagnostic messages.
    fig_size_factor : float, optional
        Scaling factor for figure size (default: 1).
    special_name : str, optional
        Additional text for plot titles (default: '').
    verbose : int, optional
        Verbosity level (default: 0).
    max_distance_fraction : float, optional
        Maximum fractional distance from reference point to retain components.
        Distance is measured relative to the image diagonal (default: 0.25).
    use_distance_filter : bool, optional
        If True, use distance-based filtering. If False, keep only the largest
        component (legacy behavior) (default: True).
    min_pix : int, optional
        Minimum pixel area for a component to be retained. If None (default),
        automatically set to half the area of the reference component. This
        filters out small artifacts while retaining larger separated structures
        such as merging systems or companions.
    
    Returns
    -------
    mask : ndarray
        Boolean mask of filtered components before dilation.
    data_mask_d : ndarray
        Boolean mask after binary dilation.
    
    Notes
    -----
    The distance and size-based filtering approach:
    1. Identifies a reference component (largest or closest to image center)
    2. Computes its centroid as the reference point
    3. Retains all components that satisfy BOTH criteria:
       - Centroids within max_distance_fraction of image diagonal from reference
       - Area >= min_pix (or >= 0.5 * reference component area if min_pix=None)
    
    This allows capturing extended structures, multiple components of a 
    complex source, or companions while filtering out spurious detections
    and small artifacts.
    """
    from skimage.measure import label, regionprops
    from photutils.segmentation import make_2dgaussian_kernel

    if isinstance(image, str) == True:
        data = load_fits_data(image)
    else:
        data = image
    if rms is None:
        std = mad_std(data[data > 0])
    else:
        std = rms

    if (dilation_size is None) or (dilation_size == '2x'):
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            if dilation_size == '2x':
                dilation_size = int(
                    2*np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            else:
                dilation_size = int(
                    np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            if verbose >= 1:
                if logger is not None:
                    logger.debug(f" ==>  Mask dilation size is {dilation_size} [px]")
                else:
                    print(f" ==>  Mask dilation size is {dilation_size} [px]")
        except:
            if dilation_size is None:
                dilation_size = 5
    
    if do_filtering:
        # Apply Gaussian smoothing to the data before thresholding
        sigma_smooth = 3.0  # Standard deviation for Gaussian kernel
        data_smooth = gaussian_filter(data, sigma=sigma_smooth)
        # Create initial masks
        mask = (data_smooth >= sigma * std)
        mask3 = (data_smooth >= 3 * std)
    else:
        mask = (data >= sigma * std)
        mask3 = (data >= 3 * std)

    # Label connected components in the mask
    labeled_mask = label(mask)
    regions = regionprops(labeled_mask)
    
    if use_distance_filter and regions:
        # Distance-based filtering approach
        ny, nx = data.shape
        image_center = np.array([nx / 2.0, ny / 2.0])
        image_diagonal = np.sqrt(nx**2 + ny**2)
        max_distance = max_distance_fraction * image_diagonal
        
        # Find reference component (largest or closest to center)
        # First try to find the component closest to center
        min_dist_to_center = np.inf
        reference_region = None
        
        for region in regions:
            centroid = np.array([region.centroid[1], region.centroid[0]])  # x, y order
            dist_to_center = np.linalg.norm(centroid - image_center)
            if dist_to_center < min_dist_to_center:
                min_dist_to_center = dist_to_center
                reference_region = region
        
        # If no component near center, use the largest component
        if reference_region is None or min_dist_to_center > 0.3 * image_diagonal:
            largest_area = 0
            for region in regions:
                if region.area > largest_area:
                    largest_area = region.area
                    reference_region = region
        
        # Get reference centroid and determine minimum size threshold
        if reference_region is not None:
            reference_centroid = np.array([reference_region.centroid[1], 
                                          reference_region.centroid[0]])
            
            # Set minimum pixel threshold
            if min_pix is None:
                min_pixel_threshold = 0.5 * reference_region.area
            else:
                min_pixel_threshold = min_pix
            
            if verbose >= 1:
                msg = f" ==>  Reference centroid at ({reference_centroid[0]:.1f}, {reference_centroid[1]:.1f})"
                if logger is not None:
                    logger.debug(msg)
                else:
                    print(msg)
                msg = f" ==>  Reference area: {reference_region.area} px, minimum threshold: {min_pixel_threshold:.0f} px"
                if logger is not None:
                    logger.debug(msg)
                else:
                    print(msg)
            
            # Filter components by distance from reference AND size
            mask = np.zeros_like(data, dtype=bool)
            n_kept = 0
            n_rejected_distance = 0
            n_rejected_size = 0
            
            for region in regions:
                centroid = np.array([region.centroid[1], region.centroid[0]])
                distance = np.linalg.norm(centroid - reference_centroid)
                
                # Check both distance and size criteria
                if distance <= max_distance:
                    if region.area >= min_pixel_threshold:
                        mask = mask | (labeled_mask == region.label)
                        n_kept += 1
                        if verbose >= 2:
                            msg = f"     Keeping component {region.label}: distance={distance:.1f} px, area={region.area} px"
                            if logger is not None:
                                logger.debug(msg)
                            else:
                                print(msg)
                    else:
                        n_rejected_size += 1
                        if verbose >= 2:
                            msg = f"     Rejecting component {region.label} (too small): distance={distance:.1f} px, area={region.area} px"
                            if logger is not None:
                                logger.debug(msg)
                            else:
                                print(msg)
                else:
                    n_rejected_distance += 1
                    if verbose >= 2:
                        msg = f"     Rejecting component {region.label} (too far): distance={distance:.1f} px, area={region.area} px"
                        if logger is not None:
                            logger.debug(msg)
                        else:
                            print(msg)
            
            if verbose >= 1:
                msg = f" ==>  Kept {n_kept} of {len(regions)} components (rejected: {n_rejected_distance} by distance, {n_rejected_size} by size)"
                if logger is not None:
                    logger.debug(msg)
                else:
                    print(msg)
        else:
            # Fallback: no regions found
            mask = np.zeros_like(data, dtype=bool)
    
    elif regions:
        # Legacy behavior: keep only largest component
        largest_area = 0
        largest_label = None
        for region in regions:
            if region.area > largest_area:
                largest_area = region.area
                largest_label = region.label
        # Create a new mask with only the largest component
        if largest_label is not None:
            mask = labeled_mask == largest_label
    
    data_mask = mask * data

    # Perform dilation
    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=disk(dilation_size),
                                            iterations=iterations).astype(mask.dtype)
    if dilation_type == 'square':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=square(dilation_size),
                                            iterations=iterations).astype(mask.dtype)

    # Apply the same filtering approach to the 3-sigma mask for plotting
    labeled_mask3 = label(mask3)
    regions3 = regionprops(labeled_mask3)
    
    if use_distance_filter and regions3:
        # Use same reference point, distance criterion, and size criterion
        if 'reference_centroid' in locals() and reference_region is not None:
            mask3 = np.zeros_like(data, dtype=bool)
            for region in regions3:
                centroid = np.array([region.centroid[1], region.centroid[0]])
                distance = np.linalg.norm(centroid - reference_centroid)
                # Apply both distance and size filters
                if distance <= max_distance and region.area >= min_pixel_threshold:
                    mask3 = mask3 | (labeled_mask3 == region.label)
        else:
            # Fallback to largest component
            if regions3:
                largest_area = 0
                largest_label = None
                for region in regions3:
                    if region.area > largest_area:
                        largest_area = region.area
                        largest_label = region.label
                if largest_label is not None:
                    mask3 = labeled_mask3 == largest_label
    elif regions3:
        # Legacy: largest component only
        largest_area = 0
        largest_label = None
        for region in regions3:
            if region.area > largest_area:
                largest_area = region.area
                largest_label = region.label
        if largest_label is not None:
            mask3 = labeled_mask3 == largest_label

    if PLOT == True:
        fig = plt.figure(figsize=(int(16*fig_size_factor), int(4*fig_size_factor)))
        ax0 = fig.add_subplot(1, 4, 1)
        ax0.imshow((mask3), origin='lower',cmap='magma')
        ax0.set_title(r'mask $>$ ' + str(3) + r'$\times\sigma_{\mathrm{rms}}$')
        ax0.axis('off')
        ax1 = fig.add_subplot(1, 4, 2)
        ax1.imshow((mask), origin='lower',cmap='magma')
        ax1.set_title(r'mask $>$ ' + str(sigma) + r'$\times\sigma_{\mathrm{rms}}$')
        ax1.axis('off')
        ax2 = fig.add_subplot(1, 4, 3)
        ax2.imshow(data_mask_d, origin='lower',cmap='magma')
        ax2.set_title(r'D(mask)'f'[{dilation_size}$|${iterations}]'f'{special_name}')
        ax2.axis('off')
        ax3 = fig.add_subplot(1, 4, 4)
        ax3 = eimshow(data * data_mask_d, ax=ax3, CM='magma',
                      rms=std,add_contours=True,
                      vmin_factor=0.01,vmax_factor=0.1)
        ax3.set_title(r'D(mask) $\times$ data')
        fig.subplots_adjust(wspace=0.01, hspace=0.01)
        ax3.axis('off')
        if show_figure == True:
            plt.show()
        plt.close(fig)
        plt.close('all')
        # else:
        #     plt.close()

    return (mask, data_mask_d)

# def t_mask_dilation(image, cell_size=None, sigma=6, rms=None,
#                   dilation_size=None, iterations=2, dilation_type='disk',
#                   do_filtering=True,
#                   PLOT=False, show_figure=False, logger=None,
#                   fig_size_factor=1,
#                   special_name='', verbose=0,
#                   max_distance_fraction=0.25, 
#                   use_distance_filter=True):
#     """
#     Mask dilation function with connected component filtering.
    
#     This function creates a mask from astronomical image data, filters connected 
#     components based on their proximity to a reference location, and applies 
#     binary dilation. The reference location is determined by either the largest 
#     component or the component closest to the image center.
    
#     Parameters
#     ----------
#     image : str or ndarray
#         Image filename or data array.
#     cell_size : float, optional
#         Cell size of the image in arcsec (currently unused).
#     sigma : float, optional
#         Sigma threshold for initial mask creation (default: 6).
#     rms : float, optional
#         RMS noise level. If None, computed using MAD statistics.
#     dilation_size : int or str, optional
#         Size of dilation structuring element. If None, computed from beam size.
#         If '2x', uses twice the beam-derived size.
#     iterations : int, optional
#         Number of dilation iterations (default: 2).
#     dilation_type : str, optional
#         Type of structuring element: 'disk' or 'square' (default: 'disk').
#     PLOT : bool, optional
#         Whether to create diagnostic plots (default: False).
#     show_figure : bool, optional
#         Whether to display plots or close them (default: False).
#     logger : logger object, optional
#         Logger for diagnostic messages.
#     fig_size_factor : float, optional
#         Scaling factor for figure size (default: 1).
#     special_name : str, optional
#         Additional text for plot titles (default: '').
#     verbose : int, optional
#         Verbosity level (default: 0).
#     max_distance_fraction : float, optional
#         Maximum fractional distance from reference point to retain components.
#         Distance is measured relative to the image diagonal (default: 0.5).
#     use_distance_filter : bool, optional
#         If True, use distance-based filtering. If False, keep only the largest
#         component (legacy behavior) (default: True).
    
#     Returns
#     -------
#     mask : ndarray
#         Boolean mask of filtered components before dilation.
#     data_mask_d : ndarray
#         Boolean mask after binary dilation.
    
#     Notes
#     -----
#     The distance-based filtering approach:
#     1. Identifies a reference component (largest or closest to image center)
#     2. Computes its centroid as the reference point
#     3. Retains all components with centroids within max_distance_fraction 
#        of the image diagonal from the reference point
    
#     This allows capturing extended structures, multiple components of a 
#     complex source, or companions while filtering out spurious detections.
#     """
#     from skimage.measure import label, regionprops
#     from photutils.segmentation import make_2dgaussian_kernel

#     if isinstance(image, str) == True:
#         data = load_fits_data(image)
#     else:
#         data = image
#     if rms is None:
#         std = mad_std(data[data > 0])
#     else:
#         std = rms

#     if (dilation_size is None) or (dilation_size == '2x'):
#         try:
#             omaj, omin, _, _, _ = beam_shape(image)
#             if dilation_size == '2x':
#                 dilation_size = int(
#                     2*np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
#             else:
#                 dilation_size = int(
#                     np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
#             if verbose >= 1:
#                 if logger is not None:
#                     logger.debug(f" ==>  Mask dilation size is {dilation_size} [px]")
#                 else:
#                     print(f" ==>  Mask dilation size is {dilation_size} [px]")
#         except:
#             if dilation_size is None:
#                 dilation_size = 5
    
#     if do_filtering:
#         # Apply Gaussian smoothing to the data before thresholding
#         sigma_smooth = 3.0  # Standard deviation for Gaussian kernel
#         data_smooth = gaussian_filter(data, sigma=sigma_smooth)
#         # kernel = make_2dgaussian_kernel(3.0, size=3)  # FWHM = 3.0

#         # # Determine padding size based on kernel
#         # pad_size = kernel.shape[0] // 2

#         # # Pad the data (using edge values to minimize artifacts)
#         # data_padded = np.pad(data, pad_size, mode='edge')

#         # # Convolve
#         # convolved_padded = scipy.signal.fftconvolve(data_padded, kernel, mode='same')

#         # # Crop back to original size
#         # data_smooth = convolved_padded[pad_size:-pad_size, pad_size:-pad_size]
#         # Create initial masks
#         mask = (data_smooth >= sigma * std)
#         mask3 = (data_smooth >= 3 * std)
#     else:
#         mask = (data >= sigma * std)
#         mask3 = (data >= 3 * std)



    
#     # Label connected components in the mask
#     labeled_mask = label(mask)
#     regions = regionprops(labeled_mask)
    
#     if use_distance_filter and regions:
#         # Distance-based filtering approach
#         ny, nx = data.shape
#         image_center = np.array([nx / 2.0, ny / 2.0])
#         image_diagonal = np.sqrt(nx**2 + ny**2)
#         max_distance = max_distance_fraction * image_diagonal
        
#         # Find reference component (largest or closest to center)
#         # First try to find the component closest to center
#         min_dist_to_center = np.inf
#         reference_region = None
        
#         for region in regions:
#             centroid = np.array([region.centroid[1], region.centroid[0]])  # x, y order
#             dist_to_center = np.linalg.norm(centroid - image_center)
#             if dist_to_center < min_dist_to_center:
#                 min_dist_to_center = dist_to_center
#                 reference_region = region
        
#         # If no component near center, use the largest component
#         if reference_region is None or min_dist_to_center > 0.3 * image_diagonal:
#             largest_area = 0
#             for region in regions:
#                 if region.area > largest_area:
#                     largest_area = region.area
#                     reference_region = region
        
#         # Get reference centroid
#         if reference_region is not None:
#             reference_centroid = np.array([reference_region.centroid[1], 
#                                           reference_region.centroid[0]])
            
#             if verbose >= 1:
#                 msg = f" ==>  Reference centroid at ({reference_centroid[0]:.1f}, {reference_centroid[1]:.1f})"
#                 if logger is not None:
#                     logger.debug(msg)
#                 else:
#                     print(msg)
            
#             # Filter components by distance from reference
#             mask = np.zeros_like(data, dtype=bool)
#             n_kept = 0
#             for region in regions:
#                 centroid = np.array([region.centroid[1], region.centroid[0]])
#                 distance = np.linalg.norm(centroid - reference_centroid)
                
#                 if distance <= max_distance:
#                     mask = mask | (labeled_mask == region.label)
#                     n_kept += 1
#                     if verbose >= 2:
#                         msg = f"     Keeping component {region.label} at distance {distance:.1f} px"
#                         if logger is not None:
#                             logger.debug(msg)
#                         else:
#                             print(msg)
            
#             if verbose >= 1:
#                 msg = f" ==>  Kept {n_kept} of {len(regions)} components within {max_distance:.1f} px"
#                 if logger is not None:
#                     logger.debug(msg)
#                 else:
#                     print(msg)
#         else:
#             # Fallback: no regions found
#             mask = np.zeros_like(data, dtype=bool)
    
#     elif regions:
#         # Legacy behavior: keep only largest component
#         largest_area = 0
#         largest_label = None
#         for region in regions:
#             if region.area > largest_area:
#                 largest_area = region.area
#                 largest_label = region.label
#         # Create a new mask with only the largest component
#         if largest_label is not None:
#             mask = labeled_mask == largest_label
    
#     data_mask = mask * data

#     # Perform dilation
#     if dilation_type == 'disk':
#         data_mask_d = ndimage.binary_dilation(mask,
#                                             structure=disk(dilation_size),
#                                             iterations=iterations).astype(mask.dtype)
#     if dilation_type == 'square':
#         data_mask_d = ndimage.binary_dilation(mask,
#                                             structure=square(dilation_size),
#                                             iterations=iterations).astype(mask.dtype)

#     # Apply the same filtering approach to the 3-sigma mask for plotting
#     labeled_mask3 = label(mask3)
#     regions3 = regionprops(labeled_mask3)
    
#     if use_distance_filter and regions3:
#         # Use same reference point and distance criterion
#         if 'reference_centroid' in locals() and reference_region is not None:
#             mask3 = np.zeros_like(data, dtype=bool)
#             for region in regions3:
#                 centroid = np.array([region.centroid[1], region.centroid[0]])
#                 distance = np.linalg.norm(centroid - reference_centroid)
#                 if distance <= max_distance:
#                     mask3 = mask3 | (labeled_mask3 == region.label)
#         else:
#             # Fallback to largest component
#             if regions3:
#                 largest_area = 0
#                 largest_label = None
#                 for region in regions3:
#                     if region.area > largest_area:
#                         largest_area = region.area
#                         largest_label = region.label
#                 if largest_label is not None:
#                     mask3 = labeled_mask3 == largest_label
#     elif regions3:
#         # Legacy: largest component only
#         largest_area = 0
#         largest_label = None
#         for region in regions3:
#             if region.area > largest_area:
#                 largest_area = region.area
#                 largest_label = region.label
#         if largest_label is not None:
#             mask3 = labeled_mask3 == largest_label

#     if PLOT == True:
#         fig = plt.figure(figsize=(int(16*fig_size_factor), int(4*fig_size_factor)))
#         ax0 = fig.add_subplot(1, 4, 1)
#         ax0.imshow((mask3), origin='lower',cmap='magma')
#         ax0.set_title(r'mask $>$ ' + str(3) + r'$\times\sigma_{\mathrm{rms}}$')
#         ax0.axis('off')
#         ax1 = fig.add_subplot(1, 4, 2)
#         #         ax1.legend(loc='lower left')
#         ax1.imshow((mask), origin='lower',cmap='magma')
#         ax1.set_title(r'mask $>$ ' + str(sigma) + r'$\times\sigma_{\mathrm{rms}}$')
#         ax1.axis('off')
#         ax2 = fig.add_subplot(1, 4, 3)
#         ax2.imshow(data_mask_d, origin='lower',cmap='magma')
#         ax2.set_title(r'D(mask)'f'[{dilation_size}$|${iterations}]'f'{special_name}')
#         ax2.axis('off')
#         ax3 = fig.add_subplot(1, 4, 4)
#         ax3 = eimshow(data * data_mask_d, ax=ax3, CM='magma',
#                       rms=std,add_contours=True,
#                       vmin_factor=0.01,vmax_factor=0.1)
#         ax3.set_title(r'D(mask) $\times$ data')
#         #         ax3.imshow(np.log(data*data_mask_d))
#         fig.subplots_adjust(wspace=0.01, hspace=0.01)
#         #         fig.tight_layout()
#         ax3.axis('off')
#         if show_figure == True:
#             plt.show()
#         else:
#             plt.close()

#     return (mask, data_mask_d)

# def t_mask_dilation(image, cell_size=None, sigma=6, rms=None,
#                   dilation_size=None, iterations=2, dilation_type='disk',
#                   PLOT=False, show_figure=False, logger=None,
#                   fig_size_factor=1,
#                   special_name='', verbose=0):
#     """
#     Mask dilation function with connected component filtering.
    
#     Similar to original function but only keeps the largest connected component,
#     typically at the center of the image.
#     """
#     from skimage.measure import label, regionprops
    
#     if isinstance(image, str) == True:
#         data = load_fits_data(image)
#     else:
#         data = image
#     if rms is None:
#         std = mad_std(data[data > 0])
#     else:
#         std = rms

#     if (dilation_size is None) or (dilation_size == '2x'):
#         try:
#             omaj, omin, _, _, _ = beam_shape(image)
#             if dilation_size == '2x':
#                 dilation_size = int(
#                     2*np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
#             else:
#                 dilation_size = int(
#                     np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
#             if verbose >= 1:
#                 if logger is not None:
#                     logger.debug(f" ==>  Mask dilation size is {dilation_size} [px]")
#                 else:
#                     print(f" ==>  Mask dilation size is {dilation_size} [px]")
#         except:
#             if dilation_size is None:
#                 dilation_size = 5

#     # Create initial masks
#     mask = (data >= sigma * std)
#     mask3 = (data >= 3 * std)
    
#     # Label connected components in the mask
#     labeled_mask = label(mask)
#     regions = regionprops(labeled_mask)
    
#     # Find the largest connected component
#     if regions:
#         largest_area = 0
#         largest_label = None
#         for region in regions:
#             if region.area > largest_area:
#                 largest_area = region.area
#                 largest_label = region.label
#         # Create a new mask with only the largest component
#         if largest_label is not None:
#             mask = labeled_mask == largest_label
    
#     data_mask = mask * data

#     # Perform dilation
#     if dilation_type == 'disk':
#         data_mask_d = ndimage.binary_dilation(mask,
#                                             structure=disk(dilation_size),
#                                             iterations=iterations).astype(mask.dtype)
#     if dilation_type == 'square':
#         data_mask_d = ndimage.binary_dilation(mask,
#                                             structure=square(dilation_size),
#                                             iterations=iterations).astype(mask.dtype)

#     # Clean up the 3-sigma mask for plotting using the same approach
#     labeled_mask3 = label(mask3)
#     regions3 = regionprops(labeled_mask3)
#     if regions3:
#         largest_area = 0
#         largest_label = None
#         for region in regions3:
#             if region.area > largest_area:
#                 largest_area = region.area
#                 largest_label = region.label
#         if largest_label is not None:
#             mask3 = labeled_mask3 == largest_label

#     if PLOT == True:
#         fig = plt.figure(figsize=(int(16*fig_size_factor), int(4*fig_size_factor)))
#         ax0 = fig.add_subplot(1, 4, 1)
#         ax0.imshow((mask3), origin='lower',cmap='magma')
#         ax0.set_title(r'mask $>$ ' + str(3) + r'$\times\sigma_{\mathrm{rms}}$')
#         ax0.axis('off')
#         ax1 = fig.add_subplot(1, 4, 2)
#         #         ax1.legend(loc='lower left')
#         ax1.imshow((mask), origin='lower',cmap='magma')
#         ax1.set_title(r'mask $>$ ' + str(sigma) + r'$\times\sigma_{\mathrm{rms}}$')
#         ax1.axis('off')
#         ax2 = fig.add_subplot(1, 4, 3)
#         ax2.imshow(data_mask_d, origin='lower',cmap='magma')
#         ax2.set_title(r'D(mask)'f'[{dilation_size}$|${iterations}]'f'{special_name}')
#         ax2.axis('off')
#         ax3 = fig.add_subplot(1, 4, 4)
#         ax3 = eimshow(data * data_mask_d, ax=ax3, CM='magma',
#                       rms=std,add_contours=True,
#                       vmin_factor=0.01,vmax_factor=0.1)
#         ax3.set_title(r'D(mask) $\times$ data')
#         #         ax3.imshow(np.log(data*data_mask_d))
#         fig.subplots_adjust(wspace=0.01, hspace=0.01)
#         #         fig.tight_layout()
#         ax3.axis('off')
#         if show_figure == True:
#             plt.show()
#         else:
#             plt.close()

#     return (mask, data_mask_d)



def mask_dilation_from_list(images, residuals, sigma=6, rms=None, dilation_size=None, iterations=2, dilation_type='disk',
                  PLOT=False, show_figure=True, logger=None, fig_size_factor=1, special_name='', verbose=0):
    """
    Mask dilation function for multi-frequency images.
    
    This function applies mask dilation across multiple images and computes a common mask
    that considers emission across all images. It avoids setting a smaller mask as the common mask.
    
    Parameters
    ----------
    images : list of str or ndarray
        List of image names or image data arrays.
    sigma : float
        Sigma level for the mask.
    rms : float
        RMS level for the mask.
    dilation_size : int
        Size of the dilation.
    iterations : int
        Number of iterations for the dilation.
    dilation_type : str
        Type of dilation. Options are 'disk' or 'square'.
    PLOT : bool
        Plot the mask dilation.
    show_figure : bool
        Show the figure.
    logger : logger
        Logger object.
    fig_size_factor : float
        Figure size factor.
    special_name : str
        Special name for the mask dilation.
    """

    # Function to compute mask for a single image
    def compute_mask(image, sigma, rms, dilation_size, iterations, dilation_type, verbose):
        if isinstance(image, str):
            data = load_fits_data(image)
        else:
            data = image
        if rms is None:
            std = mad_std(data)
        else:
            std = rms
        
        if (dilation_size is None) or (dilation_size == '2x'):
            try:
                omaj, omin, _, _, _ = beam_shape(image)
                if dilation_size == '2x':
                    dilation_size = int(2 * np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
                else:
                    dilation_size = int(np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
                if verbose >= 1:
                    if logger is not None:
                        logger.debug(f" ==> Mask dilation size is {dilation_size} [px]")
                    else:
                        print(f" ==> Mask dilation size is {dilation_size} [px]")
            except:
                dilation_size = 5

        mask = (data >= sigma * std)

        if dilation_type == 'disk':
            data_mask_d = ndimage.binary_dilation(mask, structure=disk(dilation_size), iterations=iterations).astype(mask.dtype)
        elif dilation_type == 'square':
            data_mask_d = ndimage.binary_dilation(mask, structure=square(dilation_size), iterations=iterations).astype(mask.dtype)

        return mask, data_mask_d

    # Compute individual masks for each image
    mask_list = []
    for kk in range(len(images)):
        rms = mad_std(load_fits_data(residuals[kk]))
        mask, dilated_mask = compute_mask(images[kk], sigma, rms, dilation_size, iterations, dilation_type, verbose)
        mask_list.append(dilated_mask)

    # Compute the common mask
    combined_mask = np.sum(mask_list, axis=0) >= len(mask_list) // 2  # Majority vote across masks

    # Optional plot of the common mask
    if PLOT:
        fig = plt.figure(figsize=(int(15 * fig_size_factor), int(4 * fig_size_factor)))
        ax = fig.add_subplot(1, 1, 1)
        ax.imshow(combined_mask, origin='lower', cmap='magma')
        ax.set_title('Common Mask')
        ax.axis('off')
        if show_figure:
            plt.show()
        else:
            plt.close()

    return combined_mask




def mask_dilation_from_mask(image, mask_init, cell_size=None, sigma=3,rms=None,
                  dilation_size=None,iterations=2, dilation_type='disk',
                  PLOT=True,show_figure=True):
    """
    Apply a dilation to an existing mask.
    """
    from scipy import ndimage
    from scipy.ndimage import morphology
    from skimage.morphology import disk, square
    from skimage.morphology import dilation

    if isinstance(image, str) == True:
        data = load_fits_data(image)
    else:
        data = image
    if rms is None:
        std = mad_std(data)
    else:
        std = rms

    if dilation_size is None:
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            dilation_size = int(
                np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
        except:
            if dilation_size is None:
                dilation_size = 5

    data_init = data * mask_init
    mask3 = (data >= 3 * std)
    # std = mad_std(data[mask_init])
    mask = (data_init >= sigma * std)
    data_mask = mask * data

    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(mask_init,
                                            structure=disk(dilation_size),
                                            iterations=iterations).\
                                                astype(mask_init.dtype)

    if dilation_type == 'square':
        data_mask_d = ndimage.binary_dilation(mask_init,
                                            structure=square(dilation_size),
                                            iterations=iterations).\
                                                astype(mask_init.dtype)

    if PLOT == True:
        fig = plt.figure(figsize=(15, 4))
        ax0 = fig.add_subplot(1, 4, 1)
        ax0.imshow((mask3), origin='lower')
        ax0.set_title(r'Mask above' + str(3) + '$\sigma$')
        ax0.axis('off')
        ax1 = fig.add_subplot(1, 4, 2)
        #         ax1.legend(loc='lower left')
        ax1.imshow((mask), origin='lower')
        ax1.set_title(r'Mask above' + str(sigma) + '$\sigma$')
        ax1.axis('off')
        ax2 = fig.add_subplot(1, 4, 3)
        ax2.imshow(data_mask_d, origin='lower')
        ax2.set_title(r'Dilated mask')
        ax2.axis('off')
        ax3 = fig.add_subplot(1, 4, 4)
        ax3 = fast_plot2(data * data_mask_d, ax=ax3, vmin_factor=0.1)
        ax3.set_title(r'Dilated mask $\times$ data')
        #         ax3.imshow(np.log(data*data_mask_d))
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        #         fig.tight_layout()
        ax3.axis('off')
        if show_figure == True:
            plt.show()
        else:
            plt.close()

    if cell_size is not None:
        if isinstance(image, str) == True:
            try:
                print((data * data_mask_d).sum() / beam_area2(image, cell_size))
                print((data * data_mask).sum() / beam_area2(image, cell_size))
                print((data).sum() / beam_area2(image, cell_size))
            except:
                print('Provide a cell size of the image.')
    return (mask, data_mask_d)


def grow_region_masks(data_2D, masks_deblended, ref_mask,
                      dilation_size=None, iterations=1, sigma=6.0,
                      grow_to_ref=False, max_steps=None):
    """
    Grow each deblended core into a non-overlapping aperture.

    This is the sub-region growth `structural_morphology` performs, lifted out so
    that the decomposition path can measure per-region quantities over exactly the
    same apertures -- otherwise the two sets of region fluxes would not be
    comparable.

    Two growth modes, selected by `grow_to_ref`.

    `grow_to_ref=False` (default) -- FIXED growth, the historical behaviour and
    what `structural_morphology` still does inline. Purely GEOMETRIC:
    `mask_dilation_from_mask`'s second return is a `disk(dilation_size)` dilation
    repeated `iterations` times, which does not depend on the pixel data or on
    `sigma` (`sigma` only shapes the discarded first return). Three things bound
    it, in order:

    1. `ref_mask` -- never grow outside the reference aperture;
    2. the other regions' UNDILATED cores -- never invade a neighbour's detection;
    3. `claimed_mask` -- never take a pixel an earlier (lower-index) region
       already took. Two adjacent growing masks can otherwise both claim pixels in
       the buffer between them, double-counting that flux.

    `grow_to_ref=True` -- FILL growth. The deblended cores are detected at a
    higher SNR threshold than `ref_mask` is, so a fixed number of dilation steps
    leaves most of `ref_mask` unclaimed and the region fluxes cannot sum to the
    whole-source total. Here every region instead advances one dilation step per
    round, ROUND-ROBIN, until nobody grows any further. The regions then partition
    `ref_mask` rather than nibbling at it.

    Round-robin rather than region-by-region is the whole point: running the fixed
    loop above to convergence one region at a time would let region 0 swallow the
    entire reference mask before region 1 ever moved. Advancing one step each per
    round puts the boundary between two equally sized neighbours near the midpoint.

    Non-overlap is enforced DURING growth here, not trimmed afterwards:
    `ndimage.binary_dilation`'s `mask=` argument is a constrained (geodesic)
    dilation -- only pixels True in `mask` may be switched on -- so a region can
    never leave `ref_mask` nor enter a neighbour in the first place. Ties within a
    round go to the lower region index, the same convention `claimed_mask` uses.

    Pixels of `ref_mask` not connected to any core are never claimed, by
    construction. That is deliberate: they are reported as a coverage shortfall by
    the caller rather than being forced into an arbitrary region.

    DEGENERATE CASE: if one core lies wholly inside another, the inner region is
    walled in and keeps only its single brightest pixel (see the seed rescue
    below). That is a faithful signal -- two detections with no territory that is
    unambiguously the inner one's -- and it is visible as a near-zero region row.
    The alternative, an empty aperture, is much worse: `measures()` falls through
    to `mask` and quietly reports most of the source for that region.

    Parameters
    ----------
    data_2D : np.ndarray
        Image data. Only its shape matters for the geometric dilation.
    masks_deblended : sequence of np.ndarray
        The deblended core masks, one per detected region (`SE.masks`).
    ref_mask : np.ndarray
        Reference aperture bounding all growth.
    dilation_size : int, optional
        Radius of one dilation step. Passed straight through.
    iterations : int, optional
        Number of dilation steps. Fixed mode only; ignored when `grow_to_ref`.
    sigma : float, optional
        Forwarded to `mask_dilation_from_mask`; affects only its discarded first
        return value. Kept so callers can mirror their existing call verbatim.
        Fixed mode only.
    grow_to_ref : bool, optional
        Grow until the regions stop expanding instead of for `iterations` steps.
    max_steps : int, optional
        Hard cap on the number of rounds in fill mode, so a pathological input
        cannot spin. Defaults to enough rounds to cross the image.

    Returns
    -------
    list of np.ndarray
        Boolean apertures, one per region, in the order given, guaranteed
        mutually disjoint.
    """
    n = len(masks_deblended)
    ref_bool = np.asarray(ref_mask).astype(bool)
    cores = [np.asarray(m).astype(bool) for m in masks_deblended]

    if not grow_to_ref:
        claimed_mask = np.zeros(ref_bool.shape, dtype=bool)
        grown = []
        for j in range(n):
            _, mask_new = mask_dilation_from_mask(data_2D,
                                                  masks_deblended[j],
                                                  sigma=sigma,
                                                  PLOT=False,
                                                  iterations=iterations,
                                                  dilation_size=dilation_size,
                                                  show_figure=False)
            mask_new = mask_new.astype(bool) & ref_bool
            for l in range(n):
                if l != j:
                    mask_new[cores[l]] = False
            mask_new = mask_new & ~claimed_mask
            claimed_mask = claimed_mask | mask_new
            grown.append(mask_new)
        return grown

    if dilation_size is None:
        dilation_size = 5
    dilation_size = max(int(dilation_size), 1)
    selem = disk(dilation_size)
    if max_steps is None:
        # Enough rounds for one region to cross the whole image; the loop
        # normally exits long before this on "nobody grew".
        max_steps = int(np.ceil(max(ref_bool.shape) / dilation_size)) + 2

    # A core lying partly outside ref_mask keeps only its inside part -- growth
    # is confined to ref_mask, so the outside part could never be claimed anyway.
    # The cores themselves are not guaranteed disjoint (SEP's `masks` are scaled
    # ellipses, which can overlap), and `binary_dilation` leaves input pixels on
    # even where `mask` forbids growth, so any overlap present at the start would
    # survive every round. Resolve it here, lower index first, so the returned
    # apertures are disjoint from round zero.
    grown = []
    claimed_mask = np.zeros(ref_bool.shape, dtype=bool)
    for j in range(n):
        start_j = (cores[j] & ref_bool) & ~claimed_mask
        claimed_mask = claimed_mask | start_j
        grown.append(start_j)

    # A core wholly contained in an earlier region's core is left with an empty
    # seed, and an empty mask dilates to an empty mask, so that region would
    # never grow at all and would silently measure nothing. This is not a corner
    # case: SEP's `masks` are ellipses scaled off the second moments and routinely
    # overlap each other several times over. Hand such a region back its single
    # brightest pixel, taken from whichever region currently holds it, which keeps
    # the apertures disjoint and lets the round-robin sort out the rest.
    _data = load_fits_data(data_2D) if isinstance(data_2D, str) else data_2D
    _data = np.nan_to_num(np.asarray(_data), nan=-np.inf)
    for j in range(n):
        inside_j = cores[j] & ref_bool
        if grown[j].any() or not inside_j.any():
            continue
        flat = np.where(inside_j.ravel(), _data.ravel(), -np.inf)
        seed_idx = np.unravel_index(int(np.argmax(flat)), _data.shape)
        for l in range(n):
            if l != j and grown[l][seed_idx]:
                grown[l][seed_idx] = False
        grown[j][seed_idx] = True

    for _ in range(int(max_steps)):
        changed = False
        for j in range(n):
            others = np.zeros(ref_bool.shape, dtype=bool)
            for l in range(n):
                if l != j:
                    others |= grown[l] | cores[l]
            allowed = ref_bool & ~others
            new_j = ndimage.binary_dilation(grown[j], structure=selem,
                                            iterations=1, mask=allowed)
            new_j = new_j.astype(bool)
            if int(new_j.sum()) != int(grown[j].sum()):
                changed = True
            grown[j] = new_j
        if not changed:
            break
    return grown


def convert_seg_map_to_masks(seg_map):
    """
    Convert a segmentation map to a list of masks, where each mask corresponds to a unique segment.

    Parameters:
        seg_map (2D array): Segmentation map where different segments have different integer values,
                            and masked regions have a value of zero.

    Returns:
        list: A list of boolean masks, where each mask corresponds to a segment in the segmentation map.
    """
    # Find the unique segment labels, excluding zero (the background)
    unique_segments = np.unique(seg_map)
    unique_segments = unique_segments[unique_segments != 0]

    # Create a mask for each segment
    masks = [(seg_map == segment) for segment in unique_segments]

    return masks

def convert_masks_to_segmentation(masks_regions):
    # Assuming all masks have the same shape
    height, width = masks_regions[0].shape
    
    # Create an empty array for the segmentation image
    segmentation_image = np.zeros((height, width), dtype=int)
    
    # Loop over each mask and assign a unique integer value
    for idx, mask in enumerate(masks_regions, start=1):
        # idx+1 ensures the values start from 1 instead of 0
        segmentation_image[mask == 1] = idx
    
    return segmentation_image


def split_overlapping_masks(mask1, mask2):
    """
    # Compute the overlap between the masks

    GIven two masks, find the overlapping region between the two and
    split half of the pixels in this region attributing them to one maks and the
    other half to the other mask.

    # Check for overlap between the masks and split the overlapping region
    mask1_new, mask2_new = split_overlapping_masks(mask1, mask2)

    # Show the new masks
    print(mask1_new)
    print(mask2_new)

    """

    overlap_mask = (mask1 > 0) & (mask2 > 0)

    # If there is no overlap, return the original masks
    if not np.any(overlap_mask):
        return mask1, mask2

    # Split the overlapping region in half
    half_mask1 = mask1.copy()
    half_mask2 = mask2.copy()
    overlap_indices = np.where(overlap_mask)
    num_overlapping_pixels = len(overlap_indices[0])
    half_num_overlapping_pixels = num_overlapping_pixels // 2
    for i in range(num_overlapping_pixels):
        row, col = overlap_indices[0][i], overlap_indices[1][i]
        if i < half_num_overlapping_pixels:
            half_mask1[row, col] = mask1[row, col]
            half_mask2[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = mask2[row, col]

    # Return the split masks
    return half_mask1, half_mask2




def split_overlapping_masks3(mask1, mask2, mask3):
    """
    # Compute the overlap between the masks

    The same as `split_overlapping_masks`, but for three masks.
    """
    overlap_mask1 = (mask1 > 0) & (mask2 > 0) & (mask3 == 0)
    overlap_mask2 = (mask1 > 0) & (mask3 > 0) & (mask2 == 0)
    overlap_mask3 = (mask2 > 0) & (mask3 > 0) & (mask1 == 0)

    # If there is no overlap, return the original masks
    if not np.any(overlap_mask1) and not np.any(overlap_mask2) and not np.any(
            overlap_mask3):
        return mask1, mask2, mask3

    # Split the overlapping region in half for each overlapping pair
    num_overlapping_pixels1 = np.sum(overlap_mask1)
    num_overlapping_pixels2 = np.sum(overlap_mask2)
    num_overlapping_pixels3 = np.sum(overlap_mask3)
    half_num_overlapping_pixels1 = num_overlapping_pixels1 // 2
    half_num_overlapping_pixels2 = num_overlapping_pixels2 // 2
    half_num_overlapping_pixels3 = num_overlapping_pixels3 // 2

    half_mask1 = mask1.copy()
    half_mask2 = mask2.copy()
    half_mask3 = mask3.copy()

    overlap_indices1 = np.where(overlap_mask1)
    for i in range(num_overlapping_pixels1):
        row, col = overlap_indices1[0][i], overlap_indices1[1][i]
        if i < half_num_overlapping_pixels1:
            half_mask1[row, col] = mask1[row, col]
            half_mask2[row, col] = 0
            half_mask3[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = mask2[row, col]
            half_mask3[row, col] = 0

    overlap_indices2 = np.where(overlap_mask2)
    for i in range(num_overlapping_pixels2):
        row, col = overlap_indices2[0][i], overlap_indices2[1][i]
        if i < half_num_overlapping_pixels2:
            half_mask1[row, col] = mask1[row, col]
            half_mask2[row, col] = 0
            half_mask3[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = 0
            half_mask3[row, col] = mask3[row, col]

    overlap_indices3 = np.where(overlap_mask3)
    for i in range(num_overlapping_pixels3):
        row, col = overlap_indices3[0][i], overlap_indices3[1][i]
        if i < half_num_overlapping_pixels3:
            half_mask1[row, col] = 0
            half_mask2[row, col] = mask2[row, col]
            half_mask3[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = 0
            half_mask3[row, col] = mask3[row, col]

    # Return the split masks
    return half_mask1, half_mask2, half_mask3


def trail_vector(vx, vy, v0=np.asarray([1, 0])):
    """
    Compute the trailing vector angle (in degrees) and its magnitude.
    The angle (PA) is measured counterclockwise from the positive x-axis.

    Parameters:
    vx, vy : float
        Components of the input vector.
    v0 : array-like, optional
        Reference vector (default is [1, 0] for the positive x-axis).

    Returns:
    angle_PA : float
        Position angle (PA) measured counterclockwise from the positive x-axis, in degrees.
    norm_vec : float
        Magnitude of the input vector.
    """
    from scipy.linalg import norm
    v = np.asarray([vx, vy])
    norm_vec = norm(v)
    
    # Handle zero vector case
    if norm_vec == 0:
        raise ValueError("Input vector (vx, vy) has zero magnitude.")
    
    # Normalize input vector
    v_hat = v / norm_vec
    
    # Compute the angle using arctan2 (counterclockwise convention)
    angle_radians = np.arctan2(v_hat[1], v_hat[0])  # atan2(vy, vx)
    angle_PA = np.degrees(angle_radians)  # Convert to degrees
    
    # Ensure the angle is in the range [0, 360)
    if angle_PA < 0:
        angle_PA += 360
    return angle_PA, norm_vec


def calculate_radii(image, mask):
    """
    Calculate circular radii of the emission in a 2D numpy image.

    Parameters:
        image (np.ndarray): The 2D numpy image containing the radio emission of a galaxy.
        sigma_level (float): The sigma level to use for the contour.
        background_std (float): The standard deviation of the background.

    Returns:
        float: The circular radii of the emission.
    """
    # Calculate the threshold level for the given sigma level
    #     threshold = sigma_level * background_std

    # Create a binary mask of the emission above the threshold level
    #     mask = image > threshold

    # Calculate the center of mass of the emission
    from scipy.ndimage import measurements
    com = measurements.center_of_mass(image * mask)

    # Calculate the distance of each pixel from the center of mass
    y, x = np.indices(mask.shape)
    distances = np.sqrt((x - com[1]) ** 2 + (y - com[0]) ** 2)

    # Calculate the median distance of the pixels above the threshold level
    median_distance = np.nanmedian(distances[mask])
    mean_distance = np.nanmean(distances[mask])
    std_distance = np.nanstd(distances[mask])
    # mad_std_distance = mad_std(distances[mask])

    return median_distance, mean_distance, std_distance



# Testing functions, not used anywhere.
def shannon_entropy_2d(arr):
    rows, cols = arr.shape
    result = np.zeros((rows, cols))
    for i in range(2, rows-2):
        for j in range(2, cols-2):
            box = arr[i-2:i+3, j-2:j+3].ravel()
            p = box / np.nansum(box)
            entropy = -np.nansum(p * np.log2(p))
            result[i, j] = entropy
    return result

def cvi(imname):
    try:
        os.system('casaviewer ' + imname)
    except:
        try:
            os.system('~/casaviewer ' + imname)
        except:
            pass





def get_image_statistics(imagename,cell_size=None,
                         mask_component=None,mask=None,
                         residual_name=None,region='', dic_data=None,
                         sigma_mask=6,apply_mask=True,
                         fracX=0.1, fracY=0.1):
    """
    Get some basic image statistics.



    """
    if dic_data is None:
        dic_data = {}
        dic_data['#imagename'] = os.path.basename(imagename)

    if cell_size is None:
        cell_size = get_cell_size(imagename)

    image_data = load_fits_data(imagename)
    if mask is not None:
        image_data = np.nan_to_num(image_data*mask,nan=0)
        apply_mask = False
    #     dic_data['imagename'] = imagename
    if apply_mask == True:
        omask, mask = mask_dilation(imagename, sigma=sigma_mask)
        image_data = np.nan_to_num(image_data*mask,nan=0)
    if mask_component is not None:
        image_data = np.nan_to_num(image_data * mask_component,nan=0)
        if mask is None:
            mask = mask_component
        # mask = mask_component
    if (mask_component is None) and (apply_mask == False) and (mask is None):
        unity_mask = np.ones(load_fits_data(imagename).shape)
        omask, mask = unity_mask, unity_mask
        image_data = np.nan_to_num(load_fits_data(imagename) * mask,nan=0)

    stats_im = imstat(imagename=imagename, region=region)

    box_edge, imhd = create_box(imagename, fracX=fracX, fracY=fracY)
    stats_box = imstat(imagename=imagename, box=box_edge)

    # determine the flux peak and positions of image
    flux_peak_im = stats_im['max'][0]
    flux_min_im = stats_im['min'][0]
    dic_data['max_im'] = flux_peak_im
    dic_data['min_im'] = flux_min_im
    """
    x0max,y0max = peak_center(image_data)
    dic_data['x0'], dic_data['y0'] = x0max,y0max
    # determine momentum centres.
    x0m, y0m, _, _ = momenta(image_data, PArad_0=None, q_0=None)
    dic_data['x0m'], dic_data['y0m'] = x0m, y0m

    #some geometrical measures
    # calculate PA and axis-ratio
    PA, q, x0col, y0col, PAm, qm, PAmi, qmi, PAmo, qmo,\
        x0median,y0median,\
        x0median_i,y0median_i,x0median_o,y0median_o = cal_PA_q(image_data)

    dic_data['PA'], dic_data['q'] = PA, q
    dic_data['PAm'], dic_data['qm'] = PAm, qm
    dic_data['PAm'], dic_data['qm'] = PAm, qm
    dic_data['PAmi'], dic_data['qmi'] = PAmi, qmi
    dic_data['PAmo'], dic_data['qmo'] = PAmo, qmo
    dic_data['x0m_i'], dic_data['y0m_i'] = x0median_i, y0median_i
    dic_data['x0m_o'], dic_data['y0m_o'] = x0median_o, y0median_o
    """

    # determine the rms and std of residual and of image
    rms_im = stats_im['rms'][0]
    rms_box = stats_box['rms'][0]
    sigma_im = stats_im['sigma'][0]
    sigma_box = stats_box['sigma'][0]

    dic_data['rms_im'] = rms_im
    dic_data['rms_box'] = rms_box
    dic_data['sigma_im'] = sigma_im
    dic_data['sigma_box'] = sigma_box

    # determine the image and residual flux
    flux_im = stats_im['flux'][0]
    flux_box = stats_box['flux'][0]
    dic_data['flux_im'] = flux_im
    dic_data['flux_box'] = flux_box
    sumsq_im = stats_im['sumsq'][0]
    sumsq_box = stats_box['sumsq'][0]

    q_sq = sumsq_im / sumsq_box
    q_flux = flux_im / flux_box
    # flux_ratio = flux_re/flux_im
    dic_data['q_sq'] = q_sq
    dic_data['q_flux'] = q_flux

    snr = flux_im / rms_box
    snr_im = flux_im / rms_im

    dr_e = []
    frac_ = np.linspace(0.05, 0.85, 10)
    frac_image = 0.10
    '''
    Each loop below run a sliding window,
    one to the x direction and the other to y-dircetion.
    This is to get a better estimate (in multiple regions)
    of the background rms and therefore SNR
    Each window has a fraction frac_image of the image size.
    '''
    for frac in frac_:
        box, _ = create_box(imagename, fracX=frac, fracY=frac_image)
        st = imstat(imagename, box=box)
        snr_tmp = flux_peak_im / st['rms'][0]
        dr_e.append(snr_tmp)

    dr_e2 = []
    for frac in frac_:
        box, _ = create_box(imagename, fracX=frac_image, fracY=frac)
        st = imstat(imagename, box=box)
        snr_tmp = flux_peak_im / st['rms'][0]
        dr_e2.append(snr_tmp)
    #average of the SNR -- DINAMIC RANGE
    DR_SNR_E = (np.mean(dr_e) + np.mean(dr_e2)) / 2

    dic_data['snr'] = snr
    dic_data['snr_im'] = snr_im
    dic_data['DR_SNR_E'] = DR_SNR_E

    DR_pk_rmsbox = flux_peak_im / rms_box
    DR_pk_rmsim = flux_peak_im / rms_im
    dic_data['DR_pk_rmsbox'] = DR_pk_rmsbox
    dic_data['DR_pk_rmsim'] = DR_pk_rmsim

    dic_data['bmajor'] = imhd['restoringbeam']['major']['value']
    dic_data['bminor'] = imhd['restoringbeam']['minor']['value']
    dic_data['positionangle'] = imhd['restoringbeam']['positionangle']['value']



    # if residual_name is not None:
    #     data_res = load_fits_data(residual_name)
    #     flux_res_error = 3 * np.sum(data_res * mask) \
    #                      / beam_area2(imagename, cell_size)
    #     # rms_res =imstat(residual_name)['flux'][0]
    #     flux_res = np.sum(load_fits_data(residual_name)) / beam_area2(imagename, cell_size)
    #
    #     res_error_rms =np.sqrt(
    #         np.sum((abs(data_res * mask -
    #                     np.mean(data_res * mask))) ** 2 * np.sum(mask))) / \
    #                    beam_area2(imagename,cell_size)
    #
    #     try:
    #         total_flux_tmp = dic_data['total_flux_mask']
    #     except:
    #         total_flux_tmp = flux_im
    #         total_flux_tmp = total_flux(image_data,imagename,mask=mask)
    #
    #     # print('Estimate #1 of flux error (based on sum of residual map): ')
    #     # print('Flux = ', total_flux_tmp * 1000, '+/-',
    #     #       abs(flux_res_error) * 1000, 'mJy')
    #     # print('Fractional error flux = ', flux_res_error / total_flux_tmp)
    #     print('-----------------------------------------------------------------')
    #     print('Estimate of flux error (based on rms of '
    #           'residual x area): ')
    #     print('Flux = ', total_flux_tmp * 1000, '+/-',
    #           abs(res_error_rms) * 1000, 'mJy')
    #     print('Fractional error flux = ', res_error_rms / total_flux_tmp)
    #     print('-----------------------------------------------------------------')
    #
    #     dic_data['max_residual'] = np.max(data_res * mask)
    #     dic_data['min_residual'] = np.min(data_res * mask)
    #     dic_data['flux_residual'] = flux_res
    #     dic_data['flux_error_res'] = abs(flux_res_error)
    #     dic_data['flux_error_res_2'] = abs(res_error_rms)
    #     dic_data['mad_std_residual'] = mad_std(data_res)
    #     dic_data['rms_residual'] = rms_estimate(data_res)

    #     print(' Flux=%.5f Jy/Beam' % flux_im)
    #     print(' Flux peak (image)=%.5f Jy' % flux_peak_im, 'Flux peak (residual)=%.5f Jy' % flux_peak_re)
    #     print(' flux_im/sigma_im=%.5f' % snr_im, 'flux_im/sigma_re=%.5f' % snr)
    #     print(' rms_im=%.5f' % rms_im, 'rms_re=%.5f' % rms_re)
    #     print(' flux_peak_im/rms_im=%.5f' % peak_im_rms, 'flux_peak_re/rms_re=%.5f' % peak_re_rms)
    #     print(' sumsq_im/sumsq_re=%.5f' % q)
    return (dic_data)


def level_statistics(img, cell_size=None, mask_component=None,
                    sigma=6, do_PLOT=False, crop=False,data_2D = None,data_res=None,
                    box_size=256, bkg_to_sub=None, apply_mask=True,
                    mask=None,rms=None,
                    results=None, dilation_size=None, iterations=2,
                    add_save_name='', SAVE=False, show_figure=False, ext='.jpg'):
    """
    Function old name: plot_values_std

    Slice the intensity values of an image into distinct regions.
    Then, compute information for each bin level of the emission.
    The implemented splitting is:

        1. Inner region: peak intensity -> 0.1 * peak intensity
        2. Mid-region: 0.1 * peak intensity -> 10 * rms
        3. Low-region: 10 * rms -> 6 * rms
        4. Uncertain region: 6 * rms -> 3 * rms

    Parameters
    ----------
    img : str
        Path to the image.
    cell_size : float, optional
        Cell size of the image. The default is None. In that case, the function
        get_cell_size will attempt to estimate it from the header of the image.
    mask_component : array, optional
        The default is None. This is designed to be used when an image is complex,
        and one would like to study multiple components of the emission separately,
        each one at a time.
    sigma : float, optional
        The default is 6. This is the number of standard deviations to be used during
        mask dilation.
    do_PLOT : bool, optional
        The default is False. If True, the function will plot and save the image.
    crop : bool, optional
        The default is False. If True, the function will crop the image to a square
        of size box_size.
    box_size : int, optional
        The default is 256. This is the size of the square to be used if crop=True.
    data_2D : array, optional
        The default is None. If not None, the function will use this array instead
        and consider header information from img to be used with the array data_2D.

    """
    if cell_size is None:
        cell_size = get_cell_size(img)
    if data_2D is not None:
        g_ = data_2D
    else:
        g_ = load_fits_data(img)
    g = g_.copy()
    if rms is None:
        std = mad_std(g_)
    else:
        std = rms

    if bkg_to_sub is not None:
        g = g - bkg_to_sub
    if mask_component is not None:
        g = g * mask_component

    beam_area_ = beam_area2(img)

    if mask is not None:
        g = g * mask
        apply_mask = False  # do not calculate the mask again, in case is True.
        g = g * mask
    if apply_mask == True:
        _, mask_dilated = mask_dilation(img, cell_size=cell_size, sigma=sigma,
                                        dilation_size=dilation_size, rms=rms,
                                        iterations=iterations,
                                        PLOT=False)
        g = g * mask_dilated


    g = np.nan_to_num(g,nan=0)
    
    if mask_component is not None:
        levels = np.geomspace(g.max(), (1 * std), 5)
        levels_top = np.geomspace(g.max(), g.max() * 0.1, 3)
        try:
            levels_mid = np.geomspace(g.max() * 0.1, (10 * std), 5)
        except:
            levels_mid = np.asarray([0])
        try:
            levels_low = np.geomspace(10 * std, (6.0  * std), 2)
            levels_uncertain = np.geomspace(6.0 * std, (3.0 * std), 3)
        except:
            levels_low = np.asarray([0])
            levels_uncertain = np.asarray([0])

    else:
        if apply_mask is not False:
            # print('asdasd', g.max(), std)
            levels = np.geomspace(np.nanmax(g), (1 * std), 5)
            levels_top = np.geomspace(np.nanmax(g), np.nanmax(g) * 0.1, 3)
            try:
                levels_mid = np.geomspace(np.nanmax(g) * 0.1, (10 * std), 5)
            except:
                levels_mid = np.asarray([0])
            try:
                levels_low = np.geomspace(10 * std, (6.0  * std), 2)
                levels_uncertain = np.geomspace(6 * std, (3.0 * std), 3)
            except:
                levels_low = np.asarray([0])
                levels_uncertain = np.asarray([0])
        else:
            levels = np.geomspace(np.nanmax(g), (3 * std), 5)
            levels_top = np.geomspace(np.nanmax(g), np.nanmax(g) * 0.1, 3)
            # levels_mid = np.geomspace(np.nanmax(g) * 0.1, (10 * std + dl), 5)
            try:
                levels_mid = np.geomspace(np.nanmax(g) * 0.1, (10 * std), 5)
            except:
                levels_mid = np.asarray([0])
            try:
                levels_low = np.geomspace(10 * std, (6.0  * std), 2)
                levels_uncertain = np.geomspace(3 * std, (1.0 * std), 3)
            except:
                levels_low = np.asarray([0])
                levels_uncertain = np.asarray([0])

    # pix_inner = g[g >= levels_top[-1]]
    # pix_mid = g[np.where((g < levels_top[-1]) & (g >= levels_mid[-1]))]
    # pix_low = g[np.where((g < levels_mid[-1]) & (g >= levels_low[-1]))]
    # pix_uncertain = g[np.where((g < levels_low[-1]) & (g >= levels_uncertain[-1]))]
    pix_inner_mask = (g >= levels_top[-1])
    pix_inner = g * pix_inner_mask
    pix_mid_mask = (((g < levels_top[-1]) & (g >= levels_mid[-1])))
    pix_mid = g * pix_mid_mask
    pix_low_mask = (((g < levels_mid[-1]) & (g >= levels_low[-1])))
    pix_low = g * pix_low_mask
    pix_uncertain_mask = (((g < levels_low[-1]) & (g >= levels_uncertain[-1])))
    pix_uncertain = g * pix_uncertain_mask
    inner_flux = pix_inner.sum() / beam_area_
    mid_flux = pix_mid.sum() / beam_area_
    low_flux = pix_low.sum() / beam_area_
    uncertain_flux = pix_uncertain.sum() / beam_area_

    total_flux = low_flux + mid_flux + inner_flux + uncertain_flux
    ext_flux = low_flux + mid_flux + uncertain_flux
    pix_area = len(g[g >= 3 * std])
    number_of_beams = pix_area / beam_area_
    n_beams_inner = np.sum(pix_inner_mask) / beam_area_
    n_beams_mid = np.sum(pix_mid_mask) / beam_area_
    n_beams_low = np.sum(pix_low_mask) / beam_area_
    n_beams_uncertain = np.sum(pix_uncertain_mask) / beam_area_
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(img)

    # print('Low Flux (extended) Jy                    > ', low_flux, ' >> ratio=',
    #       low_flux / total_flux)
    # print('Mid Flux (outer core + inner extended) Jy > ', mid_flux, ' >> ratio=',
    #       mid_flux / total_flux)
    # print('Inner Flux (core) Jy                      > ', inner_flux,
    #       ' >> ratio=', inner_flux / total_flux)
    # print('Uncertain Flux (<5std)                    > ', uncertain_flux,
    #       ' >> ratio=', uncertain_flux / total_flux)
    # print('Total Flux Jy                             > ', total_flux)
    # print('Total area (in # ob beams)                > ', number_of_beams)
    # print('Total inner area (in # ob beams)          > ', n_beams_inner)
    # print('Total mid area (in # ob beams)            > ', n_beams_mid)
    # print('Total low area (in # ob beams)            > ', n_beams_low)
    # print('Total uncertain area (in # ob beams)      > ', n_beams_uncertain)
    # print('Inner Flux (core) fraction                > ',
    #       inner_flux / total_flux)
    # print('Outer Flux (ext)  fraction                > ', ext_flux / total_flux)

    results['peak_of_flux'] = np.max(g)
    results['total_flux'] = total_flux
    results['inner_flux'] = inner_flux
    results['low_flux'] = low_flux
    results['mid_flux'] = mid_flux
    results['uncertain_flux'] = uncertain_flux

    results['inner_flux_f'] = inner_flux / total_flux
    results['low_flux_f'] = low_flux / total_flux
    results['mid_flux_f'] = mid_flux / total_flux
    results['uncertain_flux_f'] = uncertain_flux / total_flux

    results['number_of_beams'] = number_of_beams
    results['n_beams_inner'] = n_beams_inner
    results['n_beams_mid'] = n_beams_mid
    results['n_beams_low'] = n_beams_low
    results['n_beams_uncertain'] = n_beams_uncertain

    if do_PLOT == True:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot()
        vmin = 1 * std
        vmax = g_.max()
        # norm = visualization.simple_norm(g, stretch='log')#, max_percent=max_percent_lowlevel)
        norm = simple_norm(g, stretch='asinh', asinh_a=0.005, vmin=vmin,
                           vmax=vmax)

        if crop == True:
            try:
                xin, xen, yin, yen = do_cutout(img, box_size=box_size,
                                               center=center, return_='box')
                g = g[xin:xen, yin:yen]
            except:
                try:
                    max_x, max_y = np.where(g == g.max())
                    xin = max_x[0] - box_size
                    xen = max_x[0] + box_size
                    yin = max_y[0] - box_size
                    yen = max_y[0] + box_size
                    g = g[xin:xen, yin:yen]
                except:
                    pass
        try:
            im_plot = ax.imshow(g, cmap='magma_r', norm=norm, alpha=1.0,
                                origin='lower')
            ax.contour(g, levels=levels_top[::-1], colors='lime',
                       alpha=1.0)  # cmap='Reds', linewidths=0.75)
            ax.contour(g, levels=levels_mid[::-1], colors='yellow',
                       linewidths=0.75)
            ax.contour(g, levels=levels_low[::-1],
                       colors='#56B4E9')  # cmap='Greens', linewidths=0.75)
            ax.contour(g, levels=levels_uncertain[::-1], colors='grey',
                       linewidths=0.4)
        except:
            pass
        # im_plot.colorbar()
        #         plt.subplots_adjust(wspace=0, hspace=0)
        #         fig.tight_layout()

        if SAVE is not None:
            plt.savefig(
                img.replace('.fits', '_std_levels') + add_save_name + ext,
                dpi=300,
                bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()
    return (results)


def compute_asymetries(imagename,mask,mask_component=None,
                       bkg_to_sub=None,
                       centre=None,results=None):
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(imagename)

    if isinstance(imagename, str) == True:
        image_data = load_fits_data(imagename)
    else:
        image_data = imagename

    if bkg_to_sub is not None:
        image_data = image_data - bkg_to_sub
    if mask_component is not None:
        image_data = image_data * mask_component
    if centre is None:
        pass
    if centre is not None:
        x0A, y0A = centre
    if (mask is None) and (mask_component is not None):
        mask = mask_component
    else:
        unity_mask = np.ones(image_data.shape) == 1
        omask, mask = unity_mask, unity_mask

    try:
        BGrandom, BGmedian, BGmin, BGstd, BGx0, BGy0 \
            = background_asymmetry(image_data, mask, pre_clean=False)
    except:
        print('Error computing background assymetry.')
        BGrandom, BGmedian, BGmin, BGstd, BGx0, BGy0 = 0.0, 0.0, 0.0, 0.0, x0A, y0A

    x0A0fit, y0A0fit = fmin(assimetria0, (x0A, y0A),
                            args=(image_data, mask,), disp=0)
    x0A1fit, y0A1fit = fmin(assimetria1, (x0A, y0A),
                            args=(image_data, mask,), disp=0)
    A0 = assimetria0((x0A0fit, y0A0fit), image_data, mask) - BGmedian
    A1 = assimetria1((x0A1fit, y0A1fit), image_data, mask)
    results['A_BK_median'] = BGmedian
    results['A0'] = A0
    results['A1'] = A1
    results['x0A0fit'] = x0A0fit
    results['y0A0fit'] = y0A0fit
    results['x0A1fit'] = x0A1fit
    results['y0A1fit'] = y0A1fit
    return(results)


def convex_shape(mask):
    from scipy.spatial import ConvexHull
    indices = np.transpose(np.nonzero(mask))
    hull = ConvexHull(indices)
    convex_area = hull.area
    convex_perimeter = hull.volume
    return(convex_area,convex_perimeter)

def shape_measures(imagename, residualname, z, mask_component=None, sigma_mask=6,
             last_level=3.0, vmin_factor=1.0, plot_catalog=False,data_2D=None,
             npixels=128, fwhm=81, kernel_size=21, dilation_size=None,
             main_feature_index=0, results_final={}, iterations=2,
             fracX=0.10, fracY=0.10, deblend=False, bkg_sub=False,
             bkg_to_sub=None, rms=None,
             apply_mask=True, do_PLOT=False, SAVE=True, show_figure=True,
             mask=None,do_measurements='all',
             add_save_name=''):
    """
    Main function that perform other function calls responsible for
    all relevant calculations on the images.
    """
    cell_size = get_cell_size(imagename)
    """
    One beam area is one element resolution, so we avoind finding sub-components 
    that are smaller than the beam area. 
    """
    min_seg_pixels = beam_area2(imagename)

    if mask is not None:
        mask = mask
        apply_mask = False

    if apply_mask == True:
        _, mask_dilated = mask_dilation(imagename, cell_size=cell_size,
                                        sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations, rms=rms,
                                        PLOT=True)
        mask = mask_dilated
    else:
        mask = None

    if data_2D is not None:
        data_2D = data_2D
    else:
        data_2D = load_fits_data(imagename)

    results_final = None #start dict to store measurements.

    levels, fluxes, agrow, plt, \
        omask2, mask2, results_final = make_flux_vs_std(imagename,
                                                        cell_size=cell_size,
                                                        residual=residualname,
                                                        mask_component=mask_component,
                                                        last_level=last_level,
                                                        sigma_mask=sigma_mask,
                                                        apply_mask=False,
                                                        data_2D = data_2D,
                                                        mask=mask,
                                                        rms=rms,
                                                        vmin_factor=vmin_factor,
                                                        results=results_final,
                                                        show_figure=show_figure,
                                                        bkg_to_sub=bkg_to_sub,
                                                        add_save_name=add_save_name,
                                                        SAVE=SAVE)
    error_petro = False
    if z is not None:
        z = z
    else:
        z = 0.01

    results_final['error_petro'] = error_petro
    if mask_component is not None:
        r, ir = get_profile(load_fits_data(imagename) * mask_component)
    else:
        r, ir = get_profile(imagename)
    #     rpix = r / cell_size
    rpix = r.copy()
    r_list_arcsec = rpix * cell_size
    Rp_arcsec = results_final['C90radii'] * cell_size
    R50_arcsec = results_final['C50radii'] * cell_size
    pix_to_pc = pixsize_to_pc(z=z, cell_size=cell_size)
    r_list_pc = rpix * pix_to_pc
    Rp_pc = results_final['C90radii'] * pix_to_pc
    R50_pc = results_final['C50radii'] * pix_to_pc
    Rp_arcsec = results_final['C90radii'] * cell_size
    R50_arcsec = results_final['C50radii'] * cell_size
    r_list_arcsec = rpix * cell_size
    Rp_arcsec = results_final['C90radii'] * cell_size
    R50_arcsec = results_final['C50radii'] * cell_size


    results_final = cosmo_stats(imagename=imagename, z=z, results=results_final)

    results_final['pix_to_pc'] = pix_to_pc
    results_final['cell_size'] = cell_size
    # if error_petro == True:
    # results_final['area_beam_2Rp'] = area_beam[
    #     int(results_final['C90radii'])]
    # results_final['area_beam_R50'] = area_beam[
    #     int(results_final['C50radii'])]

    df = pd.DataFrame.from_dict(results_final, orient='index').T
    df.to_csv(imagename.replace('.fits', add_save_name + '_area_stats.csv'),
              header=True,
              index=False)

    return (results_final, mask)


def compute_flux_density(imagename, residualname, mask=None, sigma=6,
                         systematic_error_fraction=0.05):
    beam_area = beam_area2(imagename)
    image_data = load_fits_data(imagename)
    rms = mad_std(load_fits_data(residualname))
    if mask is None:
        _, mask = mask_dilation(imagename, show_figure=False,
                                PLOT=False, rms=rms, sigma=sigma,
                                      iterations=2)
    total_flux_density = np.nansum(image_data * mask) / beam_area

    # compute error in flux density
    data_res = load_fits_data(residualname)
    res_error_rms = np.sqrt(np.nansum(
        (abs(data_res * mask - np.nanmean(data_res * mask))) ** 2 * np.nansum(
            mask))) / beam_area

    total_flux_density_residual = np.nansum(data_res * mask) / beam_area

    # Calculate systematic error as a fraction of the total flux density of the image
    systematic_error = systematic_error_fraction * total_flux_density

    # Calculate total flux density error in quadrature
    total_flux_density_error = np.sqrt(
        systematic_error ** 2 + (rms / beam_area) ** 2 + total_flux_density_residual ** 2)

    # res_error_rms = 3 * np.nansum(data_res * mask)/beam_area

    print('-----------------------------------------------------------------')
    print('Estimate of flux error (based on rms of '
          'residual x area): ')
    print('Flux Density = ', total_flux_density * 1000, '+/-',
          total_flux_density_error * 1000, 'mJy')
    print('Fractional error flux = ', total_flux_density_error / total_flux_density)
    print('-----------------------------------------------------------------')
    # print(f"Flux Density = {total_flux_density*1000:.2f} mJy")
    return (total_flux_density, total_flux_density_error)

def measures(imagename, residualname, z, mask_component=None, sigma_mask=6,
             last_level=3.0, vmin_factor=1.0, plot_catalog=False,
             data_2D=None,data_res=None, flux_units = 'Jy/beam',
             flux_conversion_factor=None,
             npixels=128, fwhm=81, kernel_size=21, dilation_size=None,
             main_feature_index=0, results_final={}, iterations=2,
             fracX=0.10, fracY=0.10, deblend=False, bkg_sub=False,
             bkg_to_sub=None, rms=None,do_petro=True,
             error_map=None, rms_map=None, invvar_map=None, variance_map=None,
             weight_map=None, use_residual_as_error=False,
             crop=False, box_size=256,
             apply_mask=True, do_PLOT=False, SAVE=False, show_figure=True,
             mask=None,do_measurements='',compute_A=False,
             add_save_name='',logger=None,verbose=0):
    """
    Main function that perform other function calls responsible for
    all relevant calculations on the images.
    """
    cell_size = get_cell_size(imagename)
    """
    One beam area is one element resolution, so we avoind finding sub-components
    that are smaller than the beam area.
    """
    try:
        BEAM_AREA = beam_area2(imagename)
    except:
        BEAM_AREA = 1.0
    try:
        min_seg_pixels = beam_area2(imagename)
    except:
        min_seg_pixels = 15.0

    if mask is not None:
        mask = mask
        apply_mask = False
        omask = None
        if verbose >= 1:
            if logger is not None:
                logger.info(f"  >> Using provided mask.")
            else:
                print('     >> INFO: Using provided mask.')

    if apply_mask == True:
        # if logger is not None:
        #     logger.info(f"  CALC >> Performing mask dilation.")
        # else:
        #     print('     >> CALC: Performing mask dilation.')
        original_mask, mask_dilated = mask_dilation(imagename, cell_size=cell_size,
                                        sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations, rms=rms,
                                        PLOT=do_PLOT,verbose=verbose)
        mask = mask_dilated
        omask = original_mask
    # else:
    #     print('     >> WARN: Not using any mask.')
    #     mask = None
    # """
    # Basic background estimation.
    # """
    # # bkg_ = sep_background(crop_image,apply_mask=False,mask=mask,
    # #                       bw=11, bh=11, fw=12, fh=12)
    # bkg_ = sep_background(imagename,apply_mask=True,mask=None,
    #                       bw=11, bh=11, fw=12, fh=12)
    # bkg_to_sub = bkg_.back()

    if data_2D is not None:
        data_2D = data_2D
    else:
        data_2D = load_fits_data(imagename)

    # if verbose >= 1:
    #     if logger is not None:
    #         logger.info(f"  CALC >> Performing level statistics.")
    #     else:
    #         print('     >> CALC: Performing level statistics.')

    # results_final = level_statistics(img=imagename, cell_size=cell_size,
    #                                 mask_component=mask_component,
    #                                 mask=mask, apply_mask=False,
    #                                 data_2D=data_2D,
    #                                 sigma=sigma_mask, do_PLOT=do_PLOT,
    #                                 results=results_final, bkg_to_sub=bkg_to_sub,
    #                                 show_figure=False,
    #                                 rms=rms,
    #                                 add_save_name=add_save_name,
    #                                 SAVE=SAVE, ext='.jpg')
    if verbose >=1:
        if logger is not None:
            logger.info(f"  CALC >> Computing image properties.")
        else:
            print('     >> CALC: Computing image properties.')

    levels, fluxes, Lgrow, Lgrow_err, Lgrow_norm, Lgrow_err_norm, radii, agrow, \
        omask2, mask2, results_final = compute_image_properties(imagename,
                                                        residual=residualname,
                                                        cell_size=cell_size,
                                                        redshift=z,
                                                        mask_component=mask_component,
                                                        last_level=last_level,
                                                        sigma_mask=sigma_mask,
                                                        apply_mask=False,
                                                        crop=crop,box_size=box_size,
                                                        mask=mask,
                                                        rms=rms,
                                                        data_2D=data_2D,
                                                        data_res=data_res,
                                                        vmin_factor=vmin_factor,
                                                        # results=results_final,
                                                        show_figure=show_figure,
                                                        bkg_to_sub=bkg_to_sub,
                                                        verbose=verbose,
                                                        flux_units=flux_units,
                                                        flux_conversion_factor=flux_conversion_factor,
                                                        add_save_name=add_save_name,
                                                        SAVE=SAVE,logger=logger)
    # r_list, area_arr, area_beam, p, \
    #     flux_arr, results_final = do_petrofit(imagename, cell_size,
    #                                           mask_component=mask_component,
    #                                           PLOT=do_PLOT,
    #                                           sigma_mask=sigma_mask,
    #                                           dilation_size=dilation_size,
    #                                           npixels=npixels, fwhm=fwhm,
    #                                           kernel_size=kernel_size,
    #                                           results=results_final,
    #                                           apply_mask=apply_mask)
    omask = omask2.copy()
    error_petro = False
    if do_petro == True:
        # try:
        if verbose >= 1:
            if logger is not None:
                logger.info(f"++>> Computing Petrosian properties.")
            else:
                print('++>> Computing Petrosian properties.')
        _residual_for_error = None
        if use_residual_as_error and error_map is None and rms_map is None \
                and invvar_map is None and variance_map is None and weight_map is None:
            if data_res is not None:
                _residual_for_error = data_res
            elif residualname is not None:
                _residual_for_error = load_fits_data(residualname)

        error_arr_input, _error_source = resolve_flux_error_map(
            data_2D, error_map=error_map, rms_map=rms_map, invvar_map=invvar_map,
            variance_map=variance_map, weight_map=weight_map,
            residual_map=_residual_for_error, rms=None, verbose=(verbose > 0))

        r_list, area_arr, area_beam, p, flux_arr, error_arr, results_final, cat, \
            sorted_idx_list, segm, segm_deblend = \
            compute_petrosian_properties(data_2D, imagename,
                                            mask_component=mask_component,
                                            global_mask=mask,
                                            source_props=results_final,
                                            apply_mask=False,
                                            error=error_arr_input,
                                            sigma_level=sigma_mask,
                                            bkg_sub=bkg_sub, bkg_to_sub=bkg_to_sub,
                                            vmin=vmin_factor, plot=do_PLOT,
                                            deblend=deblend,
                                            fwhm=fwhm, kernel_size=kernel_size,
                                            show_figure=show_figure,
                                            verbose = verbose,
                                            add_save_name=add_save_name,
                                            npixels=npixels,logger=logger)
        error_petro = False
        # except Exception as e:
        #     if logger is not None:
        #         logger.warning(f"-->> ERROR when computing Petrosian properties. "
        #                         f"Will flag error_petro as True.")
        #     else:
        #         print("-->> ERROR when computing Petrosian properties. Will "
        #                 "flag error_petro as True.")
        #     error_petro = True
    else:
        error_petro = True

    results_final['error_petro'] = error_petro
    if mask_component is not None:
        centre = nd.maximum_position(data_2D * mask_component)[::-1]
        r, ir = get_profile(load_fits_data(imagename) * mask_component)
    else:
        centre = nd.maximum_position(data_2D)[::-1]
        r, ir = get_profile(data_2D,center=centre)
    #     rpix = r / cell_size
    rpix = r.copy()

    if z is not None:
        pix_to_pc = pixsize_to_pc(z=z, cell_size=cell_size)
        results_final['C95radii_pc'] = results_final['C95radii'] * pix_to_pc
        results_final['C50radii_pc'] = results_final['C50radii'] * pix_to_pc
        results_final['C95radii_asec'] = results_final['C95radii'] * cell_size
        results_final['C50radii_asec'] = results_final['C50radii'] * cell_size
        
        results_final['C95radii_pc_err'] = results_final['C95radii_err'] * pix_to_pc
        results_final['C50radii_pc_err'] = results_final['C50radii_err'] * pix_to_pc
        results_final['C95radii_asec_err'] = results_final['C95radii_err'] * cell_size
        results_final['C50radii_asec_err'] = results_final['C50radii_err'] * cell_size
        
        if error_petro == False:
            r_list_arcsec = r_list * cell_size
            pix_to_pc = pixsize_to_pc(z=z, cell_size=cell_size)
            r_list_pc = r_list * pix_to_pc
            Rp_arcsec = results_final['Rp'] * cell_size
            R50_arcsec = results_final['R50'] * cell_size
            results_final['Rp_pc'] = results_final['Rp'] * pix_to_pc
            results_final['R50_pc'] = results_final['R50'] * pix_to_pc
            results_final['Rp_asec'] = results_final['Rp'] * cell_size
            results_final['R50_asec'] = results_final['R50'] * cell_size
                        
        if error_petro == True:
            r_list_arcsec = rpix * cell_size
            pix_to_pc = pixsize_to_pc(z=z, cell_size=cell_size)
            r_list_pc = rpix * pix_to_pc
            results_final['Rp_pc'] = results_final['C95radii'] * pix_to_pc
            results_final['R50_pc'] = results_final['C50radii'] * pix_to_pc
            results_final['Rp_asec'] = results_final['C95radii'] * cell_size
            results_final['R50_asec'] = results_final['C50radii'] * cell_size


    
    if do_measurements=='all':
        try:
            x0c, y0c = results_final['x0m'], results_final['y0m']
            if compute_A == True:
                print('--==>> Computing asymetries...')
                results_final = compute_asymetries(imagename=imagename,
                                                mask=mask,
                                                bkg_to_sub=bkg_to_sub,
                                                mask_component=mask_component,
                                                centre=(x0c, y0c),
                                                results=results_final)

            # idx_R50 = np.where(flux_arr < 0.5 * results_final['flux_rp'])[0][-1]
            # idx_Rp = np.where(r_list < 2 * results_final['Rp'])[0][-1]
            # idx_Cradii = np.where(r_list < results_final['Cradii'])[0][-1]
            # idx_C50radii = np.where(r_list < results_final['C50radii'])[0][-1]

            # idx_R50 = np.where(flux_arr < 0.5 * results_final['flux_rp'])[0][-1]
            # flux_arr[idx_R50], r_list[idx_R50], results_final['R50']
            # if mask_component is None:
            # if verbose >= 1:
                # print('--==>> Computing image statistics...')
            # results_final = get_image_statistics(imagename=imagename,
            #                                     mask_component=mask_component,
            #                                     mask=mask,
            #                                     sigma_mask=sigma_mask,
            #                                     apply_mask=False,
            #                                     residual_name=residualname, fracX=fracX,
            #                                     fracY=fracY,
            #                                     dic_data=results_final,
            #                                     cell_size=cell_size)
        except:
            pass
    

    results_final = cosmo_stats(imagename=imagename, z=z, results=results_final)

    results_final['pix_to_pc'] = pix_to_pc
    results_final['cell_size'] = cell_size
    if error_petro == False:
        results_final['area_beam_2Rp'] = area_beam[int(results_final['Rp'])]
        results_final['area_beam_R50'] = area_beam[int(results_final['R50'])]
    if error_petro == True:
        results_final['area_beam_2Rp'] = results_final['A95']/BEAM_AREA
        results_final['area_beam_R50'] = results_final['A50']/BEAM_AREA

    df = pd.DataFrame.from_dict(results_final, orient='index').T
    df.to_csv(imagename.replace('.fits', add_save_name + '_stats.csv'),
              header=True,
              index=False)
    return (results_final, mask, omask)



def check_flux_duplicates(cumulative_flux, radii):
    """
    Check if the array of fluxes contain duplicated values.
    """
    cumulative_flux_norm = cumulative_flux / cumulative_flux[-1]
    
    unique_elements = np.unique(cumulative_flux_norm)
    has_duplicates = len(unique_elements) != len(cumulative_flux_norm)
    if has_duplicates:
        # raise ValueError("Cumulative flux contains duplicates.")
        # print("Cumulative flux contains duplicates.")
        # print(".... attempting removal of duplicates.")
        unique_flux, indices = np.unique(cumulative_flux_norm, return_inverse=True)
        average_radii = np.zeros_like(unique_flux, dtype=float)
        average_cumulative_flux = np.zeros_like(unique_flux, dtype=float)
        # Calculate average of y-values corresponding to each unique x-value
        for i in range(len(unique_flux)):
            average_radii[i] = np.mean(radii[indices == i])
            average_cumulative_flux[i] = np.mean(cumulative_flux[indices == i])
            
        cumulative_flux_norm = unique_flux.copy()
        radii = average_radii.copy()
        cumulative_flux = average_cumulative_flux.copy()
        
    return(cumulative_flux,radii)

# def find_fractional_radius(cumulative_flux, radii,
#                            fraction,
#                            rms=None, beam_area=None,
#                            do_plot=False):
#     """
#     Calculate the radius that encloses a given fraction of the total flux,
#     with data-driven uncertainty on the radius and enclosed area.

#     Parameters
#     ----------
#     cumulative_flux : numpy array
#         Cumulative flux densities as a function of radius.
#     radii : numpy array
#         Radii corresponding to the cumulative flux densities (pixels).
#     fraction : float
#         Fraction of total flux for which to calculate the enclosing radius.
#         Must be between 0 and 1.
#     rms : float, optional
#         Per-beam RMS noise of the image (same units as cumulative_flux).
#         Radio-specific. If provided together with beam_area, the noise-propagated
#         term is computed and compared against the profile scatter term.
#     beam_area : float, optional
#         Beam area in pixels. Radio-specific. Used together with rms.
#     do_plot : bool, optional
#         If True, plot the cumulative profile and fractional radius.

#     Returns
#     -------
#     tuple
#         radius_at_fraction : float
#             Radius enclosing the requested flux fraction (pixels).
#         fraction_flux_true : float
#             Absolute flux enclosed within radius_at_fraction.
#         area_at_fraction : float
#             Pixel area enclosed within radius_at_fraction.
#         fraction_radius_error_mean : float
#             1-sigma uncertainty on the fractional radius, derived from the
#             local scatter of the cumulative flux profile propagated via the
#             local slope dR/dF. Works for any image type. If rms and beam_area
#             are also provided, the thermal noise term is included and the
#             dominant contribution is used.
#         fraction_area_error_mean : float
#             1-sigma uncertainty on the enclosed area, propagated from
#             fraction_radius_error_mean via dA = 2*pi*R*dR.
#     """
#     from scipy.interpolate import interp1d
#     from astropy.stats import mad_std

#     if not (0 < fraction < 1):
#         raise ValueError("Fraction must be between 0 and 1.")

#     cumulative_flux, radii = check_flux_duplicates(cumulative_flux, radii)
#     cumulative_flux_norm = cumulative_flux / cumulative_flux[-1]

#     total_flux = cumulative_flux_norm[-1]
#     total_flux_true = cumulative_flux[-1]
#     fraction_flux = total_flux * fraction
#     fraction_flux_true = total_flux_true * fraction

#     spline = interp1d(cumulative_flux_norm, radii, kind='linear')
#     radius_at_fraction = float(spline(fraction_flux))
#     area_at_fraction = radii_to_area(radius_at_fraction)[0]

#     # ------------------------------------------------------------------
#     # Local slope dR/dF at the fraction point via finite difference
#     # Compact source -> steep profile -> small dR/dF -> small error
#     # Diffuse source -> flat profile -> large dR/dF -> large error
#     # ------------------------------------------------------------------
#     delta = 1e-3
#     f_lo = max(fraction - delta, float(cumulative_flux_norm[0]))
#     f_hi = min(fraction + delta, float(cumulative_flux_norm[-1]))
#     dR_dF = (float(spline(f_hi)) - float(spline(f_lo))) / (f_hi - f_lo)

#     # ------------------------------------------------------------------
#     # Data-driven sigma_F: local MAD of annular flux increments near the
#     # fraction point. Works for any image type (radio, optical, IR...).
#     # Captures profile irregularity from noise, PSF sidelobes, or sky background.
#     # ------------------------------------------------------------------
#     annular_fluxes_norm = np.diff(cumulative_flux_norm)
#     idx = int(np.searchsorted(cumulative_flux_norm, fraction))
#     half_win = max(3, len(annular_fluxes_norm) // 8)
#     lo = max(0, idx - half_win)
#     hi = min(len(annular_fluxes_norm), idx + half_win)

#     if (hi - lo) >= 3:
#         sigma_F_local = mad_std(annular_fluxes_norm[lo:hi])
#     else:
#         sigma_F_local = mad_std(annular_fluxes_norm)  # fallback to global scatter

#     fraction_radius_error_mean = 3 * abs(dR_dF) * sigma_F_local

#     # ------------------------------------------------------------------
#     # Radio-specific: add thermal noise term if rms and beam_area provided
#     # Take the dominant contribution between profile scatter and thermal noise
#     # ------------------------------------------------------------------
#     if rms is not None and beam_area is not None:
#         n_beams = max(area_at_fraction / beam_area, 1.0)
#         sigma_F_noise = (rms * np.sqrt(n_beams)) / total_flux_true
#         noise_radius_error = abs(dR_dF) * sigma_F_noise
#         fraction_radius_error_mean = max(fraction_radius_error_mean, noise_radius_error)

#     # Propagate radius uncertainty to area: dA = 2*pi*R*dR
#     fraction_area_error_mean = 2 * np.pi * radius_at_fraction * fraction_radius_error_mean

#     # ------------------------------------------------------------------
#     if do_plot:
#         plt.plot(radii, cumulative_flux_norm, 'o', label='Data Points')

#         flux_range = np.linspace(float(np.min(cumulative_flux_norm)),
#                                  float(np.max(cumulative_flux_norm)), 1000)
#         plt.plot(spline(flux_range), flux_range, '-', label='Spline Interpolation')

#         plt.plot(radius_at_fraction, fraction_flux, 'rx',
#                  markersize=10, label=f'Fractional Radius ({fraction*100:.1f}%)')
#         plt.axvline(radius_at_fraction, color='r', linestyle='--', label=f'R{fraction*100:.0f}')
#         plt.axvline(radius_at_fraction - fraction_radius_error_mean,
#                     color='gray', linestyle=':', label='$\pm 3\sigma$')
#         plt.axvline(radius_at_fraction + fraction_radius_error_mean,
#                     color='gray', linestyle=':')

#         plt.ylabel('Cumulative Flux')
#         plt.xlabel('Radius (pixels)')
#         plt.legend()
#         plt.show()

#     return (radius_at_fraction, fraction_flux_true, area_at_fraction,
#             fraction_radius_error_mean, fraction_area_error_mean)

def find_fractional_radius(cumulative_flux, radii,
                           fraction,
                           cumulative_flux_err=None,
                           do_plot=False):
    """
    Calculate the radius that encloses a given fraction of the total flux,
    with bins-independent uncertainty propagation.

    Parameters
    ----------
    cumulative_flux : numpy array
        Cumulative flux densities as a function of radius.
    radii : numpy array
        Radii corresponding to the cumulative flux densities (pixels).
    fraction : float
        Fraction of total flux. Must be between 0 and 1.
    cumulative_flux_err : numpy array, optional
        1-sigma uncertainty on the normalised cumulative flux at each radius,
        built as sqrt(cumsum(sigma_i^2)) / total_flux in compute_image_properties.
        If provided, uncertainty on the radius is propagated via the local
        slope dR/dF. This is bins-independent by construction.
        If not provided, falls back to a slope-aware +/- 5% fraction perturbation.
    do_plot : bool, optional
        Plot the cumulative profile and fractional radius with error bounds.

    Returns
    -------
    tuple
        radius_at_fraction, fraction_flux_true, area_at_fraction,
        fraction_radius_error_mean, fraction_area_error_mean
    """
    from scipy.interpolate import interp1d
    from astropy.stats import mad_std

    if not (0 < fraction < 1):
        raise ValueError("Fraction must be between 0 and 1.")

    cumulative_flux, radii = check_flux_duplicates(cumulative_flux, radii)
    cumulative_flux_norm = cumulative_flux / cumulative_flux[-1]

    total_flux = cumulative_flux_norm[-1]
    total_flux_true = cumulative_flux[-1]
    fraction_flux = total_flux * fraction
    fraction_flux_true = total_flux_true * fraction

    spline = interp1d(cumulative_flux_norm, radii, kind='linear')
    radius_at_fraction = float(spline(fraction_flux))
    area_at_fraction = radii_to_area(radius_at_fraction)[0]

    # ------------------------------------------------------------------
    # Local slope dR/dF at the fraction point via finite difference
    # ------------------------------------------------------------------
    delta = 1e-3
    f_lo = max(fraction - delta, float(cumulative_flux_norm[0]))
    f_hi = min(fraction + delta, float(cumulative_flux_norm[-1]))
    dR_dF = (float(spline(f_hi)) - float(spline(f_lo))) / (f_hi - f_lo)

    # ------------------------------------------------------------------
    # Uncertainty propagation
    # ------------------------------------------------------------------
    if cumulative_flux_err is not None:
        # Interpolate the normalised cumulative error at the fraction point.
        # cumulative_flux_err is already normalised (/ total_flux) and
        # bins-independent - it depends only on total area within R.
        err_spline = interp1d(cumulative_flux_norm, cumulative_flux_err,
                              kind='linear')
        sigma_F = float(err_spline(fraction_flux))
        fraction_radius_error_mean = abs(dR_dF) * sigma_F
    else:
        # Fallback: slope-aware +/-5% perturbation (no noise information)
        # if fraction > 0.95:
        #     dr_err_up = 0.025
        # else:
        #     dr_err_up = 0.05
        # if fraction <= 0.1:
        #     dr_err_lo = 0.025
        # else:
        #     dr_err_lo = 0.05
        dr_err_up = 0.025 if fraction > 0.95 else 0.05
        dr_err_lo = 0.025 if fraction <= 0.1 else 0.05

        lower_fraction = max(fraction - dr_err_lo, float(cumulative_flux_norm[0]))
        upper_fraction = min(fraction + dr_err_up, float(cumulative_flux_norm[-1]))
        lower_radius = float(spline(lower_fraction))
        upper_radius = float(spline(upper_fraction))
        fraction_radius_error_mean = 1*np.mean([abs(upper_radius - radius_at_fraction),
                                              abs(lower_radius - radius_at_fraction)])

    # # Pixelisation floor: cannot resolve a radius change smaller than 0.5 px.
    # fraction_radius_error_mean = max(fraction_radius_error_mean, 0.5)

    # # Propagate to area: dA = 2*pi*R*dR
    # fraction_area_error_mean = 2 * np.pi * radius_at_fraction * fraction_radius_error_mean
    
    # ------------------------------------------------------------------
    # Pixelisation floor only: a radius cannot be localised to better than
    # half a pixel.
    #
    # REMOVED (2026-08-27): a `local_dr` floor at half the local radial step
    # of the growth curve, commented as a resolution floor "reflecting beam
    # size for radio, PSF for optical". It reflects neither. `radii` is built
    # from `levels = _build_levels(peak_g, floor_g, bins)`, so the step is set
    # by `bins` (a function argument, default 64), not by any instrument
    # property. For a circular exponential it is exactly
    # h*ln(peak/floor)/(bins-1): doubling `bins` halved every radius error,
    # and the term is essentially fraction-independent, so it contributed no
    # fraction ordering -- it only flattened the errors onto a constant.
    # The spline genuinely interpolates between grid points, so the grid
    # spacing is an interpolation artefact, not a measurement limit. Where the
    # curve really is too coarse to trust, the honest signals are the
    # `growth_curve_n_levels` and `peak_to_floor_dynamic_range` diagnostics
    # exported by compute_image_properties, not an inflated error bar.
    # ------------------------------------------------------------------
    fraction_radius_error_mean = max(fraction_radius_error_mean, 0.5)

    # Propagate to area: dA = 2*pi*R*dR
    fraction_area_error_mean = 2 * np.pi * radius_at_fraction * fraction_radius_error_mean

    # ------------------------------------------------------------------
    if do_plot:
        plt.plot(radii, cumulative_flux_norm, 'o', label='Data Points')
        # plt.errorbar(radii, cumulative_flux_norm, yerr=10*cumulative_flux_err, fmt='none',
        #              ecolor='orange', alpha=0.7, label='Cumulative Flux Error')
        plt.fill_between(radii, cumulative_flux_norm - 3*cumulative_flux_err,
                         cumulative_flux_norm + 3*cumulative_flux_err, color='grey',
                         alpha=0.5, 
                         label='Cumulative Flux Uncertainty'
                         )
        flux_range = np.linspace(float(np.min(cumulative_flux_norm)),
                                 float(np.max(cumulative_flux_norm)), 1000)
        plt.plot(spline(flux_range), flux_range, '-', label='Spline Interpolation')
        plt.plot(radius_at_fraction, fraction_flux, 'rx',
                 markersize=10, label=f'R{fraction*100:.0f}%')
        plt.axvline(radius_at_fraction, color='r', linestyle='--')
        plt.axvline(radius_at_fraction - fraction_radius_error_mean,
                    color='gray', linestyle=':', label='$\pm 1\sigma$')
        plt.axvline(radius_at_fraction + fraction_radius_error_mean,
                    color='gray', linestyle=':')
        plt.ylabel('Cumulative Flux')
        plt.xlabel('Radius (pixels)')
        plt.legend()
        plt.show()

    return (radius_at_fraction, fraction_flux_true, area_at_fraction,
            fraction_radius_error_mean, fraction_area_error_mean)

def find_fractional_radius_bk(cumulative_flux, radii,
                           fraction,
                           rms=None, beam_area=None,
                           do_plot=False):
    """
    Calculate the radius that encloses a given fraction of the total flux,
    with noise-propagated uncertainty on the radius and enclosed area.

    Parameters
    ----------
    cumulative_flux : numpy array
        Cumulative flux densities as a function of radius.
    radii : numpy array
        Radii corresponding to the cumulative flux densities (pixels).
    fraction : float
        Fraction of total flux for which to calculate the enclosing radius.
        Must be between 0 and 1.
    rms : float, optional
        Per-beam RMS noise of the image (same units as cumulative_flux).
        If provided together with beam_area, enables noise-propagated uncertainty.
    beam_area : float, optional
        Beam area in pixels. Used to count independent beams within the
        aperture for noise propagation.
    do_plot : bool, optional
        If True, plot the cumulative profile and fractional radius.

    Returns
    -------
    tuple
        radius_at_fraction : float
            Radius enclosing the requested flux fraction (pixels).
        fraction_flux_true : float
            Absolute flux enclosed within radius_at_fraction.
        area_at_fraction : float
            Pixel area enclosed within radius_at_fraction.
        fraction_radius_error_mean : float
            1-sigma uncertainty on the fractional radius.
            If rms and beam_area are provided, derived from noise propagation
            via the local slope dR/dF of the cumulative profile.
            Otherwise falls back to a slope-aware +/-5% fraction perturbation.
        fraction_area_error_mean : float
            1-sigma uncertainty on the enclosed area, propagated from
            fraction_radius_error_mean.
    """
    from scipy.interpolate import interp1d

    if not (0 < fraction < 1):
        raise ValueError("Fraction must be between 0 and 1.")

    cumulative_flux, radii = check_flux_duplicates(cumulative_flux, radii)
    cumulative_flux_norm = cumulative_flux / cumulative_flux[-1]

    total_flux = cumulative_flux_norm[-1]
    total_flux_true = cumulative_flux[-1]
    fraction_flux = total_flux * fraction
    fraction_flux_true = total_flux_true * fraction

    # Use linear interpolation to find radius at the requested fraction
    spline = interp1d(cumulative_flux_norm, radii, kind='linear')

    radius_at_fraction = spline(fraction_flux)
    area_at_fraction = radii_to_area(radius_at_fraction)[0]

    # ------------------------------------------------------------------
    # Uncertainty estimation
    # ------------------------------------------------------------------
    # Local slope dR/dF at the fraction point via finite difference.
    # This is the key: a steep profile (compact source) gives small dR/dF
    # and thus small radius uncertainty; a flat diffuse profile gives large dR/dF.
    delta = 1e-3 # * total_flux  # perturbation in flux space for finite difference
    f_lo = max(fraction - delta, cumulative_flux_norm[0])
    f_hi = min(fraction + delta, cumulative_flux_norm[-1])
    dR_dF = (float(spline(f_hi)) - float(spline(f_lo))) / (f_hi - f_lo)

    if rms is not None and beam_area is not None:
        n_beams = max(area_at_fraction / beam_area, 1.0)
        sigma_F_norm = (rms * np.sqrt(n_beams)) / total_flux_true
        noise_radius_error = abs(dR_dF) * sigma_F_norm

        # slope-perturbation term (profile shape / mask boundary uncertainty)
        dr_err_up = 0.025 if fraction > 0.95 else 0.05
        dr_err_lo = 0.025 if fraction <= 0.1 else 0.05
        lower_fraction = max(fraction - dr_err_lo, float(cumulative_flux_norm[0]))
        upper_fraction = min(fraction + dr_err_up, float(cumulative_flux_norm[-1]))
        slope_radius_error = np.mean([abs(float(spline(upper_fraction)) - radius_at_fraction),
                                      abs(float(spline(lower_fraction)) - radius_at_fraction)])

        # take the dominant term
        fraction_radius_error_mean = max(noise_radius_error, slope_radius_error)
    else:
        # Fallback: slope-aware perturbation of +/-5% in fraction space.
        # Still uses the local slope so compact vs diffuse sources
        # get different errors, unlike the old hardcoded approach.
        # if fraction > 0.95:
        #     dr_err_up = 0.025
        # else:
        #     dr_err_up = 0.05
        # if fraction <= 0.1:
        #     dr_err_lo = 0.025
        # else:
        #     dr_err_lo = 0.05
        dr_err_up = 0.025 if fraction > 0.95 else 0.05
        dr_err_lo = 0.025 if fraction <= 0.1 else 0.05

        lower_fraction = max(fraction - dr_err_lo, float(cumulative_flux_norm[0]))
        upper_fraction = min(fraction + dr_err_up, float(cumulative_flux_norm[-1]))

        lower_radius = float(spline(lower_fraction))
        upper_radius = float(spline(upper_fraction))

        radius_error_range = (abs(upper_radius - radius_at_fraction),
                              abs(lower_radius - radius_at_fraction))
        fraction_radius_error_mean = np.mean(radius_error_range)

    # Propagate radius uncertainty to area: A = pi*R^2, so dA = 2*pi*R * dR
    fraction_area_error_mean = 2 * np.pi * radius_at_fraction * fraction_radius_error_mean

    # ------------------------------------------------------------------
    if do_plot:
        plt.plot(radii, cumulative_flux_norm, 'o', label='Data Points')

        flux_range = np.linspace(np.min(cumulative_flux_norm),
                                 np.max(cumulative_flux_norm), 1000)
        plt.plot(spline(flux_range), flux_range, '-', label='Spline Interpolation')

        plt.plot(radius_at_fraction, fraction_flux, 'rx',
                 markersize=10, label=f'Fractional Radius ({fraction*100:.1f}%)')
        plt.axvline(radius_at_fraction, color='r', linestyle='--')
        plt.axvline(radius_at_fraction - fraction_radius_error_mean,
                    color='gray', linestyle=':')
        plt.axvline(radius_at_fraction + fraction_radius_error_mean,
                    color='gray', linestyle=':', label='+/- 1 sigma radius')

        plt.ylabel('Cumulative Flux')
        plt.xlabel('Radius')
        plt.legend()
        plt.show()

    return (radius_at_fraction, fraction_flux_true, area_at_fraction,
            fraction_radius_error_mean, fraction_area_error_mean)


# def find_fractional_radius(cumulative_flux, radii, 
#                            fraction,do_plot=False):
#     """
#     Calculate the radius that encloses a given fraction of the total flux.
    
#     Parameters:
#     - cumulative_flux: numpy array of cumulative flux densities
#     - radii: numpy array of radii corresponding to the cumulative flux densities
#     - fraction: the fraction of the total flux for which to calculate the radius
    
#     Returns:
#     - radius_at_fraction: the radius that encloses the specified fraction of the total flux
#     """
#     from scipy.interpolate import interp1d
#     # Ensure the fraction is between 0 and 1
#     if not (0 < fraction < 1):
#         raise ValueError("Fraction must be between 0 and 1.")

#     cumulative_flux,radii = check_flux_duplicates(cumulative_flux,radii)
#     cumulative_flux_norm = cumulative_flux / cumulative_flux[-1]
#     # Calculate the total flux to find the target flux for the given fraction
#     total_flux = cumulative_flux_norm[-1]
#     total_flux_true= cumulative_flux[-1]
#     fraction_flux = total_flux * fraction
#     fraction_flux_true = total_flux_true * fraction
    
#     # Use spline interpolation to create a smooth function of cumulative flux vs radius
#     spline = interp1d(cumulative_flux_norm, radii, kind='linear')
    
#     # Use the interpolated function to find the radius for the target flux
#     radius_at_fraction = spline(fraction_flux)

#     area_at_fraction = radii_to_area(radius_at_fraction)[0]

#     # if fraction >=0.95:
#     #     dr_err_up = 0.025
#     # else:
#     #     dr_err_up = 0.05
#     # if fraction <= 0.1:
#     #     dr_err_lo = 0.025
#     # else:
#     #     dr_err_lo = 0.05

#     if fraction > 0.95:
#         dr_err_up = 0.025
#     else:
#         dr_err_up = 0.05
#     if fraction <= 0.1:
#         dr_err_lo = 0.025
#     else:
#         dr_err_lo = 0.05

    
#     # Calculate the lower and upper fractions for the error range
#     lower_fraction = max(fraction - dr_err_lo, 0)
#     upper_fraction = min(fraction + dr_err_up, 1)

#     # Calculate the flux values for the lower and upper fractions
#     lower_fraction_flux = lower_fraction * total_flux
#     upper_fraction_flux = upper_fraction * total_flux

#     # Calculate the radii at the lower and upper fractions using the spline interpolation
#     lower_radius = spline(lower_fraction_flux)
#     upper_radius = spline(upper_fraction_flux)

#     # Calculate the error range as the difference between the lower and upper radii
#     radius_error_range = (abs(upper_radius-radius_at_fraction),abs(lower_radius-radius_at_fraction))
#     # print(f"Fractional Radius with error range: {radius_at_fraction:.2f} +/- {radius_error_range[1]:.2f}/{radius_error_range[0]:.2f}")
#     fraction_radius_error_mean = np.mean(radius_error_range)
#     # Calculate the areas at the lower and upper radii
#     lower_area = radii_to_area(lower_radius)[0]
#     upper_area = radii_to_area(upper_radius)[0]

#     # Calculate the area error range as the difference between the lower and upper areas
#     area_error_range = (abs(upper_area - area_at_fraction), abs(lower_area - area_at_fraction))
#     # print(f"Fractional Area with error range: {area_at_fraction:.2f} +/- {area_error_range[1]:.2f}/{area_error_range[0]:.2f}")
#     # print(area_error_range[0],area_error_range[1],np.mean(area_error_range))
#     fraction_area_error_mean = np.mean(area_error_range)
#     if do_plot:
#         # Plot the data points
#         plt.plot(radii, cumulative_flux_norm, 'o', label='Data Points')
        
#         # Create a smooth range of flux values for plotting the interpolation
#         flux_range = np.linspace(np.min(cumulative_flux_norm), np.max(cumulative_flux_norm), 1000)
        
#         # Plot the spline interpolation
#         plt.plot(spline(flux_range), flux_range, '-', label='Spline Interpolation')
        
#         # Plot a cross at the fractional radius
#         plt.plot(radius_at_fraction, fraction_flux, 'rx', markersize=10, label=f'Fractional Radius ({fraction*100:.1f}%)')
#         plt.axvline(lower_radius)
#         plt.axvline(radius_at_fraction)
#         plt.axvline(upper_radius)
        
#         # Add labels and legend
#         plt.ylabel('Cumulative Flux')
#         plt.xlabel('Radius')
#         plt.legend()
        
#         # Show the plot
#         plt.show()
    
#     return(radius_at_fraction, fraction_flux_true, area_at_fraction, 
#            fraction_radius_error_mean, fraction_area_error_mean)

def R50_to_fwhm(R50,R50_err,q_ratio=1,scale = 1):
    """
    Calculate the beam size in arcseconds for a given R50 value.
    (Lucatelli et al., 2024)
    """
    theta_maj = 2*R50*scale
    theta_min = theta_maj*q_ratio
    
    theta_maj_err = 2*R50_err*scale
    theta_min_err = theta_maj_err*q_ratio
    
    theta = (theta_maj*theta_min)**0.5
    theta_err = (theta_maj_err*theta_min_err)**0.5
    return(theta_maj,theta_maj_err,
           theta_min,theta_min_err,
           theta,theta_err)


def concentration_index_errors(C80radii, C20radii, 
                               C90radii, C50radii, 
                               C80radii_err, C20radii_err, 
                               C90radii_err, C50radii_err):
    """
    Given fractional radii and their errors, calculate the concentration indices
    and their errors.
    
    Parameters:
    - C80radii: the radius enclosing 80% of the total flux
    - C20radii: the radius enclosing 20% of the total flux
    - C90radii: the radius enclosing 90% of the total flux
    - C50radii: the radius enclosing 50% of the total flux
    - C80radii_err: the error in the radius enclosing 80% of the total flux
    - C20radii_err: the error in the radius enclosing 20% of the total flux
    - C90radii_err: the error in the radius enclosing 90% of the total flux
    - C50radii_err: the error in the radius enclosing 50% of the total flux
    
    Returns:
    - C1: the concentration index C1
    - C1_err: the error in the concentration index C1
    - C2: the concentration index C2
    - C2_err: the error in the concentration index C2
    """
    C1 = np.log10(C80radii / C20radii)
    C2 = np.log10(C90radii / C50radii)

    dC1_dC80radii = 1 / (C80radii * np.log(10))
    dC1_dC20radii = -1 / (C20radii * np.log(10))
    C1_err = np.sqrt((dC1_dC80radii * C80radii_err)**2 + (dC1_dC20radii * C20radii_err)**2)

    dC2_dC90radii = 1 / (C90radii * np.log(10))
    dC2_dC50radii = -1 / (C50radii * np.log(10))
    C2_err = np.sqrt((dC2_dC90radii * C90radii_err)**2 + (dC2_dC50radii * C50radii_err)**2)
    # print(f"C1 = {C1:.2f} +/- {C1_err:.2f}")
    # print(f"C2 = {C2:.2f} +/- {C2_err:.2f}")
    return C1, C1_err, C2, C2_err


# def deprecated(new_func_name):
#     def decorator(old_func):
#         def wrapper(*args, **kwargs):
#             warnings.warn(f"'{old_func.__name__}' is deprecated and will be "
#                           f"removed in a future version. "
#                           f"Use '{new_func_name}' instead.",
#                           category=DeprecationWarning, stacklevel=2)
#             return old_func(*args, **kwargs)
#         return wrapper
#     return decorator



def compute_image_properties_old(img, residual,
                             cell_size=None, mask_component=None,
                             redshift = None,
                             aspect=1, last_level=3.0, mask=None, 
                             bins = 64,
                             data_2D=None,data_res=None,
                             do_petro=False,
                             dilation_size=None, iterations=2,
                             dilation_type='disk',do_fit_ellipse=False,
                             sigma_mask=6, rms=None, 
                             systematic_error_fraction = 0.05,
                             results=None, bkg_to_sub=None,
                             apply_mask=True, 
                             vmax = None, vmin_factor=3, vmax_factor=0.5,
                             crop=False, box_size=256,
                             SAVE=True, add_save_name='', plot_savemane=None,
                             show_figure=True,
                             figsize=(6, 5), image_title = None,
                             save_csv = False,
                             flux_units = 'Jy/beam',
                             ext='.jpg',logger=None,verbose=0):
    """
    Params
    ------
    img : str
        The path to the image to be analyzed.
    residual : str
        The path to the residual image associated with the image.
    cell_size : float, optional
        The default is None. The size of the pixel in arcseconds.
        If None, `get_cell_size` will be used to get from the header.
    mask_component : 2D np array, optional
        The default is None. If not None, this is the mask for a specific component
        when performing a multi-component source analysis.
    aspect : float, optional (experimental)
        The default is 1. The aspect ratio of the image for plotting.
    last_level : float, optional
        New threshold level (as multiple of sigma_mad) to be used inside the existing mask.
    bins : int, optional
        Number of bins to slice intensity levels of the image. 
    mask : 2D np array, optional
        The default is None. If not None, the function `mask_dilation` will determine it with
        default parameters.
    data_2D : 2D np array, optional
        The default is None. This can be used to pass a 2D array directly to the function.
        For example, when calculating properties from an array without reading from a file,
        e.g. a model image, you can use this. But, to obtain meaningful physical units,
        you must provide the corresponding image file where this array was derived from.
    dilation_size : int, optional
        The default is None. The size of the dilation to be used in the mask dilation.
        If None, the default value will be the size of the restoring beam.
    iterations : int, optional
        The default is 2. The number of iterations to be used in the mask dilation.
        If signs of over-dilation are present, you can set to 1.
    dilation_type : str, optional
        The default is 'disk'. The type of dilation to be used in the mask dilation.
    sigma_mask : float, optional
        The default is 6. The sigma level to be used in the mask dilation.
    rms : float, optional
        The default is None. The rms value to be used in the mask dilation.
        If None, the function `mad_std` will be used to calculate it from the residual image,
        if provided. If the residual image is not provided, the function will use the image itself.
        But, in that case, the result may not be accurate (overestimated) if the image size is
        comparable in size to the size of the source structure.
    results : dict, optional
        The default is None. A dictionary to store the results.
        You can pass an existing external dictionary, so the results will be appended to it.
    bkg_to_sub : 2D np array, optional (EXPERIMENTAL)
        The default is None. The background to be subtracted from the image.
    apply_mask : bool, optional
        The default is True. If True, the mask dilation will be calcualted from the image.
    vmin_factor : float, optional
        The default is 3. The factor (as a multiple of sigma_mad) to be used in the vmin
        calculation for the image plot.
    vmax_factor : float, optional
        The default is 0.5. The factor (as a multiple of peak brightness) to be used in the vmax
        calculation for the image plot.
    crop : bool, optional
        The default is False. If True, the image will be cropped to a box_size.
    box_size : int, optional
        The default is 256. The size of the box to be used in the image cropping.
    SAVE : bool, optional
        The default is True. If True, the image plot will be saved to a file.
    add_save_name : str, optional
        The default is ''. A string to be added to the image plot file name.
    show_figure : bool, optional
        The default is True. If True, the image plot will be shown.
    ext : str, optional
        The default is '.jpg'. The file extension to be used in the image plot file name.
    logger : logging.Logger, optional
        The default is None. A logger object to be used to log messages.

    """
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(img)

    if flux_units == 'Jy/beam':
        flux_conversion = 1e3 # Jy/beam to mJy/beam
    if flux_units == 'nanomaggies':
        flux_conversion = 3.631e-3 # nanomaggies to mJy
    if flux_units == '1nJy':
        flux_conversion = 1e-6  # nJy to mJy
    if flux_units == '10nJy':
        flux_conversion = 1e-5  # nJy to mJy
    if flux_units == 'mJy/beam':
        flux_conversion = 1.0

    from skimage.draw import disk
    if data_2D is not None:
        g_ = data_2D
    else:
        g_ = load_fits_data(img)
    if residual is not None:
        if isinstance(residual, str) == True:
            # res_ = load_fits_data(residual)
            res_ = load_fits_data(residual)
        else:
            res_ = residual
        res = np.nan_to_num(res_.copy(),nan=0)
    
    g = np.nan_to_num(g_.copy(),nan=0)
    g_original = g.copy()
    

    if bkg_to_sub is not None:
        g = g - bkg_to_sub


    scale_units = 'arcsec'
    if cell_size is None:
        try:
            cell_size = get_cell_size(img)
        except:
            scale_units = 'px'
            cell_size = 1.0

    try:
        beam_area_ = beam_area2(img, cellsize=cell_size)
        flag_freq = False
        g_hd = imhead(img)
        try:
            freq = g_hd['refval'][2] / 1e9
        except:
            freq = getfreqs([img])[0] / 1e9
        # print(freq)
        # omaj = g_hd['restoringbeam']['major']['value']
        # omin = g_hd['restoringbeam']['minor']['value']
        omaj,omin, BPA, _, _ = beam_shape(img)
        # beam_area_ = beam_area(omaj, omin, cellsize=cell_size)
    except Exception as e:
        if logger is not None:
            logger.warning(f"  -->> ERROR: {e}")
            logger.warning(f"  -->> Not a radio image?")
        else:
            print(f"     -->> ERROR: {e}")
            print(f"     -->> Not a radio image?")
        flag_freq = True
        freq = 1e15 / 1e9
        omaj = 1
        omin = 1
        beam_area_ = 1

    if rms is not None:
        std = rms
    else:
        if residual is not None:
            if isinstance(residual, str) == True:
                std = mad_std(np.nan_to_num(load_fits_data(residual), nan=0))
            else:
                std = mad_std(np.nan_to_num(residual, nan=0))
        else:
            std = mad_std(np.nan_to_num(g_.copy(),nan=0))
    if mask is not None:
        mask = mask
        omask = mask
        g = g * mask  # *(g_>3*mad_std(g_)) + 0.1*mad_std(g_)
        # res = res * mask  # *(g_>3*mad_std(g_)) + 0.1*mad_std(g_)
        # total_flux = np.sum(g_ * (g_ > 3 * std)) / beam_area_
        total_flux = np.sum(g) / beam_area_
        total_flux_nomask = np.sum(g_) / beam_area_
        apply_mask = False

    if apply_mask == True:
        if data_res is None:
            omask, mask = mask_dilation(img, sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations,
                                        dilation_type=dilation_type,
                                        show_figure=False,
                                        # show_figure=show_figure,
                                        verbose=verbose)
        else:
            omask, mask = mask_dilation(data_2D, sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations,
                                        dilation_type=dilation_type,
                                        show_figure=False,
                                        # show_figure=show_figure,
                                        verbose=verbose)
    

        # eimshow(mask)
        g = g * mask  # *(g_>3*mad_std(g_)) + 0.1*mad_std(g_)
        # res = res * mask
        # total_flux = np.sum(g_ * (g_ > 3 * std)) / beam_area_
        total_flux = np.nansum(g) / beam_area_
        total_flux_nomask = np.nansum(g_) / beam_area_
        # g = g_

    if mask_component is not None:
        #
        g = g * mask_component

        # res = res * mask_component

        total_flux_nomask = np.nansum(g) / beam_area_
        # total_flux = np.sum(g * (g > 3 * std)) / beam_area_
        total_flux = np.nansum(g) / beam_area_
        if mask is not None:
            mask = mask * mask_component
            omask = mask * mask_component
        else:
            mask = mask_component
            omask = mask_component

    # if (apply_mask is None) and  (mask is None):
    #     mask = mask_component
        # g = g_
    if (mask_component is None) and (apply_mask == False) and (mask is None):
        # g = g_
        total_flux = np.nansum(g * (g > 3 * std)) / beam_area_
        total_flux_nomask = np.nansum(g) / beam_area_

        unity_mask = np.ones(g.shape) == 1
        omask, mask = unity_mask, unity_mask

    # total_flux_nomask = np.sum(g_) / beam_area_
    # if (mask is None) and (mask_component is not None):
    #     mask = mask_component


    
    results['total_flux_nomask'] = total_flux_nomask
    results['total_flux_mask'] = total_flux
    results['peak_of_flux'] = np.nanmax(g)
    results['bmajor'] = omaj
    results['bminor'] = omin
    results['freq'] = freq*1e9
    results['flag_freq'] = flag_freq
    results['cell_size'] = cell_size
    results['flux_units'] = flux_units
    results['flux_conversion'] = flux_conversion
    

    
    # flux_mc, flux_error_m = mc_flux_error(img, g,
    #               res,
    #               num_threads=6, n_samples=1000)
    #
    # results['flux_mc'] = flux_mc
    # results['flux_error_mc'] = flux_error_m

    systematic_error = systematic_error_fraction * results['total_flux_mask']

    if residual is not None:
        if data_res is None:
            if isinstance(residual, str) == True:
                data_res = load_fits_data(residual)
            else:
                data_res = residual
        else:
            data_res = data_res
        data_res = np.nan_to_num(data_res,nan=0)
        # rms_res =imstat(residual_name)['flux'][0]
        total_flux_density_residual = np.nansum(data_res * mask) / beam_area_
        # total_flux_density_residual = np.sqrt(np.nansum((data_res * mask)**2.0)) / beam_area_
        mask_area = np.nansum(mask)
        noise_term = std * np.sqrt(mask_area / beam_area_)

        total_flux_density_error = np.sqrt(
            systematic_error**2 +
            noise_term**2 +
            total_flux_density_residual**2
        )

        results['max_residual']          = np.nanmax(data_res * mask)
        results['min_residual']          = np.nanmin(data_res * mask)
        results['flux_residual']         = total_flux_density_residual
        results['total_flux_error']      = abs(total_flux_density_error)
        results['mad_std_residual']      = mad_std(data_res,ignore_nan=True)
        results['rms_residual']          = rms_estimate(data_res)


    else:
        if verbose > 0:
            if logger is not None:
                logger.warning(f" !!>  No residual image provided."
                               f" Make sure to set the systematic_error_fraction to"
                               f" at least 10% (0.1) to account for flux calibration errors.")
            else:
                print(f" !!>  No residual image provided."
                      f" Make sure to set the systematic_error_fraction to"
                      f" at least 10% (0.1) to account for flux calibration errors.")
        noise_term = std * np.sqrt(np.nansum(mask) / beam_area_)
        total_flux_density_error = np.sqrt(systematic_error**2 + noise_term**2)
        results['total_flux_error'] = total_flux_density_error

    # print('total_flux_density_residual = ',total_flux_density_residual)
    # print('systematic_error = ', systematic_error)
    # print('(std / beam_area_) = ', (std / beam_area_))
    if verbose >= 1:
        print('-----------------------------------------------------------------')
        print('Flux Density and uncertainty (systematic + noise x sqrt(N_beams) + residual):')
        print(f"Flux Density = {results['total_flux_mask'] * flux_conversion:.3f} "
            f"+/- {results['total_flux_error'] * flux_conversion:.3f} mJy")
        print(f"Fractional error = {(results['total_flux_error'] / results['total_flux_mask']):.3f}")
        print('-----------------------------------------------------------------')

    try:
        levels = np.geomspace(np.nanmax(g), last_level * std,bins)
        levels2 = np.geomspace(np.nanmax(g), last_level * std,100)
        levels_ellipse = np.geomspace(np.nanmax(g), last_level * std,32)
    except:
        try:
            levels = np.geomspace(np.nanmax(g), last_level * np.nanstd(g),bins)
            levels2 = np.geomspace(np.nanmax(g), last_level * np.nanstd(g),100)
            levels_ellipse = np.geomspace(np.nanmax(g), last_level * np.nanstd(g),32)
        except:
            levels = np.linspace(np.nanmax(g), last_level * std,bins)
            levels2 = np.linspace(np.nanmax(g), last_level * std,100)
            levels_ellipse = np.linspace(np.nanmax(g), last_level * std,32)
            
    peak_position = np.unravel_index(np.nanargmax(g), g.shape)
    if residual is not None:
        results['peak_at_residual'] = res[peak_position]
        peak_stats = calculate_peak_error_annulus(g, res, peak_position, get_beam_size_px(img)[0],
                                 inner_factor=1.5, outer_factor=3.0,
                                 min_beams=5)
        results['peak_error'] = peak_stats['peak_error']# * np.sqrt(beam_area_/2)
        results['peak_snr'] = peak_stats['snr']
    else:
        peak_stats = None
        # results['peak_error'] = results['peak_of_flux']*0.1
        results['peak_error'] = std  # global mad_std, same units as peak
        results['peak_snr'] = results['peak_of_flux'] / std

    try:
        results['convex_error_flag'] = False
        cp = convex_morpho(g, mask, do_plot=True, n_boot=0, noise_std=6*std)
        results['PA_convex']                  = cp['PA_convex']
        results['PA_convex_err']              = cp['PA_convex_err']
        results['q_convex']                   = cp['q_convex']
        results['q_convex_err']               = cp['q_convex_err']
        results['centroid_convex']            = cp['centroid_convex']
        results['convex_major_diameter']      = cp['major_diameter']
        results['convex_major_diameter_err']  = cp['major_diameter_err']
        results['convex_minor_diameter']      = cp['minor_diameter']
        results['convex_minor_diameter_err']  = cp['minor_diameter_err']
        results['q_diam']                     = cp['q_diam']
        results['q_diam_err']                 = cp['q_diam_err']
        results['convex_axes_swapped']        = cp['axes_swapped']
    except Exception:
        results['convex_error_flag']          = True
        results['PA_convex']                  = 0.0
        results['PA_convex_err']              = float('nan')
        results['q_convex']                   = 0.0
        results['q_convex_err']               = float('nan')
        results['centroid_convex']            = [0, 0]
        results['convex_major_diameter']      = 0.0
        results['convex_major_diameter_err']  = float('nan')
        results['convex_minor_diameter']      = 0.0
        results['convex_minor_diameter_err']  = float('nan')
        results['q_diam']                     = 0.0
        results['q_diam_err']                 = float('nan')
        results['convex_axes_swapped']        = False
    if verbose > 1:
        if not results['convex_error_flag']:
            print(f"Convex hull PA: {results['PA_convex']:.2f} deg +/- {results['PA_convex_err']:.2f} deg  "
                  f"q: {results['q_convex']:.3f} +/- {results['q_convex_err']:.3f}  "
                  f"major: {results['convex_major_diameter']:.2f} +/- {results['convex_major_diameter_err']:.2f}px  "
                  f"minor: {results['convex_minor_diameter']:.2f} +/- {results['convex_minor_diameter_err']:.2f}px  "
                  f"q_diam: {results['q_diam']:.3f}")
            if results['convex_axes_swapped']:
                print("  (convex axes labelled by measured extent: the largest "
                      "extent lies along the second eigenvector, so PA_convex "
                      "is rotated 90 deg from the major-eigenvector direction)")
        else:
            print("Convex hull properties could not be calculated.")

    # Geometry
    x0max, y0max = peak_center(g * mask)
    results['x0'], results['y0'] = x0max, y0max
    # determine momentum centres.
    x0m, y0m, a_mom, b_mom, q_mom, PAdeg_mom = momenta(g * mask, 
                                                       PArad_0=None, q_0=None)
    results['x0m'], results['y0m'] = x0m, y0m
    results['a_mom'], results['b_mom'] = a_mom, b_mom
    results['q_mom'], results['PAdeg_mom'] = q_mom, PAdeg_mom
    if verbose > 1:
        print(f"Peak center: (x0, y0) = ({results['x0']:.2f}, {results['y0']:.2f})")
        print(f"Moment center: (x0m, y0m) = ({results['x0m']:.2f}, {results['y0m']:.2f})")
        print(f"Moment major axis (a_mom): {results['a_mom']:.2f} pixels")
        print(f"Moment minor axis (b_mom): {results['b_mom']:.2f} pixels")
        print(f"Moment axis ratio (q_mom): {results['q_mom']:.2f}")
        print(f"Moment position angle (PAdeg_mom): {results['PAdeg_mom']:.2f} degrees")


    fluxes = []
    fluxes_err = []  # per-annulus flux uncertainty
    areas = []
    for i in range(len(levels)):
        if i == 0:
            condition = (g >= levels[i])
            flux = np.nansum(g * condition) / beam_area_
            area = np.nansum(condition)
            fluxes.append(flux)
            areas.append(area)
            local_std = np.nanstd(g[condition]) if area > 0 else std
            fluxes_err.append(local_std * np.sqrt(max(area, 1)) / beam_area_)
            # plt.imshow(condition.astype(float), origin='lower', alpha=0.5)
            # plt.show()
        else:
            condition = ((g < levels[i - 1]) & (g >= levels[i]))
            flux = np.nansum(g * condition) / beam_area_
            area = np.nansum(condition)
            fluxes.append(flux)
            areas.append(area)
            local_std = np.nanstd(g[condition]) if area > 0 else std
            fluxes_err.append(local_std * np.sqrt(max(area, 1)) / beam_area_)
            # plt.imshow(condition.astype(float), origin='lower', alpha=0.5)
            # plt.show()
    fluxes = np.asarray(fluxes)
    fluxes_err = np.asarray(fluxes_err)
    areas = np.asarray(areas)

    agrow = areas.copy()
    results['total_flux_levels'] = np.nansum(fluxes)

    Lgrow = np.nancumsum(fluxes)
    # Cumulative error via quadrature sum - bins-independent
    Lgrow_err = np.sqrt(np.nancumsum(fluxes_err ** 2))

    _radii = []
    for i in range(0, len(levels)):
        mask_at_level_i = g >= levels[i]
        circular_radii = np.sqrt(np.nansum(mask_at_level_i) / np.pi)
        _radii.append(circular_radii)
    radii = np.asarray(_radii)

    Lgrow, radii = check_flux_duplicates(Lgrow, radii)
    # apply same deduplication to error array
    Lgrow_err = Lgrow_err[:len(Lgrow)]  # safe after dedup (duplicates merge adjacent bins)
    Lgrow_norm = Lgrow / np.nansum(fluxes)
    Lgrow_err_norm = Lgrow_err / np.nansum(fluxes)
    
    """testing
    fluxes, flux_errors, Lgrow, Lgrow_errors, areas = compute_flux_with_uncertainties(
        g, levels, beam_area_, noise_rms=std, method='noise_based'
    )
    
    agrow = areas.copy()
    _radii = []
    for i in range(0, len(levels)):
        mask_at_level_i = g >= levels[i]
        circular_radii = np.sqrt(np.nansum(mask_at_level_i) / np.pi)
        # print(circular_radii)
        _radii.append(circular_radii)
    radii = np.asarray(_radii)
    # print(Lgrow_norm)
    
    # Lgrow_errors,_ = check_flux_duplicates(Lgrow_errors,radii)
    Lgrow,radii = check_flux_duplicates(Lgrow,radii)
    Lgrow_norm = Lgrow / np.nansum(fluxes)
    results['total_flux_levels'] = np.nansum(fluxes)
    # results['fluxes'] = fluxes
    # results['Lgrow'] = Lgrow
    # results['radii'] = radii
    # results['flux_errors'] = flux_errors
    # results['Lgrow_errors'] = Lgrow_errors
    # results['total_flux_error'] = np.sqrt(np.nansum(flux_errors**2))
    """

    # rad_prof = elliptical_radial_profile(
    #     image = g,                  # already prepared (masked, bkg-subtracted)
    #     x0 = results['x0'],      # peak centre - or results['x0m'] for moment centre
    #     y0 = results['y0'],
    #     pa_deg = results['PA_convex'],      # from cal_PA_q / do_fit_ellipse block
    #     q = results['q_convex'],
    #     beam_area = beam_area_,         # already computed just above
    #     rms = std,                # already computed from mad_std or passed-in rms
    #     mask = mask,              # already computed or passed in, same shape as g
    # )
    # # -- drop-in replacements for the old level-based variables --------------
    # fluxes = rad_prof['fluxes']
    # fluxes_err = rad_prof['fluxes_err']
    # areas = rad_prof['areas']
    # agrow = rad_prof['agrow']
    # Lgrow = rad_prof['Lgrow']
    # Lgrow_err = rad_prof['Lgrow_err']
    # Lgrow_norm = rad_prof['Lgrow_norm']
    # Lgrow_err_norm = rad_prof['Lgrow_err_norm']
    # radii = rad_prof['radii'] 
    
    ####################################################################
    ## THIS NEEDS A BETTER IMPLEMENTATION USING SPLINE!!!!  ############
    ####################################################################
    mask_L20 = Lgrow_norm < 0.2
    mask_L50 = Lgrow_norm < 0.5
    mask_L80 = Lgrow_norm < 0.8
    try:
        mask_L90 = (Lgrow_norm > 0.89) & (Lgrow_norm < 0.91)
        mask_L95 = (Lgrow_norm > 0.95) & (Lgrow_norm < 0.97)
    except:
        mask_L90 = (Lgrow_norm > 0.85) & (
                    Lgrow_norm < 0.95)  # in case there not enough pixels
        mask_L95 = (Lgrow_norm > 0.92) & (
                    Lgrow_norm < 0.97)  # in case there not enough pixels

    mask_L20_idx = [i for i, x in enumerate(mask_L20) if x]
    mask_L50_idx = [i for i, x in enumerate(mask_L50) if x]
    mask_L80_idx = [i for i, x in enumerate(mask_L80) if x]
    mask_L90_idx = [i for i, x in enumerate(mask_L90) if x]
    mask_L95_idx = [i for i, x in enumerate(mask_L95) if x]
    """
    The tries and exceptions bellow are temporary solutions to avoid 
    errors when there are not enough pixels or no signal at all, when computing 
    the growth curve.    
    """
    try:
        sigma_20 = levels[mask_L20_idx[-1]]
        flag20 = False
    except:
        flag20 = True
        try:
            sigma_20 = levels[mask_L50_idx[-1]]
        except:
            sigma_20 = last_level * std

    try:
        sigma_50 = levels[mask_L50_idx[-1]]
        flag50 = False
    except:
        flag50 = True
        sigma_50 = last_level * std
    try:
        sigma_80 = levels[mask_L80_idx[-1]]
        sigma_90 = levels[mask_L90_idx[-1]]
        sigma_95 = levels[mask_L95_idx[-1]]
        flag9095 = False
    except:
        sigma_80 = last_level * std
        sigma_90 = last_level * std
        sigma_95 = last_level * std
        flag9095 = True

    # up_low_std = (np.max(gal) / (3 * mad_std(gal))) * 0.7
    # inner_shell_mask = ((gal<(up_low_std+5)*mad_std(gal)) & (gal>up_low_std*mad_std(gal)))
    # inner_shell_mask = ((g <sigma_50) & (g >sigma_50*0.95))
    # outer_shell_mask = ((g < sigma_90) & (g > sigma_90 * 0.95))

    inner_mask = (g > (sigma_50)) * mask
    outer_mask90 = (g > (sigma_90)) * mask
    inner_perimeter = perimeter_crofton(inner_mask, 4)
    outer_perimeter90 = perimeter_crofton(outer_mask90, 4)
    outer_perimeter = perimeter_crofton(mask, 4)
    


    if do_fit_ellipse:
        try:
            # some geometrical measures
            # calculate PA and axis-ratio
            # region_split = [i for i, x in enumerate(levels > sigma_50) if x][-1]
            region_split = [i for i, x in enumerate(levels_ellipse > sigma_50) if x][-1]
            
            PA, q, x0col, y0col, PAm, qm, \
                PAmi, qmi, PAmo, qmo, \
                x0median, y0median, \
                x0median_i, y0median_i, \
                x0median_o, y0median_o,profiles = cal_PA_q(g * mask, Isequence=levels_ellipse,
                                                region_split=region_split,
                                                # SAVENAME=img.replace('.fits','_ellipsefit') + ext
                                                SAVENAME=img.replace('.fits','_ellipsefit')
                                                )
        except:
                PA, q, x0col, y0col, PAm, qm = 0.0, 0.0, x0max,y0max, 0.0, 0.0
                PAmi, qmi, PAmo, qmo = 0.0, 0.0, 0.0, 0.0
                x0median, y0median = x0max,y0max
                x0median_i, y0median_i = x0max,y0max
                x0median_o, y0median_o = x0max,y0max
                profiles = None

        results['PA'], results['q'] = PA, q
        results['PAm'], results['qm'] = PAm, qm
        results['PAmi'], results['qmi'] = PAmi, qmi
        results['PAmo'], results['qmo'] = PAmo, qmo
        results['x0m_i'], results['y0m_i'] = x0median_i, y0median_i
        results['x0m_o'], results['y0m_o'] = x0median_o, y0median_o

        vx = results['x0'] - results['x0m']
        vy = results['y0'] - results['y0m']
        try:
            TvPA, Tvlength = trail_vector(vx=vx, vy=vy, v0=np.asarray([1, 0]))
        except:
            TvPA, Tvlength = 0.0, 0.0
        results['TvPA'] = TvPA
        results['Tvlength'] = Tvlength

    try:
        L20_norm = Lgrow_norm[mask_L20_idx[-1]]  # ~ 0.2
        L20 = Lgrow[mask_L20_idx[-1]]
    except:
        try:
            L20_norm = Lgrow_norm[mask_L50_idx[-1]]  # ~ 0.2
            L20 = Lgrow[mask_L50_idx[-1]]
        except:
            L20 = 0.0
            L20_norm = 0.9999
    try:
        L50_norm = Lgrow_norm[mask_L50_idx[-1]]  # ~ 0.5
        L50 = Lgrow[mask_L50_idx[-1]]
    except:
        L50_norm = 0.9999
        L50 = 0.0
    try:
        L80_norm = Lgrow_norm[mask_L80_idx[-1]]  # ~ 0.8
        L80 = Lgrow[mask_L80_idx[-1]]
    except:
        L80_norm = 0.9999
        L80 = 0.0
    try:
        """
        Not enough pixels
        """
        L90_norm = Lgrow_norm[mask_L90_idx[-1]]  # ~ 0.9
        L90 = Lgrow[mask_L90_idx[-1]]
        L95_norm = Lgrow_norm[mask_L95_idx[-1]]  # ~ 0.9
        L95 = Lgrow[mask_L95_idx[-1]]
        flagL9095 = False

    except:
        flagL9095 = True
        try:
            try:
                L90_norm = Lgrow_norm[mask_L80_idx[-1]]  # ~ 0.9
                L90 = Lgrow[mask_L80_idx[-1]]
                L95_norm = Lgrow_norm[mask_L80_idx[-1]]  # ~ 0.9
                L95 = Lgrow[mask_L80_idx[-1]]
            except:
                L90_norm = Lgrow_norm[-1]  # ~ 0.9
                L90 = Lgrow[-1]
                L95_norm = Lgrow_norm[-1]  # ~ 0.9
                L95 = Lgrow[-1]

        except:
            L90_norm = 0.9999
            L90 = fluxes.sum()
            L95_norm = 0.9999
            L95 = fluxes.sum()


    if flag_freq == False:
        try:
            TB20 = T_B(omaj, omin, freq, L20)
        except:
            TB20 = 0.0
        TB50 = T_B(omaj, omin, freq, L50)
        TB80 = T_B(omaj, omin, freq, L80)
        TB90 = T_B(omaj, omin, freq, L90)

        TB = T_B(omaj, omin, freq, total_flux)
    else:
        TB20 = 0.0
        TB50 = 0.0
        TB80 = 0.0
        TB90 = 0.0
        TB = 0.0

    levels_20 = np.asarray([sigma_20])
    levels_50 = np.asarray([sigma_50])
    levels_80 = np.asarray([sigma_80])
    levels_90 = np.asarray([sigma_90])
    levels_95 = np.asarray([sigma_95])

    levels_3sigma = np.asarray([3 * std])

    g20 = ((g * mask) > sigma_20)
    g50 = ((g * mask) > sigma_50)
    g80 = ((g * mask) > sigma_80)
    g90 = ((g * mask) > sigma_90)
    g95 = ((g * mask) > sigma_95)

    try:
        C20radii, L20, npix20, C20radii_err, npix20_err = \
            find_fractional_radius(cumulative_flux = Lgrow, 
                                   radii = radii, 
                                   fraction = 0.20,
                                   cumulative_flux_err=Lgrow_err_norm,
                                #    rms=std, beam_area=beam_area_,
                                   )
        L20_norm = 0.20
    except:
        if logger is not None:
            logger.warning(f" !!==> Not enough pixels to calculate L20-R20."
                        f"       Using flux/size of first pixel.")
        C20radii = radii[0]
        L20 = Lgrow[0]
        L20_norm = Lgrow_norm[0]
        npix20 = radii_to_area(C20radii)[0]
        C20radii_err, npix20_err = C20radii/2, npix20/2
    
    try:
        if logger is not None:
            logger.warning(f" !!==> Not enough pixels to calculate L50-R50."
                           f"       Using flux/size of first pixel.")
        C50radii, L50, npix50, C50radii_err, npix50_err = \
            find_fractional_radius(cumulative_flux = Lgrow, 
                                   radii = radii, 
                                   fraction = 0.50,
                                   cumulative_flux_err=Lgrow_err_norm,
                                #    do_plot=True
                                #    rms=std, beam_area=beam_area_,
                                   )
        L50_norm = 0.50
    except:
        C50radii = radii[0]
        L50 = Lgrow[0]
        L50_norm = Lgrow_norm[0]
        npix50 = radii_to_area(C50radii)[0]
        C50radii_err, npix50_err = C50radii/2, npix50/2
    
    try:
        C80radii, L80, npix80, C80radii_err, npix80_err = \
            find_fractional_radius(cumulative_flux = Lgrow, 
                                   radii = radii, 
                                   fraction = 0.80,
                                   cumulative_flux_err=Lgrow_err_norm,
                                #    rms=std, beam_area=beam_area_,
                                   )
        C90radii, L90, npix90, C90radii_err, npix90_err = \
            find_fractional_radius(cumulative_flux = Lgrow, 
                                   radii = radii, 
                                   fraction = 0.90,
                                   cumulative_flux_err=Lgrow_err_norm,
                                #    rms=std, beam_area=beam_area_,
                                   )
        C95radii, L95, npix95, C95radii_err, npix95_err = \
            find_fractional_radius(cumulative_flux = Lgrow, 
                                   radii = radii, 
                                   fraction = 0.95,
                                   cumulative_flux_err=Lgrow_err_norm,
                                #    rms=std, beam_area=beam_area_,
                                #    do_plot=True
                                   )
        C99radii, L99, npix99, C99radii_err, npix99_err = \
            find_fractional_radius(cumulative_flux = Lgrow, 
                                   radii = radii, 
                                   fraction = 0.99,
                                   cumulative_flux_err=Lgrow_err_norm)
    except:
        if logger is not None:
            logger.warning(f" !!==> Not enough pixels to calculate L80-R80, L90-R90, L95-R95."
                        f"       Using flux/size of first pixel.")
        C80radii = radii[-1]
        C90radii = radii[-1]
        C95radii = radii[-1]
        C99radii = radii[-1]
        L80 = Lgrow[-1]
        L90 = Lgrow[-1]
        L95 = Lgrow[-1]
        L99 = Lgrow[-1]
        npix80 = radii_to_area(C80radii)[0]
        npix90 = radii_to_area(C90radii)[0]
        npix95 = radii_to_area(C95radii)[0]
        npix99 = radii_to_area(C99radii)[0]
        C80radii_err, npix80_err = C80radii/2, npix80/2
        C90radii_err, npix90_err = C90radii/2, npix90/2
        C95radii_err, npix95_err = C95radii/2, npix95/2
        C99radii_err, npix99_err = C99radii/2, npix99/2
        L80_norm = 0.80
        L90_norm = 0.90
        L95_norm = 0.95
        L99_norm = 0.99

    # _frac_fluxes = np.asarray([0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 ,
    #                            0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95, 0.99])
    # _frac_fluxes = np.asarray([0.85, 0.9 , 0.95, 0.99])
    # frac_radii  = {}
    # frac_fluxes = {}
    # for frac_flux in _frac_fluxes:
    #     # try:
    #     radius_percent, _, _ = find_fractional_radius(Lgrow, radii, fraction = frac_flux)
    #     frac_radii[f'R{int(frac_flux*100)}'] = radius_percent
    #     frac_fluxes[f'S{int(frac_flux*100)}'] = frac_flux
        # except:
        #     pass
    # print(frac_radii.keys())
    # print(frac_fluxes)
    # print(np.sqrt(((frac_radii['R99'] - frac_radii['R95'])**2.0 + (frac_radii['R90'] - frac_radii['R95'])**2.0  + (frac_radii['R95'] - frac_radii['R85'])**2.0)/3.0))
    
    # plt.figure()
    # plt.plot(frac_fluxes,frac_radii,'o')
    # plt.show()
    # np.mean(frac_radii)

    # L80_norm = 0.80
    # L90_norm = 0.90
    # L95_norm = 0.95

    
    A20 = npix20 / beam_area_
    A50 = npix50 / beam_area_
    A80 = npix80 / beam_area_
    A90 = npix90 / beam_area_
    A95 = npix95 / beam_area_
    A99 = npix99 / beam_area_
    A20_err = npix20_err / beam_area_
    A50_err = npix50_err / beam_area_
    A80_err = npix80_err / beam_area_
    A90_err = npix90_err / beam_area_
    A95_err = npix95_err / beam_area_
    A99_err = npix99_err / beam_area_

    # A20, C20radii, npix20 = estimate_area((g > sigma_20) * mask, cell_size, omaj,
    #                                       omin)
    # A50, C50radii, npix50 = estimate_area((g > sigma_50) * mask, cell_size, omaj,
    #                                       omin)
    # A80, C80radii, npix80 = estimate_area((g > sigma_80) * mask, cell_size, omaj,
    #                                       omin)
    # A90, C90radii, npix90 = estimate_area((g > sigma_90) * mask, cell_size, omaj,
    #                                       omin)
    # A95, C95radii, npix95 = estimate_area((g > sigma_95) * mask, cell_size, omaj,
    #                                       omin)

    # if flag20 == True:
    #     if flag50 == False:
    #         A20, C20radii, npix20 = A50 / 2, C50radii / 2, npix50 / 2
    #     else:
    #         try:
    #             A20, C20radii, npix20 = A80 / 3, C80radii / 3, npix80 / 3
    #         except:
    #             A20, C20radii, npix20 = A95 / 4, C95radii / 4, npix95 / 4

    try:
        results['conv_P20'], results['conv_A20'] = convex_shape(g20)
    except:
        results['conv_P20'], results['conv_A20'] = 5,5

    try:
        results['conv_P50'], results['conv_A50'] = convex_shape(g50)
        results['conv_P80'], results['conv_A80'] = convex_shape(g80)
        results['conv_P90'], results['conv_A90'] = convex_shape(g90)
        results['conv_P95'], results['conv_A95'] = convex_shape(g95)
        results['conv_PT'], results['conv_AT'] = convex_shape(mask)
    except:
        results['conv_P50'], results['conv_A50'] = 5,5
        results['conv_P80'], results['conv_A80'] = 5,5
        results['conv_P90'], results['conv_A90'] = 5,5
        results['conv_P95'], results['conv_A95'] = 5, 5
        results['conv_PT'], results['conv_AT'] = 5,5


    # This is a more robust computation of source size.
    try:
        R20med,R20mean,R20std = calculate_radii(g, g20)
    except:
        R20med, R20mean, R20std = 1,1,1
    try:
        R50med,R50mean,R50std = calculate_radii(g, g50)
    except:
        R50med, R50mean, R50std = 1,1,1
    try:
        R80med,R80mean,R80std = calculate_radii(g, g80)
    except:
        R80med, R80mean, R80std = 1,1,1
    try:
        R90med,R90mean,R90std = calculate_radii(g, g90)
    except:
        R90med, R90mean, R90std = 1,1,1
    try:
        R95med,R95mean,R95std = calculate_radii(g, g95)
    except:
        R95med, R95mean, R95std = 1,1,1
    try:
        RTmed,RTmean,RTstd = calculate_radii(g, mask)
    except:
        RTmed, RTmean, RTstd = 1,1,1

    results['R20med'], results['R20mean'], \
        results['R20std'] = R20med,R20mean,R20std
    results['R50med'], results['R50mean'], \
        results['R50std'] = R50med,R50mean,R50std
    results['R80med'], results['R80mean'], \
        results['R80std'] = R80med,R80mean,R80std
    results['R90med'], results['R90mean'], \
        results['R90std'] = R90med,R90mean,R90std
    results['R95med'], results['R95mean'], \
        results['R95std'] = R95med,R95mean,R95std
    results['RTmed'], results['RTmean'], \
        results['RTstd'] = RTmed,RTmean,RTstd

    # # print(C20radii, C50radii, C80radii, C90radii)
    # C1 = np.log10(C80radii / C20radii)
    # C2 = np.log10(C90radii / C50radii)
    
    C1, C1_err, C2, C2_err = concentration_index_errors(C80radii, C20radii, C90radii, C50radii, 
                                                        C80radii_err, C20radii_err, C90radii_err, C50radii_err)

    AC1 = np.log10(A80 / A20)
    AC2 = np.log10(A90 / A50)

    CAC1 = np.log10(area_to_radii(results['conv_A80']) / area_to_radii(results['conv_A20']))
    CAC2 = np.log10(area_to_radii(results['conv_A90']) / area_to_radii(results['conv_A50']))



    area_total, Cradii, npix_total = estimate_area(mask, cell_size, omaj, omin)
    o_area_total, o_Cradii, o_npix_total = estimate_area(omask, cell_size, omaj,
                                                         omin)

    mask_outer = (g < sigma_50) & (g > last_level * std) * mask
    A50_100, C50_100radii, npix50_100 = estimate_area(mask_outer, cell_size,
                                                      omaj, omin)
    # gaussianity = 1 / (L50 / g.max())
    gaussianity_L50 = g.max() / L50
    gaussianity = sigma_50 / g.max()

    mask_outer_full = (g < sigma_50) * mask
    #     plt.imshow(mask_outer_full)
    A50_full, C50_full_radii, npix50_full = estimate_area(mask_outer_full,
                                                          cell_size, omaj, omin)

    radii_ratio = C50radii / C50_100radii
    radii_ratio_full = C50radii / C50_full_radii
    area_ratio = A50 / A50_100
    area_ratio_full = A50 / A50_full

    results['max'] = np.nanmax(g)
    results['std_image'] = std
    try:
        results['std_residual'] = mad_std(load_fits_data(residual))
    except:
        pass

    
    results['L50'] = L50
    results['sigma_50'] = sigma_50
    results['sigma_90'] = sigma_90
    results['sigma_95'] = sigma_95
    results['o_area_total'] = o_area_total
    results['o_Cradii'] = o_Cradii
    results['o_npix_total'] = o_npix_total
    results['area_total'] = area_total
    results['Cradii'] = Cradii
    results['npix_total'] = npix_total
    results['inner_perimeter'] = inner_perimeter
    results['outer_perimeter'] = outer_perimeter
    results['outer_perimeter90'] = outer_perimeter90

    results['TB20'] = TB20
    results['TB50'] = TB50
    results['TB80'] = TB80
    results['TB'] = TB

    #     results['nbeams_total'] = nbeams_total

    results['C1'] = C1
    results['C2'] = C2
    results['C1_err'] = C1_err
    results['C2_err'] = C2_err
    results['AC1'] = AC1
    results['AC2'] = AC2
    results['CAC1'] = CAC1
    results['CAC2'] = CAC2

    results['L20'] = L20
    results['A20'] = A20
    results['C20'] = sigma_20 / std
    results['C20radii'] = C20radii
    results['npix20'] = npix20
    results['C20radii_err'] = C20radii_err
    results['npix20_err'] = npix20_err

    results['L50'] = L50
    results['A50'] = A50
    results['C50'] = sigma_50 / std
    results['C50radii'] = C50radii
    results['C50radii_err'] = C50radii_err
    results['npix50'] = npix50
    results['npix50_err'] = npix50_err
    
    results['C50theta'] = C50radii * np.sqrt(6 * np.log(2))
    results['C50theta_err'] = C50radii_err * np.sqrt(6 * np.log(2))
    results['C50theta_dec'] = np.sqrt(results['C50theta']**2.0 - (results['bmajor']/cell_size)**2.0)
    # results['C50theta_dec_err'] = np.sqrt(results['C50theta_err']**2.0 - (results['bmajor']/cell_size)**2.0)
    results['C50theta_dec_err'] = results['C50theta_dec'] * results['C50theta_err']
    
    
    

    results['L80'] = L80
    results['A80'] = A80
    results['C80'] = sigma_80 / std
    results['C80radii'] = C80radii
    results['npix80'] = npix80
    results['C80radii_err'] = C80radii_err
    results['npix80_err'] = npix80_err

    results['L90'] = L90
    results['A90'] = A90
    results['C90'] = sigma_90 / std
    results['C90radii'] = C90radii
    results['npix90'] = npix90
    results['C90radii_err'] = C90radii_err
    results['npix90_err'] = npix90_err
    results['flag20'] = flag20
    results['flag50'] = flag50
    results['flag9095'] = flag9095
    results['flagL9095'] = flagL9095


    results['L95'] = L95
    results['A95'] = A95
    results['C95'] = sigma_95 / std
    results['C95radii'] = C95radii
    results['npix95'] = npix95
    results['C95radii_err'] = C95radii_err
    results['npix95_err'] = npix95_err

    results['L99'] = L99
    results['A99'] = A99
    results['C99radii'] = C99radii
    results['npix99'] = npix99
    results['C99radii_err'] = C99radii_err
    results['npix99_err'] = npix99_err
    
    results['A20_err'] = A20_err
    results['A50_err'] = A50_err
    results['A80_err'] = A80_err
    results['A90_err'] = A90_err
    results['A95_err'] = A95_err
    results['A99_err'] = A99_err


    results['gaussianity'] = gaussianity
    results['gaussianity_L50'] = gaussianity_L50

    results['A50_100'] = A50_100
    results['C50_100radii'] = C50_100radii
    results['npix50_100'] = npix50_100

    results['A50_full'] = A50_full
    results['C50_full_radii'] = C50_full_radii
    results['npix50_full'] = npix50_full

    results['radii_ratio'] = radii_ratio
    results['radii_ratio_full'] = radii_ratio_full
    results['area_ratio'] = area_ratio
    results['area_ratio_full'] = area_ratio_full
    results['beam_area'] = beam_area_
    
    if redshift is not None:
        pix_to_pc = pixsize_to_pc(z=redshift,cell_size = cell_size)
        scale_size = 1 * pix_to_pc
        scale_size_units = 'pc'
    else:
        scale_size = 1 * cell_size
        scale_size_units = scale_units

    if verbose > 0:
            if logger is not None:
                print_logger_header(title="Basic Source Properties",
                                    logger=logger)
                logger.debug(f" ==>  Peak of Flux="
                            f"{results['peak_of_flux']*flux_conversion:.3f} +/- "
                            f"{results['peak_error']*flux_conversion:.3f} [mJy/beam] "
                            f"(SNR={results['peak_snr']:.1f})")
                logger.debug(f" ==>  Total Flux Inside Mask="
                            f"{results['total_flux_mask']*flux_conversion:.3f} +/- "
                            f"{results['total_flux_error']*flux_conversion:.3f} [mJy]")
                logger.debug(f" ==>  Total Flux Image="
                            f"{results['total_flux_nomask'] * flux_conversion:.3f} [mJy]")
                logger.debug(f" ==>  Half-Light Radii="
                            f"{results['C50radii']:.3f} [px]")
                logger.debug(f" ==>  Total Source Size="
                            f"{results['C95radii']:.3f} [px]")
                if do_fit_ellipse:
                    logger.debug(f" ==>  Source Global Axis Ratio="
                                f"{results['qm']:.3f}")
                    logger.debug(f" ==>  Source Global PA="
                                f"{results['PAm']:.3f} [degrees]")
                    logger.debug(f" ==>  Inner Axis Ratio="
                                f"{results['qmi']:.3f}")
                    logger.debug(f" ==>  Outer Axis Ratio="
                                f"{results['qmo']:.3f}")
                    logger.debug(f" ==>  Inner PA="
                                f"{results['PAmi']:.3f} [degrees]")
                    logger.debug(f" ==>  Outer PA="
                                f"{results['PAmo']:.3f} [degrees]")
            else:
                print(f" ==>  Peak of Flux      = "
                    f"{results['peak_of_flux']*flux_conversion:.3f} +/- {results['peak_error']*flux_conversion:.3f} [mJy/beam] "
                    f"(SNR={results['peak_snr']:.1f})")
                print(f" ==>  Snu (within mask) = "
                    f"{results['total_flux_mask']*flux_conversion:.3f} +/- {results['total_flux_error'] * flux_conversion:.3f} [mJy]")
                print(f" ==>  Snu (image)       = "
                    f"{results['total_flux_nomask'] * flux_conversion:.3f} [mJy] (no independent error estimate)")
                print(f" ==>  R50               = "
                    f"{cell_size * results['C50radii']:.4f} +/- {cell_size * results['C50radii_err']:.4f} [{scale_units}]")
                print(f" ==>  R95               = "
                    f"{cell_size * results['C95radii']:.4f} +/- {cell_size * results['C95radii_err']:.4f} [{scale_units}]")
                print(f" ==>  R50               = "
                    f"{results['C50radii']:.4f} +/- {results['C50radii_err']:.4f} [px]")
                print(f" ==>  R95               = "
                    f"{results['C95radii']:.4f} +/- {results['C95radii_err']:.4f} [px]")
                if scale_size_units == 'pc':
                    A50_phy = pix_area_to_kpc_area(results['npix50'],scale_size)
                    A50_phy_err = pix_area_to_kpc_area(results['npix50_err'],scale_size)
                    A95_phy = pix_area_to_kpc_area(results['npix95'],scale_size)
                    A95_phy_err = pix_area_to_kpc_area(results['npix95_err'],scale_size)

                    print(f" ==>  A50               = "
                        f"({1e3*A50_phy:.4f} +/- {1e3*A50_phy_err:.4f}) x 1e3 [kpc^-2]")
                    print(f" ==>  A95               = "
                        f"({1e3*A95_phy:.4f} +/- {1e3*A95_phy_err:.4f}) x 1e3 [kpc^-2]")
                                
                if redshift is not None:
                    print(f" ==>  R50               = "
                        f"{scale_size * results['C50radii']:.4f} +/- {scale_size * results['C50radii_err']:.4f} [{scale_size_units}]")
                    print(f" ==>  R95               = "
                        f"{scale_size * results['C95radii']:.4f} +/- {scale_size * results['C95radii_err']:.4f} [{scale_size_units}]")
                # if do_fit_ellipse:
                #     print(f" ==>  Source Global Axis Ratio="
                #           f"{results['qm']:.2f}")
                #     print(f" ==>  Source Global PA="
                #           f"{results['PAm']:.2f} [degrees]")
                #     print(f" ==>  Inner Axis Ratio="
                #           f"{results['qmi']:.2f}")
                #     print(f" ==>  Outer Axis Ratio="
                #           f"{results['qmo']:.2f}")
                #     print(f" ==>  Inner PA="
                #           f"{results['PAmi']:.2f} [degrees]")
                #     print(f" ==>  Outer PA="
                #           f"{results['PAmo']:.2f} [degrees]")

    # omask = omask2.copy()
    error_petro = False
    if do_petro == True:
        try:
            if verbose >= 1:
                if logger is not None:
                    logger.info(f"++>> Computing Petrosian properties.")
                else:
                    print('++>> Computing Petrosian properties.')
            r_list, area_arr, area_beam, p, flux_arr, error_arr, results, cat, \
                segm, segm_deblend, sorted_idx_list = \
                compute_petrosian_properties(g, img,
                                                mask_component=mask_component,
                                                global_mask=mask,
                                                source_props=results,
                                                apply_mask=False,
                                                # error = ,
                                                sigma_level=sigma_mask,
                                                bkg_to_sub=bkg_to_sub,
                                                vmin=vmin_factor, 
                                                # plot=show_figure,
                                                # deblend=deblend,
                                                # fwhm=fwhm, kernel_size=kernel_size,
                                                show_figure=show_figure,
                                                plot=show_figure,
                                                verbose = verbose,
                                                add_save_name=add_save_name,
                                                # npixels=npixels,
                                                logger=logger)
            error_petro = False
        except Exception as e:
            if logger is not None:
                logger.warning(f"-->> ERROR when computing Petrosian properties. "
                                f"Will flag error_petro as True.")
            else:
                print("-->> ERROR when computing Petrosian properties. Will "
                        "flag error_petro as True.")
            error_petro = True
    else:
        error_petro = True

    results['error_petro'] = error_petro

    if crop == True:
        try:
            xin, xen, yin, yen = do_cutout_2D(img, 
                                              box_size=box_size, 
                                              center=None,
                                              centre_mode = 'image_centre',
                                              return_='box')
            g = g[xin:xen, yin:yen]
            g_original = g_original[xin:xen, yin:yen]
        except:
            try:
                max_x, max_y = np.where(g == np.nanmax(g))
                xin = max_x[0] - box_size
                xen = max_x[0] + box_size
                yin = max_y[0] - box_size
                yen = max_y[0] + box_size
                g = g[xin:xen, yin:yen]
                g_original = g_original[xin:xen, yin:yen]
            except:
                pass
    
    if show_figure == True:
        # Calculate image aspect ratio to determine proper space allocation
        img_height, img_width = g_original.shape
        img_aspect = img_width / img_height

        # Calculate optimal space allocation
        # Give the image subplot just enough width to display properly
        # Give the scatter plot the remaining space (but ensure it's reasonable)
        plot_aspect = figsize[0]/figsize[1]  # Desired aspect ratio for scatter plot (width/height)
        height_ref = figsize[1]  # Reference height

        # Calculate widths needed for each subplot
        img_width_needed = height_ref * img_aspect
        plot_width_needed = height_ref * plot_aspect

        # Create width ratios for gridspec
        total_width = img_width_needed + plot_width_needed
        width_ratios = [plot_width_needed/total_width, img_width_needed/total_width]

        # Set figure size based on calculated needs
        fig_width = total_width * 1.0  # Add some padding
        fig_height = height_ref
        figsize = (fig_width, fig_height)


        
        # Create figure with calculated proportions
        # fig = plt.figure(figsize=figsize)
        fig = matplotlib.figure.Figure(figsize=figsize)
        FigureCanvasAgg(fig)  # attach renderer in-place
        
        gs = gridspec.GridSpec(1, 2, figure=fig, width_ratios=width_ratios, wspace=0.1)

        # Create subplots
        ax1 = fig.add_subplot(gs[0, 0])  # Scatter plot
        ax2 = fig.add_subplot(gs[0, 1])  # Image

        # First subplot - scatter plot
        ax1.scatter(radii*cell_size, Lgrow / results['total_flux_mask'],
                    # label='Norm Masked Flux'
                    )

        ax1.axvspan(C50radii*cell_size - C50radii_err*cell_size,  # left bound
                    C50radii*cell_size + C50radii_err*cell_size,  # right bound
                    alpha=0.3, color='grey', 
                    # hatch='x.',
                    # hatch='////' ,
                    # label=f'$R_{{50}}$ uncertainty'
                    )

        ax1.axvline(C50radii*cell_size,
                    label=r"$R_{50}\sim~$"f"{C50radii*cell_size:0.3f}$\pm"f"{C50radii_err*cell_size:0.3f}''$",
                    ls='-.', color='lime',lw=4)

        ax1.axhline(L50_norm, ls='-.', color='lime',lw=4)
        # ax1.plot(C50radii*cell_size,L50_norm,'rx',markersize=10)


        ax1.axvspan(C95radii*cell_size - C95radii_err*cell_size,  # left bound
                    C95radii*cell_size + C95radii_err*cell_size,  # right bound
                    alpha=0.3, color='grey', 
                    # label=f'$R_{{50}}$ uncertainty'
                    )
        ax1.axvline(C95radii*cell_size,
                    label=r"$R_{95}\sim~$"f"{C95radii*cell_size:0.3f}$\pm"f"{C95radii_err*cell_size:0.3f}''$",
                    color='#4daf4a',lw=3)
        # ax1.plot(C95radii*cell_size,L95_norm,'rx',markersize=10)
        ax1.set_title("Integrated Flux Density \n "
                    r"$S_{\nu} =$ "
                    f"{flux_conversion*total_flux:.3f} $\pm ~ {flux_conversion*total_flux_density_error:.3f}$ [mJy]")

        ax1.set_xlabel(fr'Projected Circular Radius $R$ [{scale_units}]')
        ax1.set_ylabel(r"Normalised  ~  FGC   $~S_{\nu}(\leq R)$")

        ax1.grid(alpha=0.5)
        # ax1.legend(loc='lower right', prop={'size': 10})
        ax1.legend(
            framealpha=0.9,
            # ncol=ncol,
            # loc=loc,
            # fontsize=int(fontsize-1),
            # fontsize=12,
            handlelength=1,
            handletextpad=0.5,
            columnspacing=0.5,
            borderaxespad=0.1
            )
        ax1.set_ylim(0, 1.00)
        # ax1.set_xscale('log')
        # Second subplot - image
        vmin = vmin_factor * std
        if vmax is None:
            vmax = vmax_factor * np.nanmax(g)
        try:
            norm = simple_norm(g_original, stretch='asinh', asinh_a=0.075, vmin=vmin,
                            vmax=vmax)
            im_plot = ax2.imshow(g_original, cmap='magma_r', origin='lower', alpha=1.0,
                                norm=norm, aspect='equal')
        except: #how to avoid errors such as ValueError: vmin must be less than or equal to vmax? 
            print('Error with normalization. Plotting without normalization.')
            im_plot = ax2.imshow(g_original, cmap='magma_r', origin='lower', alpha=1.0,
                                aspect='equal')
            
        if image_title is not None:
            ax2.set_title(image_title)

        # Add contours
        try:
            ax2.contour(g, levels=levels_50, colors='lime', linewidths=2.5,
                        linestyles='-.', alpha=1.0)
            ax2.contour(g, levels=levels_95, colors='#4daf4a', linewidths=2.0,
                        alpha=1.0)
            ax2.contour(g, levels=[last_level * std], colors='cyan', linewidths=0.6,
                        alpha=1.0)
            ax2.contour(g, levels=[6.0 * std], colors='black', linewidths=1.5,
                        alpha=0.9)
            ax2.contour(g, levels=[3.0 * std], colors='brown', linewidths=1.2,
                        alpha=0.9)
        except:
            print('Not plotting contours!')
        
        ax2.axis('off')

        # Apply tight layout
        # plt.tight_layout()
        # plt.show()
        
        try:
            if SAVE is not None:
                if plot_savemane is None:
                    fig.savefig(img.replace('.fits', '_Lgrow_levels')+add_save_name + ext,
                                dpi=300,bbox_inches='tight')
                else:
                    fig.savefig(plot_savemane,
                                dpi=300,bbox_inches='tight')
            buf = io.BytesIO()
            fig.savefig(buf, format='png', dpi=100, bbox_inches='tight')
            buf.seek(0)
            display(Image(data=buf.read()))
            buf.close()
            del fig
            gc.collect()
            # del fig
            # plt.show()
            # plt.close(fig)
            # plt.close('all')
            # plt.clf()
            # plt.cla()
        except Exception as e:
            print(f"Error saving or displaying figure: {e}")


    if save_csv == True:
        import csv
        # df_temp = pd.DataFrame(results).T
        # df_temp.to_csv(img.replace('.fits', '_image_properties')+add_save_name + '.csv',
        #                header=True,index=False)

        with open(img.replace('.fits', '_image_properties')+add_save_name + '.csv',
                  'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=results.keys())
            writer.writeheader()
            writer.writerow(results)
    return (levels, fluxes, Lgrow, radii, agrow, omask, mask, results)


# ------------------------------------------------------------------------------
# Private helpers used exclusively by compute_image_properties_v2
# ------------------------------------------------------------------------------

def _build_levels(peak, floor, n):
    """
    Return n intensity levels spaced from peak down to floor.

    Tries geomspace first (log-uniform spacing), falls back to linspace when
    the floor is non-positive or geomspace raises.
    """
    if floor <= 0 or not np.isfinite(floor) or not np.isfinite(peak):
        return np.linspace(peak, max(floor, 0.0), n)
    try:
        return np.geomspace(peak, floor, n)
    except Exception:
        return np.linspace(peak, floor, n)


def _sigma_at_fraction(fraction, Lgrow_norm_pre, levels, fallback):
    """
    Return the intensity level at which the *pre-dedup* cumulative flux
    normalised growth curve first reaches `fraction`.

    Uses linear interpolation on the monotone pair (Lgrow_norm_pre, levels).
    Both arrays must have the same length and correspond 1-to-1.

    Falls back to `fallback` (e.g. last_level * std) when the fraction is
    beyond the range covered by the growth curve.
    """
    # Guard: Lgrow_norm_pre is increasing, levels is decreasing - clip fraction
    if Lgrow_norm_pre[-1] < fraction:
        return fallback
    if Lgrow_norm_pre[0] >= fraction:
        return float(levels[0])
    # np.interp requires xp increasing -> use Lgrow_norm_pre as x, levels as y
    return float(np.interp(fraction, Lgrow_norm_pre, levels))


def _fractional_radius_safe(Lgrow, radii, Lgrow_err_norm, fraction, fallback_r):
    """
    Thin wrapper around find_fractional_radius that never raises.

    Returns (r, L, npix, r_err, npix_err, ok) where ok=False signals fallback.
    """
    try:
        r, L, npix, r_err, npix_err = find_fractional_radius(
            Lgrow, radii, fraction=fraction,
            cumulative_flux_err=Lgrow_err_norm,
            # do_plot=True
        )
        return r, L, npix, r_err, npix_err, True
    except Exception:
        r      = fallback_r
        L      = float(np.interp(fraction, Lgrow / Lgrow[-1], Lgrow))
        npix   = radii_to_area(r)[0]
        r_err  = r / 2.0
        npix_err = npix / 2.0
        return r, L, npix, r_err, npix_err, False


# ------------------------------------------------------------------------------

def compute_image_properties(img, residual,
                                 cell_size=None, mask_component=None,
                                 redshift=None,
                                 aspect=1, last_level=3.0, mask=None,
                                 bins=64,
                                 data_2D=None, data_res=None,
                                 do_petro=False,
                                 dilation_size=None, iterations=2,
                                 dilation_type='disk', do_fit_ellipse=False,
                                 sigma_mask=6, rms=None,
                                 error_map=None, rms_map=None, invvar_map=None,
                                 variance_map=None, weight_map=None,
                                 use_residual_as_error=False,
                                 systematic_error_fraction=0.05,
                                 results=None, bkg_to_sub=None,
                                 apply_mask=True,
                                 vmax=None, vmin_factor=1.5, vmax_factor=0.7,
                                 crop=False, box_size=256,
                                 SAVE=True, add_save_name='', plot_savemane=None,
                                 show_figure=True,
                                 figsize=(6, 5), image_title=None,
                                 save_csv=False,
                                 flux_units='Jy/beam',
                                 flux_conversion_factor=None,
                                 psf_fwhm_px=None,
                                 ext='.jpg', logger=None, verbose=0):
    """
    Compute morpho-photometric properties of a radio (or optical) image.

    Params
    ------
    img : str
        The path to the image to be analyzed.
    residual : str
        The path to the residual image associated with the image.
    cell_size : float, optional
        The default is None. The size of the pixel in arcseconds.
        If None, `get_cell_size` will be used to get from the header.
    mask_component : 2D np array, optional
        The default is None. If not None, this is the mask for a specific component
        when performing a multi-component source analysis.
    aspect : float, optional (experimental)
        The default is 1. The aspect ratio of the image for plotting.
    last_level : float, optional
        New threshold level (as multiple of sigma_mad) to be used inside the existing mask.
    bins : int, optional
        Number of bins to slice intensity levels of the image. 
    mask : 2D np array, optional
        The default is None. If not None, the function `mask_dilation` will determine it with
        default parameters.
    data_2D : 2D np array, optional
        The default is None. This can be used to pass a 2D array directly to the function.
        For example, when calculating properties from an array without reading from a file,
        e.g. a model image, you can use this. But, to obtain meaningful physical units,
        you must provide the corresponding image file where this array was derived from.
    dilation_size : int, optional
        The default is None. The size of the dilation to be used in the mask dilation.
        If None, the default value will be the size of the restoring beam.
    iterations : int, optional
        The default is 2. The number of iterations to be used in the mask dilation.
        If signs of over-dilation are present, you can set to 1.
    dilation_type : str, optional
        The default is 'disk'. The type of dilation to be used in the mask dilation.
    sigma_mask : float, optional
        The default is 6. The sigma level to be used in the mask dilation.
    rms : float, optional
        The default is None. The rms value to be used in the mask dilation.
        If None, the function `mad_std` will be used to calculate it from the residual image,
        if provided. If the residual image is not provided, the function will use the image itself.
        But, in that case, the result may not be accurate (overestimated) if the image size is
        comparable in size to the size of the source structure.
    results : dict, optional
        The default is None. A dictionary to store the results.
        You can pass an existing external dictionary, so the results will be appended to it.
    bkg_to_sub : 2D np array, optional (EXPERIMENTAL)
        The default is None. The background to be subtracted from the image.
    apply_mask : bool, optional
        The default is True. If True, the mask dilation will be calcualted from the image.
    vmin_factor : float, optional
        The default is 3. The factor (as a multiple of sigma_mad) to be used in the vmin
        calculation for the image plot.
    vmax_factor : float, optional
        The default is 0.5. The factor (as a multiple of peak brightness) to be used in the vmax
        calculation for the image plot.
    crop : bool, optional
        The default is False. If True, the image will be cropped to a box_size.
    box_size : int, optional
        The default is 256. The size of the box to be used in the image cropping.
    SAVE : bool, optional
        The default is True. If True, the image plot will be saved to a file.
    add_save_name : str, optional
        The default is ''. A string to be added to the image plot file name.
    show_figure : bool, optional
        The default is True. If True, the image plot will be shown.
    ext : str, optional
        The default is '.jpg'. The file extension to be used in the image plot file name.
    flux_units : str, optional
        The default is 'Jy/beam'. Names the units the RAW pixel array is
        already in, which is what sets the scale of every reported flux
        density. Recognised without further input: 'Jy/beam'/'Jy'/'mJy' (array
        in Jy), 'mJy/beam' (already mJy), 'nanomaggies' (Legacy Survey/DECam),
        '1nJy'/'10nJy', and 'any' (units unknown, numbers left raw).
    flux_conversion_factor : float, optional
        The default is None. mJy per native pixel unit, overriding the lookup
        above. Required for any instrument whose native unit has no name of its
        own -- HSC/LSST zero-point counts, HST e/s, JWST PHOTFNU -- which
        `read_data` reports as flux_units='mJy/px'. `read_data` computes the
        matching value, so the two travel together:

            flux_units=input_data.flux_units,
            flux_conversion_factor=input_data.flux_conversion_factor

        Without it, an unrecognised `flux_units` leaves the fluxes as raw pixel
        sums and warns.
    logger : logging.Logger, optional
        The default is None. A logger object to be used to log messages.

    Notes
    -----
    Besides the global noise columns (`std`-derived quantities, `mad_std_residual`,
    `rms_residual`, all measured over the whole map), the returned `results` dict
    carries a per-region noise estimate:

        local_rms       `mad_std` of the residual over the region being measured
                        -- `mask` for a whole-source call, `mask_component` for a
                        deblended sub-region. Reported only: nothing else in this
                        function is computed from it.
        local_rms_npix  Number of pixels the estimate is based on.
        local_rms_from  'residual', or 'image' when no residual was given and the
                        science image was sigma-clipped instead.

    They propagate unchanged through `measures` and into the `structural_morphology`
    DataFrame, one value per row/component.

    The peak flux density carries a three-term error budget (see
    `calculate_peak_error_annulus` for the derivation and references):

        peak_error        sqrt(noise^2 + pix^2 + cal^2)
        peak_error_noise  local rms from a 1.5-3 beam annulus on the residual
        peak_error_pix    peak-sampling term, `peak_of_flux` being a pixel maximum
        peak_error_cal    systematic_error_fraction * peak -- the same knob that
                          sets the integrated-flux systematic, so the two agree

        peak_snr          UNCHANGED legacy noise-only ratio. It feeds the
                          Condon (1997)/Fomalont (1999) radius-error floor below
                          and `theta_Rn_*_err` in morphen.py, which require a
                          noise-only SNR.
        peak_snr_noise    peak / peak_error_noise (improved noise-only SNR)
        peak_snr_meas     peak / peak_error

        peak_of_flux_interp     sub-pixel peak; `peak_of_flux` stays the maximum
        peak_pixelisation_bias  mean peak loss from pixel sampling (not applied)
        peak_max_bias_sigma     sqrt(2 ln N_beams), the expected noise maximum
        peak_rms_local / peak_rms_local_clipped / peak_local_nbeams /
        peak_annulus_flag       annulus diagnostics
    """
    # --------------------------------------------------------------------------
    # Setup: units, data loading, RMS, masking
    # --------------------------------------------------------------------------
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(img)

    # mJy per native pixel unit, for the reported flux densities. Every value
    # here (and every `flux_conversion_factor` an override supplies) is in that
    # one unit, so the raw sums below are always just multiplied by it.
    flux_conversion_map = {
        'Jy/beam':    1e3,
        'Jy':         1e3,
        'mJy':        1e3,   # read_data's name for "raw array is Jy(/beam)"
        'nanomaggies': 3.631e-3,
        '1nJy':        1e-6,
        '10nJy':       1e-5,
        'mJy/beam':    1.0,
        'any':         1.0,  # units unknown -- report the raw numbers as-is
    }
    if flux_conversion_factor is not None:
        # An explicit factor wins: it is the only way to describe an instrument
        # whose native unit has no name of its own -- HSC/LSST zero-point
        # counts, HST e/s, JWST PHOTFNU -- which is exactly what read_data's
        # `flux_conversion_factor` carries. Pass the pair straight through:
        #     flux_units=input_data.flux_units,
        #     flux_conversion_factor=input_data.flux_conversion_factor
        flux_conversion = float(flux_conversion_factor)
    elif flux_units in flux_conversion_map:
        flux_conversion = flux_conversion_map[flux_units]
    else:
        # Previously a silent .get(flux_units, 1.0), which turned an unknown
        # unit into a confident, wrongly-scaled flux density. 'mJy/px' hit this
        # for every zero-point/PHOTFNU-calibrated optical image.
        flux_conversion = 1.0
        _msg = (f"flux_units={flux_units!r} has no known mJy-per-pixel-unit "
                f"conversion and no flux_conversion_factor was given -- the "
                f"reported flux densities are RAW pixel sums, not mJy. Pass "
                f"flux_conversion_factor=input_data.flux_conversion_factor "
                f"alongside flux_units to scale them.")
        print(f"!! WARNING !! {_msg}")
        if logger is not None:
            logger.warning(_msg)

    from skimage.draw import disk  # noqa: F401 – kept for v1 compatibility
    if data_2D is not None:
        g_ = data_2D
    else:
        g_ = load_fits_data(img)
    if residual is not None:
        res_ = load_fits_data(residual) if isinstance(residual, str) else residual
        res  = np.nan_to_num(res_.copy(), nan=0)

    g           = np.nan_to_num(g_.copy(), nan=0)
    g_original  = g.copy()

    if bkg_to_sub is not None:
        g = g - bkg_to_sub

    scale_units = 'arcsec'
    if cell_size is None:
        try:
            cell_size = get_cell_size(img)
        except Exception:
            scale_units = 'px'
            cell_size   = 1.0

    try:
        beam_area_  = beam_area2(img, cellsize=cell_size)
        flag_freq   = False
        g_hd        = imhead(img)
        try:
            freq = g_hd['refval'][2] / 1e9
        except Exception:
            freq = getfreqs([img])[0] / 1e9
        omaj, omin, BPA, _, _ = beam_shape(img)
    except Exception as e:
        _warn = (logger.warning if logger else print)
        _warn(f"  -->> ERROR: {e}")
        _warn(f"  -->> Not a radio image?")
        flag_freq  = True
        freq       = 1e15 / 1e9
        omaj = omin = 1
        beam_area_  = 1

    if rms is not None:
        std = rms
    elif residual is not None:
        _res_arr = np.nan_to_num(
            load_fits_data(residual) if isinstance(residual, str) else residual, nan=0
        )
        std = mad_std(_res_arr)
    else:
        std = mad_std(np.nan_to_num(g_.copy(), nan=0))

    # -- Masking ---------------------------------------------------------------
    if mask is not None:
        omask                = mask
        g                    = g * mask
        total_flux           = np.nansum(g)  / beam_area_
        total_flux_nomask    = np.nansum(g_) / beam_area_
        apply_mask           = False

    if apply_mask:
        _src = img if data_res is None else data_2D
        omask, mask = mask_dilation(_src, sigma=sigma_mask,
                                    dilation_size=dilation_size,
                                    iterations=iterations,
                                    dilation_type=dilation_type,
                                    show_figure=False,
                                    verbose=verbose)
        g                    = g * mask
        total_flux           = np.nansum(g)  / beam_area_
        total_flux_nomask    = np.nansum(g_) / beam_area_

    if mask_component is not None:
        # `total_flux_nomask` here means "flux without the mask_component
        # restriction", i.e. within the broader `mask`/`apply_mask` region
        # computed just above -- not the same as an actually-unmasked sum.
        # Falls back to the full unmasked image only if neither a mask nor
        # apply_mask was in effect before this point.
        total_flux_nomask   = total_flux if (mask is not None or apply_mask) \
                               else np.nansum(g_) / beam_area_
        g                    = g * mask_component
        total_flux           = np.nansum(g) / beam_area_
        if mask is not None:
            mask  = mask  * mask_component
            omask = omask * mask_component
        else:
            mask = omask = mask_component

    if (mask_component is None) and (not apply_mask) and (mask is None):
        total_flux           = np.nansum(g * (g > 3 * std)) / beam_area_
        total_flux_nomask    = np.nansum(g) / beam_area_
        unity_mask           = np.ones(g.shape, dtype=bool)
        omask = mask         = unity_mask

    mask_pix_area = np.nansum(mask)
    number_of_beams_in_mask = mask_pix_area / beam_area_
    
    # -- Basic scalars stored immediately -------------------------------------
    results.update({
        'total_flux_nomask': total_flux_nomask,
        'total_flux_mask':   total_flux,
        'peak_of_flux':      np.nanmax(g),
        'bmajor':            omaj,
        'bminor':            omin,
        'freq':              freq * 1e9,
        'flag_freq':         flag_freq,
        'cell_size':         cell_size,
        'flux_units':        flux_units,
        'flux_conversion':   flux_conversion,
        'number_of_beams_in_mask': number_of_beams_in_mask
    })

    systematic_error = systematic_error_fraction * results['total_flux_mask']

    if residual is not None:
        if data_res is None:
            data_res = load_fits_data(residual) if isinstance(residual, str) else residual
        data_res                             = np.nan_to_num(data_res, nan=0)
        mask_area                            = np.nansum(mask)
        total_flux_density_residual          = np.nansum(data_res * mask) / beam_area_
        noise_term                           = std * np.sqrt(mask_area / beam_area_)
        total_flux_density_error             = np.sqrt(
            systematic_error**2 + noise_term**2 + total_flux_density_residual**2
        )
        results['max_residual']          = np.nanmax(data_res * mask)
        results['min_residual']          = np.nanmin(data_res * mask)
        results['flux_residual']         = total_flux_density_residual
        results['total_flux_error']      = abs(total_flux_density_error)
        results['mad_std_residual']      = mad_std(data_res, ignore_nan=True)
        results['rms_residual']          = rms_estimate(data_res)
        # Noise local to the region actually being measured. No extra geometry is
        # needed here: by this point `mask` IS that region -- the whole-source mask
        # for a global call, or `mask * mask_component` for a deblended sub-region
        # (see the masking cascade above) -- so the same two lines give the
        # index-0 value and the per-component values.
        _lrms, _lnpix                    = compute_local_rms(data_res, mask,
                                                             use_sigma_clip=False)
        results['local_rms']             = _lrms
        results['local_rms_npix']        = _lnpix
        results['local_rms_from']        = 'residual'
    else:
        if verbose > 0:
            _warn = (logger.warning if logger else print)
            _warn(" !!>  No residual image provided. "
                  "Set systematic_error_fraction >= 0.1 to account for calibration errors.")
        noise_term                   = std * np.sqrt(np.nansum(mask) / beam_area_)
        total_flux_density_error     = np.sqrt(systematic_error**2 + noise_term**2)
        results['total_flux_error']  = total_flux_density_error
        # No residual to measure the local noise on -- fall back to the science
        # image inside the same region, sigma-clipped so the source itself does
        # not set the scale. `local_rms_from` records which map was used.
        _lrms, _lnpix                = compute_local_rms(g, mask,
                                                         use_sigma_clip=True)
        results['local_rms']         = _lrms
        results['local_rms_npix']    = _lnpix
        results['local_rms_from']    = 'image'

    if verbose >= 1:
        print('-' * 65)
        print('Flux Density (systematic + noise*sqrt(N_beams) + residual):')
        print(f"  S_nu = {results['total_flux_mask']*flux_conversion:.3f}"
              f" +/- {results['total_flux_error']*flux_conversion:.3f} mJy")
        print(f"  Fractional error = "
              f"{results['total_flux_error']/results['total_flux_mask']:.3f}")
        print('-' * 65)

    # --------------------------------------------------------------------------
    # Intensity levels - built once with the robust helper
    # --------------------------------------------------------------------------
    peak_g  = np.nanmax(g)
    floor_g = 1 * last_level * std
    
    # Diagnostic: how much dynamic range the growth curve actually has to
    # work with. Small/faint masks (e.g. deblended sub-regions) can have
    # peak_g only a few sigma above floor_g, which starves the curve of
    # usable levels and inflates the gradient-based size uncertainty below --
    # this makes that explicit instead of leaving it implicit in the errors.
    results['peak_to_floor_dynamic_range'] = (
        float(peak_g / floor_g) if floor_g > 0 else float('nan')
    )
    # print('Quick dynamic range = ', results['peak_to_floor_dynamic_range'])

    levels         = _build_levels(peak_g, floor_g, bins)
    if results['peak_to_floor_dynamic_range'] > 100:
        levels_ellipse = _build_levels(2*peak_g, floor_g, 64)
    else:
        levels_ellipse = _build_levels(peak_g, floor_g, 32)
    



    # -- Peak position & geometry ----------------------------------------------
    peak_position = np.unravel_index(np.nanargmax(g), g.shape)

    # Resolution element setting the annulus scale. `psf_fwhm_px` wins when
    # given -- without it an optical image falls through `beam_shape`'s
    # no-restoringbeam branch, which returns Omaj = Omin = cell_size, i.e. a
    # 1-pixel "beam" and a 1.5-3 *pixel* annulus.
    if psf_fwhm_px is not None:
        _peak_beam_px = float(psf_fwhm_px)
    else:
        try:
            _peak_beam_px = float(get_beam_size_px(img)[0])
        except Exception:
            _peak_beam_px = 1.0

    if residual is not None:
        results['peak_at_residual'] = res[peak_position]
        peak_stats                  = calculate_peak_error_annulus(
            g, res, peak_position, _peak_beam_px,
            inner_factor=1.5, outer_factor=3.0, min_beams=5,
            clip=False,
            systematic_error_fraction=systematic_error_fraction,
            peak_at_residual=results['peak_at_residual'],
            n_beams_in_mask=number_of_beams_in_mask,
        )
        results['peak_error']       = peak_stats['peak_error']
        results['peak_error_noise'] = peak_stats['error_noise']
        results['peak_error_pix']   = peak_stats['error_pix']
        results['peak_error_cal']   = peak_stats['error_cal']
        # peak_snr stays the legacy noise-only ratio (peak / sigma-clipped
        # annulus rms): it drives the Condon (1997)/Fomalont (1999) size-error
        # floor below and `theta_Rn_*_err` in morphen.py, both of which require
        # a noise-only SNR. Migrating them to `peak_snr_noise` is a separate,
        # deliberate decision -- not a side effect of adding a calibration term.
        results['peak_snr']            = peak_stats['snr']
        results['peak_snr_noise']      = peak_stats['snr_noise']
        results['peak_snr_meas']       = peak_stats['snr_meas']
        results['peak_of_flux_interp'] = peak_stats['peak_interp']
        results['peak_pixelisation_bias'] = peak_stats['pixelisation_bias']
        results['peak_max_bias_sigma'] = peak_stats['max_bias_sigma']
        results['peak_rms_local']         = peak_stats['rms_local']
        results['peak_rms_local_clipped'] = peak_stats['rms_local_clipped']
        results['peak_rms_local_std']     = peak_stats['rms_local_std']
        results['peak_local_nbeams']      = peak_stats['n_beams']
        results['peak_annulus_flag']      = peak_stats['flag']
    else:
        # No residual to measure a local noise on: the global `std` is the only
        # noise term available. The other two terms of the budget do not depend
        # on the residual, so they still apply.
        _k_pix = 4.0 * np.log(2.0) / max(_peak_beam_px, 1.0) ** 2
        _peak  = abs(results['peak_of_flux'])
        results['peak_error_noise'] = std
        results['peak_error_pix']   = _peak * (_k_pix / np.sqrt(90.0))
        results['peak_error_cal']   = systematic_error_fraction * _peak
        results['peak_error'] = float(np.sqrt(results['peak_error_noise']**2
                                              + results['peak_error_pix']**2
                                              + results['peak_error_cal']**2))
        results['peak_snr']       = results['peak_of_flux'] / std
        results['peak_snr_noise'] = results['peak_of_flux'] / std
        results['peak_snr_meas']  = results['peak_of_flux'] / results['peak_error']
        _pk_y, _pk_x = int(peak_position[0]), int(peak_position[1])
        _pk_interp = results['peak_of_flux']
        if 0 < _pk_y < g.shape[0] - 1 and 0 < _pk_x < g.shape[1] - 1:
            _, _px = _parabolic_peak_1d(float(g[_pk_y, _pk_x - 1]),
                                        float(results['peak_of_flux']),
                                        float(g[_pk_y, _pk_x + 1]))
            _, _py = _parabolic_peak_1d(float(g[_pk_y - 1, _pk_x]),
                                        float(results['peak_of_flux']),
                                        float(g[_pk_y + 1, _pk_x]))
            _pk_interp = results['peak_of_flux'] + (_px - results['peak_of_flux']) \
                                                 + (_py - results['peak_of_flux'])
        results['peak_of_flux_interp']    = _pk_interp
        results['peak_pixelisation_bias'] = _peak * (_k_pix / 6.0)
        results['peak_max_bias_sigma']    = (
            float(np.sqrt(2.0 * np.log(number_of_beams_in_mask)))
            if number_of_beams_in_mask > 1 else float('nan'))
        results['peak_rms_local']         = std
        results['peak_rms_local_clipped'] = std
        results['peak_rms_local_std']     = std
        results['peak_local_nbeams']      = float('nan')
        results['peak_annulus_flag']      = 'no_residual'

    try:
        results['convex_error_flag'] = False
        # do_plot=False here: the overlay is drawn later, onto ax2, only if
        # show_figure -- see the `if show_figure:` block below.
        cp = convex_morpho(g, mask, do_plot=False, n_boot=0, noise_std=6*std)
        results['PA_convex']                  = cp['PA_convex']
        results['PA_convex_err']              = cp['PA_convex_err']
        results['q_convex']                   = cp['q_convex']
        results['q_convex_err']               = cp['q_convex_err']
        results['centroid_convex']            = cp['centroid_convex']
        results['convex_major_diameter']      = cp['major_diameter']
        results['convex_major_diameter_err']  = cp['major_diameter_err']
        results['convex_minor_diameter']      = cp['minor_diameter']
        results['convex_minor_diameter_err']  = cp['minor_diameter_err']
        results['q_diam']                     = cp['q_diam']
        results['q_diam_err']                 = cp['q_diam_err']
        results['convex_axes_swapped']        = cp['axes_swapped']
    except Exception:
        results['convex_error_flag']          = True
        results['PA_convex']                  = 0.0
        results['PA_convex_err']              = float('nan')
        results['q_convex']                   = 0.0
        results['q_convex_err']               = float('nan')
        results['centroid_convex']            = [0, 0]
        results['convex_major_diameter']      = 0.0
        results['convex_major_diameter_err']  = float('nan')
        results['convex_minor_diameter']      = 0.0
        results['convex_minor_diameter_err']  = float('nan')
        results['q_diam']                     = 0.0
        results['q_diam_err']                 = float('nan')
        results['convex_axes_swapped']        = False
    if verbose > 1:
        if not results['convex_error_flag']:
            print(f"Convex hull PA: {results['PA_convex']:.2f} deg +/- {results['PA_convex_err']:.2f} deg  "
                  f"q: {results['q_convex']:.3f} +/- {results['q_convex_err']:.3f}  "
                  f"major: {results['convex_major_diameter']:.2f} +/- {results['convex_major_diameter_err']:.2f}px  "
                  f"minor: {results['convex_minor_diameter']:.2f} +/- {results['convex_minor_diameter_err']:.2f}px  "
                  f"q_diam: {results['q_diam']:.3f}")
            if results['convex_axes_swapped']:
                print("  (convex axes labelled by measured extent: the largest "
                      "extent lies along the second eigenvector, so PA_convex "
                      "is rotated 90 deg from the major-eigenvector direction)")
        else:
            print("Convex hull properties could not be calculated.")

    x0max, y0max = peak_center(g * mask)
    results['x0'], results['y0'] = x0max, y0max
    x0m, y0m, a_mom, b_mom, q_mom, PAdeg_mom = momenta(g * mask,
                                                        PArad_0=None, q_0=None)
    results['x0m'],      results['y0m']      = x0m,    y0m
    results['a_mom'],    results['b_mom']    = a_mom,  b_mom
    results['q_mom'],    results['PAdeg_mom'] = q_mom, PAdeg_mom
    if verbose > 1:
        print(f"Peak centre:   ({results['x0']:.2f}, {results['y0']:.2f})")
        print(f"Moment centre: ({results['x0m']:.2f}, {results['y0m']:.2f})")
        print(f"Moment axes:   a={results['a_mom']:.2f}px  b={results['b_mom']:.2f}px  "
              f"q={results['q_mom']:.2f}  PA={results['PAdeg_mom']:.2f} deg")

    vx = results['x0'] - results['x0m']
    vy = results['y0'] - results['y0m']
    try:
        TvPA, Tvlength = trail_vector(vx=vx, vy=vy, v0=np.asarray([1, 0]))
    except Exception:
        TvPA = Tvlength = 0.0
    results['TvPA']     = TvPA
    results['Tvlength'] = Tvlength
    
    results['peak_position_Ra_Dec'] = get_radec_from_pixel(image_input=img, 
                                                    pixel_coords = (results['x0'], results['y0']))
    results['moment_position_Ra_Dec'] = get_radec_from_pixel(image_input=img, 
                                                    pixel_coords = (results['x0m'], results['y0m']))
    if not results['convex_error_flag']:
        results['centroid_convex_position_Ra_Dec'] = get_radec_from_pixel(image_input=img, 
                                                        pixel_coords = results['centroid_convex'])
    else:
        results['centroid_convex_position_Ra_Dec'] = (0.0, 0.0)
    

    # --------------------------------------------------------------------------
    # Annular flux / area decomposition along intensity levels
    # --------------------------------------------------------------------------
    fluxes      = np.empty(len(levels))
    fluxes_err  = np.empty(len(levels))  # thermal (per-beam noise) term only
    fluxes_res  = np.empty(len(levels))  # non-thermal term: real residual/error map, or a systematic-fraction proxy
    areas       = np.empty(len(levels))

    # What spatially-resolved noise information do we actually have for the
    # non-thermal term? Real residual map (signed, matches how
    # total_flux_density_residual is defined above) takes priority; failing
    # that, whatever error/rms/variance/weight map the caller supplied
    # (resolve_flux_error_map already implements this precedence for the
    # Petrosian error path -- reused here rather than duplicated); only when
    # NEITHER is available do we fall back to a systematic_error_fraction
    # proxy per bin.
    _aux_noise_map = None
    if residual is None:
        _aux_noise_map, _ = resolve_flux_error_map(
            g, error_map=error_map, rms_map=rms_map, invvar_map=invvar_map,
            variance_map=variance_map, weight_map=weight_map, rms=None, verbose=False)
    _has_real_noise_map = (residual is not None) or (_aux_noise_map is not None)

    def _bin_aux_flux(cond, bin_flux):
        if residual is not None:
            return np.nansum(data_res * cond) / beam_area_
        if _aux_noise_map is not None:
            return np.nansum(_aux_noise_map * cond) / beam_area_
        return systematic_error_fraction * bin_flux

    # first bin: everything at or above levels[0] (the peak level)
    cond0       = g >= levels[0]
    fluxes[0]   = np.nansum(g * cond0) / beam_area_
    areas[0]    = np.nansum(cond0)
    fluxes_err[0] = std * np.sqrt(max(areas[0], 1) / beam_area_)
    fluxes_res[0] = _bin_aux_flux(cond0, fluxes[0])

    for i in range(1, len(levels)):
        cond           = (g < levels[i - 1]) & (g >= levels[i])
        fluxes[i]      = np.nansum(g * cond) / beam_area_
        areas[i]       = np.nansum(cond)
        fluxes_err[i]  = std * np.sqrt(max(areas[i], 1) / beam_area_)
        fluxes_res[i]  = _bin_aux_flux(cond, fluxes[i])

    agrow                     = areas.copy()
    results['total_flux_levels'] = np.nansum(fluxes)
    results['flux_levels_completeness'] = results['total_flux_levels'] / results['total_flux_mask']

    # -- Growth-curve uncertainty: two flavors -------------------------------
    # A globally-correlated flux-scale bias behaves very differently from
    # spatially-varying (thermal/residual) noise once it feeds a *local*
    # slope-based radius-uncertainty propagation (find_fractional_radius):
    #
    # - Lgrow_err_shape: thermal + non-thermal term, BOTH summed independently
    #   per annulus (bins-independent quadrature). This is what reflects
    #   uncertainty in WHERE a given flux fraction sits along the profile, so
    #   it is what drives C20/50/80/90radii_err. When real residual/error-map
    #   data exists, the non-thermal term is instead accumulated the same way
    #   total_flux_density_residual is (running signed sum, then abs) since
    #   it's then a genuine per-annulus flux measurement, not a proxy -- but
    #   the systematic_error_fraction proxy specifically must stay
    #   bins-independent, because cumsum(systematic_error_fraction*fluxes)
    #   collapses to exactly systematic_error_fraction*Lgrow -- i.e. the same
    #   linear-in-Lgrow term excluded below -- which is what caused the
    #   radius errors to blow up before this fix.
    # - Lgrow_err (total): additionally includes the global systematic term
    #   scaled to the cumulative flux itself. Correct for absolute-flux
    #   quantities (matches results['total_flux_error']) and the plotted
    #   uncertainty band, but not for radius uncertainty: a uniform flux-
    #   scale bias multiplies the whole profile and so does not move where a
    #   given fraction of the light falls.
    Lgrow              = np.nancumsum(fluxes)
    Lgrow_err_thermal   = np.sqrt(np.nancumsum(fluxes_err ** 2))
    if _has_real_noise_map:
        Lgrow_err_extra = np.abs(np.nancumsum(fluxes_res))
    else:
        Lgrow_err_extra = np.sqrt(np.nancumsum(fluxes_res ** 2))
    Lgrow_err_shape      = np.sqrt(Lgrow_err_thermal**2 + Lgrow_err_extra**2)
    Lgrow_err_systematic = systematic_error_fraction * Lgrow
    Lgrow_err = np.sqrt(Lgrow_err_shape**2 + Lgrow_err_systematic**2)

    radii = np.array([
        np.sqrt(np.nansum(g >= lev) / np.pi) for lev in levels
    ])

    # Save pre-dedup arrays: these share length = bins and align with `levels`.
    # radii_pre / Lgrow_pre are used for the gradient-based uncertainty below.
    Lgrow_norm_pre = np.nancumsum(fluxes) / np.nansum(fluxes)  # pre-dedup, len == bins
    Lgrow_pre      = Lgrow.copy()
    radii_pre      = radii.copy()  # pre-dedup, same indices as levels[]

    Lgrow, radii = check_flux_duplicates(Lgrow, radii)
    Lgrow_err        = Lgrow_err[:len(Lgrow)]
    Lgrow_err_shape  = Lgrow_err_shape[:len(Lgrow)]
    Lgrow_norm   = Lgrow / np.nansum(fluxes)
    # Lgrow_norm   = Lgrow / total_flux
    Lgrow_err_norm = Lgrow_err / np.nansum(fluxes)
    # Lgrow_err_norm = Lgrow_err / total_flux
    Lgrow_err_shape_norm = Lgrow_err_shape / np.nansum(fluxes)

    # -- Diagnostics for an external Monte Carlo validation of C*radii_err --
    # Purely additive (new dict keys only, no change to the positional return
    # tuple or to any existing key), so no caller can break. Exposes exactly
    # what a notebook needs to redraw per-bin flux realizations the same way
    # this function builds Lgrow_err_shape internally (thermal + non-thermal,
    # both bins-independent), then re-run check_flux_duplicates +
    # find_fractional_radius per realization and compare the empirical
    # scatter against C20/50/80/90radii_err below.
    results['mc_levels']              = levels.copy()
    results['mc_fluxes']              = fluxes.copy()          # pre-dedup, aligned with levels
    results['mc_fluxes_err_thermal']  = fluxes_err.copy()      # pre-dedup, aligned with levels
    results['mc_fluxes_err_extra']    = fluxes_res.copy()      # pre-dedup, aligned with levels
    results['mc_radii_pre']           = radii_pre.copy()       # pre-dedup, aligned with levels
    results['Lgrow_err_shape']        = Lgrow_err_shape        # post-dedup, aligned with returned Lgrow/radii
    results['Lgrow_err_shape_norm']   = Lgrow_err_shape_norm   # post-dedup, aligned with returned Lgrow/radii
    # print(Lgrow_norm[-1])
    # print((Lgrow / np.nansum(fluxes))[-1])
    # print(np.nansum(fluxes))
    # print(total_flux)

    # Number of distinct points the fractional-radius interpolation actually
    # has to work with, post-dedup. Low counts (few resolvable intensity
    # levels between peak_g and floor_g) mean the size/uncertainty estimates
    # below are based on a coarse curve -- see peak_to_floor_dynamic_range.
    results['growth_curve_n_levels'] = int(len(Lgrow))

    # --------------------------------------------------------------------------
    # Intensity thresholds at each flux fraction
    # (uses pre-dedup growth curve so indices align with the levels array)
    # --------------------------------------------------------------------------

    sigma_20 = _sigma_at_fraction(0.20, Lgrow_norm_pre, levels, floor_g)
    sigma_50 = _sigma_at_fraction(0.50, Lgrow_norm_pre, levels, floor_g)
    sigma_80 = _sigma_at_fraction(0.80, Lgrow_norm_pre, levels, floor_g)
    sigma_90 = _sigma_at_fraction(0.90, Lgrow_norm_pre, levels, floor_g)
    sigma_95 = _sigma_at_fraction(0.95, Lgrow_norm_pre, levels, floor_g)
    sigma_99 = _sigma_at_fraction(0.99, Lgrow_norm_pre, levels, floor_g)

    # backward-compat flags: True when the fraction wasn't reached
    flag20   = sigma_20 <= floor_g
    flag50   = sigma_50 <= floor_g
    flag9095 = sigma_95 <= floor_g
    flagL9095 = flag9095  # alias kept for v1 compatibility

    # Boolean masks at each sigma level
    g20 = ((g * mask) > sigma_20)
    g50 = ((g * mask) > sigma_50)
    g80 = ((g * mask) > sigma_80)
    g90 = ((g * mask) > sigma_90)
    g95 = ((g * mask) > sigma_95)
    g99 = ((g * mask) > sigma_99)

    # Flux values at each fraction level (for brightness temperature etc.)
    # Use pre-dedup Lgrow_pre: same length as Lgrow_norm_pre (= len(levels)).
    L20 = float(np.interp(0.20, Lgrow_norm_pre, Lgrow_pre))
    L50 = float(np.interp(0.50, Lgrow_norm_pre, Lgrow_pre))
    L80 = float(np.interp(0.80, Lgrow_norm_pre, Lgrow_pre))
    L90 = float(np.interp(0.90, Lgrow_norm_pre, Lgrow_pre))
    L95 = float(np.interp(0.95, Lgrow_norm_pre, Lgrow_pre))
    L99 = float(np.interp(0.99, Lgrow_norm_pre, Lgrow_pre))
    L20_norm = 0.20; L50_norm = 0.50; L80_norm = 0.80
    L90_norm = 0.90; L95_norm = 0.95

    inner_mask        = (g > sigma_50) * mask
    outer_mask90      = (g > sigma_90) * mask
    inner_perimeter   = perimeter_crofton(inner_mask, 4)
    outer_perimeter90 = perimeter_crofton(outer_mask90, 4)
    outer_perimeter   = perimeter_crofton(mask, 4)

    # -- Ellipse fitting -------------------------------------------------------
    if do_fit_ellipse:
        try:
            region_split = [i for i, x in enumerate(levels_ellipse > sigma_50) if x][-1]
            PA, q, x0col, y0col, PAm, qm, \
                PAmi, qmi, PAmo, qmo, \
                x0median, y0median, \
                x0median_i, y0median_i, \
                x0median_o, y0median_o, profiles = cal_PA_q(
                    g * mask, Isequence=levels_ellipse,
                    region_split=region_split,
                    SAVENAME=img.replace('.fits', '_ellipsefit'),
                )
        except Exception:
            PA = q = PAm = qm = PAmi = qmi = PAmo = qmo = 0.0
            x0col = y0col = x0max
            x0median = y0median = x0max
            x0median_i = y0median_i = x0max
            x0median_o = y0median_o = y0max
            profiles = None

        results.update({
            'PA': PA, 'q': q,
            'PAm': PAm, 'qm': qm,
            'PAmi': PAmi, 'qmi': qmi,
            'PAmo': PAmo, 'qmo': qmo,
            'x0m_i': x0median_i, 'y0m_i': y0median_i,
            'x0m_o': x0median_o, 'y0m_o': y0median_o,
        })



    # -- Brightness temperature ------------------------------------------------
    if not flag_freq:
        try:
            TB20 = T_B(omaj, omin, freq, L20)
        except Exception:
            TB20 = 0.0
        TB50 = T_B(omaj, omin, freq, L50)
        TB80 = T_B(omaj, omin, freq, L80)
        TB90 = T_B(omaj, omin, freq, L90)
        TB   = T_B(omaj, omin, freq, total_flux)
    else:
        TB20 = TB50 = TB80 = TB90 = TB = 0.0

    levels_20  = np.asarray([sigma_20])
    levels_50  = np.asarray([sigma_50])
    levels_80  = np.asarray([sigma_80])
    levels_90  = np.asarray([sigma_90])
    levels_95  = np.asarray([sigma_95])
    levels_99  = np.asarray([sigma_99])
    levels_3sigma = np.asarray([3 * std])

    # --------------------------------------------------------------------------
    # Fractional-flux radii - loop over fractions with safe fallback
    # --------------------------------------------------------------------------
    # Resolution/SNR floor for the radius error.
    # The thermal-noise term from find_fractional_radius is tiny for high-SNR
    # sources, so a floor is needed. The resolution element sets its scale:
    #   - Radio:   beam minor axis in pixels  (omin / cell_size)
    #   - Optical: PSF FWHM supplied via psf_fwhm_px (overrides beam when given)
    # A single formula covers both -- see the floor block inside the loop below.
    if psf_fwhm_px is not None:
        _res_pix = float(psf_fwhm_px)
    else:
        _res_pix = max(1.0, omin / cell_size) if not flag_freq else 1.0

    # Source SNR driving the floor. Deliberately a property of the *source*
    # (per sub-region, since `g` is already restricted to the region being
    # measured), not of the isophotal level -- see the floor block below.
    _snr_src = max(1.0, float(results.get('peak_snr', 1.0)))
    if not np.isfinite(_snr_src):
        _snr_src = 1.0
    _sigma_at_frac_map = {20: sigma_20, 50: sigma_50, 80: sigma_80,
                          90: sigma_90, 95: sigma_95, 99: floor_g}

    _FRACS = [(0.20, 20), (0.50, 50), (0.80, 80), (0.90, 90), (0.95, 95), (0.99, 99)]
    _C = {}  # tag -> (r, L, npix, r_err, npix_err)

    for frac, tag in _FRACS:
        fallback_r = radii[0] if frac < 0.5 else radii[-1]
        r, L_frac, npix, r_err, npix_err, ok = _fractional_radius_safe(
            Lgrow, radii, Lgrow_err_shape_norm, frac, fallback_r,
        )
        if not ok and logger is not None:
            logger.warning(f" !!==> find_fractional_radius failed for f={frac:.2f}."
                           f" Using fallback r={r:.2f}px.")
        # -- Gradient-based uncertainty ----------------------------------------
        # The Lgrow propagation (|dR/dF|*sigma_F) gives negligible values for
        # detected sources because the cumulative-flux error sigma_F is tiny
        # relative to the total flux.  The physically correct estimate is:
        #
        #   sigma_R_f = std / |dI/dR|_at_sigma_f
        #
        # "How far must the radius move to change the isophotal intensity by 1 rms?"
        # This naturally gives large errors for outer fractions (flat profile ≈ noise)
        # and small errors for inner fractions (steep profile).
        _sig_frac = _sigma_at_frac_map.get(tag, floor_g)
        _idx_g    = int(np.argmin(np.abs(levels - _sig_frac)))
        _idx_g    = int(np.clip(_idx_g, 1, len(levels) - 2))
        _dI_g     = abs(float(levels[_idx_g - 1]) - float(levels[_idx_g + 1]))
        _dR_g     = abs(float(radii_pre[_idx_g + 1]) - float(radii_pre[_idx_g - 1]))
        if _dR_g > 0.01:
            _sigma_R_grad = float(std) * _dR_g / max(_dI_g, 1e-30)
        else:
            _sigma_R_grad = 0.0   # flat profile in this bin - gradient dominates elsewhere

        # -- Resolution / SNR floor -------------------------------------------
        # Hard minimum on how well a radius can be localised, in two pieces:
        #   - 0.5 px, the pixelisation limit;
        #   - resolution element / (2 * source SNR), the standard positional /
        #     size uncertainty scaling (Condon 1997; Fomalont 1999), clamped so
        #     the floor can never exceed half a resolution element.
        # Identical for radio and optical: `_res_pix` is the beam minor axis or
        # the supplied PSF FWHM, and nothing else in the expression is
        # instrument-specific.
        #
        # This uses peak_snr rather than the isophotal ratio sigma_f/std that
        # earlier versions used. sigma_f is bounded below by
        # floor_g = last_level*std by construction, so sigma_f/std saturates at
        # ~3-5 for *every* source regardless of brightness; combined with a
        # 1/sqrt scaling that pinned C95radii_err -- and always C99radii_err,
        # whose sigma_f is exactly floor_g -- at the constant 0.5*_res_pix,
        # identical for a peak_snr=7 sub-region and a peak_snr=200 one.
        #
        # The fraction dependence of the error is carried by _sigma_R_grad
        # above, which grows naturally for outer fractions as the profile
        # flattens; it does not belong in a floor.
        _r_err_floor = max(0.5, min(0.5 * _res_pix, _res_pix / (2.0 * _snr_src)))

        r_err    = max(r_err, _sigma_R_grad, _r_err_floor)
        npix_err = 2 * np.pi * r * r_err
        _C[tag] = (r, L_frac, npix, r_err, npix_err)

    # Unpack for readability and downstream compatibility
    C20radii, L20, npix20, C20radii_err, npix20_err = _C[20]
    C50radii, L50, npix50, C50radii_err, npix50_err = _C[50]
    C80radii, L80, npix80, C80radii_err, npix80_err = _C[80]
    C90radii, L90, npix90, C90radii_err, npix90_err = _C[90]
    C95radii, L95, npix95, C95radii_err, npix95_err = _C[95]
    C99radii, L99, npix99, C99radii_err, npix99_err = _C[99]

    A20 = npix20 / beam_area_;  A20_err = npix20_err / beam_area_
    A50 = npix50 / beam_area_;  A50_err = npix50_err / beam_area_
    A80 = npix80 / beam_area_;  A80_err = npix80_err / beam_area_
    A90 = npix90 / beam_area_;  A90_err = npix90_err / beam_area_
    A95 = npix95 / beam_area_;  A95_err = npix95_err / beam_area_
    A99 = npix99 / beam_area_;  A99_err = npix99_err / beam_area_

    # --------------------------------------------------------------------------
    # Ellipticity-corrected semi-major-axis radii  (new in v2)
    #
    # The levels-based Cxx radii are equivalent circular radii of isophotal
    # areas, which for an elliptical source approximate the geometric mean of
    # the semi-axes:  Cxx ≈ a_xx * sqrt(q).
    # Dividing by sqrt(q) recovers the semi-major-axis effective radius.
    # See experiments_radial_quantities.py – Experiments 1 & 3.
    # --------------------------------------------------------------------------
    _q_convex_ok = (not results.get('convex_error_flag', True)
                    and results.get('q_convex', 0.0) > 0.05)
    q_eff    = results['q_convex'] if _q_convex_ok else results['q_mom']
    q_eff    = float(np.clip(q_eff, 0.05, 1.0))
    sqrt_q   = np.sqrt(q_eff)
    results['q_eff'] = q_eff

    # Empirical sigma(q_eff) from Exp D noise-perturbation MC (see experiments_uncertainty.py).
    # sigma(q) ≈ 0.02 / SNR_eff, calibrated at q=0.5; conservative floor of 0.005.
    snr_eff       = peak_g / max(std, 1e-30)
    sigma_q_eff   = float(np.clip(0.02 / max(snr_eff, 1.0), 0.005, 0.5))
    results['sigma_q_eff'] = sigma_q_eff

    for _tag, _r, _r_err in [
        (20, C20radii, C20radii_err),
        (50, C50radii, C50radii_err),
        (80, C80radii, C80radii_err),
        (90, C90radii, C90radii_err),
        (95, C95radii, C95radii_err),
        (99, C99radii, C99radii_err),
    ]:
        _r_major       = _r / sqrt_q
        # sigma^2(R_major) = (sigma_R/sqrt q)^2 + (R_major/(2 x q)  x  sigma_q)^2
        _term_r        = _r_err / sqrt_q
        _term_q        = _r_major / (2.0 * q_eff) * sigma_q_eff
        _r_major_err   = float(np.sqrt(_term_r**2 + _term_q**2))
        results[f'C{_tag}radii_major']     = _r_major
        results[f'C{_tag}radii_major_err'] = _r_major_err

    # -- Convex-hull shapes at each level -------------------------------------
    try:
        results['conv_P20'], results['conv_A20'] = convex_shape(g20)
    except Exception:
        results['conv_P20'], results['conv_A20'] = 5, 5

    try:
        results['conv_P50'], results['conv_A50'] = convex_shape(g50)
        results['conv_P80'], results['conv_A80'] = convex_shape(g80)
        results['conv_P90'], results['conv_A90'] = convex_shape(g90)
        results['conv_P95'], results['conv_A95'] = convex_shape(g95)
        results['conv_PT'],  results['conv_AT']  = convex_shape(mask)
    except Exception:
        for _k in ['conv_P50', 'conv_A50', 'conv_P80', 'conv_A80',
                   'conv_P90', 'conv_A90', 'conv_P95', 'conv_A95',
                   'conv_PT',  'conv_AT']:
            results[_k] = 5

    # -- Radial statistics at each level --------------------------------------
    for _label, _mask in [('20', g20), ('50', g50), ('80', g80),
                           ('90', g90), ('95', g95), ('T', mask)]:
        try:
            med, mean, std_ = calculate_radii(g, _mask)
        except Exception:
            med = mean = std_ = 1
        results[f'R{_label}med']  = med
        results[f'R{_label}mean'] = mean
        results[f'R{_label}std']  = std_

    # -- Concentration indices -------------------------------------------------
    C1, C1_err, C2, C2_err = concentration_index_errors(
        C80radii, C20radii, C90radii, C50radii,
        C80radii_err, C20radii_err, C90radii_err, C50radii_err,
    )
    AC1  = np.log10(A80 / A20)
    AC2  = np.log10(A90 / A50)
    CAC1 = np.log10(area_to_radii(results['conv_A80']) / area_to_radii(results['conv_A20']))
    CAC2 = np.log10(area_to_radii(results['conv_A90']) / area_to_radii(results['conv_A50']))

    # -- Area totals & shell quantities ----------------------------------------
    area_total,   Cradii,   npix_total   = estimate_area(mask,  cell_size, omaj, omin)
    o_area_total, o_Cradii, o_npix_total = estimate_area(omask, cell_size, omaj, omin)

    mask_outer      = (g < sigma_50) & (g > floor_g) * mask
    A50_100, C50_100radii, npix50_100 = estimate_area(mask_outer, cell_size, omaj, omin)
    mask_outer_full = (g < sigma_50) * mask
    A50_full, C50_full_radii, npix50_full = estimate_area(mask_outer_full, cell_size,
                                                           omaj, omin)

    gaussianity      = sigma_50 / g.max()
    gaussianity_L50  = g.max() / L50
    radii_ratio      = C50radii / C50_100radii
    radii_ratio_full = C50radii / C50_full_radii
    area_ratio       = A50 / A50_100
    area_ratio_full  = A50 / A50_full

    # -- Populate results dict -------------------------------------------------
    results['max']       = np.nanmax(g)
    results['std_image'] = std
    try:
        results['std_residual'] = mad_std(load_fits_data(residual))
    except Exception:
        pass

    results.update({
        'L50': L50, 'sigma_50': sigma_50,
        'sigma_90': sigma_90, 'sigma_95': sigma_95,
        'o_area_total': o_area_total, 'o_Cradii': o_Cradii, 'o_npix_total': o_npix_total,
        'area_total': area_total, 'Cradii': Cradii, 'npix_total': npix_total,
        'inner_perimeter': inner_perimeter,
        'outer_perimeter': outer_perimeter,
        'outer_perimeter90': outer_perimeter90,
        'TB20': TB20, 'TB50': TB50, 'TB80': TB80, 'TB': TB,
        'C1': C1, 'C2': C2, 'C1_err': C1_err, 'C2_err': C2_err,
        'AC1': AC1, 'AC2': AC2, 'CAC1': CAC1, 'CAC2': CAC2,
        'L20': L20, 'A20': A20, 'C20': sigma_20 / std,
        'C20radii': C20radii, 'npix20': npix20,
        'C20radii_err': C20radii_err, 'npix20_err': npix20_err,
        'L50': L50, 'A50': A50, 'C50': sigma_50 / std,
        'C50radii': C50radii, 'C50radii_err': C50radii_err,
        'npix50': npix50, 'npix50_err': npix50_err,
        'C50theta':     C50radii * np.sqrt(6 * np.log(2)),
        'C50theta_err': C50radii_err * np.sqrt(6 * np.log(2)),
        'C50theta_dec': np.sqrt(
            (C50radii * np.sqrt(6 * np.log(2)))**2
            - (results['bmajor'] / cell_size)**2
        ),
        'C50theta_dec_err': (
            C50radii * np.sqrt(6 * np.log(2))
            * C50radii_err * np.sqrt(6 * np.log(2))
        ),
        'L80': L80, 'A80': A80, 'C80': sigma_80 / std,
        'C80radii': C80radii, 'npix80': npix80,
        'C80radii_err': C80radii_err, 'npix80_err': npix80_err,
        'L90': L90, 'A90': A90, 'C90': sigma_90 / std,
        'C90radii': C90radii, 'npix90': npix90,
        'C90radii_err': C90radii_err, 'npix90_err': npix90_err,
        'flag20': flag20, 'flag50': flag50,
        'flag9095': flag9095, 'flagL9095': flagL9095,
        'L95': L95, 'A95': A95, 'C95': sigma_95 / std,
        'C95radii': C95radii, 'npix95': npix95,
        'C95radii_err': C95radii_err, 'npix95_err': npix95_err,
        'L99': L99, 'A99': A99,
        'C99radii': C99radii, 'npix99': npix99,
        'C99radii_err': C99radii_err, 'npix99_err': npix99_err,
        'A20_err': A20_err, 'A50_err': A50_err,
        'A80_err': A80_err, 'A90_err': A90_err,
        'A95_err': A95_err, 'A99_err': A99_err,
        'gaussianity': gaussianity, 'gaussianity_L50': gaussianity_L50,
        'A50_100': A50_100, 'C50_100radii': C50_100radii, 'npix50_100': npix50_100,
        'A50_full': A50_full, 'C50_full_radii': C50_full_radii, 'npix50_full': npix50_full,
        'radii_ratio': radii_ratio, 'radii_ratio_full': radii_ratio_full,
        'area_ratio': area_ratio, 'area_ratio_full': area_ratio_full,
        'beam_area': beam_area_,
    })

    if redshift is not None:
        pix_to_pc    = pixsize_to_pc(z=redshift, cell_size=cell_size)
        scale_size       = 1 * pix_to_pc
        scale_size_units = 'pc'
    else:
        scale_size       = 1 * cell_size
        scale_size_units = scale_units

    if verbose > 0:
        _log = logger.debug if logger else print
        _log(f" ==>  Peak of Flux      = "
             f"{results['peak_of_flux']*flux_conversion:.3f}"
             f" +/- {results['peak_error']*flux_conversion:.3f} [mJy/beam]"
             f" (SNR={results['peak_snr']:.1f})")
        _log(f" ==>  Snu (within mask) = "
             f"{results['total_flux_mask']*flux_conversion:.3f}"
             f" +/- {results['total_flux_error']*flux_conversion:.3f} [mJy]")
        _log(f" ==>  Snu (image)       = "
             f"{results['total_flux_nomask']*flux_conversion:.3f} (no independent error estimate)")
        _log(f" ==>  R50               = "
             f"{cell_size*C50radii:.4f} +/- {cell_size*C50radii_err:.4f} [{scale_units}]")
        _log(f" ==>  R50_major         = "
             f"{cell_size*results['C50radii_major']:.4f}"
             f" +/- {cell_size*results['C50radii_major_err']:.4f} [{scale_units}]"
             f"  (q_eff={q_eff:.2f})")
        _log(f" ==>  R95               = "
             f"{cell_size*C95radii:.4f} +/- {cell_size*C95radii_err:.4f} [{scale_units}]")
        _log(f" ==>  R95_major         = "
             f"{cell_size*results['C95radii_major']:.4f}"
             f" +/- {cell_size*results['C95radii_major_err']:.4f} [{scale_units}]"
             f"  (q_eff={q_eff:.2f})")
        if redshift is not None:
            _log(f" ==>  R50 (phys)        = "
                 f"{scale_size*C50radii:.4f}"
                 f" +/- {scale_size*C50radii_err:.4f} [{scale_size_units}]")
            _log(f" ==>  R50_major (phys)  = "
                 f"{scale_size*results['C50radii_major']:.4f}"
                 f" +/- {scale_size*results['C50radii_major_err']:.4f} [{scale_size_units}]")
            _log(f" ==>  R95 (phys)        = "
                 f"{scale_size*C95radii:.4f}"
                 f" +/- {scale_size*C95radii_err:.4f} [{scale_size_units}]")
            _log(f" ==>  R95_major (phys)  = "
                 f"{scale_size*results['C95radii_major']:.4f}"
                 f" +/- {scale_size*results['C95radii_major_err']:.4f} [{scale_size_units}]")

    # --------------------------------------------------------------------------
    # Petrosian properties (unchanged from v1)
    # --------------------------------------------------------------------------
    error_petro = False
    if do_petro:
        try:
            if verbose >= 1:
                (logger.info if logger else print)('++>> Computing Petrosian properties.')
            _residual_for_error = data_res if (use_residual_as_error and error_map is None
                                                and rms_map is None and invvar_map is None
                                                and variance_map is None and weight_map is None) else None

            error_arr_input, _error_source = resolve_flux_error_map(
                g, error_map=error_map, rms_map=rms_map, invvar_map=invvar_map,
                variance_map=variance_map, weight_map=weight_map,
                residual_map=_residual_for_error, rms=None, verbose=(verbose > 0))

            r_list, area_arr, area_beam, p, flux_arr, error_arr, results, cat, \
                sorted_idx_list, segm, segm_deblend = \
                compute_petrosian_properties(
                    g, img,
                    mask_component=mask_component,
                    global_mask=mask,
                    source_props=results,
                    apply_mask=False,
                    error=error_arr_input,
                    sigma_level=sigma_mask,
                    bkg_to_sub=bkg_to_sub,
                    vmin=vmin_factor,
                    show_figure=show_figure,
                    plot=show_figure,
                    verbose=verbose,
                    add_save_name=add_save_name,
                    logger=logger,
                )
        except Exception as e:
            (logger.warning if logger else print)(
                f"-->> ERROR in Petrosian properties: {e}")
            error_petro = True
    else:
        error_petro = True
    results['error_petro'] = error_petro

    # --------------------------------------------------------------------------
    # Optional cropping of the image to a box around the source
    # --------------------------------------------------------------------------
    # mask_for_plot tracks mask cropped the same way as g/g_original below, so
    # the convex_morpho overlay (drawn later, on the cropped ax2 image) uses a
    # mask that matches g's shape -- mask itself is left uncropped since it is
    # part of this function's return value.
    mask_for_plot = mask
    if crop:
        try:
            xin, xen, yin, yen = do_cutout_2D(img, box_size=box_size,
                                               center=None,
                                               centre_mode='image_centre',
                                               return_='box')
            g          = g[xin:xen, yin:yen]
            g_original = g_original[xin:xen, yin:yen]
            mask_for_plot = mask[xin:xen, yin:yen]
        except Exception:
            try:
                max_x, max_y = np.where(g == np.nanmax(g))
                xin = max_x[0] - box_size;  xen = max_x[0] + box_size
                yin = max_y[0] - box_size;  yen = max_y[0] + box_size
                g          = g[xin:xen, yin:yen]
                g_original = g_original[xin:xen, yin:yen]
                mask_for_plot = mask[xin:xen, yin:yen]
            except Exception:
                pass

    # --------------------------------------------------------------------------
    # Plot
    # --------------------------------------------------------------------------
    if show_figure:
        img_height, img_width = g_original.shape
        img_aspect  = img_width / img_height
        plot_aspect = figsize[0] / figsize[1]
        height_ref  = figsize[1]

        img_width_needed  = height_ref * img_aspect
        plot_width_needed = height_ref * plot_aspect
        total_width       = img_width_needed + plot_width_needed
        width_ratios      = [plot_width_needed / total_width,
                              img_width_needed  / total_width]
        figsize = (total_width * 1.0, height_ref)

        fig = matplotlib.figure.Figure(figsize=figsize)
        FigureCanvasAgg(fig)
        gs  = gridspec.GridSpec(1, 2, figure=fig, width_ratios=width_ratios, wspace=0.1)
        ax1 = fig.add_subplot(gs[0, 0])
        ax2 = fig.add_subplot(gs[0, 1])

        # ax1.scatter(radii * cell_size, Lgrow / results['total_flux_mask'])
        ax1.scatter(radii * cell_size, Lgrow_norm, color='purple', s=20, 
                    # label='Cumulative Flux'
                    )
        
        ax1.fill_between(radii * cell_size, Lgrow_norm - 1*Lgrow_err_norm,
                    Lgrow_norm + 1*Lgrow_err_norm, color='purple',
                    alpha=0.35, 
                    # label='Cumulative Flux Uncertainty'
                    )

        ax1.axvspan(C50radii * cell_size - C50radii_err * cell_size,
                    C50radii * cell_size + C50radii_err * cell_size,
                    alpha=0.3, color='grey')
        ax1.axvline(C50radii * cell_size,
                    label=(r"$R_{50}=~$"
                           f"{C50radii*cell_size:0.3f}"
                           r"$\pm$"
                           f"{C50radii_err*cell_size:0.3f}$''$"),
                    ls='-.', color='lime', lw=4)
        ax1.axhline(L50_norm, ls='-.', color='lime', lw=4)

        ax1.axvspan(C95radii * cell_size - C95radii_err * cell_size,
                    C95radii * cell_size + C95radii_err * cell_size,
                    alpha=0.3, color='grey')
        ax1.axvline(C95radii * cell_size,
                    label=(r"$R_{95}=~$"
                           f"{C95radii*cell_size:0.3f}"
                           r"$\pm$"
                           f"{C95radii_err*cell_size:0.3f}$''$"),
                    color='#4daf4a', lw=3)
        # ax1.axvline(C99radii * cell_size,
        #             label=(r"$R_{99}=~$"
        #                    f"{C99radii*cell_size:0.3f}"
        #                    r"$\pm$"
        #                    f"{C99radii_err*cell_size:0.3f}$''$"),
        #             color="olive", lw=3)
        # ax1.axvspan(C99radii * cell_size - C99radii_err * cell_size,
        #             C99radii * cell_size + C99radii_err * cell_size,
        #             alpha=0.3, color='olive')

        # ax1.set_title("Integrated Flux Density\n"
        #               r"$S_{\nu} =$ "
        ax1.set_title(r"Integrated Flux $=$ "
                      f"{flux_conversion*total_flux:.3f}"
                      r" $\pm$ "
                      f"{flux_conversion*total_flux_density_error:.3f} [mJy]")
        ax1.set_xlabel(fr'Projected Circular Radius $R$ [{scale_units}]')
        ax1.set_ylabel(r"Normalised FGC $S_{\nu}(\leq R)$")
        ax1.grid(alpha=0.5)
        ax1.legend(framealpha=0.9, handlelength=1, handletextpad=0.5,
                   columnspacing=0.5, borderaxespad=0.1)
        ax1.set_ylim(0.0, 1.0)
        # ax1.axhline(1.0,color='black',lw=2)
        # ax1.semilogx()
        # ax1.semilogy()

        vmin_plot = vmin_factor * std
        if vmax is None:
            vmax = vmax_factor * np.nanmax(g)

        # Offset coordinates centered on the image, matching eimshow's
        # "Offset [<units>]" axes -- extent must be passed consistently to
        # imshow/contour below so everything on ax2 lines up in the same
        # data space.
        ny_ax2, nx_ax2 = g_original.shape
        dx_ax2 = nx_ax2 / 2.0 * cell_size
        dy_ax2 = ny_ax2 / 2.0 * cell_size
        extent_ax2 = [-dx_ax2, dx_ax2, -dy_ax2, dy_ax2]

        try:
            norm   = simple_norm(g_original, stretch='asinh', asinh_a=0.075,
                                 vmin=vmin_plot, vmax=vmax)
            ax2.imshow(g_original, cmap='magma_r', origin='lower', norm=norm,
                       aspect='equal', extent=extent_ax2)
        except Exception:
            ax2.imshow(g_original, cmap='magma_r', origin='lower', aspect='equal',
                       extent=extent_ax2)

        if image_title is not None:
            ax2.set_title(image_title)

        try:
            ax2.contour(g, levels=levels_50,        colors='lime',   linewidths=2.5,
                        linestyles='-.', alpha=1.0, extent=extent_ax2)
            ax2.contour(g, levels=levels_95,        colors='#4daf4a', linewidths=2.0,
                        extent=extent_ax2)
            ax2.contour(g, levels=levels_99,        colors='olive', linewidths=2.0,
                        linestyles='dashed',
                        # linestyle=(0, (3, 5, 1, 5)),
                        extent=extent_ax2)
            ax2.contour(g, levels=[floor_g],        colors='cyan',   linewidths=0.6,
                        extent=extent_ax2)
            ax2.contour(g, levels=[6.0 * std],      colors='black',  linewidths=1.5,
                        alpha=0.9, extent=extent_ax2)
            ax2.contour(g, levels=[3.0 * std],      colors='brown',  linewidths=1.2,
                        alpha=0.9, extent=extent_ax2)
        except Exception:
            print('Not plotting contours!')

        # Merge the convex-morpho overlay (major/minor axes, centroid, hull)
        # into this same ax2 panel rather than popping up its own figure --
        # only convex_morpho() called standalone (ax=None) gets its own plot.
        # cell_size makes it plot in the same offset coordinates as ax2 above.
        if not results.get('convex_error_flag', True):
            try:
                convex_morpho(g, mask_for_plot, do_plot=True, ax=ax2,
                              n_boot=0, noise_std=6 * std, cell_size=cell_size)
            except Exception:
                pass

        ax2.set_xlabel(fr'Offset [{scale_units}]', labelpad=3)
        # ax2.set_ylabel(fr'Offset [{scale_units}]', labelpad=3)
        ax2.set_xlim(-dx_ax2, dx_ax2)
        ax2.set_ylim(-dy_ax2, dy_ax2)
        ax2.grid(which='both', axis='both', color='gray', linewidth=0.6, alpha=0.4)

        try:
            if SAVE:
                _save = (plot_savemane if plot_savemane is not None
                         else img.replace('.fits', '_Lgrow_levels') + add_save_name + ext)
                fig.savefig(_save, dpi=300, bbox_inches='tight')
            show_or_display(fig, dpi=100)
            del fig
            gc.collect()
        except Exception as e:
            print(f"Error saving or displaying figure: {e}")

    # --------------------------------------------------------------------------
    # Optional CSV save
    # --------------------------------------------------------------------------
    if save_csv:
        import csv
        _csvpath = img.replace('.fits', '_image_properties') + add_save_name + '.csv'
        with open(_csvpath, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=results.keys())
            writer.writeheader()
            writer.writerow(results)

    return (levels, fluxes, Lgrow, Lgrow_err, Lgrow_norm, Lgrow_err_norm, radii, agrow, omask, mask, results)


def _parabolic_peak_1d(f_minus, f_0, f_plus):
    """
    Sub-pixel vertex of the parabola through three equally spaced samples.

    Returns (offset, peak_value) with `offset` in pixels from the centre sample.
    Falls back to (0.0, f_0) when the three points are not concave (no maximum)
    or the vertex lands outside the central pixel -- both signs that the
    quadratic approximation does not hold here.
    """
    denom = f_minus - 2.0 * f_0 + f_plus
    if not np.isfinite(denom) or denom >= 0:
        return 0.0, f_0
    d = 0.5 * (f_minus - f_plus) / denom
    if not np.isfinite(d) or abs(d) > 1.0:
        return 0.0, f_0
    return d, f_0 - 0.25 * (f_minus - f_plus) * d


def calculate_peak_error_annulus(image, residual, peak_position, beam_size_pix,
                                 inner_factor=1.5, outer_factor=3.0,
                                 min_beams=5, clip=False,
                                 systematic_error_fraction=0.05,
                                 peak_at_residual=None,
                                 n_beams_in_mask=None):
    """
    Peak flux density and its full error budget.

    The error is built as

        peak_error = sqrt(sigma_local^2 + sigma_pix^2 + (f_cal * S_peak)^2)

    Each term, and why it is there:

    ``sigma_local`` -- local noise, from a `inner_factor`-`outer_factor` beam
        annulus around the peak, measured on the **residual** map. This is the
        correct noise term for a peak: Condon (1997, PASP 109, 166) eq. 41 for a
        fitted amplitude reduces exactly to sigma(A) = mu (the local rms) for an
        unresolved source, and to 1.01 mu at twice the beam. Note that the same
        formula gives *less* than mu for very extended sources (0.53 mu at
        5 theta_beam) because it describes a fit that averages noise over many
        beams -- that reduction does NOT apply here, where the peak is a single
        pixel maximum, so the per-pixel local rms is what is used.

        Three refinements over a plain `mad_std` of the annulus:
          - no sigma-clipping by default (`clip=False`). Near a bright source
            the local noise is dynamic-range-limited -- calibration sidelobes,
            not thermal noise -- and 3-sigma clipping deletes exactly those.
            `mad_std` is already robust; clipping on top of it only removes
            signal (up to ~40% low with 10% artefact pixels in the annulus).
            `rms_local_clipped` is still returned for comparison.
          - floored at |residual at the peak pixel|, a direct lower bound on the
            local model error, passed in as `peak_at_residual`.
          - inflated by sqrt(1 + 1/(2 N_beams)) for the finite-sample
            uncertainty of the rms estimate itself (~16% for the default
            annulus, which holds ~19 independent beams).

    ``sigma_pix`` -- peak sampling. `peak_of_flux` is a pixel maximum, and the
        true peak almost never sits on a pixel centre. For a Gaussian of FWHM
        theta (px) and a sub-pixel offset uniform over the central pixel, the
        fractional loss is k*(dx^2+dy^2) with k = 4 ln2 / theta^2, giving a mean
        bias k/6 and an rms k/sqrt(90) -- 1.8% and 1.1% at 5 px/beam, 5.0% and
        3.1% at 3 px/beam. The rms enters the budget; the bias is reported as
        `pixelisation_bias` but NOT applied, since `peak_of_flux` must keep its
        meaning. `peak_interp` gives the sub-pixel corrected value instead
        (3-point parabolic interpolation in x and y, the AIPS MAXFIT approach).

    ``f_cal * S_peak`` -- absolute flux-density scale. Standard practice for any
        reported flux density: 3% for NVSS (Condon et al. 1998), 5% for
        VLA-COSMOS 3 GHz (Smolcic et al. 2017), 10% for LoTSS (Shimwell et al.
        2019, 2022); the VLA scale itself is 1-3% below 12 GHz rising to ~5% at
        Ka/Q (Perley & Butler 2017). This is the term whose absence made peak
        errors come out an order of magnitude below the integrated-flux error on
        the same row. It reuses `systematic_error_fraction`, the same knob that
        already sets the integrated-flux systematic, so the two cannot diverge.

    Parameters
    ----------
    image : ndarray
        Science image.
    residual : ndarray
        Residual image (post-cleaning). The noise is measured on this.
    peak_position : tuple
        (y, x) position of the peak in pixels.
    beam_size_pix : float
        Resolution element FWHM in pixels -- the beam for radio data, or the PSF
        FWHM for optical (pass `psf_fwhm_px`). Sets the annulus scale; a value
        of 1 px means no beam was found in the header and the annulus is only a
        few pixels wide, which is flagged.
    inner_factor, outer_factor : float, optional
        Annulus radii in units of `beam_size_pix`. Defaults 1.5 and 3.0.
    min_beams : int, optional
        Minimum number of independent beams wanted in the annulus (default 5).
        Below this the rms estimate is too noisy to trust and the outer radius
        is expanded.
    clip : bool, optional
        Sigma-clip the annulus before taking `mad_std`. Default False; see
        `sigma_local` above.
    systematic_error_fraction : float, optional
        Fractional flux-scale calibration error, default 0.05.
    peak_at_residual : float, optional
        Residual value at the peak pixel, used as a floor on `sigma_local`.
    n_beams_in_mask : float, optional
        Beams spanned by the region being measured. Only used to report
        `max_bias_sigma`.

    Returns
    -------
    dict
        'peak_flux'         peak flux density (the pixel maximum)
        'peak_error'        full budget (the three terms in quadrature)
        'error_noise'       sigma_local
        'error_pix'         sigma_pix
        'error_cal'         f_cal * S_peak
        'rms_local'         annulus rms as used (unclipped unless clip=True)
        'rms_local_clipped' annulus rms with 3-sigma clipping, for comparison
        'rms_local_std'     plain std of the annulus -- diagnostic only; a large
                            rms_local_std/rms_local means the residual there is
                            artefact-dominated and the robust estimator hides it
        'snr'               peak / rms_local_clipped -- the LEGACY, noise-only
                            SNR. Deliberately unchanged: it feeds the Condon
                            (1997)/Fomalont (1999) size-error floors, which
                            require a noise-only SNR and would silently shrink
                            if a calibration term were folded in.
        'snr_noise'         peak / error_noise (the improved noise-only SNR)
        'snr_meas'          peak / peak_error
        'peak_interp'       sub-pixel interpolated peak
        'pixelisation_bias' mean peak loss from pixel sampling (not applied)
        'max_bias_sigma'    sqrt(2 ln N_beams), the expected noise maximum over
                            the region -- a positive bias on the peak of faint
                            sub-regions (Eddington-type; Condon 1997)
        'n_pixels'          pixels in the annulus
        'n_beams'           independent beams in the annulus
        'inner_radius', 'outer_radius'  radii used (pixels)
        'flag'              'ok', or '+'-joined subset of 'expanded',
                            'edge_truncated', 'few_beams', 'no_beam_scale'
    """
    from astropy.stats import mad_std, sigma_clip

    ny, nx = residual.shape
    y_c, x_c = int(peak_position[0]), int(peak_position[1])
    peak_flux = float(image[y_c, x_c])

    flags = []
    beam_size_pix = float(beam_size_pix)
    if not np.isfinite(beam_size_pix) or beam_size_pix <= 1.0:
        # No restoring beam in the header (`beam_shape` falls back to one pixel
        # per resolution element) and no psf_fwhm_px supplied. The annulus is
        # then only a few pixels wide -- usable, but say so.
        beam_size_pix = max(beam_size_pix, 1.0)
        flags.append('no_beam_scale')

    # Gaussian beam area, pi/(4 ln2) * FWHM^2. The previous version used a
    # circle of diameter FWHM (pi/4 * FWHM^2), so `min_beams=5` actually asked
    # for ~3.5 beams.
    beam_area_px = (np.pi / (4.0 * np.log(2.0))) * beam_size_pix ** 2

    inner_radius = inner_factor * beam_size_pix
    outer_radius = outer_factor * beam_size_pix

    y, x = np.ogrid[:ny, :nx]
    r = np.sqrt((x - x_c) ** 2 + (y - y_c) ** 2)

    min_pixels = int(min_beams * beam_area_px)

    annulus_mask = (r >= inner_radius) & (r <= outer_radius)
    n_pixels = int(np.sum(annulus_mask))

    # An annulus running off the edge of a cutout is silently smaller and
    # one-sided; compare against the area it should have had.
    _expected = np.pi * (outer_radius ** 2 - inner_radius ** 2)
    if n_pixels < 0.9 * _expected:
        flags.append('edge_truncated')

    if n_pixels < min_pixels:
        import warnings
        warnings.warn(
            f"Annulus has only {n_pixels} pixels (< {min_pixels}). "
            f"Expanding outer radius for better statistics."
        )
        flags.append('expanded')
        while n_pixels < min_pixels and outer_radius < min(residual.shape) / 2:
            outer_radius *= 1.2
            annulus_mask = (r >= inner_radius) & (r <= outer_radius)
            n_pixels = int(np.sum(annulus_mask))

    n_beams = n_pixels / beam_area_px if beam_area_px > 0 else 0.0
    if n_beams < min_beams:
        flags.append('few_beams')

    annulus_pixels = residual[annulus_mask]
    annulus_pixels = annulus_pixels[np.isfinite(annulus_pixels)]
    if annulus_pixels.size == 0:
        rms_raw = rms_clipped = rms_std = float('nan')
    else:
        rms_raw = float(mad_std(annulus_pixels, ignore_nan=True))
        compressed = sigma_clip(annulus_pixels, sigma=3, maxiters=5).compressed()
        rms_clipped = (float(mad_std(compressed, ignore_nan=True))
                       if len(compressed) > 0 else rms_raw)
        # Plain standard deviation, reported as a diagnostic only. `mad_std` is
        # robust by construction, so it under-reports a heavy-tailed residual:
        # with 15% strong artefact pixels in the annulus it returns ~0.36x the
        # actual scatter (dropping the sigma-clip only recovers ~14% of that).
        # rms_std / rms_local is therefore the direct test of whether the
        # residual near this source is artefact-dominated -- and if it is, the
        # robust estimator is hiding it and sigma_local is a floor, not the
        # noise. Not used in the budget: on a poorly-cleaned extended source the
        # annulus can still hold real emission, which would inflate it.
        rms_std = float(np.nanstd(annulus_pixels))

    rms_local = rms_clipped if clip else rms_raw

    # -- sigma_local: local rms, floored and inflated -------------------------
    error_noise = rms_local
    if peak_at_residual is not None and np.isfinite(peak_at_residual):
        error_noise = max(error_noise, abs(float(peak_at_residual)))
    if n_beams > 0 and np.isfinite(error_noise):
        error_noise *= np.sqrt(1.0 + 1.0 / (2.0 * n_beams))

    # -- sigma_pix: peak sampling ---------------------------------------------
    _k = 4.0 * np.log(2.0) / beam_size_pix ** 2
    pixelisation_bias = abs(peak_flux) * (_k / 6.0)
    error_pix = abs(peak_flux) * (_k / np.sqrt(90.0))

    # Sub-pixel peak, separable 3-point parabolic interpolation. Reported, not
    # substituted: `peak_of_flux` keeps its "maximum pixel" meaning.
    peak_interp = peak_flux
    if 0 < y_c < ny - 1 and 0 < x_c < nx - 1:
        _, px = _parabolic_peak_1d(float(image[y_c, x_c - 1]), peak_flux,
                                   float(image[y_c, x_c + 1]))
        _, py = _parabolic_peak_1d(float(image[y_c - 1, x_c]), peak_flux,
                                   float(image[y_c + 1, x_c]))
        peak_interp = peak_flux + (px - peak_flux) + (py - peak_flux)

    # -- calibration ----------------------------------------------------------
    error_cal = systematic_error_fraction * abs(peak_flux)

    peak_error = float(np.sqrt(np.nansum(
        np.array([error_noise, error_pix, error_cal], dtype=float) ** 2)))

    # Expected maximum of the noise over the region: the peak is a maximum over
    # many correlated pixels, so it is biased high for faint sub-regions.
    max_bias_sigma = float('nan')
    if n_beams_in_mask is not None and np.isfinite(n_beams_in_mask) \
            and n_beams_in_mask > 1:
        max_bias_sigma = float(np.sqrt(2.0 * np.log(n_beams_in_mask)))

    def _ratio(num, den):
        return num / den if (den is not None and np.isfinite(den) and den > 0) \
            else np.nan

    return {
        'peak_flux': peak_flux,
        'peak_error': peak_error,
        'error_noise': error_noise,
        'error_pix': error_pix,
        'error_cal': error_cal,
        'rms_local': rms_local,
        'rms_local_clipped': rms_clipped,
        'rms_local_std': rms_std,
        'snr': _ratio(peak_flux, rms_clipped),
        'snr_noise': _ratio(peak_flux, error_noise),
        'snr_meas': _ratio(peak_flux, peak_error),
        'peak_interp': peak_interp,
        'pixelisation_bias': pixelisation_bias,
        'max_bias_sigma': max_bias_sigma,
        'n_pixels': n_pixels,
        'n_beams': n_beams,
        'inner_radius': inner_radius,
        'outer_radius': outer_radius,
        'flag': '+'.join(flags) if flags else 'ok',
    }


def compute_local_rms(data, footprint, use_sigma_clip=False):
    """
    Robust noise level *inside* a given footprint.

    Unlike the global `std` used throughout `compute_image_properties` (which is
    the `rms` argument, else `mad_std(residual)`, else `mad_std(image)`, always
    over the whole map), this is measured only over the pixels of `footprint`.
    Fed the residual map and the mask of the region being measured, it gives the
    noise where the flux was actually summed -- which differs from region to
    region, since bright cores leave larger cleaning residuals than faint
    extended emission.

    Parameters
    ----------
    data : ndarray
        Map to measure. Normally the residual image; the science image when no
        residual is available (then set `use_sigma_clip=True`).
    footprint : ndarray
        Region to measure over. Booleans, or the float/int masks this module
        produces via `mask * mask_component` -- both are accepted.
    use_sigma_clip : bool, optional
        Sigma-clip (3 sigma, 5 iterations) before taking the MAD. Default False,
        which is correct for a residual map. Set True when `data` still contains
        source emission, so the source does not set the noise scale.

    Returns
    -------
    rms : float
        `mad_std` of the selected pixels, or NaN if the footprint is empty.
    npix : int
        Number of pixels the estimate is based on.
    """
    from astropy.stats import mad_std, sigma_clip

    fp = np.asarray(footprint, dtype=bool)
    npix = int(fp.sum())
    if npix == 0:
        return float('nan'), 0

    vals = np.asarray(data)[fp]
    if use_sigma_clip:
        clipped = sigma_clip(vals, sigma=3, maxiters=5)
        compressed = clipped.compressed()
        if len(compressed) > 0:
            vals = compressed  # otherwise keep the unclipped values

    return float(mad_std(vals, ignore_nan=True)), npix


def calculate_robust_flux_error(data, data_res, mask, total_flux, total_flux_density_residual,
                              systematic_error_fraction, beam_area_, total_pixels, res_std_mad,
                              n_bootstrap=1000, confidence_level=0.68):
    """
    Calculate robust flux density error using pre-computed statistics and maps.
    
    Parameters:
    -----------
    data : numpy.ndarray
        2D array of the science image
    data_res : numpy.ndarray
        2D array of the residual image
    mask : numpy.ndarray
        2D array boolean mask
    total_flux : float
        Pre-computed total flux density within mask
    total_flux_density_residual : float
        Pre-computed total residual flux density within mask
    systematic_error_fraction : float
        Systematic error as a fraction (e.g., 0.05 for 5%)
    beam_area_ : int
        Beam area in pixels
    total_pixels : int
        Total number of pixels within mask
    res_std_mad : float
        Pre-computed MAD standard deviation of residual map within mask
    n_bootstrap : int
        Number of bootstrap iterations
    confidence_level : float
        Confidence level for error estimation (0.68 = 1sigma)
        
    Returns:
    --------
    dict
        Dictionary containing various error estimates and statistics
    """
    # Calculate number of independent beams
    n_beams = total_pixels / beam_area_
    
    # Get masked pixels for bootstrap
    masked_data = data[mask]
    
    # Bootstrap resampling
    bootstrap_fluxes = np.zeros(n_bootstrap)
    for i in range(n_bootstrap):
        # Resample pixels with replacement
        indices = np.random.randint(0, total_pixels, size=total_pixels)
        bootstrap_fluxes[i] = np.sum(masked_data[indices]) / beam_area_
    
    # Calculate errors using different methods
    
    # 1. MAD-based error (using pre-computed MAD std)
    # mad_error = res_std_mad * np.sqrt(total_pixels) / beam_area_
    mad_error = res_std_mad * np.sqrt(total_pixels) / np.sqrt(n_beams)
    
    # 2. Bootstrap-based error
    bootstrap_std = np.std(bootstrap_fluxes)
    bootstrap_percentiles = np.percentile(bootstrap_fluxes, 
                                        [50 - confidence_level*50, 50 + confidence_level*50])
    bootstrap_error = (bootstrap_percentiles[1] - bootstrap_percentiles[0]) / 2
    
    # 3. Systematic error
    systematic_error = systematic_error_fraction * total_flux
    
    # 4. Residual-based error (using your pre-computed total residual)
    # residual_error = abs(total_flux_density_residual) / beam_area_
    residual_error = abs(total_flux_density_residual) / np.sqrt(n_beams)
    
    # Combine errors (using bootstrap as statistical component)
    total_error = np.sqrt(bootstrap_error**2 + systematic_error**2)
    
    # Alternative total error using MAD
    total_error_mad = np.sqrt(mad_error**2 + systematic_error**2)
    
    # Alternative total error using residual
    total_error_residual = np.sqrt(residual_error**2 + systematic_error**2)
    
    # Calculate signal-to-noise ratios
    snr = total_flux / total_error
    snr_mad = total_flux / total_error_mad
    snr_residual = total_flux / total_error_residual
    
    return {
        'total_flux': total_flux,
        'total_error': total_error,
        'total_error_mad': total_error_mad,
        'total_error_residual': total_error_residual,
        'bootstrap_error': bootstrap_error,
        'mad_error': mad_error,
        'residual_error': residual_error,
        'systematic_error': systematic_error,
        'snr': snr,
        'snr_mad': snr_mad,
        'snr_residual': snr_residual,
        'confidence_level': confidence_level,
        'n_beams': n_beams,
        'bootstrap_percentiles': bootstrap_percentiles
    }



def structural_morphology(imagelist, residuallist,
                          indices, masks_deblended,
                          zd, ref_mask=None, data_2D=None,sigma_mask=6.0,
                          iterations = 1,
                          iterations_subregions = 1,
                          dilation_size=None,
                          min_sigma = 3.0,
                          flux_units='Jy/beam',
                          flux_conversion_factor=None,
                          do_PLOT=False, show_figure=False,
                          do_petro = True,
                          sigma_loop_init=6.0, do_measurements='all',
                          verbose=0):
    """
    From the emission  of a given source and its deblended components,
    run in each component the morphometry analysis.

    A list of images is accepted with a common deblended mask for all.

    This was originally intended to be used with multi-resolution images.
    First, a common-representative image is processed with a source detection and
    deblending algorithm. Then, those detected/deblended regions are used to
    run a forced morphometry on all multi-resolution images.

    << Finish documentation >>

    """
    results_conc = []
    missing_data_im = []
    missing_data_re = []
    masks = []
    for i in tqdm(range(len(imagelist))):
        # try:
        crop_image = imagelist[i]
        data_2D = load_fits_data(crop_image)

        if residuallist is not None:
            crop_residual = residuallist[i]
            residual_2D = load_fits_data(crop_residual)
            std = mad_std(residual_2D)
            if verbose >= 1:
                print('Using RMS from residual')
        else:
            crop_residual = None
            residual_2D = None
            std = mad_std(data_2D)
            if verbose >= 1:
                print('Using RMS from data ')
        #         crop_residual = residuallist[i]
        cell_size = get_cell_size(crop_image)
        # std = mad_std(data_2D)
        try:
            npixels = int(2*beam_area2(crop_image))
        except:
            npixels = 50
        #             residual_2D = load_fits_data(crop_residual)



        processing_results_source = {}  # store calculations only for source
        processing_results_source['#imagename'] = os.path.basename(
            crop_image)
        # first, run the analysis for the entire source structure
        processing_results_source, mask, _ = measures(imagename=crop_image,
                                                        residualname=crop_residual,
                                                        z=zd, deblend=False,
                                                        apply_mask=True,
                                                        mask = ref_mask,
                                                        results_final=processing_results_source,
                                                        flux_units=flux_units,
                                                        flux_conversion_factor=flux_conversion_factor,
                                                        plot_catalog=False,
                                                        rms=std,
                                                        bkg_sub=False,
                                                        bkg_to_sub=None,
                                                        mask_component=None,
                                                        npixels=npixels, fwhm=121,
                                                        kernel_size=121,
                                                        sigma_mask=sigma_mask,
                                                        last_level=min_sigma,
                                                        iterations=iterations,
                                                        dilation_size=None,
                                                        do_measurements=do_measurements,
                                                        do_PLOT=do_PLOT,
                                                        do_petro = do_petro,
                                                        show_figure=show_figure,
                                                        add_save_name='',
                                                        verbose=verbose)
        flag_subcomponent = 0
        processing_results_source['freq'] = getfreqs([crop_image])[0]
        # `comp_ID` is an INTEGER everywhere: 0 is the whole source, 1..n are
        # the deblended components. It used to be written as a string here and
        # as an int in `compute_model_properties`, which is why joining the two
        # frames needed `== '0'` in one and `== 1.0` in the other. Integers also
        # sort correctly past ten components.
        processing_results_source['comp_ID'] = 0
        processing_results_source['flag_subcomponent'] = flag_subcomponent
        # ref_mask = mask*ref_mask
        if ref_mask is None:
            ref_mask = mask
        results_conc.append(processing_results_source)
        # bkg_ = sep_background(crop_image, apply_mask=True, mask=None, bw=11,
        #                       bh=11, fw=12, fh=12)

        omaj, omin, _, _, _ = beam_shape(crop_image)
        if dilation_size is None:
            # Full beam-based footprint per dilation step, matching the
            # default mask_dilation_from_mask() itself would use. Growth is
            # already bounded by ref_mask, so there's no reason to halve it.
            # Recomputed every image since beam size can vary across a
            # multi-resolution imagelist.
            dilation_size_i = int(np.sqrt(omaj * omin) / (2 * get_cell_size(crop_image)))
        else:
            dilation_size_i = dilation_size
        # print('dilation_size=', dilation_size_i)
        # masks_expanded = {}
        # Tracks pixels already assigned to an earlier (lower j) component in
        # this image. Excluding these from every subsequent mask_new keeps
        # sub-region masks strictly non-overlapping -- two adjacent growing
        # masks can otherwise both claim pixels in the buffer zone between
        # them (only the original undilated cores are mutually exclusive
        # below, not the dilated footprints), double-counting that flux.
        claimed_mask = np.zeros(ref_mask.shape, dtype=bool)

        if len(indices) > 1:
            for j in range(len(indices)):
                # ii = str(i+1)
                sigma_loop = sigma_loop_init  # reset the loop
                processing_results_components = {}  # store calculation only for individual components of the soruce
                processing_results_components['#imagename'] = os.path.basename(crop_image)
                

                mask_component = masks_deblended[j]
                data_component = mask_component*data_2D.copy()
                add_save_name = 'comp_' + str(j+1)
                # print('Component id ', processing_results_components['comp_ID'])
                try:
                    # mask_new = mask_component.copy()
                    # NOTE: `sigma` below only shapes the discarded first
                    # return value of mask_dilation_from_mask (an intensity
                    # threshold mask). The mask we actually use, `mask_new`,
                    # is its second return -- a purely geometric dilation
                    # (disk(dilation_size_i), iterations_subregions steps)
                    # that does not depend on sigma_loop or the pixel data.
                    _, mask_new = \
                        mask_dilation_from_mask(data_2D,
                                                mask_component,
                                                sigma=sigma_loop,
                                                PLOT=False,iterations=iterations_subregions,
                                                dilation_size=dilation_size_i,
                                                show_figure=False)
                    mask_new = mask_new.astype(bool) & ref_mask.astype(bool)  # avoid growing beyond the reference mask
                    # dilated masks must not overlap >> non conservation of flux
                    for l in range(len(indices)):
                        if l != j:
                            mask_new[masks_deblended[l]] = False  # never invade another component's detected core
                        else:
                            pass
                    # Also drop pixels already claimed by an earlier (lower j)
                    # component -- see claimed_mask note above.
                    mask_new = mask_new & ~claimed_mask
                    # masks_expanded[f'mask_ex_{j}'] = mask_new
                    # plt.figure()
                    # plt.imshow(mask_new*ref_mask)
                    # plt.show()
                    processing_results_components, mask, _ = \
                        measures(crop_image, crop_residual, z=zd,
                                    deblend=False, apply_mask=False,
                                    plot_catalog=False,
                                    flux_units=flux_units,
                                    flux_conversion_factor=flux_conversion_factor,
                                    bkg_sub=False,
                                    mask = ref_mask,
                                    mask_component=mask_new, rms=std,
                                    iterations=iterations_subregions, npixels=50, fwhm=121,
                                    kernel_size=121, sigma_mask=sigma_loop,
                                    last_level=min_sigma,
                                    # bkg_to_sub = bkg_.back(),
                                    dilation_size=dilation_size_i,
                                    add_save_name=add_save_name,
                                    do_measurements=do_measurements,
                                    do_PLOT=do_PLOT, show_figure=show_figure,
                                    do_petro = do_petro,
                                    verbose=verbose,
                                    results_final=processing_results_components)
                    flag_subcomponent = 0
                    # print('Component id ', processing_results_components['comp_ID'])
                    processing_results_components['freq'] = getfreqs([crop_image])[0]
                    processing_results_components['comp_ID'] = j + 1
                    processing_results_components['flag_subcomponent'] = flag_subcomponent
                    claimed_mask = claimed_mask | mask_new
                except:
                    try:
                        error_mask = True
                        while error_mask and sigma_loop > 1.0:
                            try:
                                # mask_new = mask_component.copy()
                                _, mask_new = \
                                    mask_dilation_from_mask(load_fits_data(crop_image),
                                                            mask_component,
                                                            rms=std,
                                                            sigma=sigma_loop,
                                                            iterations=3,
                                                            PLOT=False,
                                                            dilation_size=dilation_size_i,
                                                            show_figure=False)
                                mask_new = mask_new.astype(bool) & ref_mask.astype(bool)

                                # dilated masks must not overlap >> non conservation of flux
                                for l in range(len(indices)):
                                    if l != j:
                                        mask_new[masks_deblended[l]] = False
                                    else:
                                        pass
                                # see claimed_mask note above: also exclude pixels
                                # already claimed by an earlier component.
                                mask_new = mask_new & ~claimed_mask

                                # masks_expanded[f'mask_ex_{j}'] = mask_new

                                if sigma_loop >= 3.0:
                                    last_level = min_sigma
                                if sigma_loop < 3.0:
                                    last_level = sigma_loop - 0.5

                                (processing_results_components, mask,
                                    _) = measures(crop_image, crop_residual,
                                                z=zd, deblend=False,
                                                apply_mask=False,
                                                plot_catalog=False,
                                                bkg_sub=False,
                                                flux_units=flux_units,
                                                flux_conversion_factor=flux_conversion_factor,
                                                mask = ref_mask,
                                                mask_component=mask_new,
                                                rms=std, iterations=3,
                                                npixels=1000, fwhm=121,
                                                kernel_size=121,
                                                sigma_mask=sigma_loop,
                                                last_level=last_level,
                                                add_save_name=add_save_name,
                                                dilation_size=dilation_size_i,
                                                do_measurements=do_measurements,
                                                do_PLOT=do_PLOT,
                                                do_petro = do_petro,
                                                show_figure=show_figure,
                                                verbose=verbose,
                                                results_final=processing_results_components)
                                # print('Component id ', processing_results_components['comp_ID'])
                                processing_results_components['subreg_sigma'] = sigma_loop
                                error_mask = False
                                flag_subcomponent = 1
                                processing_results_components['freq'] = getfreqs([crop_image])[0]
                                processing_results_components['comp_ID'] = j + 1
                                processing_results_components['flag_subcomponent'] = flag_subcomponent
                                claimed_mask = claimed_mask | mask_new
                            except Exception as e:
                                # Handle the error, and decrease p by 0.5
                                print(
                                    f"Error occurred with sigma={sigma_loop}: {e}")
                                sigma_loop -= 0.5
                                print("Reducing sigma to=", sigma_loop)

                        if not error_mask:
                            print("Function call successful with sigma_mad=",
                                    sigma_loop)
                        else:
                            print(
                                "Unable to call function with any value of sigma_mad.")
                    except:
                        print(
                            'Last attempt to perform morphometry, '
                            'with mininum threshold allowed.')
                        processing_results_components, mask, _ = measures(
                            crop_image, crop_residual, z=zd, deblend=False,
                            apply_mask=False,mask = ref_mask,
                            plot_catalog=False, bkg_sub=False,
                            flux_units=flux_units,
                            flux_conversion_factor=flux_conversion_factor,
                            # bkg_to_sub = bkg_.back(),
                            mask_component=mask_new, rms=std,
                            dilation_size=dilation_size_i, iterations=2,
                            npixels=int(beam_area2(crop_image)),
                            fwhm=81, kernel_size=81, sigma_mask=2.0,
                            last_level=0.5,
                            add_save_name=add_save_name,
                            do_measurements=do_measurements,
                            do_PLOT=do_PLOT, show_figure=show_figure,
                            do_petro = do_petro,
                            verbose=verbose,
                            results_final=processing_results_components)
                        processing_results_components['subreg_sigma'] = 1.0
                        processing_results_components['comp_ID'] = j + 1
                        processing_results_components['freq'] = getfreqs([crop_image])[0]
                        flag_subcomponent = 1
                        processing_results_components['flag_subcomponent'] = flag_subcomponent
                        claimed_mask = claimed_mask | mask_new

                results_conc.append(processing_results_components)
                # masks.append(masks_expanded)
            processing_results_source['ncomps'] = len(indices)

            # Diagnostic: how much of the reference mask's area/flux was
            # actually claimed by the union of all component sub-masks.
            # A shortfall here is the direct explanation for
            # sum(component flux) < total_flux_mask of the whole source --
            # it's uncovered ref_mask area (typically diffuse/extended
            # emission farther from every deblended core than
            # iterations_subregions * dilation_size reaches), not noise.
            try:
                beam_area_i = beam_area2(crop_image, cellsize=cell_size)
                ref_area  = float(np.nansum(ref_mask))
                ref_flux  = processing_results_source.get('total_flux_mask', np.nan)
                covered_area = float(np.nansum(claimed_mask))
                covered_flux = float(np.nansum(data_2D * claimed_mask)) / beam_area_i
                processing_results_source['subregions_area_completeness'] = (
                    covered_area / ref_area if ref_area > 0 else float('nan')
                )
                processing_results_source['subregions_flux_completeness'] = (
                    covered_flux / ref_flux if ref_flux not in (0, None) and not np.isnan(ref_flux)
                    else float('nan')
                )
            except Exception as e:
                if verbose >= 1:
                    print(f"-->> Could not compute subregion coverage diagnostics: {e}")
                processing_results_source['subregions_area_completeness'] = float('nan')
                processing_results_source['subregions_flux_completeness'] = float('nan')
        else:
            processing_results_source['ncomps'] = 1
        # except:
        #     print(f'Some error occured for image file {os.path.basename(crop_image)}.')
        #     missing_data_im.append(os.path.basename(crop_image))
        #     missing_data_re.append(os.path.basename(crop_residual))
        #     pass
    return (pd.DataFrame(results_conc), processing_results_source, missing_data_im)



make_flux_vs_std = deprecated("make_flux_vs_std",
                              "compute_image_properties")(compute_image_properties)




"""
 __  __                  _                          _              
|  \/  | ___  _ __ _ __ | |__   ___  _ __ ___   ___| |_ _ __ _   _ 
| |\/| |/ _ \| '__| '_ \| '_ \ / _ \| '_ ` _ \ / _ \ __| '__| | | |
| |  | | (_) | |  | |_) | | | | (_) | | | | | |  __/ |_| |  | |_| |
|_|  |_|\___/|_|  | .__/|_| |_|\___/|_| |_| |_|\___|\__|_|   \__, |
                  |_|                                        |___/
                  
#Morphometry 
"""

def background_asymmetry(img, mask, pre_clean=False):
    """
    <<<Morfometryka-core part>>>
    """
    def measure_asymmetry_patch(pos, patch, img):
        (x0, y0) = pos
        rot_cell = rot180(patch, x0, y0)
        sub = patch - rot_cell
        return np.sum(abs(sub)) / np.sum(abs(mask * img))

    Mo,No = img.shape
    gridsize = Mo // 10  # 10% of the size of the image
    n_pix = gridsize ** 2
    xcells = Mo // gridsize
    ycells = No // gridsize
    asymmetry_grid = np.zeros((xcells, ycells))
    gal_area = mask.sum()

    for xi in range(xcells):
        for yi in range(ycells):
            cell_mask = mask[xi * gridsize:(xi + 1) * gridsize,
                        yi * gridsize:(yi + 1) * gridsize]

            if cell_mask.sum() > 0:
                asymmetry_grid[xi, yi] = 0
                continue

            cell_img = img[xi * gridsize:(xi + 1) * gridsize, yi * gridsize:(yi + 1) * gridsize]
            x0, y0 = fmin(measure_asymmetry_patch, (gridsize // 2, gridsize // 2), args=(cell_img, img), disp=0)
            asymmetry_grid[xi, yi] = (gal_area / n_pix) * measure_asymmetry_patch((x0, y0), cell_img, img)
            del cell_img

    linear = asymmetry_grid[np.where(asymmetry_grid != 0)].ravel()

    if len(linear) > 0:
        BGrandom = np.random.choice(linear, 1)[0]
        BGmedian = np.median(linear)
        BGmin = linear.min()
        BGstd = np.std(linear)
        position = np.where(asymmetry_grid == linear.min())
        x0 = position[1][0] * gridsize + gridsize // 2
        y0 = position[0][0] * gridsize + gridsize // 2

    elif pre_clean == False:
        # measure background asymmetry with original pre-clean image if it fails for the clean one
        return background_asymmetry(img, mask, pre_clean=True)
    else:
        '''
           This is a fallback for when something goes wrong with background asymmetry estimates.
           It should also appear as a QF.
        '''
        BGrandom = 0
        BGmedian = 0
        BGmin = 0
        BGstd = 0
        x0 = 0
        y0 = 0

    return BGrandom, BGmedian, BGmin, BGstd, x0, y0


def assimetria0(pos, img, mask, box=False):
    """
    <<<Morfometryka-core part>>>
    """
    # print(' @ - Computing Asymetry 0')
    (x0, y0) = pos
    # psfmask = psfmask(psfsigma, *img.shape, x0, y0)
    if (box):
        boxmask = np.zeros_like(img)
        try:
            radii_px = np.ceil(self.P.Rp * self.NRp)
        except:
            radii_px = np.ceil(self.P.Rp * 1.5)
        boxmask[int(x0 - radii_px):int(x0 + radii_px), int(y0 - radii_px):int(y0 + radii_px)] = 1
        imgorig = boxmask * img
        imgsub = boxmask * (img - rot180(img, x0, y0))
        A = np.sum(abs(imgsub)) / np.sum(abs(imgorig))
    else:
        imgorig = img * mask
        imgsub = (img - rot180(img, x0, y0)) * mask
        A = np.sum(abs(imgsub)) / np.sum(abs(imgorig))

    del imgorig, imgsub
    return A


def assimetria1(pos, img, mask,use_mask=True):
    """
    <<<Morfometryka-core part>>>
    """
    # print(' @ - Computing Asymetry 0')
    x0, y0 = pos
    A1img = np.abs(img - rot180(img, x0, y0)) / (np.sum(np.abs(img)))
    if use_mask==True:
        return np.sum(mask * A1img)
    else:
        AsySigma = 3.00
        A1mask = A1img > np.median(A1img) + AsySigma * mad_std(A1img)
        return np.sum(mask * A1mask * A1img)


def geo_mom(p, q, I, centered=True, normed=True, complex=False, verbose=False):
    """
    <<<Morfometryka-core part>>>
    return the central moment M_{p,q} of image I
    http://en.wikipedia.org/wiki/Image_moment
    F.Ferrari 2012, prior to 4th JPAS
    """


    M, N = I.shape
    x, y = np.meshgrid(np.arange(N), np.arange(M))

    M_00 = np.nansum(I)

    if centered:
        # centroids
        x_c = (1 / M_00) * np.nansum(x * I)
        y_c = (1 / M_00) * np.nansum(y * I)

        x = x - x_c
        y = y - y_c

        if verbose:
            print('centroid  at', x_c, y_c)

    if normed:
        NORM = M_00 ** (1 + (p + q) / 2.)
    else:
        NORM = 1.0

    if complex:
        XX = (x + y * 1j)
        YY = (x - y * 1j)
    else:
        XX = x
        YY = y

    M_pq = (1 / NORM) * np.nansum(XX ** p * YY ** q * I)

    return M_pq


def q_PA(image, sigma_scale=3.0):
    """
    <<<Morfometryka-core part>>>
    Adapted version of momenta from main mfmtk.
    """
    m00 = geo_mom(0, 0, image, centered=0, normed=0)
    m10 = geo_mom(1, 0, image, centered=0, normed=0)
    m01 = geo_mom(0, 1, image, centered=0, normed=0)
    m11 = geo_mom(1, 1, image, centered=0, normed=0)
    m20 = geo_mom(2, 0, image, centered=0, normed=0)
    m02 = geo_mom(0, 2, image, centered=0, normed=0)

    mu20 = geo_mom(2, 0, image, centered=1, normed=0)
    mu02 = geo_mom(0, 2, image, centered=1, normed=0)
    mu11 = geo_mom(1, 1, image, centered=1, normed=0)

    # centroids
    x0col = m10 / m00
    y0col = m01 / m00

    # Calculate eigenvalues (these give RMS^2)
    lam1_sq = (1 / 2.) * (mu20 + mu02 + np.sqrt((mu20 - mu02) ** 2 + 4 * mu11 ** 2)) / m00
    lam2_sq = (1 / 2.) * (mu20 + mu02 - np.sqrt((mu20 - mu02) ** 2 + 4 * mu11 ** 2)) / m00
    
    lam1_rms = np.sqrt(abs(lam1_sq))
    lam2_rms = np.sqrt(abs(lam2_sq))
    
    # Apply scaling to get physical size
    a_phys = max(lam1_rms, lam2_rms) * sigma_scale
    b_phys = min(lam1_rms, lam2_rms) * sigma_scale
    
    a_rms = max(lam1_rms, lam2_rms)
    b_rms = min(lam1_rms, lam2_rms)

    PA = (1 / 2.) * np.arctan2(2 * mu11, (mu20 - mu02))
    if PA < 0:
        PA = PA + np.pi
    PAdeg = np.rad2deg(PA)
    
    q = b_phys / a_phys
    
    return PAdeg, q, x0col, y0col, a_phys, b_phys


def peak_center(image):
    """
    <<<Morfometryka-core part>>>
    """
    y0max, x0max = nd.maximum_position((image))
    try:
        """
        For very small images of emission, this function breaks. In that
        case, just return the the max position.
        """
        # size of peak region to consider in interpolarion for x_peak
        dp = 2
        CenterOffset = 10
        No, Mo = image.shape
        peakimage = image[y0max - dp:y0max + dp, x0max - dp:x0max + dp]
        m00 = geo_mom(0, 0, peakimage, centered=0, normed=0)
        m10 = geo_mom(1, 0, peakimage, centered=0, normed=0)
        m01 = geo_mom(0, 1, peakimage, centered=0, normed=0)

        x0peak = x0max + m10 / m00 - dp
        y0peak = y0max + m01 / m00 - dp

        # check if center is galaxy center, i.e., should be near the image center
        # otherwise apply a penalty to pixel value proportional to the center distance^2
        if np.sqrt((x0peak - No / 2.) ** 2 + (
                y0peak - Mo / 2.) ** 2) > CenterOffset:
            # define a penalty as we move from the center
            xx, yy = np.meshgrid(np.arange(No) - No / 2.,
                                 np.arange(Mo) - Mo / 2.)
            rr2 = xx ** 2 + yy ** 2
            y0peak, x0peak = nd.maximum_position((image / rr2))
        return (x0peak, y0peak)
    except:
        return (x0max, y0max)


def momenta(image, PArad_0=None, q_0=None, sigma_scale=3.0):
    '''
    <<<Morfometryka-core part>>>
    Calculates center of mass, axis lengths and position angle
    '''

    m00 = geo_mom(0, 0, image, centered=0, normed=0)
    m10 = geo_mom(1, 0, image, centered=0, normed=0)
    m01 = geo_mom(0, 1, image, centered=0, normed=0)
    m11 = geo_mom(1, 1, image, centered=0, normed=0)
    m20 = geo_mom(2, 0, image, centered=0, normed=0)
    m02 = geo_mom(0, 2, image, centered=0, normed=0)

    mu20 = geo_mom(2, 0, image, centered=1, normed=0)
    mu02 = geo_mom(0, 2, image, centered=1, normed=0)
    mu11 = geo_mom(1, 1, image, centered=1, normed=0)

    # centroids
    x0col = m10 / m00
    y0col = m01 / m00

    # major, minor and axis ratio
    lam1_sq = (1 / 2.) * (mu20 + mu02 + np.sqrt((mu20 - mu02) ** 2 + 4 * mu11 ** 2)) / m00
    lam2_sq = (1 / 2.) * (mu20 + mu02 - np.sqrt((mu20 - mu02) ** 2 + 4 * mu11 ** 2)) / m00
    
    lam1_rms = np.sqrt(abs(lam1_sq))
    lam2_rms = np.sqrt(abs(lam2_sq))
    
    # Apply scaling to get physical size
    a = max(lam1_rms, lam2_rms) * sigma_scale
    b = min(lam1_rms, lam2_rms) * sigma_scale

    PA = (1 / 2.) * np.arctan2(2 * mu11, (mu20 - mu02))
    if PA < 0:
        PA = PA + 2*np.pi

    # self.PArad = PA
    # self.PAdeg = np.rad2deg(self.PArad)
    # self.q = self.b/self.a
    # self.PArad = PA
    # self.PAdeg = np.rad2deg(self.PArad)
    # self.q = self.b/self.a

    # mofified by lucatelli (2018)
    """This will force mfmtk to do photometry for the given input PA
    this can be useful when we want to study how the light profile changes
    as function of PA or q. This was indented to explore the difference between the
    profiles of elliptical and spiral galaxies, which may not change soo much for the
    former while it may for the later.
    """

    if PArad_0 is None:
        PArad = PA  # + np.pi/2
    else:
        PArad = PArad_0

    PAdeg = np.rad2deg(PArad)

    if q_0 is None:
        q = b / a
    else:
        q = q_0
    return (x0col, y0col, a, b, q, PAdeg)


def cal_PA_q(gal_image_0,Isequence = None,region_split=None,SAVENAME=None):
    '''
    <<<Morfometryka-core part>>>
    Estimates inner and outer PA nad q=(b/a)
    '''
    # mean Inner q,  mean outer q,  mean Inner PA,  mean Outer PA
    # from fitEllipse import main_test2
    
    # qmi, qmo, PAmi, PAmo, qm, PAm,\
    #     x0median,y0median,x0median_i,y0median_i,\
    #     x0median_o,y0median_o = main_test2(gal_image_0,
    #                                        Isequence = Isequence,
    #                                        region_split=region_split,
    #                                        SAVENAME=SAVENAME)
    
    from fit_ellipse import fit_ellipse_to_galaxy
    qmi, qmo, PAmi, PAmo, qm, PAm,x0median,y0median,x0median_i,y0median_i,\
        x0median_o,y0median_o,profiles = \
            fit_ellipse_to_galaxy(gal_image_0,
                                  Isequence = Isequence,
                                  region_split=region_split,
                                  fix_center=True,
                                  dx_dy=(10,10),
                                #   plot_results=True,
                                #   plot_profiles=True,
                                  save_name=SAVENAME)

    # global PA,  global q
    PA, q, x0col, y0col, a, b = q_PA(gal_image_0)

    # print("Initial PA and q = ", PA, q)
    # print("Median PA and q = ", PAm, qm)
    # print("Inner-Mean PA and q = ", PAmi, qmi)
    # print("Outer-Mean PA and q = ", PAmo, qmo)
    return (PA, q, x0col, y0col, PAm, qm, PAmi, qmi, PAmo, qmo,
            x0median,y0median,x0median_i,y0median_i,x0median_o,y0median_o,profiles)



def savitzky_golay_2d(z, window_size, order, derivative=None):
    """
    <<<Morfometryka-core part>>>
    http://nbviewer.ipython.org/github/pv/SciPy-CookBook/blob/master/ipython/SavitzkyGolay.ipynb
    """

    # number of terms in the polynomial expression
    n_terms = (order + 1) * (order + 2) / 2.0

    if window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size ** 2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [(k - n, n) for k in range(order + 1) for n in range(k + 1)]

    # coordinates of points
    ind = np.arange(-half_size, half_size + 1, dtype=np.float64)
    dx = np.repeat(ind, window_size)
    dy = np.tile(ind, [window_size, 1]).reshape(window_size ** 2, )

    # build matrix of system of equation
    A = np.empty((window_size ** 2, len(exps)))
    for i, exp in enumerate(exps):
        A[:, i] = (dx ** exp[0]) * (dy ** exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2 * half_size, z.shape[1] + 2 * half_size
    Z = np.zeros((new_shape))
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] = band - np.abs(np.flipud(z[1:half_size + 1, :]) - band)
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band + np.abs(np.flipud(z[-half_size - 1:-1, :]) - band)
    # left band
    band = np.tile(z[:, 0].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs(np.fliplr(z[:, 1:half_size + 1]) - band)
    # right band
    band = np.tile(z[:, -1].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, -half_size:] = band + np.abs(np.fliplr(z[:, -half_size - 1:-1]) - band)
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0, 0]
    Z[:half_size, :half_size] = band - np.abs(np.flipud(np.fliplr(z[1:half_size + 1, 1:half_size + 1])) - band)
    # bottom right corner
    band = z[-1, -1]
    Z[-half_size:, -half_size:] = band + np.abs(np.flipud(np.fliplr(z[-half_size - 1:-1, -half_size - 1:-1])) - band)

    # top right corner
    band = Z[half_size, -half_size:]
    Z[:half_size, -half_size:] = band - np.abs(np.flipud(Z[half_size + 1:2 * half_size + 1, -half_size:]) - band)
    # bottom left corner
    band = Z[-half_size:, half_size].reshape(-1, 1)
    Z[-half_size:, :half_size] = band - np.abs(np.fliplr(Z[-half_size:, half_size + 1:2 * half_size + 1]) - band)

    # solve system and convolve
    if derivative is None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')


def standartize(image, q, PArad, x0, y0):
    """
    <<<Morfometryka-core part>>>
    make a standard galaxy, id est, PA=0, q=1
    arguments are 'image' to be standartized and  its 'S' stamp and P phot classes
    """

    ##### rotate array
    R = np.array([[np.cos(PArad), np.sin(PArad)], [-np.sin(PArad), np.cos(PArad)]])

    ##### shear array
    S = np.diag([q, 1.])
    # SERSIC fit values
    # S = np.diag([self.Ss.qFit2D, 1.])

    # affine transform matrix, rotate then scale
    transform = np.dot(R, S)

    # where to transform about
    centro_i = (x0, y0)
    # contro_o: where to put center after
    centro_o = np.array(image.shape) / 2

    myoffset = centro_i - np.dot(transform, centro_o)
    bval = np.mean(image[-2:])
    stangal = nd.affine_transform(image, transform, offset=myoffset, order=2, cval=bval)

    return stangal


def polarim(image, origin=None, log=False):
    """
    <<<Morfometryka-core part>>>
    Reprojects a 2D numpy array ("image") into a polar coordinate system.
    "origin" is a tuple of (x0, y0) and defaults to the center of the image.
    http://stackoverflow.com/questions/3798333/image-information-along-a-polar-coordinate-system
    refactored by FF, 2013-2014 (see transpolar.py)
    """

    if origin is None:
        origin = np.array(image.shape) / 2.

    def cart2polar(x, y):
        r = np.sqrt(x ** 2 + y ** 2)
        #         max_radii = np.sqrt(x.max()**2+y.max()**2)
        #         rscale = x.max()/max_radii
        #         tscale = y.max()/(2*np.pi)
        theta = np.arctan2(y, x)
        return r, theta

    def polar2cart(r, theta):
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return x, y

    def cart2logpolar(x, y, M=1):
        alpha = 0.01
        r = np.sqrt(x ** 2 + y ** 2)
        rho = M * np.log(r + alpha)
        theta = np.arctan2(y, x)
        return rho, theta

    def logpolar2cart(rho, theta, M=1):
        x = np.exp(rho / M) * np.cos(theta)
        y = np.exp(rho / M) * np.sin(theta)
        return x, y

    ny, nx = image.shape
    if origin is None:
        x0, y0 = (nx // 2, ny // 2)
        origin = (x0, y0)
    else:
        x0, y0 = origin

    # Determine that the min and max r and theta coords will be...
    x, y = np.meshgrid(np.arange(nx) - x0, np.arange(ny) - y0)  # ,sparse=True )

    r, theta = cart2polar(x, y)

    # Make a regular (in polar space) grid based on the min and max r & theta
    r_i = np.linspace(r.min(), r.max(), nx)
    theta_i = np.linspace(theta.min(), theta.max(), ny)
    theta_grid, r_grid = np.meshgrid(theta_i, r_i)

    # Project the r and theta grid back into pixel coordinates
    xi, yi = polar2cart(r_grid, theta_grid)
    xi += origin[0]  # We need to shift the origin back to
    yi += origin[1]  # back to the lower-left corner...
    xi, yi = xi.flatten(), yi.flatten()
    coords = np.vstack((xi, yi))  # (map_coordinates requires a 2xn array)

    zi = nd.map_coordinates(image, coords, order=1)  # ,prefilter=False)
    galpolar = zi.reshape((nx, ny))

    r_polar = r_i
    theta_polar = theta_i

    return galpolar, r_polar, theta_polar, x, y


def Gradindex(image_data, Rp=None):
    """
    <<<Morfometryka-core part>>>
    Gradient Index
    Calculates an index based on the image gradient magnitude and orientation

    SGwindow and SGorder are Savitsky-Golay filter parameters
    F. Ferrari, 2014
    """

    def sigma_func(params):
        '''
        <<<Morfometryka-core part>>>
        calculates the sigma psi with different parameters
        called by the minimization routine '''
        (x0, y0, q, PA) = params
        #### creates standardized image
        # using segmentation geometric parameters
        # galnostarsstd = standartize(P.galnostars, S.q, S.PArad, S.y0col, S.x0col)
        # using Sersic geometric parameters
        galnostarsstd = standartize(image_data, q, PA, x0, y0)

        # creates polar imagem
        galpolar, r_polar, theta_polar, _, _ = polarim(galnostarsstd)
        #     print(galpolar)
        # print '%.5f %.5f %.5f %.5f' % (x0, y0, q, PA)

        if Rp is None:
            galpolarpetro = galpolar[0: -1, :]
        else:
            galpolarpetro = galpolar[0: 2 * int(Rp), :]
        #     galpolarpetro =  galpolar[ 0   : int(config.NRp * P.Rp), : ]

        # circular_area_radius
        #     if min(galpolarpetro.shape) <= SGwindow:
        #         SGorder -= 1

        try:
            # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
            dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, SGorder, 'both')
        except:
            SGwindow = 5
        # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
        dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, 1, 'both')

        mag = np.sqrt(dx ** 2 + dy ** 2)
        # mag = dxdy
        magmask = mag > (np.median(mag))
        ort = np.arctan2(dy, dx)
        # ort = np.arctan(dxdy)
        ortn = (ort + np.pi) % (np.pi)

        psi = circmean(ortn[magmask])
        sigma_psi = circstd(ortn[magmask])

        return sigma_psi

    def sigma_func_eval(params):
        '''
        <<<Morfometryka-core part>>>
        calculates the sigma psi with different parameters
        called by the minimization routine '''
        (x0, y0, q, PA) = params
        #### creates standardized image
        # using segmentation geometric parameters
        # galnostarsstd = standartize(P.galnostars, S.q, S.PArad, S.y0col, S.x0col)
        # using Sersic geometric parameters
        galnostarsstd = standartize(image_data, q, PA, x0, y0)

        # creates polar imagem
        galpolar, r_polar, theta_polar, _, _ = polarim(galnostarsstd)
        #     print(galpolar)
        # print '%.5f %.5f %.5f %.5f' % (x0, y0, q, PA)

        #     galpolarpetro =  galpolar[ 0   : int(config.NRp * P.Rp), : ]
        if Rp is None:
            galpolarpetro = galpolar[0: -1, :]
        else:
            galpolarpetro = galpolar[0: 2 * int(Rp), :]
        # circular_area_radius
        #     if min(galpolarpetro.shape) <= SGwindow:
        #         SGorder -= 1

        try:
            # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
            dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, SGorder, 'both')
        except:
            SGwindow = 5
        # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
        dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, 1, 'both')

        mag = np.sqrt(dx ** 2 + dy ** 2)
        # mag = dxdy
        magmask = mag > (np.median(mag))
        ort = np.arctan2(dy, dx)
        # ort = np.arctan(dxdy)
        ortn = (ort + np.pi) % (np.pi)

        psi = circmean(ortn[magmask])
        sigma_psi = circstd(ortn[magmask])

        return sigma_psi, ort, ortn, dx, dy, mag

    # SAVITSKY-GOLAY parameters
    # polynom order
    SGorder = 5

    # SG window is galaxy_size/10 and must be odd
    # CAN'T BE RELATIVE TO IMAGE SIZE... MUST BE RELATIVE TO GALAXY SIZE
    # SGwindow = int(S.Mo/10.)
    # SGwindow = int(P.Rp/2.)
    SGwindow = 5
    #     SGwindow = int(circular_area_radius/2)
    if SGwindow % 2 == 0:
        SGwindow = SGwindow + 1

    #     image_data=load_fits_data(imagelist[1])

    PA, q, x0col, y0col, PAm, qm, PAmi, qmi, PAmo, qmo,profiles = cal_PA_q(image_data)

    x0sigma, y0sigma, qsigma, PAsigma = \
        fmin(sigma_func, (x0col, y0col, qm, np.deg2rad(PAm)), ftol=0.1, xtol=1.0, disp=0)
    sigma_psi, ort, ortn, dx, dy, mag = sigma_func_eval((x0sigma, y0sigma, qsigma, PAsigma))

    #     x0sigma, y0sigma, qsigma, PAsigma = \
    #                 fmin(sigma_func, (Ss.x0Fit2D, Ss.y0Fit2D, Ss.qFit2D, np.deg2rad(Ss.PAFit2D)), ftol=0.1, xtol=1.0, disp=0)

    #     #if qsigma > 1:
    #     #   qsigma = 1/qsigma
    #     #   PAsigma = PAsigma - np.pi/2.

    #     x0sigma, y0sigma, qsigma, PAsigma = \
    #                 fmin(sigma_funcseg, (Ss.x0Fit2D, Ss.y0Fit2D, Ss.qFit2D, np.deg2rad(Ss.PAFit2D)), ftol=0.1, xtol=1.0, disp=0)
    return (sigma_psi, ort, ortn, dx, dy, mag, PAm, qm)



def evaluate_compactness(deconv_props, conv_props):
    """
    Determine if a model component of a radio structure is compact or extended.

    It uses the concentration index determined from the areas A20, A50, A80 and A90.

    """
    import pandas as pd
    ncomps = deconv_props['comp_ID'].shape[0]

    Spk_ratio = np.asarray(conv_props['peak_of_flux']) / np.asarray(deconv_props['peak_of_flux'])
    I50_ratio = np.asarray(deconv_props['I50']) / np.asarray(conv_props['I50'])
    class_criteria = {}
    for i in range(ncomps):
        class_criteria[f"comp_ID_{i + 1}"] = {}
        if Spk_ratio[i] < 0.5:
            class_criteria[f"comp_ID_{i + 1}"]['Spk_class'] = 'C'
        if Spk_ratio[i] >= 0.5:
            class_criteria[f"comp_ID_{i + 1}"]['Spk_class'] = 'D'

        if (I50_ratio[i] < 0.5) or (np.isnan(I50_ratio[i])): #nan because the component is too
            # small.
            class_criteria[f"comp_ID_{i + 1}"]['I50_class'] = 'C'
        if I50_ratio[i] >= 0.5:
            class_criteria[f"comp_ID_{i + 1}"]['I50_class'] = 'D'

        AC1_conv_check = conv_props['AC1'].iloc[i]
        AC2_conv_check = conv_props['AC2'].iloc[i]
        AC1_deconv_check = deconv_props['AC1'].iloc[i]
        AC2_deconv_check = deconv_props['AC2'].iloc[i]
        # print(AC1_conv_check)

        # if (AC1_conv_check >= 1.0) or (AC1_conv_check == np.inf):
        #     class_criteria[f"comp_ID_{i+1}"]['AC1_conv_class'] = 'C'
        # if (AC2_conv_check >= 1.0) or AC2_conv_check == np.inf:
        #     class_criteria[f"comp_ID_{i+1}"]['AC2_conv_class'] = 'C'
        if (AC1_deconv_check >= 1.0) or (np.isnan(AC1_deconv_check)):
            class_criteria[f"comp_ID_{i + 1}"]['AC1_deconv_class'] = 'C'
        if AC1_deconv_check < 1.0:
            class_criteria[f"comp_ID_{i + 1}"]['AC1_deconv_class'] = 'D'
        if AC2_deconv_check >= 0.75:
            class_criteria[f"comp_ID_{i + 1}"]['AC2_deconv_class'] = 'C'
        if AC2_deconv_check < 0.75:
            class_criteria[f"comp_ID_{i + 1}"]['AC2_deconv_class'] = 'D'

    class_results_df = pd.DataFrame(class_criteria)
    for i in range(ncomps):
        dessision_compact = np.sum(class_results_df[f"comp_ID_{i + 1}"] == 'C')
        # dessision_diffuse = np.sum(class_results_df[f"comp_ID_{i+1}"]=='D')
        if dessision_compact > 2:
            class_criteria[f"comp_ID_{i + 1}"]['final_class'] = 'C'
        if dessision_compact < 2:
            class_criteria[f"comp_ID_{i + 1}"]['final_class'] = 'D'
        if dessision_compact == 2:
            try:
                class_criteria[f"comp_ID_{i + 1}"]['final_class'] = \
                    class_criteria[(f"comp_ID_{i + 1}")]['Spk_class']
            except:
                class_criteria[f"comp_ID_{i + 1}"]['final_class'] = 'D'
    return (class_criteria)





def convex_morpho_old(image, mask, scale=1.0,do_plot=False):
    """
    Overlay semi-major and minor axes on an image and report properties.
    
    Parameters:
        image (2D array): The image data to be plotted.
        mask (2D array): mask data of the structure.
        scale (float): Scaling factor for the length of the axes vectors.
        do_plot (boolean): Plot or not the results.
    
    Returns:
        dict: Report containing basic morphometry.
    """
    # Extract y, x coordinates of the structure
    indices = np.transpose(np.nonzero(mask))
    y, x = indices[:, 0], indices[:, 1]
    points = np.column_stack((x, y))  # Convert to (x, y) coordinates
    
    # Compute ConvexHull
    hull = ConvexHull(points)

    # Compute centroid of the points
    centroid = np.mean(points, axis=0)

    # Covariance matrix and eigenvalues/vectors
    cov_matrix = np.cov(points, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

    # Order eigenvectors by eigenvalues (largest is major axis)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    # Semi-major and minor axes
    major_axis_vector = eigenvectors[:, 0]
    minor_axis_vector = eigenvectors[:, 1]

    # Scale eigenvectors by eigenvalues for visualization
    major_axis = scale * np.sqrt(eigenvalues[0]) * major_axis_vector
    minor_axis = scale * np.sqrt(eigenvalues[1]) * minor_axis_vector

    # Calculate position angle (anti-clockwise from x-axis)
    position_angle = np.arctan2(major_axis_vector[1], major_axis_vector[0])
    position_angle_degrees = np.degrees(position_angle)
    if position_angle_degrees < 0:
        position_angle_degrees += 360

    # Calculate axis ratio (minor/major)
    axis_ratio = np.sqrt(eigenvalues[1] / eigenvalues[0])

    if do_plot:
        # Plot image
        plt.figure(figsize=(5, 5))
        plt.imshow(image, cmap='gray', origin='lower')
        # plt.colorbar(label='Intensity')
        
        # Plot semi-major and minor axes
        plt.quiver(
            centroid[0], centroid[1], major_axis[0], major_axis[1],
            angles='xy', scale_units='xy', scale=1, color='r', 
            # label='Major Axis'
        )
        plt.quiver(
            centroid[0], centroid[1], minor_axis[0], minor_axis[1],
            angles='xy', scale_units='xy', scale=1, color='b', 
            # label='Minor Axis'
        )
    
        # Plot centroid
        plt.scatter(centroid[0], centroid[1], color='orange', 
                    # label='Centroid', 
                    zorder=5)
    
        # Add labels
        # plt.legend()
        # plt.title("Semi-Major and Minor Axes on Image")
        plt.xlabel("$x$ image coordinates")
        plt.ylabel("$y$ image coordinates")
        plt.axis('equal')
        plt.show()

    # Return report
    report = {
        "PA_convex": position_angle_degrees,
        "q_convex": axis_ratio,
        "centroid_convex": centroid
    }
    return report


# def compute_diameters_convex(hull,points):
#     """
#     Compute the major and minor diameters of a structure using its ConvexHull.

#     Parameters
#     ----------
#     hull : scipy.spatial.ConvexHull
#         A ConvexHull object representing the structure, containing points
#         that define the convex boundary of the structure.

#     Returns
#     -------
#     dict
#         A dictionary containing:
#         - "major_diameter" : float
#             The length of the major diameter (longest distance between any two 
#             points on the convex hull).
#         - "major_points" : tuple of ndarray
#             The two points defining the major diameter.
#         - "minor_diameter" : float
#             The length of the minor diameter (shortest perpendicular distance 
#             between parallel edges of the convex hull).
#         - "minor_points" : tuple of ndarray
#             The two points defining the minor diameter.

#     Examples
#     --------
#     >>> import numpy as np
#     >>> from scipy.spatial import ConvexHull
#     >>> # Example points
#     >>> points = np.array([[0, 0], [1, 1], [2, 2], [3, 0], [0, 3]])
#     >>> # Compute the convex hull
#     >>> hull = ConvexHull(points)
#     >>> # Compute diameters
#     >>> diameters = compute_diameters_convex(hull)
#     >>> print(diameters)

#     Notes
#     -----
#     - The major diameter is calculated as the maximum Euclidean distance
#       between any two points on the convex hull.
#     - The minor diameter is determined by finding the shortest perpendicular
#       width across the convex hull.
#     """
#     from itertools import combinations

#     # hull = ConvexHull(points)
#     hull_points = points[hull.vertices]  # Extract points on the convex hull

#     # Calculate Major Diameter (longest distance between hull points)
#     major_diameter = 0
#     major_points = None
#     for p1, p2 in combinations(hull_points, 2):
#         distance = np.linalg.norm(p1 - p2)
#         if distance > major_diameter:
#             major_diameter = distance
#             major_points = (p1, p2)

#     # Calculate Minor Diameter (shortest perpendicular distance between parallel edges)
#     minor_diameter = float('inf')
#     minor_points = None
#     num_hull_points = len(hull_points)

#     for i in range(num_hull_points):
#         # Get two consecutive points forming an edge
#         p1, p2 = hull_points[i], hull_points[(i + 1) % num_hull_points]
#         edge_vector = p2 - p1
#         edge_length = np.linalg.norm(edge_vector)
        
#         # Normalize the edge vector
#         if edge_length == 0:
#             continue
#         edge_normal = np.array([-edge_vector[1], edge_vector[0]]) / edge_length
        
#         # Project all hull points onto the edge normal and find the width
#         distances = np.abs(np.dot(hull_points - p1, edge_normal))
#         max_distance = distances.max()
#         if max_distance < minor_diameter:
#             minor_diameter = max_distance
#             # Points corresponding to the shortest perpendicular projection
#             projections = hull_points[np.abs(distances - max_distance) < 1e-6]
#             if len(projections) >= 2:
#                 minor_points = (projections[0], projections[1])

#     return {
#         "major_diameter": major_diameter,
#         "major_points": major_points,
#         "minor_diameter": minor_diameter,
#         "minor_points": minor_points
#     }


# def convex_morpho(image, mask, scale=1.0, do_plot=False):
#     """
#     Perform morphological analysis on a structure defined by a mask and overlay results on an image.

#     This function computes morphological properties of a structure, including the position angle,
#     axis ratio, centroid, major diameter, and minor diameter, based on the ConvexHull of the masked region.
#     Optionally, it can overlay the semi-major and semi-minor axes, centroid, and convex hull on the image.

#     Parameters
#     ----------
#     image : 2D ndarray
#         The image data where the structure is located.
#     mask : 2D ndarray
#         A binary mask defining the structure to be analyzed. Non-zero pixels are treated as part of the structure.
#     scale : float, optional
#         Scaling factor for the visualization of the axes (default is 1.0).
#     do_plot : bool, optional
#         Whether to plot the results overlaying the axes and convex hull on the image (default is False).

#     Returns
#     -------
#     dict
#         A dictionary containing the following morphometric properties:
#         - "PA_convex" : float
#             Position angle of the structure in degrees, measured counter-clockwise from the positive x-axis.
#         - "q_convex" : float
#             Axis ratio (minor/major) of the structure.
#         - "centroid_convex" : ndarray
#             Centroid coordinates of the structure as a 2-element array [x, y].
#         - "major_diameter" : float
#             Length of the major diameter (longest distance between any two points on the convex hull).
#         - "minor_diameter" : float
#             Length of the minor diameter (shortest perpendicular distance between parallel edges of the convex hull).

#     Examples
#     --------
#     >>> import numpy as np
#     >>> import matplotlib.pyplot as plt
#     >>> from scipy.spatial import ConvexHull

#     >>> # Create a sample image and mask
#     >>> image = np.random.rand(100, 100)
#     >>> mask = np.zeros_like(image, dtype=bool)
#     >>> mask[40:60, 45:65] = True  # Example structure

#     >>> # Perform analysis
#     >>> report = convex_morpho(image, mask, scale=2.0, do_plot=True)
#     >>> print(report)

#     Notes
#     -----
#     - The function uses the ConvexHull of the masked region to estimate morphological properties.
#     - The position angle is computed from the eigenvector of the largest eigenvalue of the covariance matrix.
#     - Axis ratio is calculated as the square root of the ratio of the smallest to largest eigenvalues.
#     - Major and minor diameters are determined using the ConvexHull geometry.
#     """
#     # Extract y, x coordinates of the structure
#     indices = np.transpose(np.nonzero(mask))
#     y, x = indices[:, 0], indices[:, 1]
#     points = np.column_stack((x, y))  # Convert to (x, y) coordinates
    
#     # Compute ConvexHull
#     hull = ConvexHull(points)

#     # Compute centroid of the points
#     centroid = np.mean(points, axis=0)

#     # Covariance matrix and eigenvalues/vectors
#     cov_matrix = np.cov(points, rowvar=False)
#     eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

#     # Order eigenvectors by eigenvalues (largest is major axis)
#     order = np.argsort(eigenvalues)[::-1]
#     eigenvalues = eigenvalues[order]
#     eigenvectors = eigenvectors[:, order]

#     # Semi-major and minor axes
#     major_axis_vector = eigenvectors[:, 0]
#     minor_axis_vector = eigenvectors[:, 1]

#     # Scale eigenvectors by eigenvalues for visualization
#     major_axis = scale * np.sqrt(eigenvalues[0]) * major_axis_vector
#     minor_axis = scale * np.sqrt(eigenvalues[1]) * minor_axis_vector

#     # Calculate position angle (anti-clockwise from x-axis)
#     position_angle = np.arctan2(major_axis_vector[1], major_axis_vector[0])
#     position_angle_degrees = np.degrees(position_angle)
#     if position_angle_degrees < 0:
#         position_angle_degrees += 360

#     # Calculate axis ratio (minor/major)
#     axis_ratio = np.sqrt(eigenvalues[1] / eigenvalues[0])

#     # Compute diameters
#     diameters = compute_diameters_convex(hull,points)

#     if do_plot:
#         # Plot image

#         plt.figure(figsize=(5, 5))
#         norm = simple_norm(image, stretch='sqrt', asinh_a=0.02, vmin= 3 * mad_std(image),
#                            vmax=0.2 * np.nanmax(image))
#         plt.imshow(image, cmap='gray', origin='lower',norm=norm)
        
#         # Plot semi-major and minor axes
#         plt.quiver(
#             centroid[0], centroid[1], major_axis[0], major_axis[1],
#             angles='xy', scale_units='xy', scale=1, color='limegreen'
#         )
#         plt.quiver(
#             centroid[0], centroid[1], minor_axis[0], minor_axis[1],
#             angles='xy', scale_units='xy', scale=1, color='red'
#         )
    
#         # Plot centroid
#         plt.scatter(centroid[0], centroid[1], color='limegreen', zorder=5)
    
        
#         plt.xlabel("$x$ image coordinates")
#         plt.ylabel("$y$ image coordinates")
#         plt.axis('equal')
#         # plt.legend()
#         plt.show()

#     # Return report
#     report = {
#         "PA_convex": position_angle_degrees,
#         "q_convex": axis_ratio,
#         "centroid_convex": centroid,
#         "major_diameter": diameters["major_diameter"],
#         "minor_diameter": diameters["minor_diameter"]
#     }
#     return report





def compute_diameters_convex(hull, points, intensities=None):
    """
    Compute the major and minor diameters of a structure using its ConvexHull.
    Now with optional intensity weighting for more robust measurements.

    Parameters
    ----------
    hull : scipy.spatial.ConvexHull
        A ConvexHull object representing the structure, containing points
        that define the convex boundary of the structure.
    points : ndarray
        Array of (x, y) coordinates.
    intensities : ndarray, optional
        Intensity values for each point. If provided, uses high-intensity
        regions for more robust diameter estimation.

    Returns
    -------
    dict
        A dictionary containing:
        - "major_diameter" : float
            The length of the major diameter (longest distance between any two 
            points on the convex hull).
        - "major_points" : tuple of ndarray
            The two points defining the major diameter.
        - "minor_diameter" : float
            The length of the minor diameter (shortest perpendicular distance 
            between parallel edges of the convex hull).
        - "minor_points" : tuple of ndarray
            The two points defining the minor diameter.
    """
    from itertools import combinations
    # If intensities provided, focus on high-intensity regions
    if intensities is not None:
        # Use 90th percentile of intensity for effective boundary
        percentile = 90
        sorted_indices = np.argsort(intensities)[::-1]
        n_select = max(int(len(points) * (percentile / 100.0)), 3)
        high_intensity_indices = sorted_indices[:n_select]
        high_intensity_points = points[high_intensity_indices]
        
        # Compute hull of high-intensity points for more robust diameter estimation
        if len(high_intensity_points) >= 3:
            try:
                intensity_hull = ConvexHull(high_intensity_points)
                hull_points = high_intensity_points[intensity_hull.vertices]
            except:
                # Fall back to original hull if high-intensity hull fails
                hull_points = points[hull.vertices]
        else:
            hull_points = points[hull.vertices]
    else:
        # Original behavior when no intensities provided
        hull_points = points[hull.vertices]

    # Calculate Major Diameter (longest distance between hull points)
    major_diameter = 0
    major_points = None
    for p1, p2 in combinations(hull_points, 2):
        distance = np.linalg.norm(p1 - p2)
        if distance > major_diameter:
            major_diameter = distance
            major_points = (p1, p2)

    # Calculate Minor Diameter (shortest perpendicular distance between parallel edges)
    minor_diameter = float('inf')
    minor_points = None
    num_hull_points = len(hull_points)

    for i in range(num_hull_points):
        # Get two consecutive points forming an edge
        p1, p2 = hull_points[i], hull_points[(i + 1) % num_hull_points]
        edge_vector = p2 - p1
        edge_length = np.linalg.norm(edge_vector)
        
        # Normalize the edge vector
        if edge_length == 0:
            continue
        edge_normal = np.array([-edge_vector[1], edge_vector[0]]) / edge_length
        
        # Project all hull points onto the edge normal and find the width
        distances = np.abs(np.dot(hull_points - p1, edge_normal))
        max_distance = distances.max()
        if max_distance < minor_diameter:
            minor_diameter = max_distance
            # Points corresponding to the shortest perpendicular projection
            projections = hull_points[np.abs(distances - max_distance) < 1e-6]
            if len(projections) >= 2:
                minor_points = (projections[0], projections[1])

    return {
        "major_diameter": major_diameter,
        "major_points": major_points,
        "minor_diameter": minor_diameter,
        "minor_points": minor_points
    }


def _weighted_percentile(values, weights, pcts):
    """
    Percentiles of `values` weighted by `weights`.

    Uses the standard cumulative-weight / midpoint convention, so it reduces to
    `np.percentile` when all weights are equal.

    Parameters
    ----------
    values : 1D ndarray
        Sample values.
    weights : 1D ndarray
        Non-negative weights, same length as `values`. Need not be normalised.
    pcts : float or sequence of float
        Percentile(s) in [0, 100].

    Returns
    -------
    ndarray
        The requested weighted percentile(s), same shape as `pcts`.
    """
    values  = np.asarray(values,  dtype=float)
    weights = np.asarray(weights, dtype=float)
    order   = np.argsort(values)
    v, w    = values[order], weights[order]
    wsum    = w.sum()
    if not np.isfinite(wsum) or wsum <= 0:
        return np.percentile(values, pcts)
    cum = (np.cumsum(w) - 0.5 * w) / wsum
    return np.interp(np.asarray(pcts, dtype=float) / 100.0, cum, v)


def convex_morpho(image, mask, scale=1.0, do_plot=True, weight_power=0.8,
                  diameter_percentile=99.0, n_boot=500, noise_std=None, 
                  n_sigmas_error = 3,
                  ax=None,
                  cell_size=None):
    """
    Morphological analysis of a masked structure via intensity-weighted moments.

    All quantities (PA, q, diameters) are derived from the same intensity-weighted
    covariance framework, so they are mutually consistent and robust to mask noise.
    Uncertainties are estimated via a noise bootstrap on the image pixels.

    Parameters
    ----------
    image : 2D ndarray
        Image data.
    mask : 2D ndarray
        Binary mask; non-zero pixels define the structure.
    scale : float, optional
        Kept for backward compatibility (unused).
    do_plot : bool, optional
        Overlay axes and convex hull on the image.
    ax : matplotlib Axes, optional
        If given, the overlay (arrows, centroid, hull outline) is drawn onto
        this existing axis instead of creating a standalone figure -- used
        when called from `compute_image_properties` to merge into its own
        image panel rather than popping up a second plot. The caller is
        responsible for drawing the underlying image on `ax` beforehand.
    cell_size : float, optional
        Only used when `ax` is given. Pixel scale (e.g. arcsec/pixel); when
        set, the overlay is drawn in offset coordinates centered on the image
        (`(pixel - shape/2) * cell_size`) to match an `ax` that already shows
        the image with a centered `extent` (as `compute_image_properties`
        does for its offset-coordinate `ax2` panel). Ignored when `ax` is
        None -- the standalone plot always uses raw pixel-index coordinates.
    weight_power : float, optional
        Exponent applied to clipped intensities before normalising weights
        (default 1.0). Higher values emphasise the bright core, reducing
        scatter in PA and q at the cost of shrinking the diameters.
    diameter_percentile : float, optional
        Upper percentile used for diameter estimation (default 99).
        The intensity-weighted (100-p)th and p-th percentiles of the projected
        pixel positions define the diameter, providing outlier resistance
        without hard masking. Because the percentiles are flux-weighted, this
        is a light-weighted extent, not the extent of the mask itself.
    n_boot : int, optional
        Number of noise realisations for uncertainty estimation (default 50).
        Set to 0 to skip bootstrap (uncertainties returned as NaN).

    Returns
    -------
    dict
        "PA_convex"          : position angle in degrees, CCW from +x, in [0, 180)
        "PA_convex_err"      : n_sigmas_error uncertainty on PA (circular std, degrees)
        "q_convex"           : minor/major axis ratio from eigenvalues
        "q_convex_err"       : n_sigmas_error uncertainty on q
        "centroid_convex"    : intensity-weighted centroid [x, y]
        "major_diameter"     : diameter along the major axis (pixels)
        "major_diameter_err" : n_sigmas_error uncertainty on major_diameter (pixels)
        "minor_diameter"     : diameter along the minor axis (pixels)
        "minor_diameter_err" : n_sigmas_error uncertainty on minor_diameter (pixels)
        "q_diam"             : minor_diameter / major_diameter, the axis ratio
                               implied by the measured diameters (always <= 1).
                               Tracks "q_convex" closely for elliptical sources;
                               a large discrepancy flags a non-elliptical shape.
        "q_diam_err"         : n_sigmas_error uncertainty on q_diam
        "axes_swapped"       : True when the largest measured extent lies along
                               the *second* eigenvector, so the axis labels were
                               taken from the extents rather than the eigenvalue
                               order and PA_convex is rotated 90 deg from the
                               major-eigenvector direction. Rare; indicates the
                               light distribution and the spatial extent
                               genuinely disagree about the elongation axis.

    Notes
    -----
    Diameters are intensity-weighted percentile extents along the principal
    axes, using the same weights that build the covariance, and the major/minor
    labels are assigned from those measured extents. `major_diameter >=
    minor_diameter` therefore holds by construction. `q_convex` remains the
    eigenvalue ratio sqrt(lambda_2 / lambda_1) and is unaffected by the label
    assignment; use `q_diam` when consistency with the returned diameters is
    what matters.
    """
    from scipy import ndimage as _ndi

    # -- Keep only the connected component that contains the peak -------------
    # Stray pixels or satellite sources that leaked into the mask are excluded
    # before any geometry is computed.
    labeled, n_comp = _ndi.label(mask)
    if n_comp > 1:
        peak_y, peak_x = np.unravel_index(
            np.nanargmax(np.where(mask, image, np.nan)), image.shape
        )
        clean_mask = labeled == labeled[peak_y, peak_x]
    else:
        clean_mask = mask.astype(bool)

    indices = np.transpose(np.nonzero(clean_mask))
    y, x = indices[:, 0], indices[:, 1]
    points = np.column_stack((x, y)).astype(float)
    intensities = image[clean_mask].flatten()

    def _moments(ints, pts, wp, p_lo, p_hi):
        """
        Core weighted-moment computation.

        Returns (PA, q, D_maj, D_min, centroid, major_vec, minor_vec, q_diam,
        axes_swapped).  Diameters are intensity-weighted percentile extents
        measured with the *same* weights that build the covariance, and the
        axis labels are assigned from those measured extents, so
        `D_maj >= D_min` always holds.
        """
        raw = np.clip(ints, 0.0, None) ** wp
        s = raw.sum()
        w = raw / s if s > 0 else np.full(len(pts), 1.0 / len(pts))
        cen = np.average(pts, weights=w, axis=0)
        ctr = pts - cen
        cov = (ctr * w[:, None]).T @ ctr
        evals, evecs = np.linalg.eigh(cov)
        order = np.argsort(evals)[::-1]
        evals, evecs = evals[order], evecs[:, order]
        mv = evecs[:, 0];  nv = evecs[:, 1]

        # Diameters: intensity-weighted percentile extents along each principal
        # axis.  These MUST be weighted with the same `w` as the covariance --
        # unweighted percentiles measure the (dilated, near-round) mask geometry
        # rather than the light distribution, which made D_min exceed D_maj and
        # left D_min/D_maj inconsistent with q.
        lo_m, hi_m = _weighted_percentile(ctr @ mv, w, [p_lo, p_hi])
        lo_n, hi_n = _weighted_percentile(ctr @ nv, w, [p_lo, p_hi])
        Dm = float(hi_m - lo_m)
        Dn = float(hi_n - lo_n)

        # Label the axes by measured extent, not by eigenvalue order.  For a
        # well-behaved source the two orderings agree; when they disagree (a
        # bright compact core across a fainter, longer structure) the source
        # really is more extended along the second eigenvector, so that becomes
        # the major axis and the PA rotates by 90 deg.  `mv, nv = nv, -mv`
        # preserves the handedness of the frame.
        if Dn > Dm:
            # new mv == old nv, so the major-axis projection extremes become
            # those already measured along nv.
            mv, nv = nv, -mv
            lo_m, hi_m = lo_n, hi_n
            Dm, Dn = Dn, Dm
            axes_swapped = True
        else:
            axes_swapped = False

        # Disambiguate the major-axis sign: point toward the more extended side.
        # Compare |p_lo| vs |p_hi| of the projection - the side with larger
        # absolute reach is where the source extends further (tail, lobe, etc.).
        # For symmetric sources both are equal; use y > 0 as a tiebreaker.
        if abs(lo_m) > abs(hi_m):
            mv = -mv;  nv = -nv
        elif abs(lo_m) == abs(hi_m):
            if mv[1] < 0 or (abs(mv[1]) < 1e-9 and mv[0] < 0):
                mv = -mv;  nv = -nv

        pa = float(np.degrees(np.arctan2(mv[1], mv[0])) % 180.0)
        q  = float(np.sqrt(max(evals[1], 0.0) / max(evals[0], 1e-30)))
        q_diam = float(Dn / Dm) if Dm > 0 else float('nan')
        return pa, q, Dm, Dn, cen, mv, nv, q_diam, axes_swapped

    p_lo = 100.0 - diameter_percentile
    p_hi = diameter_percentile

    # -- Central estimate ------------------------------------------------------
    position_angle_degrees, axis_ratio, major_diameter, minor_diameter, \
        centroid, major_vec, minor_vec, q_diam, axes_swapped = _moments(
            intensities, points, weight_power, p_lo, p_hi
        )

    # -- Noise bootstrap for uncertainties ------------------------------------
    # Use caller-supplied noise_std when available (preferred: the caller has
    # the real background noise before masking).  Fall back to MAD of unmasked
    # pixels; if those are all zero (masked image), use 10% of global std.
    if noise_std is not None and noise_std > 0:
        std_bg = float(noise_std)
    else:
        bg = image[~clean_mask]
        if len(bg) > 10:
            std_bg = float(np.median(np.abs(bg - np.median(bg))) * 1.4826)
        else:
            std_bg = float(np.std(image) * 0.1)

    if n_boot > 0 and std_bg > 0:
        _rng = np.random.default_rng(None)
        _pas, _qs, _dms, _dns, _qds = [], [], [], [], []
        for _ in tqdm(range(n_boot)):
            img_b    = image + _rng.normal(0.0, std_bg, image.shape)
            int_b    = img_b[clean_mask].flatten()
            pa_b, q_b, dm_b, dn_b, _cen_b, _mv_b, _nv_b, qd_b, _sw_b = _moments(
                int_b, points, weight_power, p_lo, p_hi)
            _pas.append(pa_b);  _qs.append(q_b)
            _dms.append(dm_b);  _dns.append(dn_b)
            _qds.append(qd_b)

        # Circular std for PA (period = 180 deg)
        _ang = np.radians(2.0 * np.array(_pas))
        _r   = np.sqrt(np.mean(np.cos(_ang)) ** 2 + np.mean(np.sin(_ang)) ** 2)
        pa_err  = float(n_sigmas_error*90.0 * np.degrees(np.sqrt(-2.0 * np.log(max(_r, 1e-9)))) / 180.0)
        q_err   = float(n_sigmas_error*np.std(_qs))
        dm_err  = float(n_sigmas_error*np.std(_dms))
        dn_err  = float(n_sigmas_error*np.std(_dns))
        qd_err  = float(n_sigmas_error*np.std(_qds))
    else:
        pa_err = q_err = dm_err = dn_err = qd_err = float('nan')

    major_axis = (major_diameter / 2.0) * major_vec
    minor_axis = (minor_diameter / 2.0) * minor_vec

    if do_plot:
        # ax is None: standalone use -- draw our own image + overlay, own figure.
        # ax given (e.g. from compute_image_properties): overlay only, onto the
        # caller's existing axis, which already has the image/contours on it.
        standalone = ax is None
        if standalone:
            fig, ax = plt.subplots(figsize=(5, 5))
            try:
                from astropy.visualization import simple_norm
                from astropy.stats import mad_std
                norm = simple_norm(image, stretch='sqrt', asinh_a=0.02,
                                   vmin=3 * mad_std(image), vmax=0.2 * np.nanmax(image))
                ax.imshow(image, cmap='gray', origin='lower', norm=norm)
            except Exception:
                ax.imshow(image, cmap='gray', origin='lower',
                          vmin=np.percentile(image[clean_mask], 1),
                          vmax=np.percentile(image[clean_mask], 99.5))

        # When merging onto a caller's ax that shows the image in offset
        # (centered, cell_size-scaled) coordinates, plot in that same data
        # space instead of raw pixel indices, so the overlay lines up with
        # the image/contours already drawn there.
        if (not standalone) and (cell_size is not None):
            ny_im, nx_im = image.shape
            c_plot   = np.array([(centroid[0] - nx_im / 2.0) * cell_size,
                                  (centroid[1] - ny_im / 2.0) * cell_size])
            maj_plot = major_axis * cell_size
            min_plot = minor_axis * cell_size
            pts_plot = (points - [nx_im / 2.0, ny_im / 2.0]) * cell_size
        else:
            c_plot, maj_plot, min_plot, pts_plot = centroid, major_axis, minor_axis, points

        # Explicit zorder so the overlay always sits above the image and any
        # contours already drawn on `ax` (contours ~z2, imshow ~z0/1) --
        # hull outline just above contours, arrows/centroid on top of that.
        ax.quiver(c_plot[0], c_plot[1], maj_plot[0], maj_plot[1],
                  angles='xy', scale_units='xy', scale=1, color='limegreen',
                  width=0.009, headwidth=4, headlength=5, zorder=10)
        ax.quiver(c_plot[0], c_plot[1], min_plot[0], min_plot[1],
                  angles='xy', scale_units='xy', scale=1, color='red',
                  width=0.009, headwidth=4, headlength=5, zorder=10)
        ax.scatter(c_plot[0], c_plot[1], color='limegreen', zorder=11)

        # Show convex hull outline for reference (not used in any measurement)
        if len(pts_plot) >= 3:
            try:
                hull = ConvexHull(pts_plot)
                hp = np.vstack([pts_plot[hull.vertices], pts_plot[hull.vertices[0]]])
                ax.plot(hp[:, 0], hp[:, 1], 'w--', lw=0.6, alpha=0.5, zorder=8)
            except Exception:
                pass

        if standalone:
            ax.set_xlabel("$x$ image coordinates")
            ax.set_ylabel("$y$ image coordinates")
            ax.axis('equal')
            plt.tight_layout()
            plt.show()

    return {
        "PA_convex":          position_angle_degrees,
        "PA_convex_err":      pa_err,
        "q_convex":           axis_ratio,
        "q_convex_err":       q_err,
        "centroid_convex":    centroid,
        "major_diameter":     major_diameter,
        "major_diameter_err": dm_err,
        "minor_diameter":     minor_diameter,
        "minor_diameter_err": dn_err,
        "q_diam":             q_diam,
        "q_diam_err":         qd_err,
        "axes_swapped":       axes_swapped,
    }


# import numpy as np
# from scipy.spatial import ConvexHull
# import matplotlib.pyplot as plt
# from itertools import combinations


# def compute_diameters_convex(hull, points, intensities=None):
#     """
#     Compute the major and minor diameters of a structure using its ConvexHull.
#     Now with optional intensity weighting for more robust measurements.

#     Parameters
#     ----------
#     hull : scipy.spatial.ConvexHull
#         A ConvexHull object representing the structure, containing points
#         that define the convex boundary of the structure.
#     points : ndarray
#         Array of (x, y) coordinates.
#     intensities : ndarray, optional
#         Intensity values for each point. If provided, uses high-intensity
#         regions for more robust diameter estimation.

#     Returns
#     -------
#     dict
#         A dictionary containing:
#         - "major_diameter" : float
#             The length of the major diameter (longest distance between any two 
#             points on the convex hull).
#         - "major_points" : tuple of ndarray
#             The two points defining the major diameter.
#         - "minor_diameter" : float
#             The length of the minor diameter (shortest perpendicular distance 
#             between parallel edges of the convex hull).
#         - "minor_points" : tuple of ndarray
#             The two points defining the minor diameter.
#     """
#     # If intensities provided, focus on high-intensity regions
#     if intensities is not None:
#         # Use 90th percentile of intensity for effective boundary
#         percentile = 90
#         sorted_indices = np.argsort(intensities)[::-1]
#         n_select = max(int(len(points) * (percentile / 100.0)), 3)
#         high_intensity_indices = sorted_indices[:n_select]
#         high_intensity_points = points[high_intensity_indices]
        
#         # Compute hull of high-intensity points for more robust diameter estimation
#         if len(high_intensity_points) >= 3:
#             try:
#                 intensity_hull = ConvexHull(high_intensity_points)
#                 hull_points = high_intensity_points[intensity_hull.vertices]
#             except:
#                 # Fall back to original hull if high-intensity hull fails
#                 hull_points = points[hull.vertices]
#         else:
#             hull_points = points[hull.vertices]
#     else:
#         # Original behavior when no intensities provided
#         hull_points = points[hull.vertices]

#     # Calculate Major Diameter (longest distance between hull points)
#     major_diameter = 0
#     major_points = None
#     for p1, p2 in combinations(hull_points, 2):
#         distance = np.linalg.norm(p1 - p2)
#         if distance > major_diameter:
#             major_diameter = distance
#             major_points = (p1, p2)

#     # Calculate Minor Diameter (shortest perpendicular distance between parallel edges)
#     minor_diameter = float('inf')
#     minor_points = None
#     num_hull_points = len(hull_points)

#     for i in range(num_hull_points):
#         # Get two consecutive points forming an edge
#         p1, p2 = hull_points[i], hull_points[(i + 1) % num_hull_points]
#         edge_vector = p2 - p1
#         edge_length = np.linalg.norm(edge_vector)
        
#         # Normalize the edge vector
#         if edge_length == 0:
#             continue
#         edge_normal = np.array([-edge_vector[1], edge_vector[0]]) / edge_length
        
#         # Project all hull points onto the edge normal and find the width
#         distances = np.abs(np.dot(hull_points - p1, edge_normal))
#         max_distance = distances.max()
#         if max_distance < minor_diameter:
#             minor_diameter = max_distance
#             # Points corresponding to the shortest perpendicular projection
#             projections = hull_points[np.abs(distances - max_distance) < 1e-6]
#             if len(projections) >= 2:
#                 minor_points = (projections[0], projections[1])

#     return {
#         "major_diameter": major_diameter,
#         "major_points": major_points,
#         "minor_diameter": minor_diameter,
#         "minor_points": minor_points
#     }


# def _unwrap_angle(angle, reference):
#     """
#     Unwrap angle relative to reference to handle 0/360 discontinuity.
    
#     Parameters
#     ----------
#     angle : float
#         Angle to unwrap (in degrees)
#     reference : float
#         Reference angle (in degrees)
    
#     Returns
#     -------
#     float
#         Unwrapped angle
#     """
#     diff = angle - reference
#     if diff > 180:
#         return angle - 360
#     elif diff < -180:
#         return angle + 360
#     return angle


# def _bootstrap_errors(points, intensities, weights, weight_power, n_bootstrap, 
#                       reference_pa):
#     """
#     Estimate uncertainties via bootstrap resampling.
    
#     Parameters
#     ----------
#     points : ndarray
#         Original (x, y) coordinates of shape (N, 2)
#     intensities : ndarray
#         Intensity values of shape (N,)
#     weights : ndarray
#         Normalized weights of shape (N,)
#     weight_power : float
#         Weight power parameter
#     n_bootstrap : int
#         Number of bootstrap iterations
#     reference_pa : float
#         Reference position angle for unwrapping (degrees)
    
#     Returns
#     -------
#     dict
#         Dictionary containing error estimates:
#         - "PA_convex_err" : float
#         - "q_convex_err" : float
#         - "centroid_convex_err" : ndarray
#         - "major_diameter_err" : float
#         - "minor_diameter_err" : float
#     """
#     n_points = len(points)
    
#     # Storage for bootstrap samples
#     pa_samples = []
#     q_samples = []
#     centroid_samples = []
#     major_diam_samples = []
#     minor_diam_samples = []

#     for iteration in tqdm(range(n_bootstrap)):
#         # Resample points according to their weights
#         # Ensure weights are valid (non-negative, finite, sum to 1)
#         valid_weights = np.copy(weights)
#         if not np.all(np.isfinite(valid_weights)) or np.any(valid_weights < 0):
#             # Fall back to uniform sampling if weights are invalid
#             valid_weights = np.ones(n_points) / n_points
        
#         indices = np.random.choice(n_points, size=n_points, replace=True, p=valid_weights)
#         boot_points = points[indices]
#         boot_intensities = intensities[indices]
        
#         # Recompute weights for bootstrap sample
#         # Handle negative values: shift to ensure all positive
#         boot_intensities_positive = boot_intensities - np.min(boot_intensities) + 1e-10
#         boot_weights = np.power(boot_intensities_positive, weight_power)
        
#         # Check for valid weights
#         if not np.all(np.isfinite(boot_weights)) or np.sum(boot_weights) == 0:
#             continue
            
#         boot_weights = boot_weights / np.sum(boot_weights)
        
#         try:
#             # Check if we have enough unique points for convex hull
#             unique_points = np.unique(boot_points, axis=0)
#             if len(unique_points) < 3:
#                 continue
            
#             # Compute convex hull
#             boot_hull = ConvexHull(unique_points)
            
#             # Weighted centroid
#             boot_centroid = np.average(boot_points, weights=boot_weights, axis=0)
            
#             # Weighted covariance
#             centered = boot_points - boot_centroid
#             cov = np.zeros((2, 2))
#             for i in range(len(boot_points)):
#                 p = centered[i].reshape(-1, 1)
#                 cov += boot_weights[i] * np.dot(p, p.T)
            
#             # Eigenvalues and eigenvectors
#             eigenvalues, eigenvectors = np.linalg.eig(cov)
#             order = np.argsort(eigenvalues)[::-1]
#             eigenvalues = eigenvalues[order]
#             eigenvectors = eigenvectors[:, order]
            
#             # Check for valid eigenvalues
#             if eigenvalues[0] <= 0 or eigenvalues[1] <= 0:
#                 continue
            
#             # PA and axis ratio
#             major_vec = eigenvectors[:, 0]
#             pa = np.degrees(np.arctan2(major_vec[1], major_vec[0]))
#             if pa < 0:
#                 pa += 360
            
#             # Unwrap PA relative to reference (handle 0/360 discontinuity)
#             pa = _unwrap_angle(pa, reference_pa)
            
#             q = np.sqrt(eigenvalues[1] / eigenvalues[0])
            
#             # Diameters (need to map back to original points for intensity lookup)
#             # Find which original indices correspond to unique_points
#             boot_intensities_for_hull = []
#             for up in unique_points:
#                 matching_idx = np.where((boot_points == up).all(axis=1))[0]
#                 if len(matching_idx) > 0:
#                     boot_intensities_for_hull.append(boot_intensities[matching_idx[0]])
#                 else:
#                     boot_intensities_for_hull.append(0.0)
#             boot_intensities_for_hull = np.array(boot_intensities_for_hull)
            
#             diams = compute_diameters_convex(boot_hull, unique_points, 
#                                             boot_intensities_for_hull)
            
#             # Store samples
#             pa_samples.append(pa)
#             q_samples.append(q)
#             centroid_samples.append(boot_centroid)
#             major_diam_samples.append(diams["major_diameter"])
#             minor_diam_samples.append(diams["minor_diameter"])
            
#         except Exception as e:
#             # Skip failed bootstrap iterations
#             continue
    
#     # Check if we have enough successful iterations
#     if len(pa_samples) < 50:
#         import warnings
#         warnings.warn(f"Only {len(pa_samples)} successful bootstrap iterations "
#                      f"out of {n_bootstrap}. Errors may be unreliable.")
    
#     # Convert to arrays
#     pa_samples = np.array(pa_samples)
#     q_samples = np.array(q_samples)
#     centroid_samples = np.array(centroid_samples)
#     major_diam_samples = np.array(major_diam_samples)
#     minor_diam_samples = np.array(minor_diam_samples)
    
#     # Compute standard deviations (1-sigma errors)
#     pa_err = np.std(pa_samples) if len(pa_samples) > 0 else np.nan
#     q_err = np.std(q_samples) if len(q_samples) > 0 else np.nan
#     centroid_err = np.std(centroid_samples, axis=0) if len(centroid_samples) > 0 else np.array([np.nan, np.nan])
#     major_err = np.std(major_diam_samples) if len(major_diam_samples) > 0 else np.nan
#     minor_err = np.std(minor_diam_samples) if len(minor_diam_samples) > 0 else np.nan
    
#     return {
#         "PA_convex_err": pa_err,
#         "q_convex_err": q_err,
#         "centroid_convex_err": centroid_err,
#         "major_diameter_err": major_err,
#         "minor_diameter_err": minor_err
#     }


# # def convex_morpho(image, mask, scale=1.0, do_plot=False, weight_power=1.0,
# #                   n_bootstrap=500, return_errors=True):
# #     """
# #     [Same docstring as before]
# #     """
# #     # Extract y, x coordinates of the structure
# #     indices = np.transpose(np.nonzero(mask))
# #     y, x = indices[:, 0], indices[:, 1]
# #     points = np.column_stack((x, y))  # Convert to (x, y) coordinates
    
# #     # Get intensity values at masked positions
# #     intensities = image[mask].flatten()
    
# #     # Handle negative intensities by shifting to ensure all positive
# #     # This is common in radio images due to noise and cleaning artifacts
# #     intensities_positive = intensities - np.min(intensities) + 1e-10
    
# #     # Apply weight power for emphasis on bright regions
# #     weights = np.power(intensities_positive, weight_power)
    
# #     # Check for valid weights
# #     if not np.all(np.isfinite(weights)) or np.sum(weights) == 0:
# #         # Fall back to uniform weighting if something goes wrong
# #         weights = np.ones(len(intensities)) / len(intensities)
# #     else:
# #         weights = weights / np.sum(weights)  # Normalize weights
    
# #     # Compute ConvexHull
# #     hull = ConvexHull(points)

# #     # Compute intensity-weighted centroid (center of light)
# #     centroid = np.average(points, weights=weights, axis=0)
    
# #     # Center points around weighted centroid
# #     centered_points = points - centroid
    
# #     # Compute intensity-weighted covariance matrix
# #     # This gives us the intensity-weighted second moments
# #     cov_matrix = np.zeros((2, 2))
# #     for i in range(len(centered_points)):
# #         p = centered_points[i].reshape(-1, 1)
# #         cov_matrix += weights[i] * np.dot(p, p.T)
    
# #     # Eigenvalues and eigenvectors
# #     eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

# #     # Order eigenvectors by eigenvalues (largest is major axis)
# #     order = np.argsort(eigenvalues)[::-1]
# #     eigenvalues = eigenvalues[order]
# #     eigenvectors = eigenvectors[:, order]

# #     # Semi-major and minor axes
# #     major_axis_vector = eigenvectors[:, 0]
# #     minor_axis_vector = eigenvectors[:, 1]

# #     # Calculate position angle (anti-clockwise from x-axis)
# #     position_angle = np.arctan2(major_axis_vector[1], major_axis_vector[0])
# #     position_angle_degrees = np.degrees(position_angle)
# #     if position_angle_degrees < 0:
# #         position_angle_degrees += 360

# #     # Calculate axis ratio (minor/major)
# #     axis_ratio = np.sqrt(eigenvalues[1] / eigenvalues[0])

# #     # Compute diameters with intensity weighting
# #     # Use original intensities (not shifted) for diameter calculation
# #     diameters = compute_diameters_convex(hull, points, intensities)
    
# #     # Use actual diameters for visualization vectors (radius = diameter/2)
# #     # Normalize eigenvectors and scale by half the diameter
# #     major_axis = (diameters["major_diameter"] / 2.0) * (major_axis_vector / np.linalg.norm(major_axis_vector))
# #     minor_axis = (diameters["minor_diameter"] / 2.0) * (minor_axis_vector / np.linalg.norm(minor_axis_vector))

# #     if do_plot:
# #         # Plot image
# #         plt.figure(figsize=(5, 5))
        
# #         # Try to use astropy normalization if available
# #         try:
# #             from astropy.visualization import simple_norm
# #             from astropy.stats import mad_std
# #             norm = simple_norm(image, stretch='sqrt', asinh_a=0.02, vmin=3*mad_std(image),
# #                              vmax=0.2*np.nanmax(image))
# #             plt.imshow(image, cmap='gray', origin='lower', norm=norm)
# #         except ImportError:
# #             # Fallback normalization if astropy not available
# #             plt.imshow(image, cmap='gray', origin='lower',
# #                       vmin=np.percentile(image[mask], 1),
# #                       vmax=np.percentile(image[mask], 99.5))
        
# #         # Plot semi-major and minor axes
# #         plt.quiver(
# #             centroid[0], centroid[1], major_axis[0], major_axis[1],
# #             angles='xy', scale_units='xy', scale=1, color='limegreen', 
# #             label='Major axis'
# #         )
# #         plt.quiver(
# #             centroid[0], centroid[1], minor_axis[0], minor_axis[1],
# #             angles='xy', scale_units='xy', scale=1, color='red',
# #             label='Minor axis'
# #         )
    
# #         # Plot centroid (now intensity-weighted)
# #         plt.scatter(centroid[0], centroid[1], color='limegreen', s=50, 
# #                    marker='x', zorder=5, label='Centroid')
    
# #         plt.xlabel("$x$ image coordinates")
# #         plt.ylabel("$y$ image coordinates")
# #         plt.legend(loc='upper right')
# #         plt.axis('equal')
# #         plt.tight_layout()
# #         plt.show()

# #     # Build report dictionary
# #     report = {
# #         "PA_convex": position_angle_degrees,
# #         "q_convex": axis_ratio,
# #         "centroid_convex": centroid,
# #         "major_diameter": diameters["major_diameter"],
# #         "minor_diameter": diameters["minor_diameter"]
# #     }
    
# #     # Bootstrap error estimation
# #     if return_errors and n_bootstrap > 0:
# #         errors = _bootstrap_errors(points, intensities, weights, weight_power, 
# #                                    n_bootstrap, position_angle_degrees)
# #         report.update(errors)
    
# #     return report



# from joblib import Parallel, delayed
# import numpy as np
# from scipy.spatial import ConvexHull
# import matplotlib.pyplot as plt
# from itertools import combinations
# from concurrent.futures import ProcessPoolExecutor
# import multiprocessing as mp


# def _single_bootstrap_iteration(args):
#     """
#     Perform a single bootstrap iteration.
#     Separated into its own function for parallelization.
    
#     Parameters
#     ----------
#     args : tuple
#         (seed, points, intensities, weight_power, reference_pa)
    
#     Returns
#     -------
#     dict or None
#         Dictionary with bootstrap results if successful, None otherwise
#     """
#     seed, points, intensities, weight_power, reference_pa = args
    
#     # Set random seed for this iteration
#     np.random.seed(seed)
    
#     n_points = len(points)
    
#     # Shift intensities to positive
#     intensities_positive = intensities - np.min(intensities) + 1e-10
#     weights = np.power(intensities_positive, weight_power)
    
#     if not np.all(np.isfinite(weights)) or np.sum(weights) == 0:
#         return None
    
#     weights = weights / np.sum(weights)
    
#     try:
#         # Resample
#         indices = np.random.choice(n_points, size=n_points, replace=True, p=weights)
#         boot_points = points[indices]
#         boot_intensities = intensities[indices]
        
#         # Recompute weights
#         boot_intensities_positive = boot_intensities - np.min(boot_intensities) + 1e-10
#         boot_weights = np.power(boot_intensities_positive, weight_power)
        
#         if not np.all(np.isfinite(boot_weights)) or np.sum(boot_weights) == 0:
#             return None
        
#         boot_weights = boot_weights / np.sum(boot_weights)
        
#         # Check for unique points
#         unique_points = np.unique(boot_points, axis=0)
#         if len(unique_points) < 3:
#             return None
        
#         # Convex hull
#         boot_hull = ConvexHull(unique_points)
        
#         # Weighted centroid
#         boot_centroid = np.average(boot_points, weights=boot_weights, axis=0)
        
#         # Weighted covariance
#         centered = boot_points - boot_centroid
#         cov = np.zeros((2, 2))
#         for i in range(len(boot_points)):
#             p = centered[i].reshape(-1, 1)
#             cov += boot_weights[i] * np.dot(p, p.T)
        
#         # Eigenvalues and eigenvectors
#         eigenvalues, eigenvectors = np.linalg.eig(cov)
#         order = np.argsort(eigenvalues)[::-1]
#         eigenvalues = eigenvalues[order]
#         eigenvectors = eigenvectors[:, order]
        
#         if eigenvalues[0] <= 0 or eigenvalues[1] <= 0:
#             return None
        
#         # PA and axis ratio
#         major_vec = eigenvectors[:, 0]
#         pa = np.degrees(np.arctan2(major_vec[1], major_vec[0]))
#         if pa < 0:
#             pa += 360
        
#         # Unwrap PA
#         pa = _unwrap_angle(pa, reference_pa)
#         q = np.sqrt(eigenvalues[1] / eigenvalues[0])
        
#         # Diameters
#         boot_intensities_for_hull = []
#         for up in unique_points:
#             matching_idx = np.where((boot_points == up).all(axis=1))[0]
#             if len(matching_idx) > 0:
#                 boot_intensities_for_hull.append(boot_intensities[matching_idx[0]])
#             else:
#                 boot_intensities_for_hull.append(0.0)
#         boot_intensities_for_hull = np.array(boot_intensities_for_hull)
        
#         diams = compute_diameters_convex(boot_hull, unique_points, 
#                                         boot_intensities_for_hull)
        
#         return {
#             'pa': pa,
#             'q': q,
#             'centroid': boot_centroid,
#             'major_diam': diams["major_diameter"],
#             'minor_diam': diams["minor_diameter"]
#         }
        
#     except Exception as e:
#         return None

# def _bootstrap_errors_parallel_joblib(points, intensities, weights, weight_power, 
#                                       n_bootstrap, reference_pa, n_jobs=-1, 
#                                       verbose=10):
#     """
#     Estimate uncertainties via parallel bootstrap resampling using joblib.
    
#     Parameters
#     ----------
#     n_jobs : int, optional
#         Number of parallel jobs. -1 uses all available cores (default).
#     verbose : int, optional
#         Verbosity level for joblib (0=silent, 10=progress bar).
#     """
#     # Generate random seeds
#     np.random.seed(None)
#     seeds = np.random.randint(0, 2**31 - 1, size=n_bootstrap)
    
#     # Parallel execution with progress bar
#     results = Parallel(n_jobs=n_jobs, verbose=verbose)(
#         delayed(_single_bootstrap_iteration)(
#             (seeds[i], points, intensities, weight_power, reference_pa)
#         )
#         for i in range(n_bootstrap)
#     )
    
#     # Filter out failed iterations
#     results = [r for r in results if r is not None]
    
#     # Check if we have enough successful iterations
#     if len(results) < 50:
#         import warnings
#         warnings.warn(f"Only {len(results)} successful bootstrap iterations "
#                      f"out of {n_bootstrap}. Errors may be unreliable.")
    
#     if len(results) == 0:
#         return {
#             "PA_convex_err": np.nan,
#             "q_convex_err": np.nan,
#             "centroid_convex_err": np.array([np.nan, np.nan]),
#             "major_diameter_err": np.nan,
#             "minor_diameter_err": np.nan
#         }
    
#     # Extract samples
#     pa_samples = np.array([r['pa'] for r in results])
#     q_samples = np.array([r['q'] for r in results])
#     centroid_samples = np.array([r['centroid'] for r in results])
#     major_diam_samples = np.array([r['major_diam'] for r in results])
#     minor_diam_samples = np.array([r['minor_diam'] for r in results])
    
#     return {
#         "PA_convex_err": np.std(pa_samples),
#         "q_convex_err": np.std(q_samples),
#         "centroid_convex_err": np.std(centroid_samples, axis=0),
#         "major_diameter_err": np.std(major_diam_samples),
#         "minor_diameter_err": np.std(minor_diam_samples)
#     }

# def convex_morpho(image, mask, scale=1.0, do_plot=False, weight_power=1.0,
#                   n_bootstrap=500, return_errors=True, n_jobs=-1, 
#                   parallel_backend='joblib', verbose=0):
#     """
#     Perform morphological analysis with parallel bootstrap error estimation.
    
#     Additional Parameters
#     ---------------------
#     n_jobs : int, optional
#         Number of parallel jobs for bootstrap. -1 uses all cores (default).
#         Set to 1 for serial execution.
#     parallel_backend : str, optional
#         'joblib' (default, with progress bar) or 'concurrent' (ProcessPoolExecutor)
#     verbose : int, optional
#         Verbosity level. 0=silent, 10=progress bar (joblib only).
    
#     [Rest of docstring unchanged]
#     """
#     # [Main computation code unchanged until bootstrap section]
    
#     # Extract y, x coordinates of the structure
#     indices = np.transpose(np.nonzero(mask))
#     y, x = indices[:, 0], indices[:, 1]
#     points = np.column_stack((x, y))
    
#     intensities = image[mask].flatten()
#     intensities_positive = intensities - np.min(intensities) + 1e-10
#     weights = np.power(intensities_positive, weight_power)
    
#     if not np.all(np.isfinite(weights)) or np.sum(weights) == 0:
#         weights = np.ones(len(intensities)) / len(intensities)
#     else:
#         weights = weights / np.sum(weights)
    
#     hull = ConvexHull(points)
#     centroid = np.average(points, weights=weights, axis=0)
    
#     centered_points = points - centroid
#     cov_matrix = np.zeros((2, 2))
#     for i in range(len(centered_points)):
#         p = centered_points[i].reshape(-1, 1)
#         cov_matrix += weights[i] * np.dot(p, p.T)
    
#     eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)
#     order = np.argsort(eigenvalues)[::-1]
#     eigenvalues = eigenvalues[order]
#     eigenvectors = eigenvectors[:, order]
    
#     major_axis_vector = eigenvectors[:, 0]
#     minor_axis_vector = eigenvectors[:, 1]
    
#     position_angle = np.arctan2(major_axis_vector[1], major_axis_vector[0])
#     position_angle_degrees = np.degrees(position_angle)
#     if position_angle_degrees < 0:
#         position_angle_degrees += 360
    
#     axis_ratio = np.sqrt(eigenvalues[1] / eigenvalues[0])
#     diameters = compute_diameters_convex(hull, points, intensities)
    
#     major_axis = (diameters["major_diameter"] / 2.0) * (major_axis_vector / np.linalg.norm(major_axis_vector))
#     minor_axis = (diameters["minor_diameter"] / 2.0) * (minor_axis_vector / np.linalg.norm(minor_axis_vector))
    
#     if do_plot:
#         plt.figure(figsize=(5, 5))
#         try:
#             from astropy.visualization import simple_norm
#             from astropy.stats import mad_std
#             norm = simple_norm(image, stretch='sqrt', asinh_a=0.02, 
#                              vmin=3*mad_std(image), vmax=0.2*np.nanmax(image))
#             plt.imshow(image, cmap='gray', origin='lower', norm=norm)
#         except ImportError:
#             plt.imshow(image, cmap='gray', origin='lower',
#                       vmin=np.percentile(image[mask], 1),
#                       vmax=np.percentile(image[mask], 99.5))
        
#         plt.quiver(centroid[0], centroid[1], major_axis[0], major_axis[1],
#                   angles='xy', scale_units='xy', scale=1, color='limegreen', 
#                   label='Major axis')
#         plt.quiver(centroid[0], centroid[1], minor_axis[0], minor_axis[1],
#                   angles='xy', scale_units='xy', scale=1, color='red',
#                   label='Minor axis')
#         plt.scatter(centroid[0], centroid[1], color='limegreen', s=50, 
#                    marker='x', zorder=5, label='Centroid')
#         plt.xlabel("$x$ image coordinates")
#         plt.ylabel("$y$ image coordinates")
#         plt.legend(loc='upper right')
#         plt.axis('equal')
#         plt.tight_layout()
#         plt.show()
    
#     report = {
#         "PA_convex": position_angle_degrees,
#         "q_convex": axis_ratio,
#         "centroid_convex": centroid,
#         "major_diameter": diameters["major_diameter"],
#         "minor_diameter": diameters["minor_diameter"]
#     }
    
#     # Parallel bootstrap error estimation
#     if return_errors and n_bootstrap > 0:
#         if n_jobs == 1:
#             # Serial execution (original implementation)
#             errors = _bootstrap_errors(points, intensities, weights, weight_power, 
#                                       n_bootstrap, position_angle_degrees)
#         elif parallel_backend == 'joblib':
#             errors = _bootstrap_errors_parallel_joblib(
#                 points, intensities, weights, weight_power, 
#                 n_bootstrap, position_angle_degrees, n_jobs, verbose
#             )
#         else:  # 'concurrent'
#             errors = _bootstrap_errors_parallel(
#                 points, intensities, weights, weight_power, 
#                 n_bootstrap, position_angle_degrees, n_jobs
#             )
        
#         report.update(errors)
    
#     return report
