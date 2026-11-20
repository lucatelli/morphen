def run_analysis_list(my_list,ref_residual,ref_image,z,mask_=None,rms=None,
                      sigma=6):
    results_conc_compact = []
    missing_data_im = []
    missing_data_re = []
#     z_d = {'VV705': 0.04019,'UGC5101':0.03937,'UGC8696':0.03734, 'VV250':0.03106}
    if rms is None:
        rms = mad_std(load_fits_data(ref_residual))
    if mask_ is None:
        _, mask_ = mask_dilation(ref_image,#imagelist_vla[k],
                                 sigma=sigma, iterations=2,
                                    dilation_size=None, PLOT=True)

    for i in tqdm(range(len(my_list))):
    # for i in tqdm(range(0,3)):
        crop_image = my_list[i]
        crop_residual = ref_residual#residuallist_vla[k]
        processing_results_model_compact = {} #store calculations only for source
        processing_results_model_compact['#modelname'] = os.path.basename(crop_image)

    #     processing_results_source,mask= measures(crop_image,crop_residual,z=zd,deblend=False,apply_mask=True,
    #                            results_final = processing_results_source,
    #                    plot_catalog = False,bkg_sub=False,bkg_to_sub = None,mask_component=None,
    #                    npixels=500,fwhm=121,kernel_size=121,sigma_mask=7.0,last_level=3.0,
    #                                              iterations=3,dilation_size=7,
    #                    do_PLOT=True,show_figure=True,add_save_name='')
        processing_results_model_compact, _, _ = measures(imagename=crop_image,
                                                residualname=crop_residual,
                                                z=z,rms=rms,
                                                mask = mask_,
#                                                 mask_component=mask_,
                                                do_petro=False,
                                                results_final=processing_results_model_compact,
                                                do_PLOT=True,dilation_size=None,
                                                apply_mask=False)
        results_conc_compact.append(processing_results_model_compact)
#     return(results_conc_compact)
    return(processing_results_model_compact)



r"""
 _____ _ _      
|  ___(_) | ___ 
| |_  | | |/ _ \
|  _| | | |  __/
|_|   |_|_|\___|

  ___                       _   _             
 / _ \ _ __   ___ _ __ __ _| |_(_) ___  _ __  
| | | | '_ \ / _ \ '__/ _` | __| |/ _ \| '_ \ 
| |_| | |_) |  __/ | | (_| | |_| | (_) | | | |
 \___/| .__/ \___|_|  \__,_|\__|_|\___/|_| |_|
      |_|   
#File Operations
"""


def get_list_names(root_path,prefix,which_data,source,sub_comp='',
                   version='v2',cutout_folder=''):
    path = root_path+\
           'data_analysis/LIRGI_sample/analysis_results/processing_images/'+\
           which_data+'/'+source+'/wsclean_images_'+\
           version+'/MFS_images/'+cutout_folder
    pathr = root_path+\
            'data_analysis/LIRGI_sample/analysis_results/processing_images/'+\
            which_data+'/'+source+'/wsclean_images_'+\
            version+'/MFS_residuals/'+cutout_folder

#     prefix = '*-MFS-image.fits'

    imlist = (glob.glob(path+prefix))
    imlist.sort()
    positives = []
    negatives = []
    for it in imlist:
        if '.-multiscale..-' in it:
            negatives.append(it)
#         if "taper_.-" in it:
#             negatives.append(it)
        else:
            positives.append(it)
    negatives.reverse()
#     negatives.sort()

    imagelist_sort = []
    for it in negatives:
        imagelist_sort.append(it)
    for it in positives:
        imagelist_sort.append(it)
    # em = [em[2],em[-1]]
    # i = 0
    # for image in imagelist_sort:
    #     print(i,'>>',os.path.basename(image))
    #     i=i+1

    imagelist_sort_res = []
    if cutout_folder == '':
        replacement = ['-image','-residual']
    else:
        if sub_comp =='':
            replacement = ['-image','-residual']
        else:
            replacement = ['image.cutout'+sub_comp+'.fits',
                           'residual.cutout'+sub_comp+'.fits']
    for i in range(len(imagelist_sort)):
        imagelist_sort_res.append(pathr +
                                  os.path.basename(imagelist_sort[i]).
                                  replace(replacement[0],replacement[1]))
    # i = 0
    # for image in imagelist_sort_res:
    #     print(i,'>>',os.path.basename(image))
    #     i=i+1
    return(np.asarray(imagelist_sort),np.asarray(imagelist_sort_res))




def get_fits_list_names(root_path,prefix='*.fits'):
    imlist = (glob.glob(root_path+prefix))
    imlist.sort()
    i = 0
    for image in imlist:
        print(i,'>>',os.path.basename(image))
        i=i+1
    return(imlist)

# def read_imfit_params(fileParams):
#     dlines = [line for line in open(fileParams) if len(line.strip()) > 0 and line[0] != "#"]
#     values = []
#     temp = []
#     for line in dlines:
#         if line.split()[0] == 'FUNCTION':
#             pass
#         else:
#             temp.append(float(line.split()[1]))
#         if line.split()[0] == 'r_e':
#             #         values['c1'] = {}
#             values.append(np.asarray(temp))
#             temp = []

#     if dlines[-2].split()[1] == 'FlatSky':
#         values.append(np.asarray(float(dlines[-1].split()[1])))
#     return (values)


#: Label columns every row of a decomposition table carries. They are written in
#: this order, ahead of the `measures()` columns, so the table is readable when
#: printed and so `pd.concat` across images/frequencies lines up.
DECOMP_LABEL_COLUMNS = ('imagename', 'freq', 'kind', 'comp_ID', 'region_ID',
                        'domain', 'is_compact', 'has_compact', 'n_comps_region',
                        'region_area_completeness', 'region_flux_completeness',
                        'region_overlap_fraction')

#: Allowed values of the `kind` column.
#:
#:   total          the whole source, measured on the data
#:   component      one fitted model component
#:   compact_sum    the sum of the components listed in `comp_ids`
#:   diffuse_sum    data minus `compact_sum` (DATA-driven, not the diffuse model)
#:   region_data    one detected region, measured on the data
#:   region_diffuse one detected region, minus that region's compact components
#:   decomposition  the `dec_*` scalars from `plot_decomp_results`
DECOMP_KINDS = ('total', 'component', 'compact_sum', 'diffuse_sum',
                'region_data', 'region_diffuse', 'decomposition')


def decomp_rows(props, kind, imagename=None, freq=None, domain='data',
                comp_ID=0, region_ID=0, is_compact=None, has_compact=None,
                n_comps_region=None, region_area_completeness=None,
                region_flux_completeness=None, region_overlap_fraction=None):
    """
    Tag one `measures()` result (or a list of them) with the decomposition
    labels, returning a list of flat dicts ready for `pd.DataFrame`.

    This is deliberately a labelling function and nothing more -- it does not
    measure, derive or rename anything. Every existing DataFrame keeps its own
    values; the long table is these same rows with `kind`/`comp_ID`/`domain`
    attached so they can be concatenated and filtered instead of being pulled
    apart by hand.

    Parameters
    ----------
    props : dict or sequence of dict
        One or more `measures()` property dicts.
    kind : str
        One of `DECOMP_KINDS`.
    comp_ID, region_ID : int
        1-indexed; 0 means "not a single component" / "not a single region".
    domain : str
        'data', 'conv' or 'deconv'.
    is_compact : bool, optional
        Component rows only: was this ID listed in `comp_ids`.
    has_compact, n_comps_region : optional
        Region rows only.
    region_area_completeness, region_flux_completeness : float, optional
        Region rows only, and the same for every region of one image: the
        fraction of the reference aperture's area and flux that the union of all
        region apertures covered. Below 1.0 means part of `mask_region` is
        disconnected from every deblended core, which is what explains
        sum(region flux) < total_flux_mask.
    region_overlap_fraction : float, optional
        Region rows only: how much the region apertures overlap each other, as a
        fraction of their union. 0.0 means a clean partition. Read it together
        with the completeness above -- apertures can cover the reference mask
        fully and still double-count, which shows up here and nowhere else.

    Returns
    -------
    list of dict
    """
    if isinstance(props, dict):
        props = [props]
    rows = []
    for entry in props:
        row = {'kind': kind, 'comp_ID': int(comp_ID),
               'region_ID': int(region_ID), 'domain': domain,
               'is_compact': is_compact, 'has_compact': has_compact,
               'n_comps_region': n_comps_region,
               'region_area_completeness': region_area_completeness,
               'region_flux_completeness': region_flux_completeness,
               'region_overlap_fraction': region_overlap_fraction}
        if imagename is not None:
            row['imagename'] = os.path.basename(imagename)
        if freq is not None:
            row['freq'] = freq
        # The measured columns win over anything above sharing their name, so a
        # `measures()` result that already carries e.g. `comp_ID` is not
        # silently overwritten by the label.
        row.update(entry)
        rows.append(row)
    return rows


def assemble_decomp_table(rows):
    """
    Turn accumulated `decomp_rows` output into a DataFrame with the label
    columns first and `comp_ID`/`region_ID` as integers.

    Empty input gives an empty frame carrying just the label columns, so callers
    can concatenate unconditionally.
    """
    if not rows:
        return pd.DataFrame(columns=list(DECOMP_LABEL_COLUMNS))
    df = pd.DataFrame(rows)
    for col in ('comp_ID', 'region_ID'):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)
    lead = [c for c in DECOMP_LABEL_COLUMNS if c in df.columns]
    return df[lead + [c for c in df.columns if c not in lead]]


def compute_model_properties(model_list,  # the model list of each component
                             which_model,  # `convolved` or `deconvolved`?
                             residualname,
                             rms,  # the native rms from the data itself.
                             mask_region = None,
                             z=None,
                             sigma_mask=5.0,
                             last_level = 1.5,
                             vmin_factor=1.0,
                             iterations = 2,
                             verbose=0):
    """
    Helper function function to calculate model component properties.

    For each model component fitted to a data using the sercic profile, perform morphometry on each
    component image, both deconvolved and convolved images.
    """
    model_properties = {}
    kk = 1
    # if which_model == 'conv':
    #     rms_model = rms
    #     dilation_size = 2
    # if which_model == 'deconv':
    #     rms_model = rms/len(model_list)
    #     dilation_size = 2
    if which_model == 'conv':
        dilation_size = dilation_size = get_dilation_size(model_list[0])
    if which_model == 'deconv':
        dilation_size = 2

    if verbose >= 1:
        print(f' ++==>> Dilation size = {dilation_size}')
        print(f' ++==>> Iterations = {iterations}')
        show_figure = True
    else:
        show_figure = False

    for model_component in model_list:
        try:
            print('Computing properties of model component: ', os.path.basename(model_component))
            model_component_data = load_fits_data(model_component)
            if which_model == 'conv':
                rms_model = rms
            if which_model == 'deconv':
                # rms_model = mad_std(model_component_data) + rms
                """testing"""
                rms_model = rms
            # if verbose >= 1:
            #     print(' --==>> STD RMS of model component: ', rms_model)
            #     print(' --==>> STD RMS of model bkg: ', rms)
            #     print(' --==>> Ratio rms_model/rms_bkg: ', rms_model/rms)

            _, mask_component = mask_dilation(model_component,
                                            rms=rms_model,
                                            sigma=sigma_mask, 
                                            dilation_size=dilation_size,
                                            iterations=iterations, 
                                            PLOT=True)
            # """#testing"""
            # if which_model == 'deconv':
            #     if np.nansum(mask_component) == 0 and mask_region is not None:
            #         mask_component = mask_region
            # """---"""
            
            """#testing"""
            if np.nansum(mask_component) == 0 and mask_region is not None:
                # _, mask_component = mask_dilation(model_component,
                #                                 rms=rms_model,
                #                                 sigma=3.0, 
                #                                 dilation_size=dilation_size,
                #                                 iterations=iterations, 
                #                                 PLOT=True)
                # if np.nansum(mask_component) == 0 and mask_region is not None:
                #     _, mask_component = mask_dilation(model_component,
                #                                     rms=mad_std(model_component_data),
                #                                     sigma=1.0, 
                #                                     dilation_size=dilation_size,
                #                                     iterations=iterations, 
                #                                     PLOT=True)    
                # mask_component = mask_component * mask_region
                mask_component = mask_region
            """---"""


            # print('number of pixesl in model mask = ', np.sum(mask_component))
            # print('number of pixesl in model mask * regions mask = ', np.sum(mask_component*mask_region))
            # print('number of pixesl in regions mask = ', np.sum(mask_region))

            properties, _, _ = measures(imagename=model_component,
                                            residualname=residualname,
                                            z=z,
                                            sigma_mask=sigma_mask,
                                            last_level = last_level,
                                            vmin_factor=vmin_factor,
                                            dilation_size=dilation_size,
                                            mask = mask_region,
                                            mask_component=mask_component * mask_region,
                                            show_figure = show_figure,
                                            apply_mask=False,
                                            # data_2D=load_fits_data(model_component),
                                            rms=rms_model)

            model_properties[f"model_c_{which_model}_{kk}_props"] = properties.copy()
            # model_properties[f"model_c_{which_model}_{kk}_props"]['model_file'] = model_component
            model_properties[f"model_c_{which_model}_{kk}_props"]['comp_ID'] = kk
            model_properties[f"model_c_{which_model}_{kk}_props"][
                'model_file'] = os.path.basename(model_component)
            kk = kk + 1
        except:
            empty_properties = {key: np.nan for key in model_properties[f"model_c_{which_model}_{1}_props"].keys()}
            model_properties[f"model_c_{which_model}_{kk}_props"] = empty_properties.copy()
            model_properties[f"model_c_{which_model}_{kk}_props"]['comp_ID'] = kk
            model_properties[f"model_c_{which_model}_{kk}_props"][
                'model_file'] = os.path.basename(model_component)
            print('Error computing properties of model component: ', os.path.basename(model_component))
            kk = kk + 1

    return (model_properties)



def format_nested_data(nested_data):
    """
    Format to a data frame a list of nested dictionaries.

    Parameters
    ----------
    nested_data : list of dictionaries
        The list of dictionaries to be formatted.

    Returns
    -------
    df : pandas.DataFrame
        The formatted data frame.
    """
    processed_data = []
    for item in nested_data:
        for model_name, props in item.items():
            props['model_name'] = model_name  # Add the model name to the dictionary
            processed_data.append(props)
    df = pd.DataFrame(processed_data)
    return(df)


def adjust_arrays(a, b):
    len_a, len_b = len(a), len(b)

    # If a is shorter, pad it with NaN values
    if len_a < len_b:
        a = np.pad(a, (0, len_b - len_a), 'constant', constant_values=np.nan)
    # If b is shorter, pad it with NaN values
    elif len_b < len_a:
        b = np.pad(b, (0, len_a - len_b), 'constant', constant_values=np.nan)

    return a, b


# import numpy as np
# import pandas as pd
PRINT_INFO_LEVELS = ('simple', 'full')

def _factor(f, row):
    """Resolve a conversion factor: a number is constant, a str is a column name."""
    return row[f] if isinstance(f, str) else f


def _has_err(row, ecol):
    """True if the error column exists for this row and is a finite number."""
    if ecol is None:
        return False
    e = row[ecol]
    try:
        return bool(np.isfinite(float(e)))
    except (TypeError, ValueError):
        return False


def _check_print_info(print_info):
    if print_info not in PRINT_INFO_LEVELS:
        raise ValueError(f"print_info must be one of {PRINT_INFO_LEVELS}, "
                         f"got {print_info!r}")

def _with_derived_columns(df):
    """
    Return a copy of `df` with derived columns added when they are missing:

      local_snr     = peak_of_flux / local_rms
      local_snr_err = peak_error   / local_rms   (noise on local_rms ignored)

    `local_rms` is the noise measured inside the region itself (residual map,
    or sigma-clipped image when no residual was given), so this SNR differs
    from `peak_snr`, which uses the global map std. A `local_snr` column
    already present in `df` is left untouched.
    """
    out = df.copy()
    if ('local_snr' not in out.columns
            and {'peak_of_flux', 'local_rms'} <= set(out.columns)):
        lrms = pd.to_numeric(out['local_rms'], errors='coerce')
        lrms = lrms.where(lrms > 0)                     # 0/neg rms -> NaN SNR
        out['local_snr'] = pd.to_numeric(out['peak_of_flux'], errors='coerce') / lrms
        if 'peak_error' in out.columns:
            out['local_snr_err'] = (pd.to_numeric(out['peak_error'], errors='coerce')
                                    / lrms)
    return out


# Columns shown with print_info='simple', in display order.
SIMPLE_MORPHO_COLS = [
    'total_flux_mask',
    'peak_of_flux',
    'C50radii_major',
    'C50radii_minor',
    'C95radii_major',
    'C95radii_minor',
    'local_rms',
    'local_snr',
]

def build_morpho_specs(print_info='full'):
    """
    Declarative spec for the final measurement table. One entry per printed
    column. Keys:
      'col'   : value column name (required)
      'err'   : uncertainty column name (optional)
      'units' : list of (label, factor). factor is a constant (e.g. 1000 for
                Jy->mJy) or a column name resolved per row ('cell_size' for
                pixel->arcsec). Use 1 for pass-through / already-correct units,
                and None for text columns (printed/saved as-is).
      'fmt'   : format spec (default '.3f')

    print_info : 'full'   -> every column below (the original table).
                 'simple' -> only SIMPLE_MORPHO_COLS: total flux, peak,
                             a/b at 50% and 95%, local_rms and local_snr.

    Conversions:
      cell_size : pixel -> arcsec  (radial / size quantities)
      1000      : Jy -> mJy        (fluxes, peaks, rms)
      1         : identity         (dimensionless, degrees, counts, already-arcsec)
      None      : text (no conversion, no formatting)

    Sizes:
      CXXradii        circularised radius sqrt(A/pi) = sqrt(a*b) -- NOT a semi-axis
      CXXradii_major  semi-major axis  a = R / sqrt(q_eff)
      CXXradii_minor  semi-minor axis  b = R * sqrt(q_eff)
    Shape (all errors are 1-sigma, from the beam-correlated bootstrap):
      q_eff +/- sigma_q_eff   axis ratio used for the major/minor conversion
      q_convex +/- q_convex_err, PA_convex +/- PA_convex_err
      PA_convex is CCW from the image +x axis, NOT east of north (for the usual
      CDELT1 < 0 orientation, PA_EoN = (PA_convex - 90) mod 180).
      sigma_q_method  'bootstrap_beam' / 'bootstrap_white' / 'rule'
    Noise:
      local_rms   noise inside the measured region (see compute_local_rms)
      local_snr   peak_of_flux / local_rms, derived by _with_derived_columns
    """
    _check_print_info(print_info)
    full = [
        {'col': 'total_flux_mask',        'err': 'total_flux_error',
         'units': [('mJy', 1000)],                          'fmt': '.3f'},
        {'col': 'peak_of_flux',           'err': 'peak_error',
         'units': [('mJy/bm', 1000)],                       'fmt': '.3f'},
        {'col': 'C50radii',               'err': 'C50radii_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'C50radii_major',         'err': 'C50radii_major_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'C50radii_minor',         'err': 'C50radii_minor_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'C95radii',               'err': 'C95radii_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'C95radii_major',         'err': 'C95radii_major_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'C95radii_minor',         'err': 'C95radii_minor_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'q_eff',                  'err': 'sigma_q_eff',
         'units': [('', 1)],                                'fmt': '.3f'},
        {'col': 'sigma_q_method',
         'units': [('', None)]},
        {'col': 'q_convex',               'err': 'q_convex_err',
         'units': [('', 1)],                                'fmt': '.3f'},
        {'col': 'PA_convex',              'err': 'PA_convex_err',
         'units': [('deg', 1)],                             'fmt': '.2f'},
        {'col': 'convex_major_diameter',  'err': 'convex_major_diameter_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'convex_minor_diameter',  'err': 'convex_minor_diameter_err',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'C1',                     'err': 'C1_err',
         'units': [('', 1)],                                'fmt': '.3f'},
        {'col': 'C2',                     'err': 'C2_err',
         'units': [('', 1)],                                'fmt': '.3f'},
        {'col': 'TvPA',
         'units': [('deg', 1)],                             'fmt': '.2f'},
        {'col': 'Tvlength',
         'units': [('arcsec', 'cell_size')],                'fmt': '.4f'},
        {'col': 'rms_residual',
         'units': [('mJy/bm', 1000)],                       'fmt': '.4f'},
        {'col': 'local_rms',
         'units': [('mJy/bm', 1000)],                       'fmt': '.4f'},
        {'col': 'number_of_beams_in_mask',
         'units': [('', 1)],                                'fmt': '.1f'},
        {'col': 'bmajor',
         'units': [('arcsec', 1)],                          'fmt': '.3f'},
    ]
    if print_info == 'full':
        return full

    by_col = {s['col']: s for s in full}
    by_col['local_snr'] = {'col': 'local_snr', 'err': 'local_snr_err',
                           'units': [('', 1)], 'fmt': '.1f'}
    return [by_col[c] for c in SIMPLE_MORPHO_COLS]


def _iter_output_columns(df, specs):
    """Yield (header, value_col, err_col, factor, fmt) for each output column."""
    for s in specs:
        if s['col'] not in df.columns:
            continue
        err = s.get('err') if s.get('err') in df.columns else None
        fmt = s.get('fmt', '.3f')
        for label, f in s['units']:
            header = f"{s['col']} ({label})" if label else s['col']
            yield header, s['col'], err, f, fmt

def _resolve_specs(print_info, specs, coord_specs):
    """Defaults for specs/coord_specs. Explicitly passed specs always win.
    'simple' shows no RA/Dec columns unless coord_specs is passed."""
    _check_print_info(print_info)
    if specs is None:
        specs = build_morpho_specs(print_info)
    if coord_specs is None:
        coord_specs = build_coord_specs() if print_info == 'full' else []
    return specs, coord_specs

def print_morpho_table(df, specs=None, coord_specs=None,
                       id_col='comp_ID', freq_col='freq', print_info='full'):
    """
    Pretty-print selected quantities with per-column unit conversions and
    paired uncertainties. comp_ID == 0 (whole source) shown as 'whole'.
    RA/Dec position columns (coord_specs) are appended as formatted strings.
    A NaN uncertainty (e.g. bootstrap switched off) prints the value alone.

    print_info : 'full'   -> all quantities + RA/Dec columns (original table).
                 'simple' -> total_flux_mask, peak_of_flux, C50/C95 major and
                             minor, local_rms, local_snr; no RA/Dec columns
                             (pass coord_specs explicitly to add them back).
    """
    import mlibs
    specs, coord_specs = _resolve_specs(print_info, specs, coord_specs)
    df = _with_derived_columns(df)
    coord_specs = [c for c in coord_specs if c['col'] in df.columns]

    t = mlibs.PrettyTable()
    fields = ['comp_ID', 'freq (GHz)']
    for header, _, _, _, _ in _iter_output_columns(df, specs):
        fields.append(header)
    for c in coord_specs:
        fields.append(c['header'])
    t.field_names = fields

    ids = pd.to_numeric(df[id_col], errors='coerce')
    for i, (_, row) in enumerate(df.iterrows()):
        label = 'whole' if ids.iloc[i] == 0 else f"{int(ids.iloc[i])}"
        values = [label, f"{row[freq_col] / 1e9:.4f}"]
        for _, vcol, ecol, f, fmt in _iter_output_columns(df, specs):
            if f is None:                               # text column
                values.append(str(row[vcol]))
                continue
            k = _factor(f, row)
            v = row[vcol] * k
            if _has_err(row, ecol):
                values.append(f"{v:{fmt}} ± {row[ecol] * k:{fmt}}")
            else:
                values.append(f"{v:{fmt}}")
        for c in coord_specs:
            values.append(_format_coord_cell(row[c['col']], c))
        t.add_row(values)

    print(t)
    return t


def save_morpho_csv(df, path, specs=None, coord_specs=None,
                    id_col='comp_ID', freq_col='freq', coord_split=True,
                    print_info='full'):
    """
    Save the converted quantities to CSV. Values and errors are numeric
    (converted to the target units); value and error get separate columns.
    Header labels carry the unit, matching the printed table. Text columns
    (factor None, e.g. sigma_q_method) are written as-is.

    coord_split : if True, RA/Dec columns are written as two numeric columns
                  ('<header> RA (deg)', '<header> Dec (deg)') in decimal degrees;
                  if False, one formatted string column per coordinate.
    print_info  : 'full' or 'simple', same column selection as print_morpho_table.
    """
    specs, coord_specs = _resolve_specs(print_info, specs, coord_specs)
    df = _with_derived_columns(df)
    coord_specs = [c for c in coord_specs if c['col'] in df.columns]

    out = pd.DataFrame()
    out['comp_ID'] = pd.to_numeric(df[id_col], errors='coerce').astype('Int64')
    out['freq_GHz'] = df[freq_col].values / 1e9

    for header, vcol, ecol, f, _ in _iter_output_columns(df, specs):
        if f is None:                                   # text column
            out[header] = df[vcol].values
            continue
        k = df[f].values if isinstance(f, str) else f
        out[header] = pd.to_numeric(df[vcol], errors='coerce').values * k
        if ecol is not None:
            out[f"{header} err"] = pd.to_numeric(df[ecol], errors='coerce').values * k

    for c in coord_specs:
        parsed = df[c['col']].apply(_parse_radec)
        if coord_split:
            out[f"{c['header']} RA (deg)"]  = [p[0] for p in parsed]
            out[f"{c['header']} Dec (deg)"] = [p[1] for p in parsed]
        else:
            out[c['header']] = df[c['col']].apply(lambda v: _format_coord_cell(v, c))

    out.to_csv(path, index=False)
    print(f"Saved {len(out)} rows x {out.shape[1]} cols -> {path}")
    return out


# ---------------------------------------------------------------------------
# LaTeX export (transposed: one row per component, S column per frequency)
# ---------------------------------------------------------------------------

def build_latex_specs():
    """
    Column layout for the compact LaTeX table. Two kinds of columns:

    Per-frequency (one output column per frequency present):
      {'kind': 'perfreq', 'col', 'err', 'factor', 'fmt',
       'header': callable(freq_ghz)->str}

    Per-component (single value, taken from the lowest-frequency row):
      {'kind': 'single', 'col', 'err'(opt), 'factor', 'fmt', 'header': str}

    factor: constant, or column name resolved on the source row.
    R50/R95 are circularised radii sqrt(a*b); a/b are the semi-major/minor axes.
    """
    return [
        {'kind': 'perfreq', 'col': 'total_flux_mask', 'err': 'total_flux_error',
         'factor': 1000, 'fmt': '.1f',
         'header': lambda ghz: rf'$S_{{{_band(ghz)}}}$'},
        {'kind': 'single', 'col': 'C50radii', 'err': 'C50radii_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$R_{50}$'},
        {'kind': 'single', 'col': 'C50radii_major', 'err': 'C50radii_major_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$a_{50}$'},
        {'kind': 'single', 'col': 'C50radii_minor', 'err': 'C50radii_minor_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$b_{50}$'},
        {'kind': 'single', 'col': 'C95radii', 'err': 'C95radii_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$R_{95}$'},
        {'kind': 'single', 'col': 'C95radii_major', 'err': 'C95radii_major_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$a_{95}$'},
        {'kind': 'single', 'col': 'C95radii_minor', 'err': 'C95radii_minor_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$b_{95}$'},
        {'kind': 'single', 'col': 'q_eff', 'err': 'sigma_q_eff',
         'factor': 1, 'fmt': '.2f', 'header': r'$q_{\mathrm{eff}}$'},
        {'kind': 'single', 'col': 'PA_convex', 'err': 'PA_convex_err',
         'factor': 1, 'fmt': '.1f', 'header': r'PA'},
        {'kind': 'single', 'col': 'convex_minor_diameter',
         'err': 'convex_minor_diameter_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$D_{\mathrm{min}}$'},
        {'kind': 'single', 'col': 'convex_major_diameter',
         'err': 'convex_major_diameter_err',
         'factor': 'cell_size', 'fmt': '.3f', 'header': r'$D_{\mathrm{maj}}$'},
        {'kind': 'single', 'col': 'number_of_beams_in_mask',
         'factor': 1, 'fmt': '.1f', 'header': r'$N_{\Theta}$'},
        {'kind': 'single', 'col': 'TvPA',
         'factor': 1, 'fmt': '.2f', 'header': r'$\Theta$'},
    ]

def _band(freq_ghz):
    """Rough radio band letter from frequency in GHz (for S-column subscripts)."""
    bands = [(4, 8, 'C'), (8, 12, 'X'), (12, 18, 'Ku'),
             (18, 27, 'K'), (27, 40, 'Ka'), (1, 2, 'L'), (2, 4, 'S')]
    for lo, hi, name in bands:
        if lo <= freq_ghz < hi:
            return rf'\mathrm{{{name}}}'
    return rf'{freq_ghz:.0f}'


def to_latex_table(df, latex_specs=None, id_col='comp_ID', freq_col='freq',
                   source_names=None, table_title=None, path=None):
    """
    Build a compact LaTeX tabular: one row per component, per-frequency fluxes
    spread across S-columns, frequency-independent quantities taken from the
    lowest-frequency row of each component. Columns missing from `df` are
    skipped.

    source_names : dict {comp_id: 'Label'} for the Source column. comp_id is
                   matched numerically (0, 1, 2, ...).
    table_title  : optional \\multicolumn banner row.
    path         : if given, write the .tex; otherwise just return the string.
    """
    if latex_specs is None:
        latex_specs = build_latex_specs()
    latex_specs = [s for s in latex_specs if s['col'] in df.columns]

    ids = pd.to_numeric(df[id_col], errors='coerce')
    freqs_ghz = np.sort((df[freq_col] / 1e9).unique())
    comp_ids = sorted(ids.dropna().unique())

    perfreq = [s for s in latex_specs if s['kind'] == 'perfreq']
    singles = [s for s in latex_specs if s['kind'] == 'single']

    # ---- header row ----
    heads = ['Source']
    for s in perfreq:
        for ghz in freqs_ghz:
            heads.append(s['header'](ghz))
    for s in singles:
        heads.append(s['header'])
    ncols = len(heads)

    def fac(f, row):
        return row[f] if isinstance(f, str) else f

    def fmt_pm(row, s):
        k = fac(s['factor'], row)
        v = row[s['col']] * k
        if s.get('err') and s['err'] in df.columns and _has_err(row, s['err']):
            return rf"{v:{s['fmt']}} $\pm$ {row[s['err']] * k:{s['fmt']}}"
        return f"{v:{s['fmt']}}"

    body = []
    for cid in comp_ids:
        sub = df[ids == cid]
        sub_ghz = (sub[freq_col] / 1e9)
        low_row = sub.loc[sub_ghz.idxmin()]     # lowest-freq row for singles

        label = None
        if source_names is not None:
            label = source_names.get(cid, source_names.get(int(cid),
                                     source_names.get(str(int(cid)))))
        if label is None:
            label = 'whole' if cid == 0 else f'{int(cid)}'

        cells = [label]
        # per-frequency flux columns
        for s in perfreq:
            for ghz in freqs_ghz:
                match = sub[np.isclose(sub_ghz.values, ghz)]
                if len(match):
                    r = match.iloc[0]
                    k = fac(s['factor'], r)
                    v = r[s['col']] * k
                    prefix = rf"({ghz:.1f})~"
                    if s.get('err') and s['err'] in df.columns and _has_err(r, s['err']):
                        cells.append(rf"{prefix}{v:{s['fmt']}} $\pm$ {r[s['err']] * k:{s['fmt']}}")
                    else:
                        cells.append(rf"{prefix}{v:{s['fmt']}}")
                else:
                    cells.append('--')
        # frequency-independent columns (from lowest-freq row)
        for s in singles:
            cells.append(fmt_pm(low_row, s))
        body.append(' & '.join(cells) + r' \\')

    # ---- assemble ----
    lines = [r'\begin{tabular}{' + 'c' + 'r' * (ncols - 1) + '}',
             r'\hline\hline']
    if table_title:
        lines.append(rf'\multicolumn{{{ncols}}}{{c}}{{\it {table_title}}} \\ \hline')
    lines.append(' & '.join(heads) + r' \\ \hline')
    lines.extend(body)
    lines.append(r'\hline')
    lines.append(r'\end{tabular}')
    tex = '\n'.join(lines)

    if path:
        with open(path, 'w') as fh:
            fh.write(tex + '\n')
        print(f"Wrote LaTeX table -> {path}")
    return tex


# ---------------------------------------------------------------------------
# RA/Dec coordinate columns  (values like (233.73842514, 23.50319123) in deg)
# ---------------------------------------------------------------------------

def _parse_radec(val):
    """Return (ra_deg, dec_deg) from a tuple/list/np.array or a '(ra, dec)' string."""
    if val is None or (np.isscalar(val) and pd.isna(val)):
        return (np.nan, np.nan)
    if isinstance(val, str):
        import ast
        try:
            val = ast.literal_eval(val.strip())
        except (ValueError, SyntaxError):
            nums = [float(x) for x in
                    val.replace('(', ' ').replace(')', ' ').replace(',', ' ').split()]
            val = nums
    ra, dec = float(val[0]), float(val[1])
    return (ra, dec)


def _radec_hmsdms(ra_deg, dec_deg, sep='hmsdms', precision=2):
    """Format decimal-degree RA/Dec as sexagesimal. Uses astropy if present,
    else a pure-python fallback."""
    if not np.isfinite(ra_deg) or not np.isfinite(dec_deg):
        return '--'
    try:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        c = SkyCoord(ra=ra_deg * u.deg, dec=dec_deg * u.deg, frame='icrs')
        return c.to_string('hmsdms', precision=precision)
    except ImportError:
        # RA in hours:min:sec, Dec in deg:arcmin:arcsec
        rh = ra_deg / 15.0
        rh_h = int(rh); rh_m = int((rh - rh_h) * 60)
        rh_s = (rh - rh_h - rh_m / 60) * 3600
        sign = '-' if dec_deg < 0 else '+'
        d = abs(dec_deg)
        dd = int(d); dm = int((d - dd) * 60); ds = (d - dd - dm / 60) * 3600
        return (f"{rh_h:02d}h{rh_m:02d}m{rh_s:0{3+precision}.{precision}f}s "
                f"{sign}{dd:02d}d{dm:02d}m{ds:0{3+precision}.{precision}f}s")


def build_coord_specs():
    """
    RA/Dec position columns. 'style': 'sexagesimal' (hms/dms) or 'decimal'.
    Printed/saved as a single 'RA Dec' string per column (no unit conversion,
    no uncertainty pairing).
    """
    return [
        {'col': 'peak_position_Ra_Dec',            'header': 'Peak RA/Dec',
         'style': 'sexagesimal', 'precision': 2},
        {'col': 'centroid_convex_position_Ra_Dec', 'header': 'Centroid RA/Dec',
         'style': 'sexagesimal', 'precision': 2},
    ]


def _format_coord_cell(val, spec):
    ra, dec = _parse_radec(val)
    if spec.get('style', 'sexagesimal') == 'decimal':
        if not (np.isfinite(ra) and np.isfinite(dec)):
            return '--'
        p = spec.get('precision', 6)
        return f"{ra:.{p}f} {dec:+.{p}f}"
    return _radec_hmsdms(ra, dec, sep='hmsdms', precision=spec.get('precision', 2))

def append_sed_to_df(df, comp_id, result_sed=None, params_sed=None,
                     id_col='comp_ID', prefix_lmfit='sed_', prefix_mcmc='sed_mcmc_'):
    """
    Append one component's SED fit onto its rows (all frequencies) in df.
    Call once per component with the results from your external
    do_fit_spec_RC_linear call.
    """
    mask = df[id_col] == comp_id
    if not mask.any():
        raise ValueError(f"No rows with {id_col} == {comp_id!r}. "
                         f"Found: {df[id_col].unique()}")

    def setcol(col, val):
        if col not in df.columns:
            df[col] = np.nan
        df.loc[mask, col] = val

    if result_sed is not None:
        for name, p in result_sed.params.items():
            setcol(f'{prefix_lmfit}{name}',     p.value)
            setcol(f'{prefix_lmfit}{name}_err', p.stderr)

    if params_sed is not None:
        for name, d in params_sed.items():
            setcol(f'{prefix_mcmc}{name}',       d['best'])
            setcol(f'{prefix_mcmc}{name}_lower', d['lower'])
            setcol(f'{prefix_mcmc}{name}_upper', d['upper'])

    return df

