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



"""
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
