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


def compute_model_properties(model_list,  # the model list of each component
                             which_model,  # `convolved` or `deconvolved`?
                             residualname,
                             rms,  # the native rms from the data itself.
                             mask_region = None,
                             z=None,
                             sigma_mask=5.0,
                             last_level = 1.0,
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
                rms_model = mad_std(model_component_data) + rms
            # if verbose >= 1:
            #     print(' --==>> STD RMS of model component: ', rms_model)
            #     print(' --==>> STD RMS of model bkg: ', rms)
            #     print(' --==>> Ratio rms_model/rms_bkg: ', rms_model/rms)

            _, mask_component = mask_dilation(model_component,
                                            rms=rms_model,
                                            sigma=sigma_mask, 
                                            dilation_size=dilation_size,
                                            iterations=iterations, 
                                            PLOT=show_figure)
            
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
                                            mask_component=mask_component,
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
            empty_properties = {key: np.nan for key in model_properties[f"model_c_{which_model}_{kk-1}_props"].keys()}
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
