
"""
 ____                           
/ ___|  ___  _   _ _ __ ___ ___ 
\___ \ / _ \| | | | '__/ __/ _ \
 ___) | (_) | |_| | | | (_|  __/
|____/ \___/ \__,_|_|  \___\___|

 _____      _                  _   _             
| ____|_  _| |_ _ __ __ _  ___| |_(_) ___  _ __  
|  _| \ \/ / __| '__/ _` |/ __| __| |/ _ \| '_ \ 
| |___ >  <| |_| | | (_| | (__| |_| | (_) | | | |
|_____/_/\_\\__|_|  \__,_|\___|\__|_|\___/|_| |_|


"""


def sep_background(imagename,mask=None,apply_mask=False,show_map=False,
                   bw=64, bh=64, fw=5, fh=5, use_beam_fraction=False,
                   b_factor=2,f_factor=2):
    """
    Use SEP to estimate the background of an image.

    Parameters
    ----------
    imagename : str
        Path to the image.
    mask : array
        Mask to be applied to the image.
    apply_mask : bool
        If True, calculate the dilated mask from the image.
    show_map : bool
        If True, show the background map.
    bw : int
        Box width for the background estimation.
    bh : int
        Box height for the background estimation.
    fw : int
        Filter width for the background estimation.
    fh : int
        Filter height for the background estimation.
    use_beam_fraction : bool
        If True, use the beam fraction sizes for the sizes of the boxes and
        filters (bw, bh, fw, fh).
    bfactor : int (optional)
        Factor to multiply the box sizes (bw, bh) by.
    factor : int (optional)
        Factor to multiply the filter sizes (fw, fh) by.

    Returns
    -------
    bkg : sep.Background
        Background object.
    """
    '''
    If using astropy.io.fits, you get an error (see bug on sep`s page).
    '''
    _data_2D = fitsio.read(imagename)
    if len(_data_2D.shape) == 4:
        data_2D = _data_2D[0][0]
    else:
        data_2D = _data_2D

    if use_beam_fraction:
        bspx = get_beam_size_px(imagename)
        bspx_x, bspx_y = int(bspx[1]), int(bspx[2])
        bspx_avg = int(bspx[0])
        print(f"Beam Size in px=({bspx_x},{bspx_y})")
        print(f"Average beam Size in px=({bspx_avg})")
        bw, bh = int(bspx_x*b_factor), int(bspx_y*b_factor)
        fw, fh = int(bspx_avg*f_factor), int(bspx_avg*f_factor)

    if (mask is None) and (apply_mask==True):
        _, mask = mask_dilation(imagename, PLOT=False,
                                sigma=3, iterations=2, dilation_size=10)
        bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)

    else:
        bkg = sep.Background(data_2D,bw=bw, bh=bh, fw=fw, fh=fh)
    bkg_rms = bkg.rms()
    bkg_image = bkg.back()
    if show_map == True:
        plt.imshow(bkg_image,origin='lower')
        plt.title(f"bkg/data=({(np.nanmax(bkg_image)/np.nanmax(data_2D-bkg_image)):.3f},{(np.nansum(bkg_image)/np.nansum(data_2D-bkg_image)):.3f})")
        plt.colorbar()
        # plt.clf()
        # plt.close()
    return(bkg)

def distances_from_reference(x_positions, y_positions, reference_coordinate):
    # Reference coordinate
    x_ref, y_ref = reference_coordinate

    # Calculate distances to the reference coordinate
    distances = np.sqrt((x_positions - x_ref)**2 + (y_positions - y_ref)**2)

    # Return the distances sorted
#     sorted_distances = np.sort(distances)

    return distances

def sep_source_ext(imagename, residualname=None, 
                   sigma=6.0, iterations=2, dilation_size=None,
                   deblend_nthresh=100, deblend_cont=0.005, maskthresh=0.0,
                   gain=1.0, filter_kernel=None, mask=None,
                   segmentation_map=False, clean_param=1.0, clean=True,
                   minarea=20, filter_type='matched', 
                   sort_by='distance',
                   bw=64, bh=64, fw=3, fh=3, ell_size_factor=None,
                   apply_mask=False, sigma_mask=6, minarea_factor=1.0,
                   npixels=None,
                   show_bkg_map=True, show_detection=False):
    """
    Simple source extraction algorithm (using SEP https://sep.readthedocs.io/en/v1.1.x/).

    Parameters
    ----------
    imagename : str
        Path to the image.
    sigma : float
        Sigma level for detection.
    iterations : int
        Number of iterations for the mask dilation.
    dilation_size : int
        Size of the dilation kernel.
    deblend_nthresh : int
        Number of thresholds for deblending.
    deblend_cont : float
        Minimum contrast ratio for deblending.
    maskthresh : float
        Threshold for the mask.
    gain : float
        Gain of the image.
    filter_kernel : array
        Filter kernel for the convolution.
    mask : array
        Mask to be applied to the image.
    segmentation_map : bool
        If True, returns the segmentation map.
    clean_param : float
        Cleaning parameter.
    clean : bool
        If True, clean the image.
    minarea : int
        Minimum area for detection.
    filter_type : str
        Type of filter to be used.
    sort_by : str
        Sort the output by flux or area.
    bw : int
        Box width for the background estimation.
    bh : int
        Box height for the background estimation.
    fw : int
        Filter width for the background estimation.
    fh : int
        Filter height for the background estimation.
    ell_size_factor : int
        Size of the ellipse to be plotted.
    apply_mask : bool
        If True, apply the mask to the image.
    sigma_mask : float
        Sigma level for the mask.
    minarea_factor : float

    """

    # filter_kernel_5x5 = np.array([
    #     [1, 1, 1, 1, 1],
    #     [1, 2, 2, 2, 1],
    #     [1, 2, 3, 2, 1],
    #     [1, 2, 2, 2, 1],
    #     [1, 1, 1, 1, 1]])

    data_2D = fitsio.read(imagename)
    if len(data_2D.shape) == 4:
        data_2D = data_2D[0][0]
        
    if residualname is not None:
        residual_2D = fitsio.read(residualname)
        if len(residual_2D.shape) == 4:
            residual_2D = residual_2D[0][0]
    
    if residualname is not None:
        m, s = np.mean(data_2D), mad_std(residual_2D)
    else:
        m, s = np.mean(data_2D), mad_std(data_2D)

    if apply_mask and mask is None:
        _, mask = mask_dilation(data_2D, sigma=sigma_mask, iterations=iterations,
                                rms=s,
                                PLOT=True,show_figure=True,
                                dilation_size=dilation_size)
        # data_2D = data_2D_ * mask

    if npixels is None:
        npixels = int(minarea * minarea_factor)
    
    print('+++++++++++++++++++++++')
    print('SEP Filter sizes:')
    print('    bw,bh=(', int(bw), int(bh),')')
    print('    fw,fh=(', int(fw), int(fh),')')
    print(f'   rms={s}')
    print(f'   npixels={npixels}')
    print(f'   deblend_nthresh={deblend_nthresh}')
    print(f'   deblend_cont={deblend_cont}')
    print(f'   clean={clean}')
    print(f'   clean_param={clean_param}')
    print('+++++++++++++++++++++++')
    bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
    # print(bkg.globalback)
    # print(bkg.globalrms)
    bkg_image = bkg.back()
    bkg_rms = bkg.rms()
    
    data_sub = data_2D  - bkg_image

    if show_bkg_map == True:
        plt.figure()
        # display bkg map.
        plt.imshow(data_sub, interpolation='nearest', cmap='gray', vmin=1*s,
                    vmax=0.8*np.max(data_sub), origin='lower')
        plt.colorbar()
        plt.show()
        plt.figure()
        plt.imshow(bkg_image)
        plt.show()
        # fast_plot2(bkg_rms)

    if mask is not None:
        data_sub = data_sub * mask
    else:
        data_sub = data_sub.copy()
    
    
    if segmentation_map == True:
        objects, seg_maps = sep.extract(data_sub, thresh=sigma * s,
                                        minarea=npixels, filter_type=filter_type,
                                        deblend_nthresh=deblend_nthresh,
                                        deblend_cont=deblend_cont,
                                        filter_kernel=filter_kernel,
                                        maskthresh=maskthresh, gain=gain,
                                        clean=clean, clean_param=clean_param,
                                        segmentation_map=segmentation_map,
                                        err=None, mask=None)
        print('++==>> INFO: Total number of Sources/Structures (deblended) = ', len(objects))
    else:
        objects = sep.extract(data_sub, thresh=sigma * s,
                              minarea=npixels, filter_type=filter_type,
                              deblend_nthresh=deblend_nthresh,
                              deblend_cont=deblend_cont, filter_kernel=filter_kernel,
                              maskthresh=maskthresh, gain=gain,
                              clean=clean, clean_param=clean_param,
                              segmentation_map=segmentation_map,
                              err=None, mask=None)

    # len(objects)
    from matplotlib.patches import Ellipse
    from skimage.draw import ellipse

    # m, s = np.mean(data_sub), np.std(data_sub)
    if show_detection == True:
        fig, ax = plt.subplots(figsize=(4, 4))
        norm = simple_norm(data_sub, stretch='asinh', asinh_a=0.02, vmin=s,
                    vmax=0.25*np.nanmax(data_sub))
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       norm=norm,
                    #    vmin=s, vmax=0.2*np.nanmax(data_sub), 
                       origin='lower')

    masks_regions = []

    if ell_size_factor is None:
        if mask is not None:
            # ell_size_factor = np.sqrt(np.sum(mask) / (np.pi))/cat[0].equivalent_radius.value
            ell_size_factor = 0.1*np.sqrt(np.nansum(mask) / (np.pi))
        else:
            ell_size_factor = 0.5
    
    y, x = np.indices(data_2D.shape[:2])
    for i in range(len(objects)):
        e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
                    width=1 * ell_size_factor * objects['a'][i],
                    height=1 * ell_size_factor * objects['b'][i],
                    angle=objects['theta'][i] * 180. / np.pi)

        xc = objects['x'][i]
        yc = objects['y'][i]
        a = ell_size_factor * objects['a'][i]
        b = ell_size_factor * objects['b'][i]
        theta = objects['theta'][i]
        rx = (x - xc) * np.cos(theta) + (y - yc) * np.sin(theta)
        ry = (y - yc) * np.cos(theta) - (x - xc) * np.sin(theta)

        inside = ((rx / a) ** 2 + (ry / b) ** 2) <= 1
        mask_ell = np.zeros_like(data_2D)
        mask_ell[inside] = True
        if show_detection == True:
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        masks_regions.append(mask_ell)

    #         plt.savefig('components_SEP.pdf',dpi=300, bbox_inches='tight')
    flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
                                         3.0, err=bkg.globalrms, gain=1.0)
    # for i in range(len(objects)):
    #     print("object {:d}: flux = {:f} +/- {:f}".format(i, flux[i], fluxerr[i]))
    # objects['b'] / objects['a'], np.rad2deg(objects['theta'])
    # print(objects)
    # sort regions from largest size to smallest size.
    # print(objects['x'], objects['y'])
    _masks_regions = [mask == 1.0 for mask in masks_regions]
    masks_regions = _masks_regions
    mask_areas = []
    mask_fluxes = []
    for mask_comp in masks_regions:
        area_mask = np.sum(mask_comp)
        sum_mask = np.sum(mask_comp * data_2D)
        mask_areas.append(area_mask)
        mask_fluxes.append(sum_mask)
    mask_areas = np.asarray(mask_areas)
    mask_fluxes = np.asarray(mask_fluxes)
    if sort_by == 'area':
        sorted_indices_desc = np.argsort(mask_areas)[::-1]
        sorted_arr_desc = mask_areas[sorted_indices_desc]
    if sort_by == 'flux':
        sorted_indices_desc = np.argsort(mask_fluxes)[::-1]
        sorted_arr_desc = mask_fluxes[sorted_indices_desc]
    if sort_by == 'distance':
        ref_centre = data_2D.shape[0] / 2, data_2D.shape[1] / 2
        distances = distances_from_reference(objects['x'],
                                             objects['y'],
                                             ref_centre
                                             )
        sorted_indices_desc = np.argsort(distances)
        sorted_arr_desc = distances[sorted_indices_desc]
        
    
    
    objects_sorted = {}
    objects_sorted['xc'] = np.asarray([1] * len(objects))
    objects_sorted['yc'] = np.asarray([1] * len(objects))
    for i in range(len(objects)):
        objects_sorted['xc'][sorted_indices_desc[i]] = objects['x'][sorted_indices_desc[i]]
        objects_sorted['yc'][sorted_indices_desc[i]] = objects['y'][sorted_indices_desc[i]]

    if show_detection == True:
        for i in range(len(objects)):
            xc = objects['x'][sorted_indices_desc[i]]
            yc = objects['y'][sorted_indices_desc[i]]
            label = str('ID' + str(i + 1))
            label_x = xc + 3 * ell_size_factor / 2
            label_y = yc + 20 * ell_size_factor / 2
            line_end_y = label_y - 5 
            
            text = Text(label_x, label_y, label, ha='center', va='center', color='red',fontsize=10)
            ax.add_artist(text)
            ax.plot([xc, label_x], [yc, line_end_y], color='red', alpha=0.4)

        plt.axis('off')
        # plt.show()
        plt.savefig(imagename + '_SE_reg.jpg', dpi=300, bbox_inches='tight')
        plt.show()

    if segmentation_map == True:
        return (masks_regions, sorted_indices_desc, bkg_image, seg_maps, objects_sorted)
    else:
        return (masks_regions, sorted_indices_desc, bkg_image,  objects_sorted)

# def sep_source_ext(imagename, sigma=10.0, iterations=2, dilation_size=None,
#                    deblend_nthresh=100, deblend_cont=0.005, maskthresh=0.0,
#                    gain=1, filter_kernel=None, mask=None,
#                    segmentation_map=False, clean_param=1.0, clean=True,
#                    minarea=20, filter_type='matched', sort_by='flux',
#                    bw=64, bh=64, fw=3, fh=3, ell_size_factor=2,
#                    apply_mask=False,sigma_mask=6,minarea_factor=1.0,
#                    show_bkg_map=False, show_detection=False):
#     """
#     Simple source extraction algorithm (using SEP https://sep.readthedocs.io/en/v1.1.x/).
#
#
#     """
#     import sep
#     import fitsio
#     import matplotlib.pyplot as plt
#     from matplotlib.text import Text
#     from matplotlib import rcParams
#
#     data_2D = fitsio.read(imagename)
#     if len(data_2D.shape) == 4:
#         data_2D = data_2D[0][0]
#     m, s = np.mean(data_2D), mad_std(data_2D)
#     bkg = sep.Background(data_2D)
#
#     if apply_mask:
#         if mask is not None:
#             data_2D = data_2D * mask
#         else:
#             _, mask = mask_dilation(data_2D, sigma=sigma_mask, iterations=iterations,
#                                     dilation_size=dilation_size)
#             data_2D = data_2D * mask
#
#     # else:
#     #     mask = None
#     bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
#     # print(bkg.globalback)
#     # print(bkg.globalrms)
#     bkg_image = bkg.back()
#     bkg_rms = bkg.rms()
#
#     if show_bkg_map == True:
#         plt.figure()
#         #display bkg map.
#         plt.imshow(data_2D, interpolation='nearest', cmap='gray', vmin=m - s,
#                    vmax=m + s, origin='lower')
#         plt.colorbar()
#         plt.close()
#         plt.figure()
#         plt.imshow(bkg_image)
#         plt.close()
#     # fast_plot2(bkg_rms)
#     data_sub = data_2D - bkg
#     if segmentation_map == True:
#         npixels = int(minarea * minarea_factor)
#         objects, seg_maps = sep.extract(data_sub, thresh=sigma,
#                                         minarea=npixels, filter_type=filter_type,
#                                         deblend_nthresh=deblend_nthresh,
#                                         deblend_cont=deblend_cont, filter_kernel=filter_kernel,
#                                         maskthresh=maskthresh, gain=gain,
#                                         clean=clean, clean_param=clean_param,
#                                         segmentation_map=segmentation_map,
#                                         err=bkg.globalrms, mask=mask)
#     else:
#         npixels = int(minarea * minarea_factor)
#         objects = sep.extract(data_sub, thresh=sigma,
#                               minarea=npixels, filter_type=filter_type,
#                               deblend_nthresh=deblend_nthresh,
#                               deblend_cont=deblend_cont, filter_kernel=filter_kernel,
#                               maskthresh=maskthresh, gain=gain,
#                               clean=clean, clean_param=clean_param,
#                               segmentation_map=segmentation_map,
#                               err=bkg.globalrms, mask=mask)
#
#     # len(objects)
#     from matplotlib.patches import Ellipse
#     from skimage.draw import ellipse
#
#     m, s = np.mean(data_sub), np.std(data_sub)
#     if show_detection == True:
#         fig, ax = plt.subplots()
#         im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
#                        vmin=m - s, vmax=m + s, origin='lower')
#
#     masks_regions = []
#
#     y, x = np.indices(data_2D.shape[:2])
#     for i in range(len(objects)):
#         e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
#                     width=2 * ell_size_factor * objects['a'][i],
#                     height=2 * ell_size_factor * objects['b'][i],
#                     angle=objects['theta'][i] * 180. / np.pi)
#
#         xc = objects['x'][i]
#         yc = objects['y'][i]
#         a = ell_size_factor * objects['a'][i]
#         b = ell_size_factor * objects['b'][i]
#         theta = objects['theta'][i]
#         rx = (x - xc) * np.cos(theta) + (y - yc) * np.sin(theta)
#         ry = (y - yc) * np.cos(theta) - (x - xc) * np.sin(theta)
#
#         inside = ((rx / a) ** 2 + (ry / b) ** 2) <= 1
#         mask_ell = np.zeros_like(data_2D)
#         mask_ell[inside] = True
#         if show_detection == True:
#             e.set_facecolor('none')
#             e.set_edgecolor('red')
#             ax.add_artist(e)
#         masks_regions.append(mask_ell)
#
#     #         plt.savefig('components_SEP.pdf',dpi=300, bbox_inches='tight')
#     flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
#                                          3.0, err=bkg.globalrms, gain=1.0)
#     for i in range(len(objects)):
#         print("object {:d}: flux = {:f} +/- {:f}".format(i, flux[i], fluxerr[i]))
#     objects['b'] / objects['a'], np.rad2deg(objects['theta'])
#
#     # sort regions from largest size to smallest size.
#     mask_areas = []
#     mask_fluxes = []
#     for mask_comp in masks_regions:
#         area_mask = np.sum(mask_comp)
#         sum_mask = np.sum(mask_comp * data_2D)
#         mask_areas.append(area_mask)
#         mask_fluxes.append(sum_mask)
#     mask_areas = np.asarray(mask_areas)
#     mask_fluxes = np.asarray(mask_fluxes)
#     if sort_by == 'area':
#         sorted_indices_desc = np.argsort(mask_areas)[::-1]
#         sorted_arr_desc = mask_areas[sorted_indices_desc]
#     if sort_by == 'flux':
#         sorted_indices_desc = np.argsort(mask_fluxes)[::-1]
#         sorted_arr_desc = mask_fluxes[sorted_indices_desc]
#
#     objects_sorted = {}
#     objects_sorted['xc'] = np.asarray([1] * len(objects))
#     objects_sorted['yc'] = np.asarray([1] * len(objects))
#     for i in range(len(objects)):
#         objects_sorted['xc'][i] = objects['x'][sorted_indices_desc[i]]
#         objects_sorted['yc'][i] = objects['y'][sorted_indices_desc[i]]
#
#     if show_detection == True:
#         for i in range(len(objects)):
#             xc = objects['x'][sorted_indices_desc[i]]
#             yc = objects['y'][sorted_indices_desc[i]]
#             label = str('ID' + str(i + 1))
#             text = Text(xc + 10 * ell_size_factor, yc + 3 * ell_size_factor, label, ha='center', va='center', color='red')
#             ax.add_artist(text)
#
#         plt.axis('off')
#         # plt.show()
#         plt.savefig(imagename + '_SE_reg.jpg', dpi=300, bbox_inches='tight')
#         plt.show()
#
#     if segmentation_map == True:
#         return (masks_regions, sorted_indices_desc, seg_maps, objects_sorted)
#     else:
#         return (masks_regions, sorted_indices_desc, objects_sorted)


def phot_source_ext(imagename, residual=None, sigma=1.0, iterations=2, dilation_size=None,
                    deblend_nthresh=5, deblend_cont=1e-6, maskthresh=0.0,
                    gain=1, filter_kernel=None, mask=None,
                    segmentation_map=False, clean_param=1.0, clean=True,
                    minarea=100, minarea_factor=1, npixels = None,
                    psf_data = None,
                    filter_type='matched', 
                    sort_by='distance', 
                    first_ID_only=False,
                    threshold_mode='sigma',
                    bw=64, bh=64, fw=3, fh=3, ell_size_factor=None,
                    apply_mask=False, sigma_mask=6,
                    show_bkg_map=False, show_detection=False,
                    SE_ref=None):
    """
    Simple source extraction algorithm (using SEP https://sep.readthedocs.io/en/v1.1.x/).


    """

    data_2D = load_fits_data(imagename)
    if len(data_2D.shape) == 4:
        data_2D = data_2D[0][0]
    # m, s = np.mean(data_2D), np.std(data_2D)
    if residual is not None:
        m, s = np.mean(data_2D), mad_std(residual)
    else:
        m, s = np.mean(data_2D), mad_std(data_2D)
    # bkg = 0.0
    if apply_mask and mask is None:
        _, mask = mask_dilation(data_2D, sigma=sigma_mask, iterations=iterations,
                                rms=s,
                                PLOT=True,show_figure=True,
                                dilation_size=dilation_size)


    # bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
    bkg_image, bkg_rms, _ = multiscale_segmentation_background(data_2D,
                                                   combine_method='finest_valid',
                                                   sigma_sigma_clip=2.5, maxiters=15,
                                                   beam_size_px=3, n_scales=8)
    
    if psf_data is not None:
        psf_fwhm = 1.0*psf_params(psf_data)
    else:
        psf_fwhm = 5.0
    bkg_image, bkg_rms, _ = multiscale_segmentation_background(data_2D,
                                                               combine_method='finest_valid',
                                                               #    combine_method='minimum',
                                                               sigma_sigma_clip=2.0, maxiters=15,
                                                               beam_size_px=int(psf_fwhm*1+1), 
                                                               #    n_scales=8,
                                                            #    do_plot=True,
                                                            #    profile_type='azimuthal'
                                                               )
    # print(bkg.globalback)
    # print(bkg.globalrms)

    # bkg_image = bkg.back()
    # bkg_rms = bkg.rms()
    
    data_sub = data_2D - bkg_image

    if show_bkg_map == True:
        # plt.figure()
        # # display bkg map.
        # plt.imshow(data_2D, interpolation='nearest', cmap='gray', vmin=3*s,
        #            vmax=0.2*np.max(data_2D), origin='lower')
        # plt.colorbar()
        # plt.close()
        # plt.show()
        # plt.figure()
        # plt.imshow(bkg_image)
        # plt.colorbar()
        # plt.show()
        # plt.close()
        
        vmax = np.percentile(data_sub,99.9)
        vmin = np.percentile(data_sub,0.1)
        norm = simple_norm(data_sub, 
                                vmax = vmax,
                                vmin = vmin,
                                stretch='asinh',
                                asinh_a=0.05)
        fig = plt.figure(figsize=(18, 6))
        ax0 = fig.add_subplot(1, 3, 2)
        ax0 = plt.imshow(data_sub, origin='lower', norm=norm, cmap='Greys_r')
        plt.title(f'data - bkg')
        
        ax1 = fig.add_subplot(1, 3, 1)
        ax1 = plt.imshow(data_2D, origin='lower', norm=norm, cmap='Greys_r')
        plt.title(f'data')
        
        ax2 = fig.add_subplot(1, 3, 3)
        ax2 = plt.imshow(bkg_image, origin='lower', 
                        #  norm=norm, 
                         cmap='Greys_r')
        plt.title(f'bkg')
        
        # fig = plt.figure(figsize=(18, 6))
        # ax0 = fig.add_subplot(1, 3, 2)
        # ax0 = eimshow(data_sub,ax=ax0,fig=fig,
        #                     plot_colorbar=True,
        #                     rms = 3*np.nanmedian(bkg_rms),
        #                     # cbar_orientation = 'horizontal',
        #                     show_axis='off',
        #                     vmin_factor=0.5,
        #                     cbar_n_points = 3,
        #                     add_contours=True)
        # ax0.set_title(f'data - bkg')
        
        # ax1 = fig.add_subplot(1, 3, 1)
        # ax1 = eimshow(data_2D,ax=ax1,fig=fig,
        #                     plot_colorbar=True,
        #                     rms = 3*np.nanmedian(bkg_rms),
        #                     # cbar_orientation = 'horizontal',
        #                     show_axis='off',
        #                     vmin_factor=0.5,
        #                     cbar_n_points=3,
        #                     plot_title = f'data',
        #                     add_contours=True)
        # ax2 = fig.add_subplot(1, 3, 3)
        # ax2 = eimshow(bkg_image,ax=ax2,fig=fig,
        #                     plot_colorbar=True,
        #                     rms = 3*np.nanmedian(bkg_rms),
        #                     # cbar_orientation = 'horizontal',
        #                     show_axis='off',
        #                     # vmin_factor=0.5,
        #                     # vmax_factor=1.0,
        #                     cbar_n_points=3,
        #                     plot_title = f'bkg',
        #                     add_contours=False)
        plt.show()
        plt.clf()
        plt.close()

    
    if mask is not None:
        data_sub = data_sub * mask
    else:
        data_sub = data_sub
        
    if (threshold_mode == 'residual') and residual is not None:
        threshold = sigma * residual
    elif threshold_mode == 'sigma':
        threshold = sigma * s
            
    # else:
    #     mask = None
    # print(data_sub)
    if npixels is None:
        npixels = int(minarea * minarea_factor)
    # print(' INFO: Uinsg min number of pixels of :', npixels)
    cat, segm, seg_maps = make_catalog(image=data_sub,
                                       threshold=threshold,
                                       deblend=True, 
                                       contrast=deblend_cont,
                                       nlevels=deblend_nthresh,
                                       npixels=npixels,
                                       figsize=(8, 8),
                                       plot=show_detection, vmin=1.0 * s)
    # cat, segm, seg_maps = make_catalog(image=data_sub,
    #                                    threshold=sigma * s,
    #                                    deblend=True, 
    #                                    contrast=deblend_cont,
    #                                    nlevels=deblend_nthresh,
    #                                    npixels=npixels,
    #                                    figsize=(20, 20),
    #                                    plot=show_detection, vmin=1.0 * s)
    if sort_by == 'flux':
        indices = list(order_cat(cat, key='segment_flux', reverse=True))
    if sort_by == 'area':
        indices = list(order_cat(cat, key='area', reverse=True))
    if sort_by == 'distance':
        ref_centre = data_2D.shape[0] / 2, data_2D.shape[1] / 2
        if SE_ref is not None:
            coords, distances, indices = \
                sorted_detected_coordinates(SE_ref.objects['xc'],
                                            SE_ref.objects['yc'],
                                            cat.xcentroid, 
                                            cat.ycentroid, 
                                            ref_centre)
            print(indices)
        else:
            
            distances = distances_from_reference(cat.xcentroid, 
                                            cat.ycentroid,
                                            ref_centre)
            indices = np.argsort(distances)

    if first_ID_only and sort_by == 'distance':
        """
        Force photometry on the central detected source only.
        """
        indices = [indices[0]]
        cat = cat[indices]

        
    masks_deblended = []
    for k in range(len(indices)):
        # print(k)
        masks_deblended.append(seg_maps == seg_maps.labels[indices[k]])


    # len(objects)
    from matplotlib.patches import Ellipse
    from skimage.draw import ellipse

    # m, s = np.mean(data_sub), np.std(data_sub)
    if show_detection == True:
        fig, ax = plt.subplots(figsize=(6, 6))
        norm = simple_norm(data_sub, stretch='sqrt', asinh_a=0.02, vmin=s,
                    vmax=0.2*np.nanmax(data_sub))
        im = ax.imshow(data_sub, interpolation='nearest', 
                    #    cmap='RdBu',
                       cmap='gray_r',
                       norm=norm, origin='lower')
        # im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
        #                vmin=s, vmax=0.2*np.nanmax(data_sub), origin='lower')

    masks_regions = []

    if ell_size_factor is None:
        if mask is not None:
            # ell_size_factor = np.sqrt(np.sum(mask) / (np.pi))/cat[0].equivalent_radius.value
            ell_size_factor = 0.05*np.sqrt(np.nansum(mask) / (np.pi))
        else:
            ell_size_factor = 0.5
    
    y, x = np.indices(data_2D.shape[:2])
    
    # objects = cat
    for i in range(len(cat)):
        source = cat[i]
        seg_mask = (seg_maps.data == i + 1)
        e = Ellipse(xy=(source.centroid[0], source.centroid[1]),
                    width=1 * ell_size_factor * source.equivalent_radius.value,
                    height=1 * ell_size_factor * (
                                1 - source.ellipticity.value) * source.equivalent_radius.value,
                    # angle=source.orientation.value * 180 / np.pi
                    angle=source.orientation.value
                    )

        xc = source.centroid[0]
        yc = source.centroid[1]
        a = ell_size_factor * source.equivalent_radius.value
        b = ell_size_factor * (
                    1 - source.ellipticity.value) * source.equivalent_radius.value
        theta = source.orientation.value
        rx = (x - xc) * np.cos(theta) + (y - yc) * np.sin(theta)
        ry = (y - yc) * np.cos(theta) - (x - xc) * np.sin(theta)

        inside = ((rx / a) ** 2 + (ry / b) ** 2) <= 1
        mask_ell = np.zeros_like(data_2D)
        mask_ell[inside] = True
        if show_detection == True:
            e.set_facecolor('none')
            # e.set_edgecolor('#009E73')
            e.set_edgecolor('orange')
            e.set_linewidth(2)
            ax.add_artist(e)
        masks_regions.append(seg_mask)

    #         plt.savefig('components_SEP.pdf',dpi=300, bbox_inches='tight')
    # flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
    #                                      3.0, err=bkg.globalrms, gain=1.0)
    # for i in range(len(objects)):
    #     print("object {:d}: flux = {:f} +/- {:f}".format(i, flux[i], fluxerr[i]))
    # objects['b'] / objects['a'], np.rad2deg(objects['theta'])

    # sort regions from largest size to smallest size.
    mask_areas = []
    mask_fluxes = []
    for mask_comp in masks_regions:
        area_mask = np.sum(mask_comp)
        sum_mask = np.sum(mask_comp * data_2D)
        mask_areas.append(area_mask)
        mask_fluxes.append(sum_mask)
    mask_areas = np.asarray(mask_areas)
    mask_fluxes = np.asarray(mask_fluxes)
    if sort_by == 'area':
        sorted_indices_desc = np.argsort(mask_areas)[::-1]
        sorted_arr_desc = mask_areas[sorted_indices_desc]
    if sort_by == 'flux':
        sorted_indices_desc = np.argsort(mask_fluxes)[::-1]
        sorted_arr_desc = mask_fluxes[sorted_indices_desc]
    if sort_by == 'distance':
        ref_centre = data_2D.shape[0] / 2, data_2D.shape[1] / 2
        if SE_ref is not None:
            coords, distances, sorted_indices_desc = \
                sorted_detected_coordinates(SE_ref.objects['xc'],
                                            SE_ref.objects['yc'],
                                            cat.xcentroid, 
                                            cat.ycentroid, 
                                            ref_centre)
            # sorted_indices_desc = np.argsort(distances)
            # sorted_arr_desc = distances[sorted_indices_desc]
                
        else:
            distances = distances_from_reference(cat.xcentroid, 
                                                cat.ycentroid,
                                                ref_centre)
            sorted_indices_desc = np.argsort(distances)
            # sorted_arr_desc = distances[sorted_indices_desc]

    objects_sorted = {}
    objects_sorted['xc'] = np.asarray([1] * len(cat))
    objects_sorted['yc'] = np.asarray([1] * len(cat))
    for i in range(len(cat)):
        source = cat[sorted_indices_desc[i]]
        objects_sorted['xc'][i] = source.centroid[0]
        objects_sorted['yc'][i] = source.centroid[1]

    if show_detection == True:
        for i in range(len(cat)):
            source = cat[sorted_indices_desc[i]]
            xc = source.centroid[0]
            yc = source.centroid[1]
            label = str('ID' + str(i + 1))
            
            label_x = xc + 3 * ell_size_factor
            label_y = yc + 20 * ell_size_factor
            
            line_end_y = label_y - 5 
            
            # text = Text(label_x, label_y, label, ha='center', va='center', color='purple')
            text = ax.text(label_x, label_y, label, ha='center', va='center', color='purple', fontweight='bold')
            text.set_fontweight('bold')
            ax.add_artist(text)
            
            ax.plot([xc, label_x], [yc, line_end_y], color='purple', alpha=0.4)
            

        plt.axis('off')
        plt.savefig(imagename + '_SE_reg.jpg', dpi=300, bbox_inches='tight')
        plt.show()



    if segmentation_map == True:
        return (masks_deblended, sorted_indices_desc, bkg_image,
                seg_maps, objects_sorted, cat)
    else:
        return (masks_deblended, sorted_indices_desc, bkg_image,
                None, objects_sorted, cat)


def astphot_source_ext(imagename, residual=None, sigma=1.0, iterations=2, dilation_size=None,
                    deblend_nthresh=5, deblend_cont=1e-6, maskthresh=0.0,
                    gain=1, filter_kernel=None, mask=None,
                    segmentation_map=False, clean_param=1.0, clean=True,
                    minarea=100, minarea_factor=1, npixels = None,
                    psf_data = None,
                    filter_type='matched', 
                    sort_by='distance', 
                    first_ID_only=False,
                    threshold_mode='sigma',
                    bw=64, bh=64, fw=3, fh=3, ell_size_factor=None,
                    apply_mask=False, sigma_mask=6,
                    show_bkg_map=False, show_detection=False,
                    SE_ref=None):
    """
    Simple source extraction algorithm using astropy/photutils.
    
    This function performs source detection, deblending, and cataloging using
    photutils detection and source catalog functionality. It provides direct
    access to source properties through the SourceCatalog object.
    """
    from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
    from matplotlib.patches import Ellipse
    from matplotlib.text import Text
    
    data_2D = load_fits_data(imagename)
    if len(data_2D.shape) == 4:
        data_2D = data_2D[0][0]
    
    if residual is not None:
        m, s = np.mean(data_2D), mad_std(residual)
    else:
        m, s = np.mean(data_2D), mad_std(data_2D)
    
    if apply_mask and mask is None:
        _, mask = mask_dilation(data_2D, sigma=sigma_mask, iterations=iterations,
                                rms=s,
                                PLOT=True, show_figure=True,
                                dilation_size=dilation_size)

    # Background estimation using multiscale segmentation
    if psf_data is not None:
        psf_fwhm = 1.0 * psf_params(psf_data)
    else:
        psf_fwhm = 5.0
    
    bkg_image, bkg_rms, _ = multiscale_segmentation_background(
        data_2D,
        combine_method='finest_valid',
        sigma_sigma_clip=2.0, 
        maxiters=15,
        beam_size_px=int(psf_fwhm * 1 + 1)
    )
    
    data_sub = data_2D - bkg_image

    if show_bkg_map == True:
        vmax = np.percentile(data_sub, 99.9)
        vmin = np.percentile(data_sub, 0.1)
        norm = simple_norm(data_sub, 
                          vmax=vmax,
                          vmin=vmin,
                          stretch='asinh',
                          asinh_a=0.05)
        fig = plt.figure(figsize=(18, 6))
        ax0 = fig.add_subplot(1, 3, 2)
        ax0 = plt.imshow(data_sub, origin='lower', norm=norm, cmap='Greys_r')
        plt.title(f'data - bkg')
        
        ax1 = fig.add_subplot(1, 3, 1)
        ax1 = plt.imshow(data_2D, origin='lower', norm=norm, cmap='Greys_r')
        plt.title(f'data')
        
        ax2 = fig.add_subplot(1, 3, 3)
        ax2 = plt.imshow(bkg_image, origin='lower', cmap='Greys_r')
        plt.title(f'bkg')
        
        plt.show()
        plt.clf()
        plt.close()

    if mask is not None:
        data_sub = data_sub * mask
    else:
        data_sub = data_sub
    
    # Determine threshold based on mode
    if (threshold_mode == 'residual') and residual is not None:
        threshold = sigma * residual
    elif threshold_mode == 'sigma':
        threshold = sigma * s
    
    if npixels is None:
        npixels = int(minarea * minarea_factor)
    
    # Source detection using photutils
    segmap = detect_sources(data_sub, threshold, npixels=npixels)
    
    # Handle case where no sources are detected
    if segmap is None:
        print("WARNING: No sources detected, returning empty results")
        empty_dict = {'xc': np.array([]), 'yc': np.array([])}
        if segmentation_map:
            return ([], [], bkg_image, None, empty_dict)
        else:
            return ([], [], bkg_image, empty_dict)
    
    # Deblend sources using photutils
    if deblend_nthresh > 0:
        seg_maps = deblend_sources(
            data_sub, 
            segmap, 
            npixels=npixels,
            nlevels=deblend_nthresh,
            contrast=deblend_cont,
            mode='exponential',
            progress_bar=False
        )
    else:
        seg_maps = segmap
    
    # Create source catalog using photutils SourceCatalog
    cat = SourceCatalog(data_2D, seg_maps, convolved_data=data_sub)
    
    if len(cat) == 0:
        print("WARNING: No sources in catalog, returning empty results")
        empty_dict = {'xc': np.array([]), 'yc': np.array([])}
        if segmentation_map:
            return ([], [], bkg_image, seg_maps, empty_dict)
        else:
            return ([], [], bkg_image, empty_dict)
    
    # Sort sources based on specified criterion
    if sort_by == 'flux':
        indices = np.argsort([src.segment_flux for src in cat])[::-1]
    elif sort_by == 'area':
        indices = np.argsort([src.area.value for src in cat])[::-1]
    elif sort_by == 'distance':
        ref_centre = data_2D.shape[0] / 2, data_2D.shape[1] / 2
        if SE_ref is not None:
            coords, distances, indices = sorted_detected_coordinates(
                SE_ref.objects['xc'],
                SE_ref.objects['yc'],
                cat.xcentroid, 
                cat.ycentroid, 
                ref_centre
            )
            print(indices)
        else:
            distances = distances_from_reference(
                cat.xcentroid, 
                cat.ycentroid,
                ref_centre
            )
            indices = np.argsort(distances)
    else:
        indices = np.arange(len(cat))

    # Handle first_ID_only option
    if first_ID_only and sort_by == 'distance':
        indices = [indices[0]]
        cat = cat[indices]

    # Create deblended masks for each source
    masks_deblended = []
    for k in range(len(indices)):
        masks_deblended.append(seg_maps.data == seg_maps.labels[indices[k]])

    # Visualization if requested
    if show_detection == True:
        fig, ax = plt.subplots(figsize=(7, 7))
        norm = simple_norm(data_sub, stretch='sqrt', asinh_a=0.02, vmin=s,
                          vmax=0.2 * np.nanmax(data_sub))
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                      norm=norm, origin='lower')

    # Calculate ellipse size factor
    if ell_size_factor is None:
        if mask is not None:
            ell_size_factor = 0.05 * np.sqrt(np.nansum(mask) / (np.pi))
        else:
            ell_size_factor = 0.5
    
    y, x = np.indices(data_2D.shape[:2])
    
    # Create region masks for each source
    masks_regions = []
    for i in range(len(cat)):
        source = cat[i]
        seg_mask = (seg_maps.data == seg_maps.labels[i])
        
        if show_detection == True:
            e = Ellipse(
                xy=(source.xcentroid, source.ycentroid),
                width=1 * ell_size_factor * source.equivalent_radius.value,
                height=1 * ell_size_factor * (1 - source.ellipticity.value) * source.equivalent_radius.value,
                angle=source.orientation.value
            )
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        
        masks_regions.append(seg_mask)

    # Calculate mask areas and fluxes for sorting validation
    mask_areas = []
    mask_fluxes = []
    for mask_comp in masks_regions:
        area_mask = np.sum(mask_comp)
        sum_mask = np.sum(mask_comp * data_2D)
        mask_areas.append(area_mask)
        mask_fluxes.append(sum_mask)
    mask_areas = np.asarray(mask_areas)
    mask_fluxes = np.asarray(mask_fluxes)
    
    # Re-sort based on specified criterion (validation)
    if sort_by == 'area':
        sorted_indices_desc = np.argsort(mask_areas)[::-1]
    elif sort_by == 'flux':
        sorted_indices_desc = np.argsort(mask_fluxes)[::-1]
    elif sort_by == 'distance':
        ref_centre = data_2D.shape[0] / 2, data_2D.shape[1] / 2
        if SE_ref is not None:
            coords, distances, sorted_indices_desc = sorted_detected_coordinates(
                SE_ref.objects['xc'],
                SE_ref.objects['yc'],
                cat.xcentroid, 
                cat.ycentroid, 
                ref_centre
            )
        else:
            distances = distances_from_reference(
                cat.xcentroid, 
                cat.ycentroid,
                ref_centre
            )
            sorted_indices_desc = np.argsort(distances)
    else:
        sorted_indices_desc = indices

    # Create objects_sorted dictionary for compatibility
    objects_sorted = {}
    objects_sorted['xc'] = np.asarray([cat[sorted_indices_desc[i]].xcentroid for i in range(len(cat))])
    objects_sorted['yc'] = np.asarray([cat[sorted_indices_desc[i]].ycentroid for i in range(len(cat))])

    # Add source labels to visualization
    if show_detection == True:
        for i in range(len(cat)):
            source = cat[sorted_indices_desc[i]]
            xc = source.xcentroid
            yc = source.ycentroid
            label = str('ID' + str(i + 1))
            
            label_x = xc + 3 * ell_size_factor
            label_y = yc + 20 * ell_size_factor
            line_end_y = label_y - 5 
            
            text = Text(label_x, label_y, label, ha='center', va='center', color='red')
            ax.add_artist(text)
            ax.plot([xc, label_x], [yc, line_end_y], color='red', alpha=0.4)

        plt.axis('off')
        plt.savefig(imagename + '_SE_reg.jpg', dpi=300, bbox_inches='tight')
        plt.show()

    # Return results with preserved interface
    if segmentation_map == True:
        return (masks_deblended, sorted_indices_desc, bkg_image, seg_maps, objects_sorted, cat)
    else:
        return (masks_deblended, sorted_indices_desc, bkg_image, None, objects_sorted, cat)

"""
Robust Source Extraction Implementation
----------------------------------------
A sophisticated, adaptive source extraction algorithm combining modern computer 
vision techniques with astronomical source detection best practices.

Key Features:
- Multi-scale adaptive background estimation
- Hierarchical multi-threshold detection
- Marker-based watershed deblending
- Extended structure recovery
- Adaptive parameter selection based on data characteristics

Author: Developed for morphen library
"""


def robust_source_ext(imagename, residualname=None, 
                      sigma=6.0, iterations=2, dilation_size=None,
                      deblend_nthresh=100, deblend_cont=0.005, maskthresh=0.0,
                      gain=1.0, filter_kernel=None, mask=None,
                      segmentation_map=False, clean_param=1.0, clean=True,
                      minarea=20, filter_type='matched', 
                      sort_by='distance',
                      bw=64, bh=64, fw=3, fh=3, ell_size_factor=None,
                      apply_mask=False, sigma_mask=6, minarea_factor=1.0,
                      npixels=None,
                      # New robust-specific parameters
                      multiscale_levels=3,  # Number of scales for analysis
                      extended_threshold_factor=0.5,  # Lower threshold for extended emission
                      watershed_connectivity=2,  # Watershed connectivity (1, 2, or 3)
                      merge_threshold=0.3,  # Threshold for merging nearby detections
                      adaptive_background=True,  # Use adaptive background estimation
                      preserve_extended=True,  # Preserve extended low-surface-brightness features
                      min_separation=None,  # Minimum separation between sources (in pixels)
                      show_bkg_map=False, show_detection=False):
    """
    Robust source extraction algorithm combining multi-scale detection with 
    intelligent deblending.
    
    This algorithm uses:
    1. Multi-scale adaptive background estimation
    2. Hierarchical threshold detection (high threshold for cores, low for extended)
    3. Marker-based watershed for controlled deblending
    4. Extended structure recovery at multiple scales
    5. Perceptual grouping to avoid over-deblending
    
    Parameters
    ----------
    imagename : str
        Path to the FITS image file.
    residualname : str, optional
        Path to the residual image for RMS estimation.
    sigma : float, default=6.0
        Detection threshold in units of RMS for compact sources.
    multiscale_levels : int, default=3
        Number of scales for multi-scale decomposition.
    extended_threshold_factor : float, default=0.5
        Factor multiplied by sigma for detecting extended emission.
        Lower values detect fainter extended structures.
    watershed_connectivity : int, default=2
        Connectivity for watershed algorithm (1, 2, or 3).
    merge_threshold : float, default=0.3
        Distance threshold (as fraction of source size) for merging nearby detections.
    adaptive_background : bool, default=True
        Use locally adaptive background estimation.
    preserve_extended : bool, default=True
        Preserve extended low-surface-brightness features.
    min_separation : float, optional
        Minimum separation between source peaks (in pixels).
        If None, automatically determined from beam size.
    Other parameters match standard source extraction functions.
    
    Returns
    -------
    masks_regions : list
        List of boolean masks for each detected source.
    sorted_indices_desc : array
        Indices for sorting sources by specified criterion.
    bkg_image : array
        Background image.
    seg_maps : array
        Segmentation map with labeled regions.
    objects_sorted : dict
        Dictionary with sorted source coordinates ('xc', 'yc').
    """
    
    # ========================
    # 1. DATA LOADING AND PREPARATION
    # ========================
    print("\n" + "="*60)
    print("ROBUST SOURCE EXTRACTION")
    print("="*60)
    
    # Load data
    _data_2D = fitsio.read(imagename)
    if len(_data_2D.shape) == 4:
        data_2D = _data_2D[0][0]
    else:
        data_2D = _data_2D
    
    # Get residual for noise estimation if available
    if residualname is not None:
        residual_data = fitsio.read(residualname)
        if len(residual_data.shape) == 4:
            residual_data = residual_data[0][0]
        s = mad_std(residual_data)
    else:
        s = mad_std(data_2D)
    
    print(f"\n[1/7] Data loaded: shape={data_2D.shape}, RMS={s:.6e}")
    
    # Set minimum area if npixels not provided
    if npixels is None:
        try:
            npixels = int(minarea * minarea_factor)
        except:
            npixels = int(data_2D.shape[0] / 30)
    
    # Determine minimum separation automatically from beam if not provided
    if min_separation is None:
        try:
            beam_px = get_beam_size_px(imagename)
            min_separation = beam_px[0] * 0.5  # Half beam as minimum separation
        except:
            min_separation = 3.0  # Default fallback
    
    print(f"      Minimum detection area: {npixels} pixels")
    print(f"      Minimum source separation: {min_separation:.1f} pixels")
    
    # ========================
    # 2. ADAPTIVE MULTI-SCALE BACKGROUND ESTIMATION
    # ========================
    print(f"\n[2/7] Computing adaptive multi-scale background...")
    
    if adaptive_background:
        # Use multi-scale background estimation
        bkg_scales = []
        for scale_factor in [1, 2, 4]:
            scale_bw = int(bw * scale_factor)
            scale_bh = int(bh * scale_factor)
            scale_fw = max(3, int(fw * scale_factor))
            scale_fh = max(3, int(fh * scale_factor))
            
            bkg_temp = sep.Background(data_2D, mask=mask, 
                                     bw=scale_bw, bh=scale_bh, 
                                     fw=scale_fw, fh=scale_fh)
            bkg_scales.append(bkg_temp.back())
        
        # Combine scales - use minimum to preserve sources
        bkg_image = np.min(bkg_scales, axis=0)
        
        # Final background object for RMS
        bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
        
        print(f"      Multi-scale background computed at {len(bkg_scales)} scales")
    else:
        # Standard background estimation
        bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
        bkg_image = bkg.back()
    
    # Background-subtracted data
    data_sub = data_2D - bkg_image
    
    if mask is not None:
        data_sub = data_sub * mask
    
    bkg_rms = bkg.rms()
    
    if show_bkg_map:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        axes[0].imshow(data_2D, origin='lower', cmap='gray')
        axes[0].set_title('Original Data')
        axes[0].axis('off')
        
        axes[1].imshow(bkg_image, origin='lower', cmap='gray')
        axes[1].set_title('Background')
        axes[1].axis('off')
        
        norm = simple_norm(data_sub, stretch='asinh', asinh_a=0.02, 
                          vmin=s, vmax=0.25*np.nanmax(data_sub))
        axes[2].imshow(data_sub, origin='lower', cmap='gray', norm=norm)
        axes[2].set_title('Background-Subtracted')
        axes[2].axis('off')
        
        plt.tight_layout()
        plt.savefig(imagename + '_robust_background.jpg', dpi=150, bbox_inches='tight')
        plt.show()
    
    # ========================
    # 3. MULTI-SCALE FEATURE DETECTION
    # ========================
    print(f"\n[3/7] Performing multi-scale source detection...")
    
    # High-threshold detection for compact sources (cores)
    thresh_high = sigma * s
    objects_core, seg_core = sep.extract(
        data_sub, thresh=thresh_high,
        minarea=int(npixels * 0.7),  # Slightly smaller for cores
        filter_type=filter_type,
        deblend_nthresh=deblend_nthresh * 2,  # More aggressive for cores
        deblend_cont=deblend_cont,
        clean=clean, clean_param=clean_param,
        segmentation_map=True,
        err=None, mask=None
    )
    
    print(f"      High-threshold detection (sigma={sigma:.1f}): {len(objects_core)} cores found")
    
    # Low-threshold detection for extended emission
    if preserve_extended:
        thresh_low = sigma * extended_threshold_factor * s
        objects_ext, seg_ext = sep.extract(
            data_sub, thresh=thresh_low,
            minarea=npixels * 2,  # Larger minimum area for extended
            filter_type='matched',
            deblend_nthresh=int(deblend_nthresh * 0.5),  # Less aggressive
            deblend_cont=deblend_cont * 10,  # Higher contrast needed
            clean=False,  # Don't clean extended structures
            clean_param=clean_param,
            segmentation_map=True,
            err=None, mask=None
        )
        
        print(f"      Low-threshold detection (sigma={sigma*extended_threshold_factor:.1f}): {len(objects_ext)} extended structures found")
    else:
        objects_ext = objects_core
        seg_ext = seg_core
    
    # ========================
    # 4. INTELLIGENT MARKER-BASED WATERSHED DEBLENDING
    # ========================
    print(f"\n[4/7] Applying marker-based watershed deblending...")
    
    # Create markers from core detections
    markers = np.zeros_like(data_sub, dtype=np.int32)
    
    # Use core detections as seeds
    for i in range(len(objects_core)):
        x_int = int(np.round(objects_core['x'][i]))
        y_int = int(np.round(objects_core['y'][i]))
        
        # Check bounds
        if 0 <= y_int < markers.shape[0] and 0 <= x_int < markers.shape[1]:
            markers[y_int, x_int] = i + 1
    
    # Apply distance transform to spread markers slightly
    from scipy.ndimage import distance_transform_edt, binary_dilation
    
    # Dilate markers slightly based on min_separation
    structure = disk(int(min_separation))
    markers_dilated = np.zeros_like(markers)
    for i in range(1, len(objects_core) + 1):
        marker_i = (markers == i)
        markers_dilated += (binary_dilation(marker_i, structure) * i).astype(np.int32)
    
    # Create watershed image - use inverted intensity so watershed flows to peaks
    watershed_image = -data_sub
    watershed_image[watershed_image < 0] = 0  # Remove negative values
    
    # Apply watershed
    from skimage.segmentation import watershed
    from skimage.filters import gaussian
    
    # Smooth slightly to reduce over-segmentation
    watershed_image_smooth = gaussian(watershed_image, sigma=0.5)
    
    # Generate mask for watershed (only process significant regions)
    watershed_mask = data_sub > (extended_threshold_factor * sigma * s)
    
    # Apply watershed
    seg_watershed = watershed(watershed_image_smooth, markers_dilated, 
                             mask=watershed_mask, 
                             connectivity=watershed_connectivity)
    
    print(f"      Watershed segmentation: {len(np.unique(seg_watershed)) - 1} regions")
    
    # ========================
    # 5. EXTENDED STRUCTURE RECOVERY AND MERGING
    # ========================
    print(f"\n[5/7] Recovering extended structures...")
    
    # Combine watershed result with low-threshold extended detections
    seg_maps = seg_watershed.copy()
    
    if preserve_extended:
        # Add extended regions that don't overlap with cores
        extended_only = (seg_ext > 0) & (seg_watershed == 0)
        
        # Label these extended regions
        from scipy.ndimage import label as nd_label
        extended_labeled, n_extended = nd_label(extended_only)
        
        # Add to segmentation map with new labels
        max_label = seg_maps.max()
        for i in range(1, n_extended + 1):
            extended_region = (extended_labeled == i)
            # Check if large enough
            if np.sum(extended_region) >= npixels:
                max_label += 1
                seg_maps[extended_region] = max_label
        
        print(f"      Added {n_extended} extended-only regions")
    
    # ========================
    # 6. SOURCE PROPERTY EXTRACTION AND MERGING
    # ========================
    print(f"\n[6/7] Computing source properties and applying merging...")
    
    # Get unique labels
    unique_labels = np.unique(seg_maps)
    unique_labels = unique_labels[unique_labels > 0]  # Remove background
    
    # Extract properties for each region
    objects_list = []
    masks_regions = []
    
    for label_id in unique_labels:
        region_mask = (seg_maps == label_id)
        region_data = data_sub * region_mask
        
        # Compute centroid
        y_coords, x_coords = np.where(region_mask)
        weights = region_data[region_mask]
        
        if np.sum(weights) > 0:
            xc = np.average(x_coords, weights=weights)
            yc = np.average(y_coords, weights=weights)
        else:
            xc = np.mean(x_coords)
            yc = np.mean(y_coords)
        
        # Compute moments for shape
        y_centered = y_coords - yc
        x_centered = x_coords - xc
        
        # Second moments
        if len(weights) > 0 and np.sum(weights) > 0:
            mxx = np.average(x_centered**2, weights=weights)
            myy = np.average(y_centered**2, weights=weights)
            mxy = np.average(x_centered * y_centered, weights=weights)
            
            # Eigenvalues for ellipse parameters
            trace = mxx + myy
            det = mxx * myy - mxy**2
            if det > 0:
                eig_sum = trace
                eig_diff = np.sqrt(trace**2 - 4*det)
                a_sq = 0.5 * (eig_sum + eig_diff)
                b_sq = 0.5 * (eig_sum - eig_diff)
                
                a = np.sqrt(np.maximum(a_sq, 0)) * 2  # Semi-major axis
                b = np.sqrt(np.maximum(b_sq, 0)) * 2  # Semi-minor axis
                
                # Position angle
                if mxx != myy:
                    theta = 0.5 * np.arctan2(2 * mxy, mxx - myy)
                else:
                    theta = 0.0
            else:
                a = 3.0
                b = 3.0
                theta = 0.0
        else:
            a = 3.0
            b = 3.0
            theta = 0.0
        
        # Total flux
        flux = np.sum(region_data)
        
        # Store properties
        objects_list.append({
            'x': xc,
            'y': yc,
            'a': a,
            'b': b,
            'theta': theta,
            'flux': flux,
            'area': np.sum(region_mask),
            'label': label_id
        })
        
        masks_regions.append(region_mask)
    
    print(f"      Total regions before merging: {len(objects_list)}")
    
    # Optional: Merge nearby detections that are likely part of the same source
    # (This helps avoid over-deblending)
    if merge_threshold > 0 and len(objects_list) > 1:
        from scipy.spatial.distance import pdist, squareform
        
        # Compute pairwise distances
        coords = np.array([[obj['x'], obj['y']] for obj in objects_list])
        sizes = np.array([np.mean([obj['a'], obj['b']]) for obj in objects_list])
        
        distances = squareform(pdist(coords))
        
        # Merge sources that are closer than merge_threshold * average_size
        merge_distance = merge_threshold * np.mean(sizes)
        
        # Simple greedy merging
        merged_labels = list(range(len(objects_list)))
        for i in range(len(objects_list)):
            for j in range(i+1, len(objects_list)):
                if distances[i, j] < merge_distance:
                    # Merge j into i
                    merged_labels[j] = merged_labels[i]
        
        # Regroup based on merges
        unique_merge_labels = list(set(merged_labels))
        objects_merged = []
        masks_merged = []
        
        for merge_id in unique_merge_labels:
            indices = [i for i, x in enumerate(merged_labels) if x == merge_id]
            
            if len(indices) == 1:
                # No merge needed
                objects_merged.append(objects_list[indices[0]])
                masks_merged.append(masks_regions[indices[0]])
            else:
                # Merge multiple regions
                merged_mask = np.zeros_like(masks_regions[0])
                for idx in indices:
                    merged_mask |= masks_regions[idx]
                
                # Recompute properties
                region_data = data_sub * merged_mask
                y_coords, x_coords = np.where(merged_mask)
                weights = region_data[merged_mask]
                
                if np.sum(weights) > 0:
                    xc = np.average(x_coords, weights=weights)
                    yc = np.average(y_coords, weights=weights)
                else:
                    xc = np.mean(x_coords)
                    yc = np.mean(y_coords)
                
                # Use maximum extents
                a = max([objects_list[idx]['a'] for idx in indices])
                b = max([objects_list[idx]['b'] for idx in indices])
                theta = objects_list[indices[0]]['theta']  # Use first
                flux = np.sum(region_data)
                
                objects_merged.append({
                    'x': xc,
                    'y': yc,
                    'a': a,
                    'b': b,
                    'theta': theta,
                    'flux': flux,
                    'area': np.sum(merged_mask),
                    'label': objects_list[indices[0]]['label']
                })
                
                masks_merged.append(merged_mask)
        
        objects_list = objects_merged
        masks_regions = masks_merged
        
        print(f"      Total regions after merging: {len(objects_list)}")
    
    # ========================
    # 7. SORTING AND OUTPUT PREPARATION
    # ========================
    print(f"\n[7/7] Sorting and preparing output...")
    
    # Convert objects_list to arrays for sorting
    if len(objects_list) == 0:
        print("\n*** WARNING: No sources detected! ***")
        # Return empty results matching expected format
        empty_objects = {'xc': np.array([]), 'yc': np.array([])}
        empty_seg = np.zeros_like(data_2D, dtype=np.int32)
        if segmentation_map:
            return [], np.array([]), bkg_image, empty_seg, empty_objects
        else:
            return [], np.array([]), bkg_image, empty_objects
    
    x_positions = np.array([obj['x'] for obj in objects_list])
    y_positions = np.array([obj['y'] for obj in objects_list])
    fluxes = np.array([obj['flux'] for obj in objects_list])
    areas = np.array([obj['area'] for obj in objects_list])
    
    # Sort according to specified criterion
    if sort_by == 'flux':
        sorted_indices_desc = np.argsort(fluxes)[::-1]
    elif sort_by == 'area':
        sorted_indices_desc = np.argsort(areas)[::-1]
    elif sort_by == 'distance':
        ref_centre = (data_2D.shape[0] / 2, data_2D.shape[1] / 2)
        distances = distances_from_reference(x_positions, y_positions, ref_centre)
        sorted_indices_desc = np.argsort(distances)
    else:
        sorted_indices_desc = np.arange(len(objects_list))
    
    # Create objects_sorted dictionary
    objects_sorted = {
        'xc': x_positions[sorted_indices_desc],
        'yc': y_positions[sorted_indices_desc]
    }
    
    # Sort masks
    masks_regions_sorted = [masks_regions[i] for i in sorted_indices_desc]
    
    # Create final segmentation map with sorted labels
    seg_maps_final = np.zeros_like(data_2D, dtype=np.int32)
    for i, idx in enumerate(sorted_indices_desc):
        seg_maps_final[masks_regions[idx]] = i + 1
    
    print(f"\n      Final source count: {len(objects_list)}")
    print(f"      Sorting criterion: {sort_by}")
    
    # ========================
    # 8. VISUALIZATION (if requested)
    # ========================
    if show_detection:
        fig, ax = plt.subplots(figsize=(10, 10))
        
        # Display background-subtracted image
        norm = simple_norm(data_sub, stretch='asinh', asinh_a=0.02, 
                          vmin=s, vmax=0.25*np.nanmax(data_sub))
        ax.imshow(data_sub, origin='lower', cmap='gray', norm=norm)
        
        # Determine ellipse size factor
        if ell_size_factor is None:
            if mask is not None:
                ell_size_factor = 0.1 * np.sqrt(np.nansum(mask) / np.pi)
            else:
                ell_size_factor = 1.0
        
        # Plot detected sources
        for i, idx in enumerate(sorted_indices_desc):
            obj = objects_list[idx]
            
            # Draw ellipse
            e = Ellipse(
                xy=(obj['x'], obj['y']),
                width=2 * ell_size_factor * obj['a'],
                height=2 * ell_size_factor * obj['b'],
                angle=obj['theta'] * 180. / np.pi
            )
            e.set_facecolor('none')
            e.set_edgecolor('red')
            e.set_linewidth(1.5)
            ax.add_artist(e)
            
            # Add label
            label = f"ID{i + 1}"
            label_x = obj['x'] + 3 * ell_size_factor
            label_y = obj['y'] + 20 * ell_size_factor
            line_end_y = label_y - 5
            
            text = Text(label_x, label_y, label, ha='center', va='center', 
                       color='red', fontsize=9, weight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                alpha=0.7, edgecolor='red'))
            ax.add_artist(text)
            ax.plot([obj['x'], label_x], [obj['y'], line_end_y], 
                   color='red', alpha=0.4, linewidth=1)
        
        ax.axis('off')
        plt.title(f'Robust Source Extraction: {len(objects_list)} sources detected', 
                 fontsize=14, weight='bold', pad=20)
        plt.tight_layout()
        plt.savefig(imagename + '_robust_detection.jpg', dpi=300, bbox_inches='tight')
        plt.show()
    
    print("\n" + "="*60)
    print("ROBUST EXTRACTION COMPLETE")
    print("="*60 + "\n")
    
    # Return results in expected format
    if segmentation_map:
        return (masks_regions_sorted, sorted_indices_desc, bkg_image, 
                seg_maps_final, objects_sorted)
    else:
        return (masks_regions_sorted, sorted_indices_desc, bkg_image, 
                objects_sorted)
    
"""
Robust Source Extraction - Optical Optimized Version
-----------------------------------------------------
A variant of the robust source extraction algorithm specifically optimized for 
optical imaging data (SDSS, EFIGI, HST, etc.) where:
- No beam size information is available
- Sources are compact with extended light profiles
- Background is typically smoother than radio
- Need to separate nearby sources more aggressively
"""


def robust_optical_source_ext(imagename, residualname=None, 
                              sigma=3.0, iterations=2, dilation_size=None,
                              deblend_nthresh=100, deblend_cont=0.005, maskthresh=0.0,
                              gain=1.0, filter_kernel=None, mask=None,
                              segmentation_map=False, clean_param=1.0, clean=True,
                              minarea=20, filter_type='matched', 
                              sort_by='distance',
                              bw=64, bh=64, fw=3, fh=3, ell_size_factor=None,
                              apply_mask=False, sigma_mask=6, minarea_factor=1.0,
                              npixels=None,
                              # Optical-optimized parameters
                              multiscale_levels=2,  # Fewer scales for optical
                              extended_threshold_factor=0.7,  # Less aggressive on extended
                              watershed_connectivity=1,  # More aggressive separation
                              merge_threshold=0.15,  # Less merging
                              adaptive_background=True,
                              preserve_extended=False,  # Don't add extended-only regions
                              min_separation_factor=5.0,  # Minimum separation as fraction of npixels
                              compact_core_emphasis=True,  # Focus on compact cores
                              show_bkg_map=False, show_detection=False):
    """
    Optical-optimized robust source extraction algorithm.
    
    Key differences from radio version:
    1. More aggressive source separation (lower merge_threshold, connectivity=1)
    2. Less emphasis on extended emission recovery
    3. Automatic minimum area scaling based on image size
    4. Compact core emphasis to avoid over-large regions
    5. Better heuristics when no beam size available
    
    Parameters
    ----------
    imagename : str
        Path to the FITS image file.
    sigma : float, default=3.0
        Detection threshold (lower for optical than radio).
    extended_threshold_factor : float, default=0.7
        Higher than radio version - less aggressive extended detection.
    watershed_connectivity : int, default=1
        More aggressive separation (4-connectivity).
    merge_threshold : float, default=0.15
        Much lower than radio - less merging of nearby sources.
    preserve_extended : bool, default=False
        Disabled for optical - focus on compact sources.
    min_separation_factor : float, default=5.0
        Minimum separation as multiple of sqrt(npixels).
    compact_core_emphasis : bool, default=True
        Weight toward compact cores rather than extended light.
    
    Returns
    -------
    Same format as standard robust_source_ext.
    """
    
    # ========================
    # 1. DATA LOADING AND PREPARATION
    # ========================
    print("\n" + "="*60)
    print("ROBUST OPTICAL SOURCE EXTRACTION")
    print("="*60)
    
    # Load data
    _data_2D = fitsio.read(imagename)
    if len(_data_2D.shape) == 4:
        data_2D = _data_2D[0][0]
    else:
        data_2D = _data_2D
    
    # Get residual for noise estimation if available
    if residualname is not None:
        residual_data = fitsio.read(residualname)
        if len(residual_data.shape) == 4:
            residual_data = residual_data[0][0]
        s = mad_std(residual_data)
    else:
        s = mad_std(data_2D)
    
    print(f"\n[1/7] Data loaded: shape={data_2D.shape}, RMS={s:.6e}")
    
    # Optical-specific minimum area heuristic
    if npixels is None:
        # For optical: use image-size dependent heuristic
        # Typical optical sources: ~10-100 pixels depending on resolution
        image_diagonal = np.sqrt(data_2D.shape[0]**2 + data_2D.shape[1]**2)
        npixels = max(10, int(image_diagonal / 80))  # Scale with image size
        npixels = int(npixels * minarea_factor)
    
    # Minimum separation for optical
    min_separation = min_separation_factor * np.sqrt(npixels)
    
    print(f"      Minimum detection area: {npixels} pixels")
    print(f"      Minimum source separation: {min_separation:.1f} pixels")
    
    # ========================
    # 2. ADAPTIVE BACKGROUND ESTIMATION
    # ========================
    print(f"\n[2/7] Computing adaptive background...")
    
    if adaptive_background:
        # Optical typically has smoother backgrounds - use fewer scales
        bkg_scales = []
        for scale_factor in [1, 2]:  # Only 2 scales for optical
            scale_bw = int(bw * scale_factor)
            scale_bh = int(bh * scale_factor)
            scale_fw = max(3, int(fw * scale_factor))
            scale_fh = max(3, int(fh * scale_factor))
            
            bkg_temp = sep.Background(data_2D, mask=mask, 
                                     bw=scale_bw, bh=scale_bh, 
                                     fw=scale_fw, fh=scale_fh)
            bkg_scales.append(bkg_temp.back())
        
        bkg_image = np.min(bkg_scales, axis=0)
        bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
    else:
        bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
        bkg_image = bkg.back()
    
    # Background-subtracted data
    data_sub = data_2D - bkg_image
    
    if mask is not None:
        data_sub = data_sub * mask
    
    bkg_rms = bkg.rms()
    
    if show_bkg_map:
        fig, axes = plt.subplots(1, 3, figsize=(15, 4))
        
        axes[0].imshow(data_2D, origin='lower', cmap='gray')
        axes[0].set_title('Original Data')
        axes[0].axis('off')
        
        axes[1].imshow(bkg_image, origin='lower', cmap='gray')
        axes[1].set_title('Background')
        axes[1].axis('off')
        
        norm = simple_norm(data_sub, stretch='asinh', asinh_a=0.02, 
                          vmin=s, vmax=0.25*np.nanmax(data_sub))
        axes[2].imshow(data_sub, origin='lower', cmap='gray', norm=norm)
        axes[2].set_title('Background-Subtracted')
        axes[2].axis('off')
        
        plt.tight_layout()
        plt.savefig(imagename + '_robust_opt_background.jpg', dpi=150, bbox_inches='tight')
        plt.show()
    
    # ========================
    # 3. COMPACT CORE DETECTION
    # ========================
    print(f"\n[3/7] Detecting compact source cores...")
    
    # High-threshold detection for cores - more aggressive for optical
    thresh_high = sigma * s
    
    if compact_core_emphasis:
        # Use smaller minimum area for core detection
        core_npixels = max(5, int(npixels * 0.5))
    else:
        core_npixels = npixels
    
    objects_core, seg_core = sep.extract(
        data_sub, thresh=thresh_high,
        minarea=core_npixels,
        filter_type=filter_type,
        deblend_nthresh=deblend_nthresh,
        deblend_cont=deblend_cont,
        clean=clean, clean_param=clean_param,
        segmentation_map=True,
        err=None, mask=None
    )
    
    print(f"      Compact core detection (sigma={sigma:.1f}): {len(objects_core)} cores found")
    
    # ========================
    # 4. AGGRESSIVE WATERSHED SEPARATION
    # ========================
    print(f"\n[4/7] Applying aggressive watershed separation...")
    
    # Create markers from core detections
    markers = np.zeros_like(data_sub, dtype=np.int32)
    
    for i in range(len(objects_core)):
        x_int = int(np.round(objects_core['x'][i]))
        y_int = int(np.round(objects_core['y'][i]))
        
        if 0 <= y_int < markers.shape[0] and 0 <= x_int < markers.shape[1]:
            markers[y_int, x_int] = i + 1
    
    # Dilate markers based on min_separation - but keep small for optical
    from scipy.ndimage import binary_dilation
    structure = disk(max(2, int(min_separation / 4)))  # Smaller dilation for optical
    markers_dilated = np.zeros_like(markers)
    for i in range(1, len(objects_core) + 1):
        marker_i = (markers == i)
        markers_dilated += (binary_dilation(marker_i, structure) * i).astype(np.int32)
    
    # Watershed with optical-appropriate settings
    from skimage.segmentation import watershed
    from skimage.filters import gaussian
    
    watershed_image = -data_sub
    watershed_image[watershed_image < 0] = 0
    
    # Less smoothing for optical
    watershed_image_smooth = gaussian(watershed_image, sigma=0.3)
    
    # More restrictive mask for optical
    watershed_mask = data_sub > (extended_threshold_factor * sigma * s)
    
    # Apply watershed with aggressive connectivity
    seg_watershed = watershed(watershed_image_smooth, markers_dilated, 
                             mask=watershed_mask, 
                             connectivity=watershed_connectivity)
    
    print(f"      Watershed segmentation: {len(np.unique(seg_watershed)) - 1} regions")
    
    # ========================
    # 5. OPTICAL-SPECIFIC REGION REFINEMENT
    # ========================
    print(f"\n[5/7] Refining regions for optical sources...")
    
    seg_maps = seg_watershed.copy()
    
    # For optical, do NOT add extended-only regions
    # Focus on compact sources only
    
    # ========================
    # 6. SOURCE PROPERTY EXTRACTION
    # ========================
    print(f"\n[6/7] Computing source properties...")
    
    unique_labels = np.unique(seg_maps)
    unique_labels = unique_labels[unique_labels > 0]
    
    objects_list = []
    masks_regions = []
    
    for label_id in unique_labels:
        region_mask = (seg_maps == label_id)
        region_data = data_sub * region_mask
        
        # Compute centroid
        y_coords, x_coords = np.where(region_mask)
        weights = region_data[region_mask]
        
        if np.sum(weights) > 0:
            xc = np.average(x_coords, weights=weights)
            yc = np.average(y_coords, weights=weights)
        else:
            xc = np.mean(x_coords)
            yc = np.mean(y_coords)
        
        # Compute moments for shape
        y_centered = y_coords - yc
        x_centered = x_coords - xc
        
        if len(weights) > 0 and np.sum(weights) > 0:
            mxx = np.average(x_centered**2, weights=weights)
            myy = np.average(y_centered**2, weights=weights)
            mxy = np.average(x_centered * y_centered, weights=weights)
            
            trace = mxx + myy
            det = mxx * myy - mxy**2
            if det > 0:
                eig_sum = trace
                eig_diff = np.sqrt(trace**2 - 4*det)
                a_sq = 0.5 * (eig_sum + eig_diff)
                b_sq = 0.5 * (eig_sum - eig_diff)
                
                a = np.sqrt(np.maximum(a_sq, 0)) * 2
                b = np.sqrt(np.maximum(b_sq, 0)) * 2
                
                if mxx != myy:
                    theta = 0.5 * np.arctan2(2 * mxy, mxx - myy)
                else:
                    theta = 0.0
            else:
                a = 3.0
                b = 3.0
                theta = 0.0
        else:
            a = 3.0
            b = 3.0
            theta = 0.0
        
        flux = np.sum(region_data)
        
        objects_list.append({
            'x': xc,
            'y': yc,
            'a': a,
            'b': b,
            'theta': theta,
            'flux': flux,
            'area': np.sum(region_mask),
            'label': label_id
        })
        
        masks_regions.append(region_mask)
    
    print(f"      Total regions before merging: {len(objects_list)}")
    
    # Very conservative merging for optical
    if merge_threshold > 0 and len(objects_list) > 1:
        from scipy.spatial.distance import pdist, squareform
        
        coords = np.array([[obj['x'], obj['y']] for obj in objects_list])
        sizes = np.array([np.mean([obj['a'], obj['b']]) for obj in objects_list])
        
        distances = squareform(pdist(coords))
        
        # Much more conservative merging for optical
        merge_distance = merge_threshold * np.mean(sizes)
        
        merged_labels = list(range(len(objects_list)))
        for i in range(len(objects_list)):
            for j in range(i+1, len(objects_list)):
                if distances[i, j] < merge_distance:
                    merged_labels[j] = merged_labels[i]
        
        unique_merge_labels = list(set(merged_labels))
        objects_merged = []
        masks_merged = []
        
        for merge_id in unique_merge_labels:
            indices = [i for i, x in enumerate(merged_labels) if x == merge_id]
            
            if len(indices) == 1:
                objects_merged.append(objects_list[indices[0]])
                masks_merged.append(masks_regions[indices[0]])
            else:
                merged_mask = np.zeros_like(masks_regions[0])
                for idx in indices:
                    merged_mask |= masks_regions[idx]
                
                region_data = data_sub * merged_mask
                y_coords, x_coords = np.where(merged_mask)
                weights = region_data[merged_mask]
                
                if np.sum(weights) > 0:
                    xc = np.average(x_coords, weights=weights)
                    yc = np.average(y_coords, weights=weights)
                else:
                    xc = np.mean(x_coords)
                    yc = np.mean(y_coords)
                
                a = max([objects_list[idx]['a'] for idx in indices])
                b = max([objects_list[idx]['b'] for idx in indices])
                theta = objects_list[indices[0]]['theta']
                flux = np.sum(region_data)
                
                objects_merged.append({
                    'x': xc,
                    'y': yc,
                    'a': a,
                    'b': b,
                    'theta': theta,
                    'flux': flux,
                    'area': np.sum(merged_mask),
                    'label': objects_list[indices[0]]['label']
                })
                
                masks_merged.append(merged_mask)
        
        objects_list = objects_merged
        masks_regions = masks_merged
        
        print(f"      Total regions after merging: {len(objects_list)}")
    
    # ========================
    # 7. SORTING AND OUTPUT
    # ========================
    print(f"\n[7/7] Sorting and preparing output...")
    
    if len(objects_list) == 0:
        print("\n*** WARNING: No sources detected! ***")
        empty_objects = {'xc': np.array([]), 'yc': np.array([])}
        empty_seg = np.zeros_like(data_2D, dtype=np.int32)
        if segmentation_map:
            return [], np.array([]), bkg_image, empty_seg, empty_objects
        else:
            return [], np.array([]), bkg_image, empty_objects
    
    x_positions = np.array([obj['x'] for obj in objects_list])
    y_positions = np.array([obj['y'] for obj in objects_list])
    fluxes = np.array([obj['flux'] for obj in objects_list])
    areas = np.array([obj['area'] for obj in objects_list])
    
    if sort_by == 'flux':
        sorted_indices_desc = np.argsort(fluxes)[::-1]
    elif sort_by == 'area':
        sorted_indices_desc = np.argsort(areas)[::-1]
    elif sort_by == 'distance':
        ref_centre = (data_2D.shape[0] / 2, data_2D.shape[1] / 2)
        distances = distances_from_reference(x_positions, y_positions, ref_centre)
        sorted_indices_desc = np.argsort(distances)
    else:
        sorted_indices_desc = np.arange(len(objects_list))
    
    objects_sorted = {
        'xc': x_positions[sorted_indices_desc],
        'yc': y_positions[sorted_indices_desc]
    }
    
    masks_regions_sorted = [masks_regions[i] for i in sorted_indices_desc]
    
    seg_maps_final = np.zeros_like(data_2D, dtype=np.int32)
    for i, idx in enumerate(sorted_indices_desc):
        seg_maps_final[masks_regions[idx]] = i + 1
    
    print(f"\n      Final source count: {len(objects_list)}")
    print(f"      Sorting criterion: {sort_by}")
    
    # ========================
    # 8. VISUALIZATION
    # ========================
    if show_detection:
        fig, ax = plt.subplots(figsize=(10, 10))
        
        norm = simple_norm(data_sub, stretch='asinh', asinh_a=0.02, 
                          vmin=s, vmax=0.25*np.nanmax(data_sub))
        ax.imshow(data_sub, origin='lower', cmap='gray', norm=norm)
        
        # Optical-appropriate ellipse sizing
        if ell_size_factor is None:
            # For optical, use smaller ellipses relative to detected size
            ell_size_factor = 0.8  # Smaller than radio default
        
        for i, idx in enumerate(sorted_indices_desc):
            obj = objects_list[idx]
            
            e = Ellipse(
                xy=(obj['x'], obj['y']),
                width=2 * ell_size_factor * obj['a'],
                height=2 * ell_size_factor * obj['b'],
                angle=obj['theta'] * 180. / np.pi
            )
            e.set_facecolor('none')
            e.set_edgecolor('red')
            e.set_linewidth(1.5)
            ax.add_artist(e)
            
            label = f"ID{i + 1}"
            label_x = obj['x'] + 2 * ell_size_factor * obj['a']
            label_y = obj['y'] + 2 * ell_size_factor * obj['b']
            
            text = Text(label_x, label_y, label, ha='center', va='center', 
                       color='red', fontsize=9, weight='bold',
                       bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                                alpha=0.7, edgecolor='red'))
            ax.add_artist(text)
        
        ax.axis('off')
        plt.title(f'Robust Optical Extraction: {len(objects_list)} sources detected', 
                 fontsize=14, weight='bold', pad=20)
        plt.tight_layout()
        plt.savefig(imagename + '_robust_opt_detection.jpg', dpi=300, bbox_inches='tight')
        plt.show()
    
    print("\n" + "="*60)
    print("ROBUST OPTICAL EXTRACTION COMPLETE")
    print("="*60 + "\n")
    
    if segmentation_map:
        return (masks_regions_sorted, sorted_indices_desc, bkg_image, 
                seg_maps_final, objects_sorted)
    else:
        return (masks_regions_sorted, sorted_indices_desc, bkg_image, 
                objects_sorted)


# def plot_bkg_info(
#     background,
#     data,
#     mad_std_func=None,
#     figsize=(14, 5),
#     cmap='viridis',
#     show=True,
#     save_path=None
# ):
#     """
#     Plot background information including 2D map and 1D profiles.
    
#     Parameters
#     ----------
#     background : 2D array
#         Background map to visualize
#     data : 2D array
#         Original data for comparison
#     mad_std_func : callable, optional
#         Function to compute MAD standard deviation. If None, uses np.std
#     figsize : tuple, optional
#         Figure size (width, height). Default is (14, 5)
#     cmap : str, optional
#         Colormap for the background image. Default is 'viridis'
#     show : bool, optional
#         Whether to display the plot. Default is True
#     save_path : str, optional
#         Path to save the figure. If None, figure is not saved
    
#     Returns
#     -------
#     fig : matplotlib.figure.Figure
#         The figure object containing the plots
#     axes : array of matplotlib.axes.Axes
#         Array containing the two subplot axes
#     """
    
#     # Use provided MAD function or fall back to standard deviation
#     if mad_std_func is None:
#         mad_std_func = np.std
    
#     # Calculate statistics
#     bkg_mad_std = mad_std(background,ignore_nan=True)
#     bkg_min = np.nanmin(background)
#     bkg_max = np.nanmax(background)
#     data_max = np.nanmax(data)
#     bkg_data_ratio = bkg_max / data_max if data_max != 0 else np.nan
    
#     # Create figure with two subplots
#     fig, axes = plt.subplots(1, 2, figsize=figsize)
    
#     # Left panel: 2D background map
#     ax1 = axes[0]
#     im = ax1.imshow(background, origin='lower', cmap=cmap, interpolation='nearest')
#     cbar = plt.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
#     cbar.set_label('Background Value', fontsize=10)
#     ax1.set_xlabel('X Pixel', fontsize=11)
#     ax1.set_ylabel('Y Pixel', fontsize=11)
#     ax1.set_title('Background Map', fontsize=12, fontweight='bold')
    
#     # Add statistics annotation
#     stats_text = (
#         f'MAD STD: {bkg_mad_std:.3f}\n'
#         f'Min: {bkg_min:.3f}\n'
#         f'Max: {bkg_max:.3f}\n'
#         f'Data Max: {data_max:.3f}\n'
#         f'Bkg/Data: {bkg_data_ratio:.6f}'
#     )
#     ax1.annotate(
#         stats_text,
#         xy=(0.02, 0.98),
#         xycoords='axes fraction',
#         fontsize=9,
#         color='black',
#         verticalalignment='top',
#         bbox=dict(boxstyle="round,pad=0.5", facecolor='wheat', 
#                   edgecolor='black', linewidth=1.5, alpha=0.9)
#     )
    
#     # Right panel: 1D profiles (mean along columns)
#     ax2 = axes[1]
#     y_coords = np.arange(data.shape[0])
    
#     ax2.plot(y_coords, np.nanmean(data, axis=1), 
#              label='Data', linewidth=2, alpha=0.8)
#     ax2.plot(y_coords, np.nanmean(background, axis=1), 
#              label='Background', linewidth=2, alpha=0.8)
#     ax2.plot(y_coords, np.nanmean(data - background, axis=1), 
#              label='Data - Background', linewidth=2, alpha=0.8, linestyle='--')
    
#     ax2.set_xlabel('Y Pixel', fontsize=11)
#     ax2.set_ylabel('Mean Value', fontsize=11)
#     ax2.set_title('Mean Profiles (averaged along X)', fontsize=12, fontweight='bold')
#     ax2.legend(loc='best', framealpha=0.9, fontsize=9)
#     ax2.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
#     ax2.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
    
#     plt.tight_layout()
    
#     # Save if requested
#     if save_path is not None:
#         plt.savefig(save_path, dpi=150, bbox_inches='tight')
#         print(f"Figure saved to: {save_path}")
    
#     # Show if requested
#     if show:
#         plt.show()
    
#     return fig, axes

def plot_bkg_info(
    background,
    data,
    mad_std_func=None,
    profile_type='mean',
    center=None,
    binsize=1,
    figsize=(14, 5),
    cmap='viridis',
    show=True,
    save_path=None
):
    """
    Plot background information including 2D map and 1D profiles.
    
    Parameters
    ----------
    background : 2D array
        Background map to visualize
    data : 2D array
        Original data for comparison
    mad_std_func : callable, optional
        Function to compute MAD standard deviation. If None, uses np.std
    profile_type : str, optional
        Type of profile to plot: 'mean' for mean along X axis, 
        'azimuthal' for radial profiles. Default is 'mean'
    center : tuple of float, optional
        Center (x, y) for azimuthal profile. If None and profile_type='azimuthal',
        uses image center
    binsize : float, optional
        Bin size for azimuthal profile in pixels. Default is 1
    figsize : tuple, optional
        Figure size (width, height). Default is (14, 5)
    cmap : str, optional
        Colormap for the background image. Default is 'viridis'
    show : bool, optional
        Whether to display the plot. Default is True
    save_path : str, optional
        Path to save the figure. If None, figure is not saved
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the plots
    axes : array of matplotlib.axes.Axes
        Array containing the two subplot axes
    """

    # Use provided MAD function or fall back to standard deviation
    if mad_std_func is None:
        mad_std_func = np.std
    
    # Calculate statistics
    bkg_mad_std = mad_std(background,ignore_nan=True)
    bkg_min = np.nanmin(background)
    bkg_max = np.nanmax(background)
    data_max = np.nanmax(data)
    bkg_data_ratio = bkg_max / data_max if data_max != 0 else np.nan
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    
    # Left panel: 2D background map
    ax1 = axes[0]
    im = ax1.imshow(background, origin='lower', cmap=cmap, interpolation='nearest')
    cbar = plt.colorbar(im, ax=ax1, fraction=0.046, pad=0.04)
    cbar.set_label('Background Value', fontsize=10)
    ax1.set_xlabel('X Pixel', fontsize=11)
    ax1.set_ylabel('Y Pixel', fontsize=11)
    ax1.set_title('Background Map', fontsize=12, fontweight='bold')
    
    # Add statistics annotation
    stats_text = (
        f'MAD STD: {bkg_mad_std:.3f}\n'
        f'Min: {bkg_min:.3f}\n'
        f'Max: {bkg_max:.3f}\n'
        f'Data Max: {data_max:.3f}\n'
        f'Bkg/Data: {bkg_data_ratio:.6f}'
    )
    ax1.annotate(
        stats_text,
        xy=(0.02, 0.98),
        xycoords='axes fraction',
        fontsize=9,
        color='black',
        verticalalignment='top',
        bbox=dict(boxstyle="round,pad=0.5", facecolor='wheat', 
                  edgecolor='black', linewidth=1.5, alpha=0.9)
    )
    
    # Right panel: Profiles
    ax2 = axes[1]
    
    if profile_type.lower() == 'mean':
        # Mean profiles along X axis
        y_coords = np.arange(data.shape[0])
        
        ax2.plot(y_coords, np.nanmean(data, axis=1), 
                 label='Data', linewidth=2, alpha=0.8)
        ax2.plot(y_coords, np.nanmean(background, axis=1), 
                 label='Background', linewidth=2, alpha=0.8)
        ax2.plot(y_coords, np.nanmean(data - background, axis=1), 
                 label='Data - Background', linewidth=2, alpha=0.8, linestyle='--')
        
        ax2.set_xlabel('Y Pixel', fontsize=11)
        ax2.set_ylabel('Mean Value', fontsize=11)
        ax2.set_title('Mean Profiles (averaged along X)', fontsize=12, fontweight='bold')
        
    elif profile_type.lower() == 'azimuthal':
        center = nd.maximum_position(data)[::-1]
        # Azimuthal profiles
        radius_data, profile_data = get_profile(
            data, center=center, binsize=binsize, 
            interpnan=True, stddev=False, return_nr=False
        )
        radius_bkg, profile_bkg = get_profile(
            background, center=center, binsize=binsize,
            interpnan=True, stddev=False, return_nr=False
        )
        radius_sub, profile_sub = get_profile(
            data - background, center=center, binsize=binsize,
            interpnan=True, stddev=False, return_nr=False
        )
        profile_sub = profile_data - profile_bkg
        # ax2.plot(radius_data, profile_data, 
        #          label='Data', linewidth=2, alpha=0.8)
        # ax2.plot(radius_bkg, profile_bkg, 
        #          label='Background', linewidth=2, alpha=0.8)
        # ax2.plot(radius_sub, profile_sub, 
        #          label='Data - Background', linewidth=2, alpha=0.8, linestyle='--')
        ax2.set_yscale('symlog', linthresh=1e-3)  # linthresh sets linear region near zero

        ax2.plot(radius_data, profile_data, 
                label='Data', linewidth=2, alpha=0.8)
        ax2.plot(radius_bkg, profile_bkg, 
                label='Background', linewidth=2, alpha=0.8)
        ax2.plot(radius_sub, profile_sub, 
                label='Data - Background', linewidth=2, alpha=0.8, linestyle='--')
        # ax2.set_xscale('log')  # log scale for radius

        ax2.set_xlabel('Radius (pixels)', fontsize=11)
        ax2.set_ylabel('Azimuthal Average', fontsize=11)
        
        # Add center marker on left panel
        if center is None:
            center = (background.shape[1] / 2, background.shape[0] / 2)
        ax1.plot(center[0], center[1], 'r+', markersize=15, markeredgewidth=2, 
                 label='Profile Center')
        ax1.legend(loc='upper right', fontsize=8)
        
        ax2.set_title('Azimuthal Profiles', fontsize=12, fontweight='bold')
        
    else:
        raise ValueError(f"profile_type must be 'mean' or 'azimuthal', got '{profile_type}'")
    
    ax2.legend(loc='best', framealpha=0.9, fontsize=9)
    ax2.grid(True, alpha=0.3, linestyle=':', linewidth=0.5)
    ax2.axhline(0, color='k', linestyle='-', linewidth=0.5, alpha=0.3)
    
    plt.tight_layout()
    
    # Save if requested
    if save_path is not None:
        plt.savefig(save_path, dpi=150, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")
    
    # Show if requested
    if show:
        plt.show()
    
    return fig, axes

def estimate_RMS_map(data, box_size = (16, 16), 
                     filter_size=(3, 3),
                     sigma=3.0, exclude_percentile = 80):
    """
    Estimate the RMS background map using Background2D from photutils.

    Parameters
    ----------
    data : 2D array
        Input image data
    box_size : tuple of int
        Size of the box to use for background estimation
    filter_size : tuple of int
        Size of the filter to smooth the background
    sigma : float
        Sigma value for sigma clipping
    exclude_percentile : float
        Percentile of pixel values to exclude when estimating background
    Returns
    -------
    background_rms : 2D array
        Estimated RMS background map
    """

    sigma_clip = SigmaClip(sigma=sigma)
    bkg_estimator = MMMBackground()
    bkg = Background2D(data = data, box_size = box_size, 
                       filter_size=filter_size,
                       exclude_percentile=exclude_percentile,
                       sigma_clip=sigma_clip, 
                       bkg_estimator=bkg_estimator)
    return bkg.background_rms


def estimate_RMS_map_bi(data, box_size = (16, 16), 
                     filter_size=(3, 3),
                     sigma=3.0, exclude_percentile = 80):
    """
    Estimate the RMS background map using Background2D from photutils.

    Parameters
    ----------
    data : 2D array
        Input image data
    box_size : tuple of int
        Size of the box to use for background estimation
    filter_size : tuple of int
        Size of the filter to smooth the background
    sigma : float
        Sigma value for sigma clipping
    exclude_percentile : float
        Percentile of pixel values to exclude when estimating background
    Returns
    -------
    background_rms : 2D array
        Estimated RMS background map
    """
    from photutils.background import BiweightScaleBackgroundRMS

    sigma_clip = SigmaClip(sigma=sigma)
    bkgrms = BiweightScaleBackgroundRMS(sigma_clip=sigma_clip)
    bkgrms_value = bkgrms.calc_background_rms(data)
    

    return bkg.background_rms

# def estimate_RMS_map(data, box_size=(50, 50), 
#                      filter_size=(25, 25),
#                      sigma=3.0, exclude_percentile=90,
#                      mask=None, coverage_map=None):
#     """
#     Estimate the RMS background map using Background2D from photutils.
    
#     Parameters
#     ----------
#     data : 2D array
#         Input image data
#     box_size : tuple of int
#         Size of boxes for background estimation (default: 50x50 pixels)
#     filter_size : tuple of int  
#         Size of median filter for smoothing (default: 50x50 pixels)
#     sigma : float
#         Sigma value for sigma clipping (default: 3.0)
#     exclude_percentile : float
#         Percentile of brightest pixels to exclude (default: 10)
#     mask : 2D boolean array, optional
#         Mask of pixels to exclude (True = exclude)
#     coverage_map : 2D array, optional
#         Coverage or sensitivity map; low-coverage regions will be masked
        
#     Returns
#     -------
#     background_rms : 2D array
#         Estimated RMS background map
#     """
#     # Build combined mask
#     combined_mask = mask.copy() if mask is not None else np.zeros_like(data, dtype=bool)
    
#     # Mask low-coverage regions if coverage map provided
#     if coverage_map is not None:
#         combined_mask |= (coverage_map < 0.1 * np.nanmedian(coverage_map))
    
#     # Mask edge regions with zero or very low signal
#     combined_mask |= (np.abs(data) < 1e-10)
    
#     sigma_clip = SigmaClip(sigma=sigma)
#     bkg_estimator = MADStdBackgroundRMS()
    
#     bkg = Background2D(data=data, 
#                        box_size=box_size,
#                        filter_size=filter_size,
#                        exclude_percentile=exclude_percentile,
#                        sigma_clip=sigma_clip,
#                        bkg_estimator=bkg_estimator,
#                        mask=combined_mask)
    
#     return bkg.background_rms

def s_segmentation_background(data, 
                              beam_size_px=3, 
                              n_scales=4,
                              profile_type='mean',
                              do_plot=False):
    """
    Simple multi-scale background estimation using segmentation masking.
    """
    n_scales = 4
    box_sizes = [64,32,16,8]
    nsigmas = [1.0, 2.0, 3.0, 4.0]
    fiters = [21, 15, 9, 5]
    # fiters = [5,9,15,21]

    
    bkg_estimator = MMMBackground()
    final_bkg = np.zeros_like(data)
    final_rms = np.zeros_like(data)



    for scale in range(n_scales):
        try:
            sigma_clip = SigmaClip(sigma=nsigmas[scale])
            bkg_estimator = MMMBackground(sigma_clip=sigma_clip)
            threshold = detect_threshold(data, nsigma=nsigmas[scale], sigma_clip=sigma_clip)
            kernel = Gaussian2DKernel(x_stddev=beam_size_px * (fiters[scale]/10.0), y_stddev=beam_size_px * (fiters[scale]/10.0))
            convolved_data = convolve(data, kernel)
            segm = detect_sources(convolved_data, threshold, npixels=10)

            source_mask = segm.data.astype(bool)

            bkg = Background2D(data, (box_sizes[scale], box_sizes[scale]), 
                            filter_size=(fiters[scale], fiters[scale]),
                            sigma_clip=sigma_clip, 
                            bkg_estimator=bkg_estimator,
                            mask=source_mask,
                            exclude_percentile=50.0)
            
            final_bkg += bkg.background
            final_rms += bkg.background_rms**2
        except Exception as e:
            print(f"Scale {scale+1}: Background estimation failed with error: {e}")
            continue

    final_rms = np.sqrt(final_rms) / n_scales
    final_bkg /= n_scales
    if do_plot:
        plot_bkg_info(final_bkg, data, show=True,profile_type=profile_type)

    return final_bkg, final_rms, None

def simple_background_estimation(data,nsigma=2.0,sigma_sigma_clip=2.0,
                                 box_size=7,filter_size=(7,7),
                                 x_stddev=1.0,y_stddev=1.0, npixels=10,
                                 exclude_percentile=10.0,
                                 profile_type='mean',
                                 do_plot=False):
    """
    Simple background estimation using segmentation to mask sources.
    """
    sigma_clip = SigmaClip(sigma=sigma_sigma_clip)
    bkg_estimator = MMMBackground()
    # First pass: detect bright sources
    threshold = detect_threshold(data, nsigma=sigma_sigma_clip, sigma_clip=sigma_clip)
    kernel = Gaussian2DKernel(x_stddev=x_stddev,y_stddev=y_stddev)
    convolved_data = convolve(data, kernel)
    segm = detect_sources(convolved_data, threshold, npixels=npixels)

    # Create mask from segmentation map
    source_mask = segm.data.astype(bool)

    # Now estimate background excluding masked regions
    bkg = Background2D(data, (box_size, box_size), filter_size=filter_size,
                    sigma_clip=sigma_clip, 
                    bkg_estimator=bkg_estimator,
                    mask=source_mask,
                    exclude_percentile=exclude_percentile)
    if do_plot:
        plot_bkg_info(bkg.background, data, show=True,profile_type=profile_type)
    return bkg.background, bkg.background_rms, None


def multiscale_segmentation_background_old(data, beam_size_px=10, 
                                      n_scales=5,
                                      detection_sigmas=None,
                                      box_sizes_in_beams=None,
                                      kernel_sizes=None,
                                      exclude_percentile=50.0,
                                      sigma_sigma_clip=2.0,
                                      maxiters=10,
                                      combine_method='finest_valid',
                                      profile_type='mean',
                                      do_plot=False):
    """
    Multi-scale background estimation using iterative source detection and masking.
    
    This function builds upon the photutils segmentation approach by applying
    source detection at multiple threshold levels and spatial scales, creating
    progressively refined background estimates that capture granular structure
    while minimising source contribution.
    
    Parameters
    ----------
    data : ndarray
        Input 2D image array
    beam_size_px : float
        Beam size in pixels, used as reference scale
    n_scales : int
        Number of detection scales to employ
    detection_sigmas : list of float, optional
        Detection thresholds in sigma units for each scale
        Default: [3.0, 2.0, 1.5, 1.0] for 4 scales
    box_sizes_in_beams : list of float, optional
        Background2D box sizes as multiples of beam size for each scale
        Default: [4, 3, 2, 1.5] for 4 scales
    kernel_sizes : list of float, optional
        Gaussian kernel sizes for source detection smoothing at each scale
        Default: [2.0, 1.5, 1.0, 0.5] for 4 scales
    combine_method : str
        Method for combining multi-scale backgrounds:
        'finest_valid' - use finest scale with sufficient unmasked area
        'weighted_average' - combine scales with adaptive weighting
        'minimum' - take minimum across scales
    exclude_percentile : float
        Percentile of pixel values to exclude when estimating background
        
    Returns
    -------
    background : ndarray
        Final background estimate
    bkg_rms : float
        Background RMS computed from clean regions
    diagnostics : dict
        Dictionary containing intermediate products for validation
    """
    
    # Set default parameters based on number of scales
    # if detection_sigmas is None:
    #     detection_sigmas = np.linspace(6.0, 1.0, n_scales)
    
    # if box_sizes_in_beams is None:
    #     box_sizes_in_beams = np.linspace(32.0, 3.0, n_scales)
    
    # if kernel_sizes is None:
    #     kernel_sizes = np.linspace(5.0, 0.5, n_scales)
    if detection_sigmas is None:
        detection_sigmas = np.linspace(5.0, 1.0, n_scales)
    
    if box_sizes_in_beams is None:
        box_sizes_in_beams = np.linspace(5.0, 1.0, n_scales)
    
    if kernel_sizes is None:
        kernel_sizes = np.linspace(5.0, 1.0, n_scales)

    # Ensure parameter lists match n_scales
    assert len(detection_sigmas) == n_scales
    assert len(box_sizes_in_beams) == n_scales
    assert len(kernel_sizes) == n_scales
    
    # Initialize storage for multi-scale results
    backgrounds = []
    masks = []
    background_rmss = []
    coverage_fractions = []
    
    # Progressive multi-scale detection and background estimation
    cumulative_mask = np.zeros_like(data, dtype=bool)
    
    for scale_idx in range(n_scales):
        nsigma = detection_sigmas[scale_idx]
        kernel_size = kernel_sizes[scale_idx]
        box_size = max(int(beam_size_px * box_sizes_in_beams[scale_idx]), 8)
        
        # Ensure box size is appropriate for image dimensions
        box_size = min(box_size, min(data.shape) // 4)
        
        # Detect sources at this threshold level
        sigma_clip = SigmaClip(sigma=sigma_sigma_clip, maxiters=maxiters)
        threshold = detect_threshold(data, nsigma=nsigma, sigma_clip=sigma_clip)
        
        # Smooth data for detection using scale-appropriate kernel
        kernel = Gaussian2DKernel(x_stddev=kernel_size)
        convolved_data = convolve(data, kernel)
        
        # Perform source detection
        segm = detect_sources(convolved_data, threshold, npixels=max(5, int(beam_size_px/2)))
        
        if segm is not None:
            # Create mask from segmentation map
            source_mask = segm.data.astype(bool)
            
            # Combine with cumulative mask from previous scales
            cumulative_mask = cumulative_mask | source_mask
        
        # Estimate background using cumulative mask
        # This excludes all sources detected at this and previous scales
        try:
            # bkg_estimator = MedianBackground(sigma=3.0, sigma_lower=2.0, sigma_upper=6.0, maxiters=10, cenfunc='median', stdfunc='std', grow=False)
            # bkg_estimator = MADStdBackgroundRMS()
            bkg_estimator = MMMBackground()
            # bkg_estimator = MedianBackground()
            
            # Filter size should be smaller for finer scales
            filter_size = max(3, int(box_size / 4))
            # print(f"Scale {scale_idx+1}: box_size={box_size}, filter_size={filter_size}, nsigma={nsigma}")
            if filter_size % 2 == 0:
                filter_size += 1  # Ensure odd
            
            bkg = Background2D(
                data, 
                (box_size, box_size), 
                filter_size=(filter_size, filter_size),
                sigma_clip=sigma_clip, 
                bkg_estimator=bkg_estimator,
                mask=cumulative_mask,
                exclude_percentile=exclude_percentile,
            )
            
            backgrounds.append(bkg.background)
            background_rmss.append(bkg.background_rms_median)
            
            # Compute fraction of image that is unmasked (coverage)
            coverage_fraction = np.sum(~cumulative_mask) / cumulative_mask.size
            coverage_fractions.append(coverage_fraction)
            
        except Exception as e:
            print(f"Warning: Background estimation failed at scale {scale_idx+1}: {e}")
            # Use previous background if available, otherwise use median
            if backgrounds:
                backgrounds.append(backgrounds[-1].copy())
                background_rmss.append(background_rmss[-1])
            else:
                backgrounds.append(np.full_like(data, np.median(data)))
                background_rmss.append(np.std(data))
            coverage_fractions.append(coverage_fraction if 'coverage_fraction' in locals() else 0.5)
        
        masks.append(cumulative_mask.copy())
    
    # Combine multi-scale backgrounds according to specified method
    if combine_method == 'finest_valid':
        # Use the finest scale that has sufficient coverage
        # Start from finest scale and work toward coarser
        final_background = backgrounds[-1].copy()
        final_rms = background_rmss[-1]
        
        for i in range(len(backgrounds)-1, -1, -1):
            if coverage_fractions[i] > 0.3:  # At least 30% unmasked
                final_background = backgrounds[i]
                final_rms = background_rmss[i]
                break
    
    elif combine_method == 'weighted_average':
        # Weight each scale by its coverage and inverse box size
        # Finer scales get higher weight where they have good coverage
        weights = np.zeros((len(backgrounds), data.shape[0], data.shape[1]))
        
        for i in range(len(backgrounds)):
            # Coverage-based weight
            coverage_weight = coverage_fractions[i]
            
            # Scale preference weight (favor finer scales)
            scale_weight = (i + 1) / len(backgrounds)
            
            # Combined weight
            weights[i] = coverage_weight * scale_weight
        
        # Normalize weights
        weight_sum = np.sum(weights, axis=0)
        weight_sum[weight_sum == 0] = 1.0  # Avoid division by zero
        
        final_background = np.zeros_like(data)
        for i in range(len(backgrounds)):
            final_background += (weights[i] / weight_sum) * backgrounds[i]
        
        final_rms = np.median(background_rmss)
    
    elif combine_method == 'minimum':
        # Take minimum across all scales
        # This is most conservative for source removal
        final_background = np.minimum.reduce(backgrounds)
        final_rms = np.min(background_rmss)
    
    else:
        raise ValueError(f"Unknown combine_method: {combine_method}")
    
    # Apply final smoothing at beam scale for physical consistency
    final_background = gaussian_filter(final_background, sigma=beam_size_px/4)
    
    # Compute final RMS from cleanest regions
    final_residual = data - final_background
    final_mask = masks[-1]  # Most complete mask
    
    if np.sum(~final_mask) > 100:
        clean_residual = final_residual[~final_mask]
        _, _, final_rms = sigma_clipped_stats(clean_residual, sigma=sigma_sigma_clip)
    
    # Prepare comprehensive diagnostics
    diagnostics = {
        'backgrounds_by_scale': backgrounds,
        'masks_by_scale': masks,
        'background_rms_by_scale': background_rmss,
        'coverage_fractions': coverage_fractions,
        'detection_sigmas_used': detection_sigmas,
        'box_sizes_used': [int(beam_size_px * bs) for bs in box_sizes_in_beams],
        'kernel_sizes_used': kernel_sizes,
        'final_residual': final_residual,
        'final_mask': final_mask,
        'combine_method': combine_method
    }
    if do_plot:
        plot_bkg_info(final_background, data, show=True,profile_type=profile_type)
    return final_background, final_rms, diagnostics


def multiscale_segmentation_background(data, beam_size_px=10, 
                                      n_scales=5,
                                      detection_sigmas=None,
                                      box_sizes_in_beams=None,
                                      kernel_sizes=None,
                                      exclude_percentile=50.0,
                                      sigma_sigma_clip=2.0,
                                      maxiters=10,
                                      combine_method='finest_valid',
                                      profile_type='mean',
                                      do_plot=False):
    """
    Multi-scale background estimation using iterative source detection and masking.
    
    This function builds upon the photutils segmentation approach by applying
    source detection at multiple threshold levels and spatial scales, creating
    progressively refined background estimates that capture granular structure
    while minimising source contribution.
    
    Parameters
    ----------
    data : ndarray
        Input 2D image array
    beam_size_px : float
        Beam size in pixels, used as reference scale
    n_scales : int
        Number of detection scales to employ
    detection_sigmas : list of float, optional
        Detection thresholds in sigma units for each scale
        Default: [3.0, 2.0, 1.5, 1.0] for 4 scales
    box_sizes_in_beams : list of float, optional
        Background2D box sizes as multiples of beam size for each scale
        Default: [4, 3, 2, 1.5] for 4 scales
    kernel_sizes : list of float, optional
        Gaussian kernel sizes for source detection smoothing at each scale
        Default: [2.0, 1.5, 1.0, 0.5] for 4 scales
    combine_method : str
        Method for combining multi-scale backgrounds:
        'finest_valid' - use finest scale with sufficient unmasked area
        'weighted_average' - combine scales with adaptive weighting
        'minimum' - take minimum across scales
    exclude_percentile : float
        Percentile of pixel values to exclude when estimating background
        
    Returns
    -------
    background : ndarray
        Final background estimate
    bkg_rms : ndarray
        Background RMS map (2D array)
    diagnostics : dict
        Dictionary containing intermediate products for validation
    """
    
    # Set default parameters based on number of scales
    if detection_sigmas is None:
        detection_sigmas = np.linspace(3.0, 1.0, n_scales)
    
    if box_sizes_in_beams is None:
        box_sizes_in_beams = np.linspace(3.0, 1.0, n_scales)
    
    if kernel_sizes is None:
        kernel_sizes = np.linspace(3.0, 1.0, n_scales)

    # Ensure parameter lists match n_scales
    assert len(detection_sigmas) == n_scales
    assert len(box_sizes_in_beams) == n_scales
    assert len(kernel_sizes) == n_scales
    
    # Initialize storage for multi-scale results
    backgrounds = []
    masks = []
    background_rmss = []  # Will store 2D RMS maps
    coverage_fractions = []
    # bkg_estimator = MMMBackground()
    # bkg_estimator = LocalBackground(inner_radius=8, outer_radius=25)
    bkg_estimator = SExtractorBackground()

    # Progressive multi-scale detection and background estimation
    cumulative_mask = np.zeros_like(data, dtype=bool)
    
    for scale_idx in range(n_scales):
        nsigma = detection_sigmas[scale_idx]
        kernel_size = kernel_sizes[scale_idx]
        box_size = max(int(beam_size_px * box_sizes_in_beams[scale_idx]), 8)
        
        # Ensure box size is appropriate for image dimensions
        box_size = min(box_size, min(data.shape) // 4)
        
        # Detect sources at this threshold level
        sigma_clip = SigmaClip(sigma=sigma_sigma_clip, maxiters=maxiters)
        threshold = detect_threshold(data, nsigma=nsigma, sigma_clip=sigma_clip)
        
        # Smooth data for detection using scale-appropriate kernel
        kernel = Gaussian2DKernel(x_stddev=kernel_size)
        convolved_data = convolve(data, kernel)
        
        # Perform source detection
        segm = detect_sources(convolved_data, threshold, npixels=max(5, int(beam_size_px/2)))
        
        if segm is not None:
            # Create mask from segmentation map
            source_mask = segm.data.astype(bool)
            
            # Combine with cumulative mask from previous scales
            cumulative_mask = cumulative_mask | source_mask
        
        # Estimate background using cumulative mask
        # This excludes all sources detected at this and previous scales
        try:
            
            # Filter size should be smaller for finer scales
            filter_size = max(3, int(box_size / 4))
            if filter_size % 2 == 0:
                filter_size += 1  # Ensure odd
            
            bkg = Background2D(
                data, 
                (box_size, box_size), 
                filter_size=(filter_size, filter_size),
                sigma_clip=sigma_clip, 
                bkg_estimator=bkg_estimator,
                mask=cumulative_mask,
                exclude_percentile=exclude_percentile,
            )
            
            backgrounds.append(bkg.background)
            background_rmss.append(bkg.background_rms)  # Store 2D RMS map
            
            # Compute fraction of image that is unmasked (coverage)
            coverage_fraction = np.sum(~cumulative_mask) / cumulative_mask.size
            coverage_fractions.append(coverage_fraction)
            
        except Exception as e:
            print(f"Warning: Background estimation failed at scale {scale_idx+1}: {e}")
            # Use previous background if available, otherwise use median
            if backgrounds:
                backgrounds.append(backgrounds[-1].copy())
                background_rmss.append(background_rmss[-1].copy())  # Copy 2D RMS map
            else:
                backgrounds.append(np.full_like(data, np.median(data)))
                background_rmss.append(np.full_like(data, np.std(data)))  # Create 2D RMS map
            coverage_fractions.append(coverage_fraction if 'coverage_fraction' in locals() else 0.5)
        
        masks.append(cumulative_mask.copy())
    
    # Combine multi-scale backgrounds according to specified method
    if combine_method == 'finest_valid':
        # Use the finest scale that has sufficient coverage
        # Start from finest scale and work toward coarser
        final_background = backgrounds[-1].copy()
        final_rms = background_rmss[-1].copy()  # Use 2D RMS map
        
        for i in range(len(backgrounds)-1, -1, -1):
            if coverage_fractions[i] > 0.3:  # At least 30% unmasked
                final_background = backgrounds[i]
                final_rms = background_rmss[i].copy()  # Use 2D RMS map
                break
    
    elif combine_method == 'weighted_average':
        # Weight each scale by its coverage and scale preference
        # Use a more balanced approach to prevent coarse scales from dominating
        weights = np.zeros((len(backgrounds), data.shape[0], data.shape[1]))
        
        for i in range(len(backgrounds)):
            # Apply power law to coverage to reduce advantage of high-coverage scales
            # This prevents coarse scales from dominating just due to better coverage
            coverage_weight = coverage_fractions[i] ** 0.5  # Square root dampens coverage differences
            
            # Scale preference weight (favor finer scales more aggressively)
            # Use exponential rather than linear to give finer scales more influence
            scale_weight = np.exp((i - len(backgrounds)/2) / len(backgrounds))
            
            # Combined weight
            weights[i] = coverage_weight * scale_weight
        
        # Normalize weights
        weight_sum = np.sum(weights, axis=0)
        weight_sum[weight_sum == 0] = 1.0  # Avoid division by zero
        
        final_background = np.zeros_like(data)
        final_rms = np.zeros_like(data)  # Initialize 2D RMS map
        
        for i in range(len(backgrounds)):
            final_background += (weights[i] / weight_sum) * backgrounds[i]
            # Combine RMS maps with same weights (proper error propagation)
            final_rms += (weights[i] / weight_sum)**2 * background_rmss[i]**2
        
        final_rms = np.sqrt(final_rms)  # Take square root for proper RMS combination
    
    elif combine_method == 'minimum':
        # Take minimum across all scales
        # This is most conservative for source removal
        final_background = np.minimum.reduce(backgrounds)
        # For RMS, take the corresponding RMS at each pixel from the minimum background
        final_rms = np.zeros_like(data)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                min_idx = np.argmin([bg[i, j] for bg in backgrounds])
                final_rms[i, j] = background_rmss[min_idx][i, j]
    
    else:
        raise ValueError(f"Unknown combine_method: {combine_method}")
    
    # Apply final smoothing at beam scale for physical consistency
    final_background = gaussian_filter(final_background, sigma=beam_size_px/4)
    final_rms = gaussian_filter(final_rms, sigma=beam_size_px/4)  # Also smooth RMS map
    
    # Prepare comprehensive diagnostics
    diagnostics = {
        'backgrounds_by_scale': backgrounds,
        'masks_by_scale': masks,
        'background_rms_by_scale': background_rmss,
        'coverage_fractions': coverage_fractions,
        'detection_sigmas_used': detection_sigmas,
        'box_sizes_used': [int(beam_size_px * bs) for bs in box_sizes_in_beams],
        'kernel_sizes_used': kernel_sizes,
        'final_residual': data - final_background,
        'final_mask': masks[-1],
        'combine_method': combine_method
    }
    if do_plot:
        plot_bkg_info(final_background, data, show=True, profile_type=profile_type)

    return final_background, final_rms, diagnostics


def _sigma_from_noise_inputs(rms_map=None, invvar_map=None, variance_map=None,
                             weight_map=None, verbose=True):
    """
    Reduce whatever noise description is available to a single sigma map.

    Surveys hand out noise in four different conventions and they are easy to
    confuse, so the conversion lives in one place:

        rms / sigma      -> used as-is
        inverse variance -> sigma = 1 / sqrt(invvar)
        variance         -> sigma = sqrt(variance)
        weight           -> sigma = 1 / sqrt(weight)

    Precedence when more than one is supplied is most-direct-first: an explicit
    rms map beats an inverse-variance map beats a variance map beats a weight
    map. Nothing is combined, because a survey's weight and variance maps are
    usually two views of the same numbers and averaging them would understate the
    noise.

    Legacy Survey and HSC "weight" images ARE inverse variance despite the name;
    if yours is, pass it as `invvar_map`. The two happen to convert identically
    here, so the distinction only matters for documentation.

    Returns
    -------
    sigma : ndarray or None
        Sigma map in data units, or None if nothing was supplied.
    source : str
        Which input it came from, for the diagnostics dict.
    """
    supplied = [(rms_map, 'rms_map'), (invvar_map, 'invvar_map'),
                (variance_map, 'variance_map'), (weight_map, 'weight_map')]
    supplied = [(v, n) for v, n in supplied if v is not None]

    if not supplied:
        return None, 'none'

    if len(supplied) > 1 and verbose:
        print(f"   Noise inputs supplied: {[n for _, n in supplied]}. "
              f"Using '{supplied[0][1]}' (precedence: rms > invvar > variance > "
              f"weight); the others are ignored.")

    value, name = supplied[0]
    value = np.asarray(value, dtype=float)

    if name == 'rms_map':
        sigma = value
    elif name in ('invvar_map', 'weight_map'):
        with np.errstate(divide='ignore', invalid='ignore'):
            sigma = 1.0 / np.sqrt(value)
    else:  # variance_map
        sigma = np.sqrt(value)

    # Zero weight / zero invvar means "no data here", which becomes inf sigma.
    # Leave it finite but huge so downstream weighting drives those pixels to
    # zero influence instead of producing NaNs in the residual.
    bad = ~np.isfinite(sigma)
    if np.any(bad):
        good = sigma[np.isfinite(sigma)]
        fill = (np.nanmedian(good) * 1e6) if good.size else 1.0
        sigma = np.where(bad, fill, sigma)
        if verbose:
            print(f"   {int(bad.sum())} pixel(s) had zero/invalid noise; set to "
                  f"{fill:.3g} (effectively unweighted).")

    return sigma, name


def compute_bkg_rms_maps(data, psf_fwhm_px=None, mask=None,
                         source_mask=None, detect_sigma=2.0, dilate_beams=3.0,
                         box_size_beams=12.0, filter_size=3,
                         exclude_percentile=10.0,
                         smooth_scale_beams=None,
                         adaptive_mask=True, mask_iterations=3,
                         faint_detect_sigma=1.0,
                         box_from_source_extent=False,
                         rms_box_size_beams=3.0, smooth_rms=False,
                         weight_map=None, variance_map=None, invvar_map=None,
                         rms_map=None,
                         min_free_fraction=0.35,
                         sigma_clip_sigma=3.0, maxiters=10,
                         do_plot=False, show_bkg_map=None, plot_name=None,
                         return_diagnostics=False, verbose=True):
    """
    Estimate a background map and an RMS map that do not absorb source flux.

    Why this exists alongside `multiscale_segmentation_background`
    -------------------------------------------------------------
    The multi-scale estimator combines several detection passes and then, with
    `combine_method='finest_valid'`, keeps the FINEST scale that still has 30%
    unmasked area. Finest means the smallest box -- around one beam. A box that
    small sitting under a galaxy many beams across gets filled by interpolation
    from neighbouring boxes that are themselves galaxy, so the "background" ends
    up tracking the source. On the J0014 test field the resulting map is 2.2x
    higher under the galaxy than at the frame edge; it should be flat.

    This function inverts that strategy: fix the MASK, not the box.

      1. mask sources iteratively -- detect, dilate, estimate a provisional
         background, subtract it, then re-detect at a fainter threshold on the
         subtracted image and grow the mask again. Two or three passes reach the
         faint wings that a single dilation misses;
      2. estimate the background on a modest box (a few PSF), which is now safe
         because its neighbours are genuinely off-source;
      3. measure the rms in SMALLER boxes, with the same mask, so the noise map
         retains real spatial structure;
      4. interpolate across the masked region and smooth the background only.

    Step 1 is the important one. An earlier version instead forced the box up to
    the source extent, which did stop the leakage (a box larger than the source
    cannot be filled from contaminated neighbours) but left only a handful of
    boxes across the frame -- a 3x3 grid interpolated up to full resolution, with
    no local structure left in either map. `box_from_source_extent=True` restores
    that behaviour.

    Parameters
    ----------
    data : ndarray
        2D image.
    psf_fwhm_px : float, optional
        PSF/beam FWHM in pixels; sets the physical scale for `dilate_beams`,
        `box_size_beams` and `smooth_scale_beams`. Defaults to 3 px if unknown.
    mask : ndarray, optional
        Pixels to ignore entirely (bad pixels, chip edges). True = ignore.
    source_mask : ndarray, optional
        Supply your own source mask and skip detection. True = source.
    detect_sigma : float, optional
        Detection threshold for the source mask, in sigma. Deliberately low --
        this mask exists to protect the background estimate, so over-masking is
        much cheaper than under-masking.
    dilate_beams : float, optional
        Dilate the source mask by this many PSF FWHM.
    box_size_beams : float, optional
        Background box size, in PSF FWHM. With `adaptive_mask` doing the work of
        keeping source light out, this can stay small enough for the background
        to follow real large-scale structure.
    exclude_percentile : float, optional
        Percentage of masked pixels above which a box is dropped and
        interpolated over. Low is good: a box with much of the source in it
        should not be measured from the few wing pixels that survived the mask.
    filter_size : int, optional
        Median filter size (in boxes) applied to the low-resolution background.
    smooth_scale_beams : float, optional
        Final Gaussian smoothing of the BACKGROUND, in PSF FWHM. Defaults to a
        quarter of the box actually used. The rms is not smoothed unless
        `smooth_rms` is set.
    adaptive_mask : bool, optional
        Grow the source mask iteratively (detect -> provisional background ->
        subtract -> re-detect fainter -> re-dilate). This is what makes a small
        box safe. Set False for a single detection pass.
    mask_iterations : int, optional
        Number of growth passes when `adaptive_mask` is on. Two or three is
        plenty; each pass costs one extra `Background2D`.
    faint_detect_sigma : float, optional
        Threshold, in sigma, for the re-detection passes. Lower than
        `detect_sigma` on purpose -- the point is to catch the wings.
    box_from_source_extent : bool, optional
        Restore the previous behaviour: force the box up to the largest masked
        blob's extent. Produces a very flat background and an essentially
        featureless rms map.
    rms_box_size_beams : float, optional
        Box size for the rms pass, in PSF FWHM. Smaller than the background box
        so the noise map can actually vary across the frame.
    smooth_rms : bool, optional
        Smooth the rms map with the same kernel as the background. Off by
        default: smoothing a noise map on the background's scale is what made it
        look modelled rather than measured.
    do_plot : bool, optional
        Draw a three-panel diagnostic: background map with the source mask
        outlined, a horizontal cut through the background at the source centre
        (flat reads as flat, leakage reads as a central bump), and the rms map.
        `show_bkg_map` is accepted as an alias.
    plot_name : str, optional
        Save the diagnostic figure here instead of only showing it.
    weight_map, variance_map, invvar_map, rms_map : ndarray, optional
        Externally supplied noise. See `_sigma_from_noise_inputs` for the
        conversions and the precedence. When any is given it REPLACES the
        measured RMS, since a survey's own noise model is better than anything
        measurable from a single cutout.
    min_free_fraction : float, optional
        Warn when less than this fraction of the frame is left unmasked. Below
        it the background is being interpolated over most of the image and is
        probably overestimated -- the usual cause is a cutout too small for the
        galaxy in it.
    return_diagnostics : bool, optional
        Also return a dict of intermediate products.

    Returns
    -------
    bkg_map : ndarray
        Background map, data units, same grid as `data`.
    rms_map : ndarray
        Sigma map, data units, same grid as `data`. This is the convention the
        fitter expects: a per-pixel standard deviation, NOT a variance and NOT a
        weight.
    diagnostics : dict, optional
    """
    data = np.asarray(data, dtype=float)
    if psf_fwhm_px is None or not np.isfinite(psf_fwhm_px) or psf_fwhm_px <= 0:
        psf_fwhm_px = 3.0
        if verbose:
            print("   No PSF size given; assuming 3 px for the background scales.")

    bad_mask = np.zeros(data.shape, dtype=bool) if mask is None else np.asarray(mask, bool)

    sigma_clip = SigmaClip(sigma=sigma_clip_sigma, maxiters=maxiters)
    kernel = Gaussian2DKernel(x_stddev=psf_fwhm_px / 2.355)
    dilation_px = int(np.ceil(dilate_beams * psf_fwhm_px))

    def _detect(image, nsigma, ignore):
        """Detect above `nsigma`, smooth first, then dilate. True = source."""
        try:
            threshold = detect_threshold(image, nsigma=nsigma,
                                         sigma_clip=sigma_clip, mask=ignore)
            segm = detect_sources(convolve(image, kernel), threshold,
                                  npixels=max(5, int(psf_fwhm_px)))
        except Exception:
            return np.zeros(image.shape, dtype=bool)
        if segm is None:
            return np.zeros(image.shape, dtype=bool)
        found = segm.data.astype(bool)
        if dilation_px > 0:
            found = nd.binary_dilation(found, iterations=dilation_px)
        return found

    # ---------------------------------------------------------------- source mask
    user_source_mask = source_mask is not None
    if user_source_mask:
        source_mask = np.asarray(source_mask, bool)
        if dilation_px > 0:
            source_mask = nd.binary_dilation(source_mask,
                                             iterations=dilation_px)
    else:
        source_mask = _detect(data, detect_sigma, bad_mask)

    # ------------------------------------------------------- adaptive mask growth
    # A single dilation only reaches `dilate_beams` past the last pixel above
    # `detect_sigma`. On J0014 that is 14 px, while the galaxy's wings run much
    # further, so a small background box still landed on real light -- which is
    # what forced the box up to the source extent before. Growing the mask
    # instead: estimate a provisional background, subtract it, and re-detect at a
    # fainter threshold on the flattened image, where the wings stand out.
    def _blob_extent(m):
        """Pixel extent of the largest connected blob in `m`."""
        lab, n = nd.label(m)
        if not n:
            return 0
        big = int(np.argmax(nd.sum(m, lab, range(1, n + 1)))) + 1
        ys_, xs_ = np.where(lab == big)
        return int(max(xs_.max() - xs_.min(), ys_.max() - ys_.min()))

    mask_growth = [int(source_mask.sum())]
    if adaptive_mask and not user_source_mask and mask_iterations > 0:
        for _ in range(int(mask_iterations)):
            # The provisional background is deliberately estimated on a box AS
            # LARGE AS THE SOURCE. That map is flat and featureless -- useless as
            # a final product, but exactly right as a scaffold: subtracting it
            # leaves the wings standing proud instead of absorbing them, which is
            # what lets the fainter re-detection actually find them. A small
            # probe box does the opposite; it tracks the source, flattens the
            # wings away, and the mask stops growing.
            _probe_box = max(int(round(box_size_beams * psf_fwhm_px)), 8,
                             _blob_extent(source_mask))
            _probe_box = min(_probe_box, max(8, min(data.shape) // 3))
            try:
                _probe = Background2D(
                    data, (_probe_box, _probe_box),
                    filter_size=(3, 3),
                    sigma_clip=sigma_clip,
                    bkg_estimator=SExtractorBackground(),
                    mask=(source_mask | bad_mask),
                    exclude_percentile=90.0,
                )
                flattened = data - np.asarray(_probe.background, dtype=float)
            except Exception:
                break
            grown = source_mask | _detect(flattened, faint_detect_sigma, bad_mask)
            gained = int(grown.sum()) - int(source_mask.sum())
            source_mask = grown
            mask_growth.append(int(source_mask.sum()))
            # Converged: nothing meaningful left to add.
            if gained <= 0.01 * source_mask.size:
                break
            # Runaway: the mask is eating the frame, so stop and let the
            # free-fraction warning below report it.
            if source_mask.sum() > (1.0 - min_free_fraction) * source_mask.size:
                break

    exclude = source_mask | bad_mask
    free_fraction = float(np.sum(~exclude)) / exclude.size

    # Extent of the largest masked blob: this, not the PSF, sets the box size.
    labels, nlab = nd.label(source_mask)
    source_extent_px = 0
    if nlab:
        sizes = nd.sum(source_mask, labels, range(1, nlab + 1))
        biggest = int(np.argmax(sizes)) + 1
        ys, xs = np.where(labels == biggest)
        source_extent_px = int(max(xs.max() - xs.min(), ys.max() - ys.min()))

    small_frame_warning = free_fraction < min_free_fraction
    if small_frame_warning and verbose:
        print(f"   WARNING: only {100*free_fraction:.1f}% of the frame is "
              f"source-free (threshold {100*min_free_fraction:.0f}%). The "
              f"background is being interpolated over most of the image and is "
              f"probably overestimated. Use a larger cutout, or reduce "
              f"dilate_beams / detect_sigma.")

    # ------------------------------------------------------------------ background
    # Box size and source protection are SEPARATE jobs. What determines whether
    # the source leaks into the background is whether the mask covers the wings;
    # what the box size determines is how much real structure survives. Measured
    # on J0014 r-band (384 px frame, galaxy 222 px across, PSF 3.76 px), as the
    # ratio of the background under the source to the background off it:
    #
    #     box     single 12 px dilation      adaptive mask (33% of frame)
    #     23 px            2.21                        1.59
    #     45 px             --                         1.04
    #     68 px             --                         0.96
    #     90 px             --                         1.00
    #    128 px            1.01                         --
    #
    # Forcing the box to 128 px did buy flatness, but at the price of a 3x3 grid
    # interpolated across the whole frame -- and the rms map inherited it, with a
    # spatial spread of 1.6% of its own median (i.e. constant). With the adaptive
    # mask, 12 beams is flat to within 4% and the rms keeps 12% spread.
    box_px = max(int(round(box_size_beams * psf_fwhm_px)), 8)
    if box_from_source_extent:
        box_px = max(box_px, source_extent_px)
    # Leave enough boxes across the frame to interpolate between.
    box_px = min(box_px, max(8, min(data.shape) // 3))
    fsize = int(filter_size)
    if fsize % 2 == 0:
        fsize += 1

    def _try_background(box, excl):
        bkg = Background2D(
            data, (box, box),
            filter_size=(fsize, fsize),
            sigma_clip=SigmaClip(sigma=sigma_clip_sigma, maxiters=maxiters),
            bkg_estimator=SExtractorBackground(),
            mask=exclude,
            exclude_percentile=excl,
        )
        return (np.asarray(bkg.background, dtype=float),
                np.asarray(bkg.background_rms, dtype=float))

    # A LOW exclude_percentile is what we want: any box with a meaningful amount
    # of source in it should be dropped and interpolated over, rather than
    # measured from the handful of wing pixels that survived the mask. Raising it
    # is only a fallback for when that leaves no usable boxes at all.
    bkg_map = measured_rms = None
    attempts = ([(box_px, p) for p in (exclude_percentile, 50.0, 90.0)]
                + [(max(8, box_px // 2), p) for p in (exclude_percentile, 50.0, 90.0)])
    for attempt_box, attempt_excl in attempts:
        try:
            bkg_map, measured_rms = _try_background(attempt_box, attempt_excl)
            if (attempt_box, attempt_excl) != (box_px, exclude_percentile) and verbose:
                print(f"   Relaxed to box={attempt_box} px, "
                      f"exclude_percentile={attempt_excl} to get a usable "
                      f"background (source covers a large part of the frame).")
            box_px, exclude_percentile = attempt_box, attempt_excl
            break
        except Exception:
            continue

    if bkg_map is None:
        if verbose:
            print("   Background2D found no usable box; falling back to a "
                  "constant sky from the source-free pixels. Treat this "
                  "background as unreliable.")
        free_pixels = data[~exclude]
        if free_pixels.size == 0:
            free_pixels = data.ravel()
        bkg_map = np.full(data.shape, float(np.nanmedian(free_pixels)))
        measured_rms = np.full(data.shape, float(mad_std(free_pixels)))
        small_frame_warning = True

    # ------------------------------------------------------------- rms, measured
    # A separate, smaller-box pass with the same mask. Sharing the background's
    # box meant the rms was the same coarse grid, so it carried no local noise
    # information at all -- it looked modelled because effectively it was.
    rms_box_px = max(int(round(rms_box_size_beams * psf_fwhm_px)), 8)
    rms_box_px = min(rms_box_px, box_px, max(8, min(data.shape) // 4))
    if rms_box_px < box_px:
        for _excl in (exclude_percentile, 50.0, 90.0):
            try:
                _, measured_rms = _try_background(rms_box_px, _excl)
                break
            except Exception:
                continue
    else:
        rms_box_px = box_px

    if smooth_scale_beams is None:
        smooth_px = box_px / 4.0
    else:
        smooth_px = smooth_scale_beams * psf_fwhm_px
    if smooth_px > 0:
        bkg_map = gaussian_filter(bkg_map, sigma=smooth_px)
        if smooth_rms:
            measured_rms = gaussian_filter(measured_rms,
                                           sigma=smooth_px * rms_box_px / box_px)

    # ----------------------------------------------------------------------- rms
    external_sigma, sigma_source = _sigma_from_noise_inputs(
        rms_map=rms_map, invvar_map=invvar_map, variance_map=variance_map,
        weight_map=weight_map, verbose=verbose)

    if external_sigma is not None:
        rms_out = external_sigma
    else:
        rms_out = measured_rms
        sigma_source = 'measured'

    # Two independent quality checks, because they fail in different ways.
    rms_level = float(np.nanmedian(rms_out))

    # (a) Flatness: is the background still tracking the source? `free_fraction`
    #     only sees the DETECTED source, so a low surface brightness halo that
    #     never got masked pulls the background up under the galaxy without
    #     tripping it. This compares the level under the source against the level
    #     off it.
    bkg_excess = np.nan
    if np.any(source_mask) and np.any(~exclude):
        bkg_excess = float(np.nanmedian(bkg_map[source_mask])
                           - np.nanmedian(bkg_map[~exclude]))
        if rms_level > 0 and abs(bkg_excess) > rms_level:
            small_frame_warning = True
            if verbose:
                print(f"   WARNING: the background is {bkg_excess:+.3g} "
                      f"({abs(bkg_excess)/rms_level:.1f} sigma) higher under the "
                      f"source than off it, so it is still tracking the source. "
                      f"Raise mask_iterations or dilate_beams, lower "
                      f"faint_detect_sigma, increase box_size_beams, or use a "
                      f"larger cutout.")

    # (b) Edge gradient: does the source reach the frame border? When it does,
    #     there is no sky anywhere in the cutout, the background comes out FLAT
    #     but at the source's own level, and check (a) cannot see it -- the map
    #     really is flat, just wrong. The signature is that the image is still
    #     falling off at the border instead of being noise.
    ny, nx = data.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    rr = np.hypot(xx - (nx - 1) / 2.0, yy - (ny - 1) / 2.0)
    r_out = 0.5 * min(ny, nx)
    inner_ring = (rr > 0.75 * r_out) & (rr <= 0.85 * r_out)
    outer_ring = (rr > 0.90 * r_out)
    edge_gradient = np.nan
    edge_sigma = np.nan
    if np.any(inner_ring) and np.any(outer_ring):
        edge_gradient = float(np.nanmedian(data[inner_ring])
                              - np.nanmedian(data[outer_ring]))
        # The noise has to come from a high-pass view of the image. Measuring it
        # directly in the outer ring does not work: when the source reaches the
        # border, the ring's spread is dominated by the very gradient this is
        # trying to detect, so the "noise" inflates and hides it.
        #
        # Differences between neighbouring pixels give that view. (Subtracting a
        # 3x3 median does not: the median of a 3x3 patch is frequently the centre
        # pixel itself, so more than half the differences are exactly zero and
        # mad_std returns exactly 0.) Adjacent-pixel differences of independent
        # noise have variance 2*sigma^2, hence the sqrt(2).
        edge_noise = float(mad_std(np.diff(data, axis=0))) / np.sqrt(2.0)
        # Compare against the standard error of the two MEDIANS, not against the
        # per-pixel noise. Each ring holds tens of thousands of pixels, so a
        # gradient far below one pixel's sigma is still overwhelming evidence.
        # Testing it per pixel is why this guard stayed silent on J0014, where
        # the profile is still declining at the corners: the drop is 0.0038
        # against a 0.0064 pixel noise (0.6 sigma per pixel) but 84 standard
        # errors of the median.
        n_in = int(np.sum(inner_ring))
        n_out = int(np.sum(outer_ring))
        sem = (1.2533 * edge_noise * np.sqrt(1.0 / max(n_in, 1)
                                             + 1.0 / max(n_out, 1)))
        edge_sigma = edge_gradient / sem if sem > 0 else 0.0
        # Two conditions, because either alone misfires: significance alone flags
        # harmless sub-percent tilts in deep images, and size alone flags noise.
        bkg_level_for_edge = abs(float(np.nanmedian(bkg_map)))
        if (edge_sigma > 5.0
                and edge_gradient > 0.25 * max(bkg_level_for_edge, rms_level)):
            small_frame_warning = True
            if verbose:
                print(f"   WARNING: the image is still falling off at the frame "
                      f"border ({edge_gradient:+.3g}, {edge_sigma:.0f} standard "
                      f"errors, {100*edge_gradient/max(bkg_level_for_edge, 1e-30):.0f}% "
                      f"of the estimated background). The source fills the "
                      f"cutout, so there is no sky in it to measure and this "
                      f"background is a source level, not a sky level. Use a "
                      f"larger cutout.")

    # How much real structure the rms map retained. A flat number here means the
    # map is an interpolated constant and the weights it produces carry no
    # spatial information.
    rms_structure = float(np.nanstd(rms_out) / rms_level) if rms_level else np.nan

    if verbose:
        residual_free = (data - bkg_map)[~exclude]
        print(f"   bkg: box={box_px} px, rms box={rms_box_px} px "
              f"(source extent {source_extent_px} px), "
              f"dilation={dilation_px} px, source-free={100*free_fraction:.1f}%")
        if len(mask_growth) > 1:
            print(f"   adaptive mask: "
                  f"{' -> '.join(f'{100.0*g/source_mask.size:.1f}%' for g in mask_growth)}"
                  f" of the frame over {len(mask_growth)-1} pass(es)")
        print(f"   bkg level={np.nanmedian(bkg_map):.5g}, "
              f"excess under source={bkg_excess:+.3g} "
              f"({abs(bkg_excess)/rms_level if rms_level else float('nan'):.2f} "
              f"sigma, should be ~0)")
        print(f"   residual median off-source={np.nanmedian(residual_free):+.3g}")
        print(f"   rms from '{sigma_source}': median={rms_level:.5g}, "
              f"spatial spread={100*rms_structure:.1f}% of the median")

    if show_bkg_map is not None:
        do_plot = bool(show_bkg_map)
    if do_plot:
        _plot_bkg_rms_diagnostic(data, bkg_map, rms_out, source_mask,
                                 box_px=box_px, rms_box_px=rms_box_px,
                                 plot_name=plot_name)

    if not return_diagnostics:
        return bkg_map, rms_out

    diagnostics = {
        'source_mask': source_mask,
        'exclude_mask': exclude,
        'free_fraction': free_fraction,
        'small_frame_warning': small_frame_warning,
        'box_size_px': box_px,
        'rms_box_size_px': rms_box_px,
        'source_extent_px': source_extent_px,
        'exclude_percentile': exclude_percentile,
        'dilation_px': dilation_px,
        'mask_growth': mask_growth,
        'smooth_px': smooth_px,
        'measured_rms': measured_rms,
        'rms_source': sigma_source,
        'rms_structure': rms_structure,
        'residual': data - bkg_map,
        'bkg_excess_under_source': bkg_excess,
        'edge_gradient': edge_gradient,
        'edge_gradient_sigma': edge_sigma,
    }
    return bkg_map, rms_out, diagnostics


def _plot_bkg_rms_diagnostic(data, bkg_map, rms_out, source_mask,
                             box_px=None, rms_box_px=None, plot_name=None):
    """
    Three-panel check on a background/rms pair produced by
    `compute_bkg_rms_maps`.

    Left    background map, with the source mask outlined.
    Middle  a horizontal cut through the background at the source centre, with
            the source extent shaded. This is the panel that matters: a correct
            background is flat across the shaded region, and one that is still
            tracking the source shows a bump there.
    Right   the rms map. If it is featureless, the noise was not measured
            locally and the weights it produces carry no spatial information.
    """
    ny, nx = data.shape
    if np.any(source_mask):
        yc, xc = nd.center_of_mass(source_mask)
        ys, xs = np.where(source_mask)
        x_lo, x_hi = int(xs.min()), int(xs.max())
    else:
        yc, xc = (ny - 1) / 2.0, (nx - 1) / 2.0
        x_lo, x_hi = 0, nx - 1
    row = int(np.clip(round(yc), 0, ny - 1))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))

    lo, hi = np.nanpercentile(bkg_map, [1, 99])
    im0 = axes[0].imshow(bkg_map, origin='lower', cmap='viridis',
                         vmin=lo, vmax=hi)
    axes[0].contour(source_mask.astype(float), levels=[0.5], colors='white',
                    linewidths=1.0)
    axes[0].axhline(row, color='red', ls=':', lw=1.0)
    axes[0].set_title(f'background (box {box_px} px)', fontsize=10)
    fig.colorbar(im0, ax=axes[0], fraction=0.046)

    axes[1].plot(bkg_map[row, :], color='C0', lw=1.5, label='background')
    axes[1].axvspan(x_lo, x_hi, color='gray', alpha=0.2, label='source extent')
    off = ~source_mask
    if np.any(off):
        axes[1].axhline(float(np.nanmedian(bkg_map[off])), color='k', ls='--',
                        lw=1.0, label='off-source median')
    axes[1].set_xlabel('x [px]')
    axes[1].set_ylabel('background')
    axes[1].set_title(f'cut at y={row} (flat is correct)', fontsize=10)
    axes[1].legend(fontsize=8)
    axes[1].grid(True, alpha=0.4, ls=':')

    lo, hi = np.nanpercentile(rms_out, [1, 99])
    im2 = axes[2].imshow(rms_out, origin='lower', cmap='magma',
                         vmin=lo, vmax=hi)
    axes[2].set_title(f'rms (box {rms_box_px} px)', fontsize=10)
    fig.colorbar(im2, ax=axes[2], fraction=0.046)

    fig.tight_layout()
    if plot_name is not None:
        fig.savefig(plot_name, dpi=150, bbox_inches='tight')
    return fig