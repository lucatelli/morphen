
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
        plt.title(f"max(bkg)/max(data)={(bkg_image.max()/data_2D.max()):.6f}")
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
        fig, ax = plt.subplots(figsize=(10, 10))
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       vmin=s, vmax=0.2*np.nanmax(data_sub), origin='lower')

    masks_regions = []

    if ell_size_factor is None:
        if mask is not None:
            # ell_size_factor = np.sqrt(np.sum(mask) / (np.pi))/cat[0].equivalent_radius.value
            ell_size_factor = 0.1*np.sqrt(np.sum(mask) / (np.pi))
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
            label_x = xc + 3 * ell_size_factor
            label_y = yc + 20 * ell_size_factor
            line_end_y = label_y - 5 
            
            text = Text(label_x, label_y, label, ha='center', va='center', color='red')
            ax.add_artist(text)
            ax.plot([xc, label_x], [yc, line_end_y], color='red', alpha=0.4)

        plt.axis('off')
        # plt.show()
        plt.savefig(imagename + '_SEP.jpg', dpi=300, bbox_inches='tight')
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
#         plt.savefig(imagename + '_SEP.jpg', dpi=300, bbox_inches='tight')
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
                    filter_type='matched', 
                    sort_by='distance', threshold_mode='sigma',
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
    if apply_mask:
        _, mask = mask_dilation(data_2D, sigma=sigma_mask, iterations=iterations,
                                rms=s,
                                PLOT=True,show_figure=True,
                                dilation_size=dilation_size)


    bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
    # print(bkg.globalback)
    # print(bkg.globalrms)

    bkg_image = bkg.back()
    bkg_rms = bkg.rms()
    
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
        
        fig = plt.figure(figsize=(18, 6))
        ax0 = fig.add_subplot(1, 3, 2)
        ax0 = eimshow(data_sub,ax=ax0,fig=fig,
                            plot_colorbar=True,
                            # cbar_orientation = 'horizontal',
                            show_axis='off',
                            vmin_factor=0.5,
                            cbar_n_points = 3,
                            add_contours=True)
        ax0.set_title(f'data - bkg')
        
        ax1 = fig.add_subplot(1, 3, 1)
        ax1 = eimshow(data_2D,ax=ax1,fig=fig,
                            plot_colorbar=True,
                            # cbar_orientation = 'horizontal',
                            show_axis='off',
                            vmin_factor=0.5,
                            cbar_n_points=3,
                            plot_title = f'data',
                            add_contours=True)
        ax2 = fig.add_subplot(1, 3, 3)
        ax2 = eimshow(bkg_image,ax=ax2,fig=fig,
                            plot_colorbar=True,
                            # cbar_orientation = 'horizontal',
                            show_axis='off',
                            # vmin_factor=0.5,
                            # vmax_factor=1.0,
                            cbar_n_points=3,
                            plot_title = f'bkg',
                            add_contours=False)
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
                                       figsize=(20, 20),
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
        
    masks_deblended = []
    for k in range(len(indices)):
        # print(k)
        masks_deblended.append(seg_maps == seg_maps.labels[indices[k]])


    # len(objects)
    from matplotlib.patches import Ellipse
    from skimage.draw import ellipse

    # m, s = np.mean(data_sub), np.std(data_sub)
    if show_detection == True:
        fig, ax = plt.subplots(figsize=(10, 10))
        norm = simple_norm(data_sub, stretch='sqrt', asinh_a=0.02, min_cut=s,
                    max_cut=0.2*np.nanmax(data_sub))
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       norm=norm, origin='lower')
        # im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
        #                vmin=s, vmax=0.2*np.nanmax(data_sub), origin='lower')

    masks_regions = []

    if ell_size_factor is None:
        if mask is not None:
            # ell_size_factor = np.sqrt(np.sum(mask) / (np.pi))/cat[0].equivalent_radius.value
            ell_size_factor = 0.05*np.sqrt(np.sum(mask) / (np.pi))
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
                    # angle=source.orientation.value * 180. / np.pi
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
            e.set_edgecolor('red')
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
            
            text = Text(label_x, label_y, label, ha='center', va='center', color='red')
            ax.add_artist(text)
            
            ax.plot([xc, label_x], [yc, line_end_y], color='red', alpha=0.4)

        plt.axis('off')
        plt.savefig(imagename + '_SEP_phot.jpg', dpi=300, bbox_inches='tight')
        plt.show()



    if segmentation_map == True:
        return (masks_deblended, sorted_indices_desc, bkg,
                seg_maps, objects_sorted)
    else:
        return (masks_deblended, sorted_indices_desc, bkg,
                objects_sorted)