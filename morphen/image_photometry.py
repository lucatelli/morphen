"""
 ____  _           _                       _              
|  _ \| |__   ___ | |_ ___  _ __ ___   ___| |_ _ __ _   _ 
| |_) | '_ \ / _ \| __/ _ \| '_ ` _ \ / _ \ __| '__| | | |
|  __/| | | | (_) | || (_) | | | | | |  __/ |_| |  | |_| |
|_|   |_| |_|\___/ \__\___/|_| |_| |_|\___|\__|_|   \__, |
                                                    |___/
                                                    
 ____      _                 _
|  _ \ ___| |_ _ __ ___  ___(_) __ _ _ __
| |_) / _ \ __| '__/ _ \/ __| |/ _` | '_ \
|  __/  __/ |_| | | (_) \__ \ | (_| | | | |
|_|   \___|\__|_|  \___/|___/_|\__,_|_| |_|
#Petrosian
"""


def do_petrofit(image, cell_size, mask_component=None, fwhm=8, kernel_size=5, npixels=32,
                main_feature_index=0, sigma_mask=7, dilation_size=10,deblend=False,
                apply_mask=True, PLOT=True, show_figure = True, results=None):
    # from petrofit.photometry import order_cat
    # from petrofit.photometry import make_radius_list

    # from petrofit source_photometry
    # from petrofit import make_catalog, plot_segments
    # from petrofit import plot_segment_residual
    # from petrofit import order_cat

    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(image)

    # sigma = fwhm * gaussian_fwhm_to_sigma
    # kernel = Gaussian2DKernel(sigma, x_size=kernel_size, y_size=kernel_size)
    data_2D_ = load_fits_data(image)
    #mad std must be computed on the original data, not on the masked array.
    std = np.std(data_2D_)

    if mask_component is not None:
        data_2D_ = data_2D_ * mask_component
    if apply_mask == True:
        omask, mask = mask_dilation(image, cell_size=cell_size,
                                    sigma=sigma_mask, dilation_size=dilation_size,
                                    PLOT=False)

        # data_2D = data_2D*(data_2D>=3*mad_std(data_2D))
        data_2D = data_2D_ * mask
        # std = np.std(data_2D)
    else:
        data_2D = data_2D_
        mask = np.ones(data_2D.shape)
        # std = mad_std(data_2D)

    plt.figure()

    cat, segm, segm_deblend = make_catalog(
        data_2D,
        threshold=1.0 * std,
        deblend=deblend,
        # kernel_size=kernel_size,
        # fwhm=fwhm,
        npixels=npixels,
        plot=PLOT, vmax=data_2D.max(), vmin=3 * std
    )

    # Display source properties
    print("Num of Targets:", len(cat))

    # Convert to table
    cat_table = cat.to_table()

    vmax = data_2D.max()
    vmin = 3 * std

    if PLOT == True:
        plt.figure()
        plot_segments(segm, image=data_2D, vmax=vmax, vmin=vmin)
        try:
            plt.figure()
            plot_segment_residual(segm, data_2D, vmax=vmax * 0.01)
            plt.figure()
            plot_segments(segm_deblend, image=data_2D, vmax=vmax, vmin=vmin)
        except:
            pass

    # Sort and get the largest object in the catalog
    sorted_idx_list = order_cat(cat, key='area', reverse=True)
    idx = sorted_idx_list[main_feature_index]  # index 0 is largest
    source = cat[idx]  # get source from the catalog

    #     try:
    #         r,ir = get_profile(image)
    #     except:
    # radii_last = 3*estimate_circular_aperture(crop_image,cellsize=0.008,std=3)
    # radii_last
    r,ir = get_profile(image)
    # radii_last = r[-1]/cell_size
    # radii_last = r[-1]
    # radii_last = 2*results['outer_perimeter']/(2.0*np.pi)
    if mask_component is not None:
        radii_last = int(3 * np.sqrt((np.sum(mask_component) / np.pi)))
    else:
        radii_last = np.sqrt(data_2D.shape[0]**2.0 + data_2D.shape[1]**2.0)/2
        # radii_last = 2 * estimate_area_nbeam(image, mask, cellsize=cell_size)
    print('Creating radii list with max value of =', radii_last)
    #     r_list = np.arange(4, radii_last, 4)
    #     r_list = make_radius_list(
    #         max_pix=radii_last, # Max pixel to go up to
    #         n=int(len(r)) # the number of radii to produce
    #     )

    r_list = make_radius_list(
        max_pix=radii_last,  # Max pixel to go up to
        n=int(radii_last)  # the number of radii to produce
    )

    #     max_pix = 3*estimate_area_nbeam(image,mask,cellsize=cell_size)
    #     # cell_size = 0.008
    #     true_resolution = 0.05
    #     n_points = max_pix * (cell_size/true_resolution)

    #     r_list = make_radius_list(
    #         max_pix=max_pix, # Max pixel to go up to
    #         n=int(n_points) # the number of radii to produce
    #     )
    # print(repr(r_list))

    flux_arr, area_arr, error_arr = source_photometry(
        # Inputs
        source,  # Source (`photutils.segmentation.catalog.SourceCatalog`)
        data_2D,  # Image as 2D array
        #     segm_deblend, # Deblended segmentation map of image
        segm,  # Deblended segmentation map of image
        r_list,  # list of aperture radii
        # Options
        cutout_size=2 * max(r_list),  # Cutout out size, set to double the max radius
        bg_sub=False,  # Subtract background
        sigma=1, sigma_type='clip',  # Fit a 2D plane to pixels within 3 sigma of the mean
        plot=PLOT, vmax=0.3 * data_2D.max(), vmin=3 * std,  # Show plot with max and min defined above
    )

    beam_A = beam_area2(image, cellsize=cell_size)
    S_flux = flux_arr / beam_A
    area_beam = area_arr / beam_A

    from petrofit.petrosian import Petrosian
    p = Petrosian(r_list, area_arr, flux_arr)

    from copy import copy

    p_copy = copy(p)
    p_copy.eta = 0.2
    p_copy.epsilon = 2

    # print('eta =', p_copy.eta)
    # print('epsilon =', p_copy.epsilon)
    # print('r_half_light (old vs new) = {:0.2f} vs {:0.2f}'.format(p.r_half_light, p_copy.r_half_light))
    # print('r_total_flux (old vs new) = {:0.2f} vs {:0.2f}'.format(p.r_total_flux, p_copy.r_total_flux))



    # print('R50 = ', p_copy.r_half_light)
    # print('Rp=', p_copy.r_petrosian)

    results['R50'] = p_copy.r_half_light
    results['Rp'] = p_copy.r_petrosian
    results['R50_2'] = p.r_half_light
    results['Rp_2'] = p.r_petrosian
    results['flux_rp'] = S_flux.max()
    results['r_total_flux_2'] = p.r_total_flux
    results['total_flux_rp_2'] = p.total_flux/beam_A

    results['R20p'] = p_copy.fraction_flux_to_r(fraction=0.2)
    results['R50p'] = p_copy.r_half_light
    results['R80p'] = p_copy.fraction_flux_to_r(fraction=0.8)
    results['R90p'] = p_copy.fraction_flux_to_r(fraction=0.9)
    results['R20p_2'] = p.fraction_flux_to_r(fraction=0.2)
    results['R50p_2'] = p.r_half_light
    results['R80p_2'] = p.fraction_flux_to_r(fraction=0.8)
    results['R90p_2'] = p.fraction_flux_to_r(fraction=0.9)

    C1p = np.log10(results['R80p'] / results['R20p'])
    C2p = np.log10(results['R90p'] / results['R50p'])

    results['C1p'] = C1p
    results['C2p'] = C2p

    results['r_total_flux'] = p_copy.r_total_flux
    results['total_flux_rp'] = p_copy.total_flux/beam_A
    results['r_total_flux_2'] = p.r_total_flux
    results['total_flux_rp_2'] = p.total_flux/beam_A

    if PLOT == True:
        plt.figure()
        #         plt.plot(r_list,S_flux/S_flux[-1])
        #         plt.xlim(0,2*p.r_petrosian)
        #         plt.semilogx()
        p_copy.plot(plot_r=True)
        plt.savefig(image.replace('.fits', '_Rpetro_flux.jpg'), dpi=300, bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()

    return (r_list, area_arr, area_beam, p_copy, S_flux, results)


def petrosian_metrics(source, data_2D, segm, mask_source,global_mask=None,
                 i='1', petro_properties={},sigma_type='clip',eta_value=0.2,
                 rlast=None, sigma=3, vmin=3, bkg_sub=False,error=None,
                 plot=False):
    if rlast is None:
        if mask_source is not None:
            _, area_convex_mask = convex_shape(mask_source)
            # rlast = int(2 * np.sqrt((np.sum(mask_source) / np.pi)))
            rlast = int(1.5*area_to_radii(area_convex_mask))
        else:
            if global_mask is not None:
                _, area_convex_mask = convex_shape(global_mask)
                # rlast = int(2 * np.sqrt((np.sum(global_mask) / np.pi)))
                rlast = int(1.5*area_to_radii(area_convex_mask))
            else:
                rlast = np.sqrt(data_2D.shape[0]**2.0 + data_2D.shape[1]**2.0)/2
    else:
        rlast = rlast

    r_list = make_radius_list(max_pix=rlast,  # Max pixel to go up to
                              n=int(rlast)  # the number of radii to produce
                             )
    cutout_size = 2 * max(r_list)

    # if (bkg_to_sub is not None) and (bkg_sub == False):
    #     data_2D = data_2D - bkg_to_sub

    flux_arr, area_arr, error_arr = source_photometry(source=source,
                                                      image=data_2D,
                                                      segm_deblend=segm,
                                                      r_list=r_list,
                                                      error=error,
                                                      cutout_size=cutout_size,
                                                      bg_sub=bkg_sub, sigma=sigma,
                                                      sigma_type=sigma_type,
                                                      plot=plot, vmax=0.3 * data_2D.max(),
                                                      vmin=vmin * mad_std(data_2D)
                                                      )
    #     fast_plot2(mask_source * data_2D)
    p = Petrosian(r_list, area_arr, flux_arr)
    from copy import copy
    p_015 = copy(p)
    if eta_value is None:
        p_015.eta = 0.2
    else:
        p_015.eta = eta_value

    p_015.epsilon = 2
    R50 = p_015.r_half_light
    Snu = p_015.total_flux
    Rp = p_015.r_petrosian
    try:
        Rpidx = int(2 * Rp)
    except:
        Rpidx = int(r_list[-1])
    petro_properties['R50'] = R50
    petro_properties['Snu'] = Snu
    petro_properties['Rp'] = Rp
    petro_properties['Rpidx'] = Rpidx
    petro_properties['rlast'] = rlast
    return (petro_properties, flux_arr, area_arr, error_arr, p_015, r_list)

def compute_petrosian_properties(data_2D, imagename, mask_component=None,
                                 global_mask=None,
                                 i=0, source_props=None,apply_mask = False,
                                 sigma_level=3, bkg_sub=False,error=None,
                                 vmin=1.0, plot=False, deblend=False,
                                 show_figure = True,plot_catalog=False,
                                 segm_reg= 'mask',vmax=0.1,bkg_to_sub=None,
                                 fwhm=121, kernel_size=81, npixels=128,
                                 verbose=0,
                                 add_save_name='',logger=None):
    # if mask:
    if source_props is None:
        source_props = {}
        source_props['#imagename'] = os.path.basename(imagename)
    if imagename is not None:
        beam_area_ = beam_area2(imagename, cellsize=None)
    if npixels is None:
        npixels = int(beam_area_)

    ii = str(i + 1)
    std = mad_std(data_2D)
    data_component = data_2D.copy()
    # if apply_mask == True:
    #     omask, mask = mask_dilation(image, cell_size=cell_size,
    #                                 sigma=sigma_mask, dilation_size=dilation_size,
    #                                 PLOT=False)
    if (bkg_to_sub is not None) and (bkg_sub == False):
        data_component = data_component  - bkg_to_sub
    if global_mask is not None:
        '''
        Global mask:
            mask that describe the full extension of the source.
        '''
        data_component = data_component * global_mask
        # sigma_level = 1.0
    if mask_component is not None:
        """
        mask_component: 2D bolean array
            If source has multiple components, this mask is the mask of
            only one component, or more, like
            mask_component = mask_component1 + mask_component2.
        """
        # sigma_level = 1.0
        data_component = data_component * mask_component


    # data_component = data_2D
    # eimshow(data_component)

    try:
        cat, segm, segm_deblend = make_catalog(image=data_component,
                                               threshold=3 * std,
                                               deblend=deblend,
                                               # kernel_size=kernel_size,
                                               # fwhm=fwhm,
                                               npixels=npixels,
                                               # because we already deblended it!
                                               plot=plot_catalog,
                                               vmax=vmax*data_component.max(),
                                               vmin=vmin * std)
    except:
        try:
            cat, segm, segm_deblend = make_catalog(image=data_component,
                                                   threshold=0.5 * std,
                                                   deblend=deblend,
                                                   # kernel_size=kernel_size,
                                                   # fwhm=fwhm,
                                                   npixels=npixels,
                                                   # because we already deblended it!
                                                   plot=plot_catalog,
                                                   vmax=vmax*data_component.max(),
                                                   vmin=vmin * std)
        except:
            cat, segm, segm_deblend = make_catalog(image=data_component,
                                                   threshold=0.01 * std,
                                                   deblend=deblend,
                                                   # kernel_size=kernel_size,
                                                   # fwhm=fwhm,
                                                   npixels=npixels,
                                                   # because we already deblended it!
                                                   plot=plot_catalog,
                                                   vmax=vmax * data_component.max(),
                                                   vmin=vmin * std)

    sorted_idx_list = order_cat(cat, key='area', reverse=True)
    idx = sorted_idx_list[0]  # index 0 is largest
    source = cat[idx]



    # if plot == True:
    #     # plt.figure()
    #     plot_segments(segm, image=data_component, vmax=vmax, vmin=vmin)
    #     try:
    #         # plt.figure()
    #         plot_segment_residual(segm, image=data_component, vmax=vmax * 0.01)
    #         # plt.figure()
    #         # plot_segments(data_component, image=data_2D, vmax=vmax, vmin=vmin)
    #     except:
    #         pass

    # source = cat[0]

    source_props['PA'] = source.orientation.value
    source_props['q'] = 1 - source.ellipticity.value
    source_props['area'] = source.area.value
    source_props['Re'] = source.equivalent_radius.value
    source_props['x0c'] = source.xcentroid
    source_props['y0c'] = source.ycentroid

    if segm_reg == 'deblended':
        segm_mask = segm_deblend
    if segm_reg == 'mask':
        segm_mask = segm




    # help function to be used if iteration required.
    source_props, flux_arr, area_arr, error_arr, p, r_list = \
        petrosian_metrics(source=source,
                          data_2D=data_component,
                          segm=segm_mask,global_mask=global_mask,
                          mask_source=mask_component,i=ii,
                          petro_properties=source_props,
                          rlast=None, sigma=sigma_level,
                          vmin=vmin, bkg_sub=bkg_sub,error=error,
                          plot=plot)

    """
    Check if Rp is larger than last element (rlast) of the R_list. If yes,
    we need to run petro_params again, with a larger rlast, at least r_last>=Rp.
    If not, R50 will be np.nan as well Snu.
    """

    if (source_props['rlast'] < 2 * source_props['Rp']) or \
            (np.isnan(p.r_total_flux)):
        if verbose>0:
            print('WARNING: Number of pixels for petro region is to small. '
                  'Looping over until good condition is satisfied.')
            print(f"Rlast        = {source_props['rlast']}")
            print(f"2*Rp           =  {2*source_props['Rp']}")
            print(f"r_total_flux = {p.r_total_flux}")


        if (mask_component is not None):
            _, area_convex_mask = convex_shape(mask_component)
            # rlast = int(2 * np.sqrt((np.sum(mask_source) / np.pi)))
            Rlast_new = int(3.0 * area_to_radii(area_convex_mask))
        else:
            if global_mask is not None:
                _, area_convex_mask = convex_shape(global_mask)
                # rlast = int(2 * np.sqrt((np.sum(mask_source) / np.pi)))
                Rlast_new = int(3.0 * area_to_radii(area_convex_mask))
            else:
                Rlast_new = 2 * source_props['Rp'] + (2 * source_props['Rp']-source_props['rlast']+1)
        # print(f'Rlast (new) = {Rlast_new}')
        source_props, flux_arr, area_arr, error_arr, p, r_list = \
            petrosian_metrics(source=source,data_2D=data_component,
                              segm=segm_mask,global_mask=global_mask,
                              mask_source=mask_component,eta_value=0.20,
                              i=ii,petro_properties=source_props,
                              rlast=Rlast_new,sigma=sigma_level,
                              vmin=vmin,bkg_sub=bkg_sub,error=error,
                              plot=plot)

    if (source_props['rlast'] < 2 * source_props['Rp']) or \
            (np.isnan(p.r_total_flux)):
        if verbose>0:
            print('WARNING: Number of pixels for petro region is to small. '
                  'Looping over until good condition is satisfied.')
            print(f"Rlast        = {source_props['rlast']}")
            print(f"2*Rp         =  {2*source_props['Rp']}")
            print(f"r_total_flux = {p.r_total_flux}")
        Rlast_new = 2 * source_props['Rp'] + (2 * source_props['Rp']-source_props['rlast']+1)
        # print(f'Rlast (new) = {Rlast_new}')
        source_props, flux_arr, area_arr, error_arr, p, r_list = \
            petrosian_metrics(source=source,data_2D=data_component,
                              segm=segm_mask,global_mask=global_mask,
                              mask_source=mask_component,eta_value=0.25,
                              i=ii,petro_properties=source_props,
                              rlast=Rlast_new,sigma=sigma_level,
                              vmin=vmin,bkg_sub=bkg_sub,error=error,
                              plot=plot)


    beam_A = beam_area2(imagename)
    S_flux = flux_arr / beam_A
    area_beam = area_arr / beam_A

    # from copy import copy
    # p_copy = copy(p)
    # p_copy.eta = 0.2
    # p_copy.epsilon = 2

    # print('eta =', p.eta)
    # print('epsilon =', p.epsilon)
    # print('r_half_light (old vs new) = {:0.2f}'.
    #       format(p.r_half_light))
    # print('r_total_flux (old vs new) = {:0.2f}'.
    #       format(p.r_total_flux))


    """
    Now, estimate the effective intensity.
    """
    r, ir = get_profile(data_component, binsize=1.0)
    try:
        I50 = ir[int(source_props['R50'])]
    except:
        source_props['R50'] = source_props['Re'] / 2
        I50 = ir[int(source_props['R50'])]
        # I50 = ir[0]*0.1
    source_props['I50'] = I50

    # print('R50 = ', p.r_half_light)
    # print('Rp=', p.r_petrosian)

    source_props['R50'] = p.r_half_light
    source_props['Rp'] = p.r_petrosian
    # source_props['R50_2'] = p.r_half_light
    # source_props['Rp_2'] = p.r_petrosian
    source_props['flux_rp'] = S_flux.max()
    source_props['r_total_flux'] = p.r_total_flux
    source_props['total_flux_rp'] = p.total_flux/beam_A

    # source_props['R20p'] = p_copy.fraction_flux_to_r(fraction=0.2)
    # source_props['R50p'] = p_copy.r_half_light
    # source_props['R80p'] = p_copy.fraction_flux_to_r(fraction=0.8)
    # source_props['R90p'] = p_copy.fraction_flux_to_r(fraction=0.9)
    source_props['R20p'] = p.fraction_flux_to_r(fraction=0.2)
    source_props['R50p'] = p.r_half_light
    source_props['R80p'] = p.fraction_flux_to_r(fraction=0.8)
    source_props['R90p'] = p.fraction_flux_to_r(fraction=0.9)

    C1p = np.log10(source_props['R80p'] / source_props['R20p'])
    C2p = np.log10(source_props['R90p'] / source_props['R50p'])

    source_props['C1p'] = C1p
    source_props['C2p'] = C2p

    # source_props['r_total_flux'] = p.r_total_flux
    # source_props['total_flux_rp'] = p.total_flux/beam_A
    # source_props['r_total_flux_2'] = p.r_total_flux
    # source_props['total_flux_rp_2'] = p.total_flux/beam_A

    if plot == True:
        plot_flux_petro(imagename, flux_arr, r_list, add_save_name)
        plt.figure()
        #         plt.plot(r_list,S_flux/S_flux[-1])
        #         plt.xlim(0,2*p.r_petrosian)
        #         plt.semilogx()
        p.plot(plot_r=True)
        plt.savefig(imagename.replace('.fits', '_Rpetro_flux'+
                                      add_save_name+'.jpg'),
                    dpi=300, bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()

        plt.figure()
        p.plot_cog()
        plt.savefig(imagename.replace('.fits', '_cog_flux'+
                                      add_save_name+'.jpg'),
                    dpi=300, bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()

    # if (global_mask is not None) and (imagename is not None):
    #     data_comp_mask = data_component * global_mask
    #     total_comp_flux = np.sum(data_comp_mask) / beam_area_
    #     source_props['c' + ii + '_total_flux'] = total_comp_flux
    return(r_list, area_arr, area_beam, p, flux_arr, error_arr, source_props,
           cat, sorted_idx_list, segm, segm_deblend)

def compute_petro_source(data_2D, mask_component=None, global_mask=None,
                         obs_type = 'radio',
                         npixels=None, nlevels=1, contrast=1,
                         imagename=None, i=0, source_props={},positions=None,
                         sigma_level=3, bkg_sub=False,
                         vmin=1, plot=False, deblend=False):
    """
    Perform petrosian photometry of a source. It can be used for the full structure of the source
    or for individual regions (incl. deblended).

    Parameters:
    -----------
    data_2D: 2D numpy array
        Image data.
    mask_component: 2D numpy array
        Mask of a given region (e.g. used when analysing multiple distinc regions of a source.
    global_mask: 2D numpy array
        Mask of the full source.
    imagename: str
        Name of the image. Used when one wants to point the original .fits header (wcs) of the
        image where the data_2D comes from.
    i: int
        Index of the component (when analysing multiple regions).
    source_props: dict
        Dictionary to store the properties of the source.
    positions: list
        List of positions of the sources.
    sigma_level: float
        Level of the sigma to be used for the thresholding.
    bkg_sub: bool
        If True, subtract the background.
    vmin: float
        Minimum value of the image.
    plot: bool
        If True, plot the results.
    deblend: bool
        If True, deblend the sources.
    """
    verbose = 0
    # if mask:
    if imagename is not None:
        try:
            beam_area_ = beam_area2(imagename, cellsize=None)
        except:
            beam_area_ = 1.0

    if isinstance(i, str) == True:
        # print(' ++==>> Computing global source properties.')
        ii = i
    else:
        ii = str(i + 1)
    
    # if isinstance(data_2D, str) == True:
    #     std = mad_std(load_fits_data(data_2D))
    # else:
    if mad_std(data_2D) == 0:
        std = mad_std(data_2D[(data_2D>0)])
    else:
        std = mad_std(data_2D)
    
    if mask_component is not None:
        data_component = data_2D * mask_component
        if npixels is None:
            npixels = int(np.nansum(mask_component)/50)
    else:
        mask_component = np.ones(data_2D.shape)
        data_component = data_2D
        if npixels is None:
            npixels = int(data_2D.shape[0]/50)
    
    # print(' ++==>> Making catalog...')
    cat, segm, segm_deblend = make_catalog(image=data_component,
                                           threshold=sigma_level * std,
                                           deblend=False,# because we already deblended it!
                                        #    npixels=npixels,
                                           npixels=npixels, nlevels=nlevels, contrast=contrast,
                                           plot=plot, vmax=np.nanmax(data_component),
                                           vmin=vmin * std)

    sorted_idx_list = order_cat(cat, key='area', reverse=True)
    # print(len(sorted_idx_list))
    source = cat[sorted_idx_list[0]]

    do_fit_ellipse = True
    if do_fit_ellipse:
        levels_ellipse = np.geomspace(np.nanmax(data_component), sigma_level * std, 32)
        try:
            sigma_50 = levels_ellipse[levels_ellipse < 0.2*np.nanmax(levels_ellipse)][0]
        except:
            sigma_50 = levels_ellipse[levels_ellipse < 0.8*np.nanmax(levels_ellipse)][0]
        region_split = [i for i, x in enumerate(levels_ellipse > sigma_50) if x][-1]
        
        PA, q, x0col, y0col, PAm, qm, \
            PAmi, qmi, PAmo, qmo, \
            x0median, y0median, \
            x0median_i, y0median_i, \
            x0median_o, y0median_o,profiles = cal_PA_q(data_component, 
                                                       Isequence=levels_ellipse,
                                                       region_split=region_split,
                                                       SAVENAME=f'{imagename}_{i}'.replace('.fits','_ellipsefit')
                                                    #    SAVENAME=img.replace('.fits','_ellipsefit') + ext
                                            )

    source_props['c' + ii + '_PA'] = source.orientation.value
    source_props['c' + ii + '_q'] = 1 - source.ellipticity.value
    source_props['c' + ii + '_area'] = source.area.value
    source_props['c' + ii + '_Re'] = source.equivalent_radius.value
    source_props['c' + ii + '_x0c'] = source.xcentroid
    source_props['c' + ii + '_y0c'] = source.ycentroid
    source_props['c' + ii + '_label'] = source.label

    # help function to be used if iteration required.
    # print(' ++==>> Computing petrosian parameters...')
    source_props, p = petro_params(source=source, data_2D=data_component, segm=segm,
                                mask_source=mask_component,positions=positions,
                                i=ii, petro_properties=source_props,
                                rlast=None, sigma=sigma_level,
                                vmin=vmin, bkg_sub=bkg_sub,
                                plot=plot)
    
    # print('Rlast', source_props['c' + ii + '_rlast'])
    # print('2Rp', 2 * source_props['c' + ii + '_Rp'])
    """
    Check if Rp is larger than last element (rlast) of the R_list. If yes,
    we need to run petro_params again, with a larger rlast, at least r_last>=Rp.
    If not, R50 will be np.nan.
    """
    if (source_props['c' + ii + '_rlast'] < 2 * source_props['c' + ii + '_Rp']) \
            or (np.isnan(p.r_total_flux)):
        if verbose>0:
            print('WARNING: Number of pixels for petro region is to small. '
                  'Looping over until good condition is satisfied.')
        # Rlast_new = 2 * source_props['c' + ii + '_Rp'] + 3

        _, area_convex_mask = convex_shape(mask_component)
        # rlast = int(2 * np.sqrt((np.sum(mask_source) / np.pi)))
        Rlast_new = int(3.0 * area_to_radii(area_convex_mask))
        # print('Rlast', source_props['c' + ii + '_rlast'])
        # print('2Rp', 2 * source_props['c' + ii + '_Rp'])
        # print('Rlast_new',Rlast_new)
        # print(' ++==>> Re-computing petrosian parameters...')
        source_props, p = petro_params(source=source, data_2D=data_component,
                                    segm=segm, mask_source=mask_component,
                                    i=ii, petro_properties=source_props,
                                    rlast=Rlast_new, sigma=sigma_level,
                                    vmin=vmin,
                                    bkg_sub=bkg_sub, plot=plot)

        if (source_props['c' + ii + '_rlast'] < 2 * source_props['c' + ii + '_Rp']) \
                or (np.isnan(p.r_total_flux)):
            if verbose>0:
                print('WARNING: Number of pixels for petro region is to small. '
                    'Looping over until good condition is satisfied.')
            # Rlast_new = 2 * source_props['Rp'] + 3
            # print(' ++==>> Re-computing petrosian parameters...')

            Rlast_new = 2 * source_props['c' + ii + '_Rp'] + 3
            # print('Rlast', source_props['c' + ii + '_rlast'])
            # print('2Rp', 2 * source_props['c' + ii + '_Rp'])
            # print('Rlast_new',Rlast_new)
            
            source_props, p = petro_params(source=source, data_2D=data_2D,
                                        segm=segm,mask_source=mask_component,
                                        i=ii, petro_properties=source_props,
                                        rlast=Rlast_new, sigma=sigma_level,
                                        vmin=vmin,eta_value=0.25,
                                        bkg_sub=bkg_sub, plot=plot)



    """
    Now, estimate the effective intensity.
    """
    r, ir = get_profile(data_component, binsize=1.0)
    # if obs_type == 'radio':
    #     I50 = np.nanmax(ir)
    # else:
    try:
        I50 = ir[int(source_props['c' + ii + '_R50'])]
    except:
        source_props['c' + ii + '_R50'] = source_props['c' + ii + '_Re'] / 2
        I50 = ir[int(source_props['c' + ii + '_R50'])]
        # I50 = ir[0]*0.1

    source_props['c' + ii + '_I50'] = I50

    if (global_mask is not None) and (imagename is not None):
        data_comp_mask = data_component * global_mask
        total_comp_flux = np.sum(data_comp_mask) / beam_area_
        source_props['c' + ii + '_total_flux'] = total_comp_flux
    return (source_props)

def petro_cat(data_2D, fwhm=24, npixels=128, kernel_size=15,
              nlevels=30, contrast=0.001,bkg_sub=False,
              sigma_level=20, vmin=5,
              deblend=True, plot=False):
    """
    Use PetroFit class to create catalogues.
    """
    cat, segm, segm_deblend = make_catalog(
        image=data_2D,
        threshold=sigma_level * mad_std(data_2D),
        # kernel_size=kernel_size,
        # fwhm=fwhm,
        nlevels=nlevels,
        deblend=deblend,
        npixels=npixels,contrast=contrast,
        plot=plot, vmax=data_2D.max(), vmin=vmin * mad_std(data_2D)
    )

    sorted_idx_list = order_cat(cat, key='area', reverse=True)
    #     idx = sorted_idx_list[main_feature_index]  # index 0 is largest
    #     source = cat[idx]  # get source from the catalog
    return (cat, segm, sorted_idx_list)


def petro_params(source, data_2D, segm, mask_source, positions=None,
                 i='1', petro_properties={},sigma_type='clip',eta_value=None,
                 rlast=None, sigma=3, vmin=3, bkg_sub=True, plot=False):
    if rlast is None:
        rlast = int(np.sqrt((np.sum(mask_source) / np.pi)))
        # _, area_convex_mask = convex_shape(mask_source)
        # Rlast_convex = int(1.0 * area_to_radii(area_convex_mask))        
        # print('Estimate for rlast =', rlast)
        # print('Estimate for rlast (convex) =', Rlast_convex)

    else:
        rlast = rlast
        

    r_list = make_radius_list(max_pix=rlast,  # Max pixel to go up to
                              n=int(rlast)  # the number of radii to produce
                             )
    cutout_size = 2 * max(r_list)
    flux_arr, area_arr, error_arr = source_photometry(source, data_2D, segm,
                                                      r_list, cutout_size=cutout_size,
                                                      # position2=positions,
                                                      bg_sub=bkg_sub, sigma=sigma,
                                                      sigma_type=sigma_type,
                                                      plot=plot, vmax=0.3 * data_2D.max(),
                                                      vmin=vmin * mad_std(data_2D)
                                                      )
    #     fast_plot2(mask_source * data_2D)
    p = Petrosian(r_list, area_arr, flux_arr)

    if eta_value is None:
        R50 = p.r_half_light
        R20 = p.fraction_flux_to_r(fraction=0.2)
        R80 = p.fraction_flux_to_r(fraction=0.8)
        C1p = np.log10(R80 / R20)
        Snu = p.total_flux
        Rp = p.r_petrosian
        p_return = p
        if plot == True:
            plt.figure()
            p.plot(plot_r=plot)
            plt.figure()
            p.plot_cog()
        #     print('    R50 =', R50)
        #     print('     Rp =', Rp)


    if eta_value is not None:
        from copy import copy
        p_new = copy(p)
        p_new.eta = eta_value
        R50 = p_new.r_half_light
        R20 = p_new.fraction_flux_to_r(fraction=0.2)
        R80 = p_new.fraction_flux_to_r(fraction=0.8)
        C1p = np.log10(R80 / R20)
        Snu = p_new.total_flux
        Rp = p_new.r_petrosian
        p_return = p_new
        if plot == True:
            plt.figure()
            p_new.plot(plot_r=plot)
            p_new.plot_cog()
        #     print('    R50 =', R50)
        #     print('     Rp =', Rp)

    try:
        Rpidx = int(2 * Rp)
    except:
        Rpidx = int(r_list[-1])
    petro_properties['c' + i + '_R50'] = R50
    petro_properties['c' + i + '_R20'] = R20
    petro_properties['c' + i + '_R80'] = R80
    petro_properties['c' + i + '_C1'] = C1p
    petro_properties['c' + i + '_Snu'] = Snu
    petro_properties['c' + i + '_Rp'] = Rp
    petro_properties['c' + i + '_Rpidx'] = Rpidx
    petro_properties['c' + i + '_rlast'] = rlast

    return (petro_properties, p_return)


def source_props(data_2D, source_props={},sigma_mask = 5,
                 fwhm=24, npixels=128, kernel_size=15, nlevels=30,
                 contrast=0.001,sigma_level=20, vmin=5,bkg_sub=False,
                 deblend=True,PLOT=False,apply_mask=False):
    '''
    From a 2D image array, perform simple source extraction, and calculate basic petrosian
    properties.
    '''
    verbose = 0
    if apply_mask:
        _, mask = mask_dilation(data_2D, PLOT=False,
                                sigma=sigma_mask, iterations=2, dilation_size=10)
        data_2D = data_2D*mask

    cat, segm, sorted_idx_list = petro_cat(data_2D, fwhm=fwhm, npixels=npixels,
                                           kernel_size=kernel_size,bkg_sub=bkg_sub,
                                           nlevels=nlevels, contrast=contrast,
                                           sigma_level=sigma_level, vmin=vmin,
                                           deblend=deblend, plot=PLOT)
    #     i = 0
    for i in range(len(sorted_idx_list)):
        ii = str(i + 1)
        seg_image = cat[sorted_idx_list[i]]._segment_img.data
        # seg_image = np.logical_not(cat[sorted_idx_list[i]].segment_ma.mask)
        source = cat[sorted_idx_list[i]]
        source_props['c' + ii + '_PA'] = source.orientation.value
        source_props['c' + ii + '_q'] = 1 - source.ellipticity.value
        source_props['c' + ii + '_area'] = source.area.value
        source_props['c' + ii + '_Re'] = source.equivalent_radius.value
        source_props['c' + ii + '_x0c'] = source.xcentroid
        source_props['c' + ii + '_y0c'] = source.ycentroid
        source_props['c' + ii + '_label'] = source.label

        label_source = source.label
        # plt.imshow(seg_image==label_source)
        mask_source = seg_image == label_source
        # mask_source = seg_image
        source_props, p = petro_params(source=source, data_2D=data_2D, segm=segm,
                                    mask_source=mask_source,
                                    i=ii, petro_properties=source_props,
                                    rlast=None, sigma=sigma_level,
                                    vmin=vmin, bkg_sub=bkg_sub,
                                    plot=PLOT)

        #         print(Rp_props['rlast'],2*Rp_props['Rp'])
        if ((source_props['c' + ii + '_rlast'] < 2 * source_props['c' + ii + '_Rp'])) \
                or (np.isnan(p.r_total_flux)):
            if verbose>0:
                print('WARNING: Number of pixels for petro region is to small. '
                    'Looping over until good condition is satisfied.')
            # Rlast_new = 2 * source_props['c' + ii + '_Rp'] + 3

            _, area_convex_mask = convex_shape(mask_source)
            # rlast = int(2 * np.sqrt((np.sum(mask_source) / np.pi)))
            Rlast_new = int(3.0 * area_to_radii(area_convex_mask))

            source_props, p = petro_params(source=source, data_2D=data_2D,
                                        segm=segm,mask_source=mask_source,
                                        i=ii, petro_properties=source_props,
                                        rlast=Rlast_new, sigma=sigma_level,
                                        vmin=vmin,
                                        bkg_sub=bkg_sub, plot=PLOT)

        if (source_props['c' + ii + '_rlast'] < 2 * source_props['c' + ii + '_Rp']) \
                or (np.isnan(p.r_total_flux)):
            if verbose>0:
                print('WARNING: Number of pixels for petro region is to small. '
                    'Looping over until good condition is satisfied.')
            # Rlast_new = 2 * source_props['Rp'] + 3
            Rlast_new = 2 * source_props['c' + ii + '_Rp'] + 3
            source_props, p = petro_params(source=source, data_2D=data_2D,
                                        segm=segm,mask_source=mask_source,
                                        i=ii, petro_properties=source_props,
                                        rlast=Rlast_new, sigma=sigma_level,
                                        vmin=vmin,eta_value=0.25,
                                        bkg_sub=bkg_sub, plot=PLOT)

        r, ir = get_profile(data_2D * mask_source, binsize=1.0)
        try:
            I50 = ir[int(source_props['c' + ii + '_R50'])]
        except:
            source_props['c' + ii + '_R50'] = source_props['c' + ii + '_Re']/2
            I50 = ir[int(source_props['c' + ii + '_R50'])]
            # I50 = ir[0]*0.1

        source_props['c' + ii + '_I50'] = I50

    source_props['ncomps'] = len(sorted_idx_list)
    return (source_props,cat, segm)