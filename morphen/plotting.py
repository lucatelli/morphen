"""
 ____  _       _   _   _
|  _ \| | ___ | |_| |_(_)_ __   __ _
| |_) | |/ _ \| __| __| | '_ \ / _` |
|  __/| | (_) | |_| |_| | | | | (_| |
|_|   |_|\___/ \__|\__|_|_| |_|\__, |
                               |___/
Plotting Functions
"""

class CustomFormatter(mticker.ScalarFormatter):
    def __init__(self, factor=1, **kwargs):
        self.factor = factor
        mticker.ScalarFormatter.__init__(self, **kwargs)

    def __call__(self, x, pos=None):
        x = x * self.factor
        if x == 0:
            return "0.00"
        return "{:.2f}".format(x)


def make_scalebar(ax, left_side, length, color='w', linestyle='-', label='',
                  fontsize=12, text_offset=0.1*u.arcsec):
    axlims = ax.axis()
    lines = ax.plot(u.Quantity([left_side.ra, left_side.ra-length]),
                    u.Quantity([left_side.dec]*2),
                    color=color, linestyle=linestyle, marker=None,
                    transform=ax.get_transform('fk5'),
                   )
    txt = ax.text((left_side.ra-length/2).to(u.deg).value,
                  (left_side.dec+text_offset).to(u.deg).value,
                  label,
                  verticalalignment='bottom',
                  horizontalalignment='center',
                  transform=ax.get_transform('icrs'),
                  color=color,
                  fontsize=fontsize,
                 )
    ax.axis(axlims)
    return lines,txt


def plot_radial_profile(imagedatas, refimage=None,
                        ax=None, centre=None, labels=None,
                        line_styles=None,
                        figsize=(5, 5)):
    if ax is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(1, 1, 1)
    else:
        pass

    if labels is None:
        _labels = [''] * len(imagedatas)
    else:
        _labels = labels

    if refimage != None:
        cell_size = get_cell_size(refimage)
        xaxis_units = '[arcsec]'
    else:
        cell_size = 1.0
        xaxis_units = '[px]'

    for i in range(len(imagedatas)):
        radius, intensity = get_profile(imagedatas[i], center=centre)
        try:
            ax.plot(radius * cell_size, intensity, 
                    label=_labels[i],
                    color=line_styles['color'][i],
                    linestyle=line_styles['linestyle'][i],
                    linewidth=line_styles['linewidth'][i],
                    )
        except:
            ax.plot(radius * cell_size, abs(intensity), label=_labels[i])
        

#     plt.semilogy()
    plt.xlabel(f"Radius {xaxis_units}")
    plt.ylabel(f"Radial Intensity [mJy/beam]")
    # plt.xlim(0,cell_size*radius[int(len(radius)/2)])
    # plt.semilogx()
    if labels != None:
        plt.legend()
    return (ax)


def azimuthal_average_profile_with_shaded_errors(image, rms_image, center, sigma=1, bin_size=1.0,
                                                log_scale=True,cell_size=1.0,
                                                figsize=(5, 5),
                                                ylabel='Azimuthal Average Intensity $I(R)$ [Jy/Beam]',
                                                xlabel='Radius [pixels]',
                                                title=None
                                                ):
    """
    Calculate the azimuthal average profile of an image and plot the errors using shaded regions.
    
    Parameters:
    - image: 2D numpy array, the intensity image.
    - rms_image: 2D numpy array, the rms (error) image.
    - center: tuple (x, y), the coordinates of the center (in pixels).
    - sigma: float, the number of standard deviations for the shaded error regions.
    - bin_size: float, size of the radial bins for averaging.

    Returns:
    - radii: numpy array, the radial distances.
    - profile: numpy array, the azimuthal average of intensity at each radius.
    - error_rms: numpy array, the error at each radius based on the rms image.
    - error_std: numpy array, the error at each radius based on the standard deviation.
    """
    
    y, x = np.indices(image.shape)
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)

    # Define bin edges and centers
    r_max = np.max(r)
    bin_edges = np.arange(0, r_max + bin_size, bin_size)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    # Create arrays to hold the azimuthal average and errors
    profile = np.zeros_like(bin_centers)
    error_rms = np.zeros_like(bin_centers)
    error_std = np.zeros_like(bin_centers)

    for i in range(len(bin_centers)):
        mask = (r >= bin_edges[i]) & (r < bin_edges[i + 1])

        if np.any(mask):
            profile[i] = np.nanmean(image[mask])
            error_rms[i] = np.sqrt(np.nanmean(rms_image[mask]**2)) # Use RMS image for errors
            error_std[i] = np.nanstd(image[mask])                  # Use std deviation of pixel intensities        


    nans_mask = ~np.isnan(profile)
    bin_centers = bin_centers[nans_mask]
    profile = profile[nans_mask]
    error_rms = error_rms[nans_mask]
    error_std = error_std[nans_mask]
    
    
    # Plotting the results with shaded error regions
    plt.figure(figsize=figsize)
    
    # Plot the profile
    plt.plot(bin_centers*cell_size, profile, 
            #  label='Profile', 
             color='black', linestyle='-', linewidth=2)

    # Plot shaded error regions
    plt.fill_between(bin_centers*cell_size, profile - sigma * error_rms, profile + sigma * error_rms, 
                     color='gray', alpha=0.3, 
                    #  label=f'RMS Errors ($\pm{sigma}\sigma$)'
                    )
    

    if log_scale:
        plt.yscale('log')
    plt.xlabel(fr'{xlabel}', fontsize=14)
    plt.ylabel(fr'{ylabel}', fontsize=14)
    plt.title(title)
    # plt.legend(loc='upper right')
    # plt.grid(True, which='both', linestyle='--', alpha=0.6)
    plt.tight_layout()
    # plt.show()

    return bin_centers, profile, error_rms, error_std


def plot_azimuthal_profile(image, rms_image, center, sigma=3, bin_size=1.0,
                           log_scale=True, cell_size=1.0,
                           figsize=(5, 5),
                           ylabel='Azimuthal Average Intensity $I(R)$ [Jy/Beam]',
                           xlabel='Radius [pixels]',
                           title=None,
                           which_error='rms',
                           min_points_per_bin=5,  # Minimum points required per bin
                           weight_by_points=True  # Weight averages by number of points
                           ):
    """
    Calculate the azimuthal average profile with improved handling of asymmetric masks.
    
    Additional Parameters:
    - min_points_per_bin: int, minimum number of valid points required in a bin
    - weight_by_points: bool, whether to weight the errors by number of points
    
    Returns:
    - radii: numpy array, the radial distances
    - profile: numpy array, the azimuthal average of intensity
    - error_rms: numpy array, the weighted RMS errors
    - error_std: numpy array, the weighted standard deviation
    - points_per_bin: numpy array, number of valid points in each bin
    """
    
    y, x = np.indices(image.shape)
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)

    r_max = np.nanmax(r)
    bin_edges = np.arange(0, r_max + bin_size, bin_size)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

    profile = np.zeros_like(bin_centers)
    error_rms = np.zeros_like(bin_centers)
    error_std = np.zeros_like(bin_centers)
    points_per_bin = np.zeros_like(bin_centers)

    for i in range(len(bin_centers)):
        mask = (r >= bin_edges[i]) & (r < bin_edges[i + 1])
        valid_pixels = ~np.isnan(image[mask])
        
        if np.nansum(valid_pixels) >= min_points_per_bin:
            points_per_bin[i] = np.nansum(valid_pixels)
            values = image[mask][valid_pixels]
            rms_values = rms_image[mask][valid_pixels]
            
            profile[i] = np.nanmean(values)
            
            # Weight errors by sqrt(N) if requested
            if weight_by_points:
                n_points = len(values)
                error_rms[i] = np.sqrt(np.nanmean(rms_values**2)) / np.sqrt(n_points)
                error_std[i] = np.std(values) / np.sqrt(n_points)
            else:
                error_rms[i] = np.sqrt(np.nanmean(rms_values**2))
                error_std[i] = np.nanstd(values)
        else:
            profile[i] = np.nan
            error_rms[i] = np.nan
            error_std[i] = np.nan

    # Remove bins with insufficient points
    valid_bins = ~np.isnan(profile)
    bin_centers = bin_centers[valid_bins]
    profile = profile[valid_bins]
    error_rms = error_rms[valid_bins]
    error_std = error_std[valid_bins]
    points_per_bin = points_per_bin[valid_bins]

    # Plotting
    plt.figure(figsize=figsize)
    
    # Plot the profile
    plt.plot(bin_centers*cell_size, profile, color='black', 
             marker='.',linestyle='-.',
             linewidth=1)

    if weight_by_points:
        sigma = 3
    else:
        sigma = 3
    # Plot shaded error regions
    
    if which_error == 'std':
        plt.fill_between(bin_centers*cell_size, 
                        profile - sigma * error_std, 
                        profile + sigma * error_std,
                        color='gray', alpha=0.3)
    if which_error == 'rms':
        plt.fill_between(bin_centers*cell_size, 
                        profile - sigma * error_rms, 
                        profile + sigma * error_rms,
                        color='gray', alpha=0.3)
    if which_error == 'both':
        plt.fill_between(bin_centers*cell_size, 
                        profile - sigma * error_rms, 
                        profile + sigma * error_rms,
                        label = 'rms error',
                        color='gray', alpha=0.5)
        plt.fill_between(bin_centers*cell_size,
                        profile - sigma * error_std, 
                        profile + sigma * error_std,
                        label = 'std error',
                        color='orange', alpha = 0.3)
        plt.legend()

    if log_scale:
        plt.yscale('log')
    plt.xlabel(fr'{xlabel}', fontsize=14)
    plt.ylabel(fr'{ylabel}', fontsize=14)
    plt.title(title)
    plt.tight_layout()

    return bin_centers, profile, error_rms, error_std, points_per_bin


def fast_plot2(imagename, crop=False, box_size=128, center=None, with_wcs=True,vmax_factor=0.5,
               vmin_factor=1, plot_colorbar=True, figsize=(5, 5), aspect=1, ax=None):
    """
    Fast plotting of an astronomical image with/or without a wcs header.

    imagename:
        str or 2d array.
        If str (the image file name), it will attempt to read the wcs and plot the coordinates axes.

        If 2darray, will plot the data with generic axes.

        support functions:
            load_fits_data() -> casa to numpy: A function designed mainly to read CASA fits images,
                     but can be used to open any fits images.

                     However, it does not read header/wcs.
                     Note: THis function only works inside CASA environment.




    """
    if ax is None:
        fig = plt.figure(figsize=figsize)
    #     ax = fig.add_subplot(1,1,1)
    #     try:
    if isinstance(imagename, str) == True:
        if with_wcs == True:
            hdu = pf.open(imagename)
            #         hdu=pf.open(img)
            ww = WCS(hdu[0].header, naxis=2)
            try:
                if len(np.shape(hdu[0].data) == 2):
                    g = hdu[0].data[0][0]
                else:
                    g = hdu[0].data
            except:
                g = load_fits_data(imagename)
        if with_wcs == False:
            g = load_fits_data(imagename)

        if crop == True:
            xin, xen, yin, yen = do_cutout(imagename, box_size=box_size, center=center, return_='box')
            g = g[xin:xen, yin:yen]

    else:
        g = imagename

    if crop == True:
        max_x, max_y = np.where(g == g.max())
        xin = max_x[0] - box_size
        xen = max_x[0] + box_size
        yin = max_y[0] - box_size
        yen = max_y[0] + box_size
        g = g[xin:xen, yin:yen]

    if mad_std(g) == 0:
        std = g.std()
    else:
        std = mad_std(g)
    if ax is None:
        if with_wcs == True and isinstance(imagename, str) == True:
            ax = fig.add_subplot(projection=ww.celestial)
            ax.set_xlabel('RA')
            ax.set_ylabel('DEC')
        else:
            ax = fig.add_subplot()
            ax.set_xlabel('x pix')
            ax.set_ylabel('y pix')

    vmin = vmin_factor * std

    #     print(g)
    vmax = vmax_factor * g.max()

    norm = simple_norm(g, stretch='sqrt', asinh_a=0.02, min_cut=vmin, max_cut=vmax)

    im_plot = ax.imshow((g), cmap='magma_r', origin='lower', alpha=1.0, norm=norm,
                        aspect=aspect)  # ,vmax=vmax, vmin=vmin)#norm=norm
    #     ax.set_title('Image')
    try:
        levels_g = np.geomspace(5.0 * g.max(), 0.1 * g.max(), 7)
        #     x = np.geomspace(1.5*mad_std(g),10*mad_std(g),4)
        levels_black = np.geomspace(3 * (mad_std(g) + 0.00001), 0.1 * g.max(), 7)
    except:
        try:
            levels_g = np.geomspace(5.0 * g.max(), 3 * (mad_std(g), 7))
            levels_black = np.asarray([0])
        except:
            levels_g = np.asarray([0])
            levels_black = np.asarray([0])
    #     xneg = np.geomspace(5*mad_std(g),vmin_factor*mad_std(g),2)
    #     y = -xneg[::-1]
    #     levels_black = np.append(y,x)

    #     levels_white = np.geomspace(g.max(),10*mad_std(g),7)
    # levels_white = np.geomspace(g.max(), 0.1 * g.max(), 5)

    #     cg.show_contour(data, colors='black',levels=levels_black,linestyle='.-',linewidths=0.2,alpha=1.0)
    #     cg.show_contour(data, colors='#009E73',levels=levels_white[::-1],linewidths=0.2,alpha=1.0)
    try:
        ax.contour(g, levels=levels_black, colors='black', linewidths=0.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)
        #     ax.contour(g, levels=levels_white[::-1],colors='#009E73',linewidths=0.2,alpha=1.0)#cmap='Reds', linewidths=0.75)
        ax.contour(g, levels=levels_g[::-1], colors='white', linewidths=0.6, alpha=1.0)  # cmap='Reds', linewidths=0.75)
    except:
        print('Not plotting contours!')
    #     cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.0,0.38,0.02,0.23]))#,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')
    try:
        cb = plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([0.91, 0.08, 0.05, 0.84]))
        cb.set_label(r"Flux [Jy/Beam]")
    except:
        pass
    # if ax==None:
    #     if plot_colorbar==True:
    #         # cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([0.91,0.08,0.05,0.82]))
    #         cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([0.91,0.08,0.05,0.32]))
    #         # cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.0,0.38,0.02,0.23]))
    #         cb.set_label(r"Flux [Jy/Beam]")
    return (ax)


def plot_flux_petro(imagename, flux_arr, r_list,
                    savefig=True, show_figure=True,
                    add_save_name = ''):
    plt.figure(figsize=(4, 3))
    cell_size = get_cell_size(imagename)
    plt.plot(r_list * cell_size, 1000 * flux_arr / beam_area2(imagename),
             color='black', lw='3')
    idx_lim = int(np.where(flux_arr / np.max(flux_arr) > 0.95)[0][0] * 1.5)
    plt.grid()
    plt.xlabel('Aperture Radii [arcsec]')
    plt.ylabel(r'$S_\nu$ [mJy]')
    plt.title('Curve of Growth')
    try:
        plt.xlim(0, r_list[idx_lim] * cell_size)
    except:
        plt.xlim(0, r_list[-1] * cell_size)
    if savefig is True:
        plt.savefig(
            imagename.replace('.fits', '_flux_aperture_'+add_save_name+'.jpg'),
            dpi=300, bbox_inches='tight')
    if show_figure == True:
        plt.show()
    else:
        plt.close()

def make_cl(image):
    std = mad_std(image)
    levels = np.geomspace(image.max() * 5, 7 * std, 10)
    return (levels[::-1])


def plot_slices_fig(data_2D, show_figure=True, label='',color=None,FIG=None,linestyle='--.'):
    plot_slice = np.arange(0, data_2D.shape[0])

    if FIG is None:
        fig = plt.figure(figsize=(6, 6))
        ax1 = fig.add_subplot(2, 1, 1)
        ax2 = fig.add_subplot(2, 1, 2)
    else:
        fig,ax1,ax2 = FIG
    ax1.plot(plot_slice, np.mean(data_2D, axis=0), linestyle=linestyle, color=color, ms=14,
             label=label)
    ax1.legend(fontsize=11)
    ax1.grid()
    ax1.set_ylabel('mean $x$ direction')
    ax1.set_xlim(data_2D.shape[0] / 2 - 0.25 * data_2D.shape[0],
                 data_2D.shape[0] / 2 + 0.25 * data_2D.shape[0])

    ax2.plot(plot_slice, np.mean(data_2D, axis=1), linestyle=linestyle, color=color, ms=14,
             label=label)
    ax2.set_xlabel('Image Slice [px]')
    ax2.set_ylabel('mean $y$ direction')
    ax2.set_xlim(data_2D.shape[0] / 2 - 0.25 * data_2D.shape[0],
                 data_2D.shape[0] / 2 + 0.25 * data_2D.shape[0])
    ax2.grid()
    # plt.semilogx()
    # plt.xlim(300,600)
    # if image_results_conv is not None:
    #     plt.savefig(
    #         image_results_conv.replace('.fits', 'result_lmfit_slices.pdf'),
    #         dpi=300, bbox_inches='tight')
    #     if show_figure == True:
    #         plt.show()
    #     else:
    #         plt.close()
    return(fig,ax1,ax2)

def plot_slices(data_2D, residual_2D, model_dict, image_results_conv=None,
                Rp_props=None, show_figure=True):
    plot_slice = np.arange(0, data_2D.shape[0])
    if Rp_props is not None:
        plotlim = Rp_props['c' + str(1) + '_rlast']
        # plotlim = 0
        # for i in range(Rp_props['ncomps']):
        #     plotlim = plotlim + Rp_props['c' + str(i + 1) + '_rlast']

    fig = plt.figure(figsize=(6, 6))
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)

    norm_plot = np.max(np.mean(data_2D, axis=0))
    ax1.plot(plot_slice, np.mean(data_2D, axis=0)/norm_plot, '--.', color='purple', ms=14,
             label='DATA')
    ax1.plot(plot_slice, np.mean(model_dict['model_total_conv']/norm_plot, axis=0), '.-',
             color='limegreen', linewidth=4, label='MODEL')
    ax1.plot(plot_slice, np.mean(model_dict['best_residual_conv']/norm_plot, axis=0), '.-',
             color='black', linewidth=4, label='RESIDUAL')
    try:
        ax1.plot(plot_slice, np.mean(residual_2D, axis=0)/norm_plot, '.-', color='grey',
                linewidth=4, label='MAP RESIDUAL')
    except:
        pass
    #     ax1.set_xlabel('$x$-slice')
    #     ax1.set_xaxis('off')
    #     ax1.set_xticks([])
    ax1.legend(fontsize=11)
    ax1.grid()
    ax1.set_ylabel('fractional mean $x$ direction')
    if Rp_props is not None:
        ax1.set_xlim(Rp_props['c1_x0c'] - plotlim, Rp_props['c1_x0c'] + plotlim)
    #     ax1.set_title('asd')
    # plt.plot(np.mean(shuffled_image,axis=0),color='red')

    norm_plot = np.max(np.mean(data_2D, axis=1))
    ax2.plot(plot_slice, np.mean(data_2D, axis=1)/norm_plot, '--.', color='purple', ms=14,
             label='DATA')
    ax2.plot(plot_slice, np.mean(model_dict['model_total_conv']/norm_plot, axis=1), '.-',
             color='limegreen', linewidth=4, label='MODEL')
    ax2.plot(plot_slice, np.mean(model_dict['best_residual_conv']/norm_plot, axis=1), '.-',
             color='black', linewidth=4, label='RESIDUAL')
    try:
        ax2.plot(plot_slice, np.mean(residual_2D, axis=1)/norm_plot, '.-', color='grey',
                linewidth=4, label='MAP RESIDUAL')
    except:
        pass
    ax2.set_xlabel('Image Slice [px]')
    ax2.set_ylabel('fractional mean $y$ direction')
    if Rp_props is not None:
        ax2.set_xlim(Rp_props['c1_y0c'] - plotlim, Rp_props['c1_y0c'] + plotlim)
    ax2.grid()
    # plt.semilogx()
    # plt.xlim(300,600)
    if image_results_conv is not None:
        plt.savefig(
            image_results_conv.replace('.fits', 'result_lmfit_slices.pdf'),
            dpi=300, bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()


def plot_fit_results(imagename, model_dict, image_results_conv,
                     sources_photometies,bkg_image=None,vmax_factor=0.1,data_2D_=None,
                     mask = None,
                     vmin_factor=3, obs_type = 'radio',
                     show_figure=True,crop=False,box_size=100):
    if data_2D_ is not None:
        data_2D = data_2D_
    else:
        data_2D = load_fits_data(imagename)

    plot_data_model_res(data_2D, modelname=model_dict['model_total_conv'],
               residualname=model_dict['best_residual_conv'],
               reference_image=imagename,
               NAME=image_results_conv[-2].replace('.fits',
                                                   '_data_model_res'),
               crop=crop, vmin_factor=vmin_factor,obs_type = obs_type,
               box_size=box_size)



    ncomponents = sources_photometies['ncomps']
    if sources_photometies is not None:
        plotlim =  4.0 * sources_photometies['c'+str(int(ncomponents))+'_Rp']
        if plotlim > data_2D.shape[0]/2:
            plotlim = data_2D.shape[0]/2
        # plotlim = 0
        # for i in range(ncomponents):
        #     plotlim = plotlim + sources_photometies['c' + str(i + 1) + '_rlast']

    model_name = image_results_conv[-2]
    residual_name = image_results_conv[-1]
    cell_size = get_cell_size(imagename)
    profile_data = {}
    # center = get_peak_pos(imagename)
    center = nd.maximum_position(data_2D*mask)[::-1]
    for i in range(ncomponents):
        component_name = image_results_conv[
            i]  # crop_image.replace('.fits','')+"_"+str(ncomponents)+"C_model_component_"+str(i+1)+special_name+'_IMFIT_opt.fits'
        Ir_r = get_profile(component_name,
                           center=center)
        profile_data['r' + str(i + 1)], profile_data['Ir' + str(i + 1)], \
        profile_data['c' + str(i + 1) + '_name'] = Ir_r[0], Ir_r[
            1], component_name

    r, ir = get_profile(data_2D, center=center)
    rmodel, irmodel = get_profile(model_name, center=center)
    rre, irre = get_profile(residual_name, center=center)
    if bkg_image is not None:
        r_bkg, ir_bkg = get_profile(bkg_image,center=center)

    # plt.plot(radiis[0],profiles[0])
    # plt.plot(radiis[1],profiles[1])
    # plt.plot(radiis[2],np.log(profiles[2]))
    # colors = ['black','purple','gray','red']
    colors = ['red', 'blue', 'teal', 'brown', 'cyan','orange','forestgreen',
              'pink', 'slategrey','darkseagreen','peru','royalblue','darkorange']

    plt.figure(figsize=(5, 5))
    plt.plot(r * cell_size, abs(ir), '--.', ms=10, color='purple', alpha=1.0,
             label='DATA')
    for i in range(ncomponents):
        #     try:
        #         plt.plot(profile_data['r'+str(i+1)],abs(profile_data['Ir'+str(i+1)])[0:r.shape[0]],'--',label='comp'+str(i+1),color=colors[i])
        plt.plot(profile_data['r' + str(i + 1)] * cell_size,
                 abs(profile_data['Ir' + str(i + 1)]), '--',
                 label='COMP_' + str(i + 1), color=colors[i])
    #     except:
    #         pass
    if bkg_image is not None:
        ir_model_data =  irmodel# + ir_bkg
        plt.plot(r_bkg * cell_size, abs(ir_bkg), '--', label='bkg', color='brown')
    else:
        ir_model_data = irmodel

    plt.plot(r * cell_size, abs(irre), '.-', label='RESIDUAL', color='black')
    plt.plot(r * cell_size, abs(ir_model_data), '--', color='limegreen', label='MODEL',
             linewidth=4)
    plt.semilogy()
    if obs_type == 'radio':
        plt.xlabel(r'Projected Radius $R$ [arcsec]')
        plt.ylabel(r'Radial Intensity $I(R)$ [Jy/beam]')
    else:
        plt.xlabel(r'Projected Radius $R$ [px]')
        plt.ylabel(r'Radial Intensity $I(R)$')
        # plt.semilogx()
    plt.legend(fontsize=11)
    plt.ylim(1e-7, -0.05 * np.log(ir[0]))
    # plt.xlim(0,3.0)
    plt.grid()
    if sources_photometies is not None:
        plt.xlim(0, plotlim * cell_size)
        idRp_main = int(sources_photometies['c1_Rp'])
        plt.axvline(r[idRp_main] * cell_size)
    plt.savefig(image_results_conv[-2].replace('.fits', 'result_lmfit_IR.pdf'),
                dpi=300, bbox_inches='tight')
    if show_figure == True:
        plt.show()
        # return(plt)
    else:
        plt.close()




# plt.savefig(config_file.replace('params_imfit.csv','result_lmfit_py_IR.pdf'),dpi=300, bbox_inches='tight')

def total_flux(data2D,image,mask=None,BA=None,
               sigma=6,iterations=3,dilation_size=7,PLOT=False,
               silent=True):

    if BA is None:
        try:
            BA = beam_area2(image)
        except:
            print('WARNING: Beam area not found, setting to 1.')
            BA = 1
    else:
        BA = BA
    if mask is None:
        _,mask = mask_dilation(data2D,sigma=sigma,iterations=iterations,
                               dilation_size=dilation_size,PLOT=PLOT)
    else:
        mask = mask
#     tf = np.sum((data2D> sigma*mad_std(data2D))*data2D)/beam_area2(image)
    blank_sum = np.sum(data2D)/BA
    sum3S = np.sum(data2D*(data2D> 3.0*mad_std(data2D)))/BA
    summask = np.sum(data2D*mask)/BA
    if silent==False:
        print('Blank Sum   = ',blank_sum)
        print('Sum  3sigma = ',sum3S)
        print('Sum mask    = ',summask)
    return(summask)

def total_flux_faster(data2D,mask):
#     tf = np.sum((data2D> sigma*mad_std(data2D))*data2D)/beam_area2(image)
    summask = np.sum(data2D*mask)
    return(summask)


def plot_decomp_results(imagename,compact,extended_model,data_2D_=None,
                        vmax_factor=0.5,vmin_factor=3,rms=None,
                        figsize=(13,13),nfunctions=None,
                        obs_type = 'radio',
                        special_name=''):

    decomp_results = {}
    if rms is None:
        rms = mad_std(load_fits_data(imagename))
    else:
        rms = rms

    max_factor =  load_fits_data(imagename).max()
#     compact = model_dict['model_c1_conv']
    if data_2D_ is not None:
        data_2D = data_2D_
    else:
        data_2D = load_fits_data(imagename)

    # rms_std_data = mad_std(data_2D)
    extended = data_2D  - compact
    if nfunctions == 1:
        residual_modeling = data_2D - (compact)
    else:
        residual_modeling = data_2D - (compact + extended_model)
    fig = plt.figure(figsize=figsize)
    ax1 = fig.add_subplot(3, 3, 1)
    # ax.yaxis.set_ticks([])
    ax1 = eimshow(imagename,ax=ax1,fig=fig,rms=rms,plot_title='Total Emission',
                 vmax_factor=vmax_factor,vmin_factor=vmin_factor)
    # cb = plt.colorbar(mappable=plt.gca().images[0],
    #                           cax=fig.add_axes([0.9, 0.65, 0.02, 0.2]))
    # ax.yaxis.set_ticks([])
    ax1.axis('off')
    ax2 = fig.add_subplot(3, 3, 2)
    # ax.yaxis.set_ticks([])
    ax2 = eimshow(compact,ax=ax2,fig=fig,rms=rms,
                 plot_title='Compact Emission',
                 vmax_factor=vmax_factor,vmin_factor=vmin_factor)
    # ax.yaxis.set_ticks([])
    ax2.axis('off')

    ax3 = fig.add_subplot(3, 3, 3)
    # ax.yaxis.set_ticks([])
    ax3 = eimshow(extended,ax=ax3,fig=fig,rms=rms,
                 plot_title='Diffuse Emission',vmax_factor=vmax_factor,
                 vmin_factor=vmin_factor)
    # ax.yaxis.set_ticks([])
    # cb = plt.colorbar(mappable=plt.gca().images[0],
    #                           cax=fig.add_axes([0.00, 0.65, 0.02, 0.2]))
    ax3.axis('off')

    ax4 = fig.add_subplot(3, 3, 4)
    slice_ext = np.sqrt(np.mean(extended,axis=0)**2.0 + np.mean(extended,axis=1)**2.0)
    if nfunctions == 1:
        slice_ext_model = np.sqrt(
            np.mean(residual_modeling, axis=0) ** 2.0 + np.mean(residual_modeling,
                                                             axis=1) ** 2.0)
    else:
        slice_ext_model = np.sqrt(
            np.mean(extended_model, axis=0) ** 2.0 + np.mean(extended_model,
                                                             axis=1) ** 2.0)
    slice_data = np.sqrt(np.mean(data_2D,axis=0)**2.0 + np.mean(data_2D,axis=1)**2.0)
    ax4.plot(slice_ext,label='COMPACT SUB')
    ax4.plot(slice_data,label='DATA')
    ax4.plot(slice_ext_model, label='EXTENDED MODEL')
    ax4.legend(prop={'size': 11})
    xlimit = [data_2D.shape[0] / 2 - 0.15 * data_2D.shape[0],
              data_2D.shape[0] / 2 + 0.15 * data_2D.shape[0]]
    ax4.set_xlim(xlimit[0],xlimit[1])
    # ax.semilogx()

    ax5 = fig.add_subplot(3, 3, 5)
    ax5.axis('off')

    try:
        omaj, omin, _, _, _ = beam_shape(imagename)
        dilation_size = int(
            np.sqrt(omaj * omin) / (2 * get_cell_size(imagename)))
    except:
        dilation_size = 10

    _, mask_model_rms_self_compact = mask_dilation(compact,
                                           sigma=1, dilation_size=dilation_size,
                                           iterations=2,PLOT=False)
    _, mask_data = mask_dilation(data_2D, rms=rms,
                                           sigma=6, dilation_size=dilation_size,
                                           iterations=2,PLOT=False)
    _, mask_model_rms_self_extended = mask_dilation(extended,
                                           sigma=6, dilation_size=dilation_size,
                                           iterations=2,PLOT=False)

    _, mask_model_rms_image_compact = mask_dilation(compact,
                                            rms=rms,
                                            sigma=1, dilation_size=dilation_size,
                                            iterations=2,PLOT=False)
    _, mask_model_rms_image_extended = mask_dilation(extended,
                                            rms=rms,
                                            sigma=6, dilation_size=dilation_size,
                                            iterations=2,PLOT=False)
    if nfunctions == 1:
        _, mask_model_rms_image_extended_model = mask_dilation(residual_modeling,
                                                rms=rms,
                                                sigma=1, dilation_size=dilation_size,
                                                iterations=2,PLOT=False)
    else:
        _, mask_model_rms_image_extended_model = mask_dilation(extended_model,
                                                rms=rms,
                                                sigma=1, dilation_size=dilation_size,
                                                iterations=2,PLOT=False)

    try:
        beam_area_px = beam_area2(imagename)
    except:
        beam_area_px = 1
    if obs_type == 'radio':
        flux_scale_factor = 1000
    else:
        flux_scale_factor = 1
    # print('Flux on compact (self rms) = ',
    #       1000*np.sum(compact*mask_model_rms_self_compact)/beam_area_px)
    # print('Flux on compact (data rms) = ',
    #       1000 * np.sum(compact * mask_model_rms_image_compact) / beam_area_px)
    flux_density_compact = flux_scale_factor*np.sum(
        compact*mask_model_rms_image_compact)/beam_area_px
    if nfunctions == 1:
        flux_density_extended_model = flux_scale_factor * np.sum(
            residual_modeling * mask_data) / beam_area_px
    else:
        flux_density_extended_model = flux_scale_factor * np.sum(
            extended_model * mask_data) / beam_area_px

    flux_density_ext_old = flux_scale_factor*total_flux(extended,imagename,BA=beam_area_px,
                                       mask = mask_model_rms_image_extended)
    flux_density_ext = flux_scale_factor*np.sum(
        extended*mask_data)/beam_area_px

    flux_data = flux_scale_factor*total_flux(data_2D,imagename,BA=beam_area_px,
                                       mask = mask_data)
    flux_density_ext_self_rms = flux_scale_factor*total_flux(extended,imagename,BA=beam_area_px,
                                       mask = mask_model_rms_self_extended)

    if nfunctions==1:
        flux_res = flux_data - (flux_density_compact)
    else:
        flux_res = flux_data - (
                    flux_density_extended_model + flux_density_compact)

    # print('Flux on extended (self rms) = ',flux_density_ext_self_rms)
    # print('Flux on extended (data rms) = ',flux_density_ext_old)
    # print('Flux on extended2 (data rms) = ', flux_density_ext)
    # print('Flux on extended model (data rms) = ', flux_density_extended_model)
    # print('Flux on data = ', flux_data)
    # print('Flux on residual = ', flux_res)

    decomp_results['flux_data'] = flux_data
    decomp_results['flux_density_ext_old'] = flux_density_ext_old
    decomp_results['flux_density_ext'] = flux_density_ext
    decomp_results['flux_density_ext_self_rms'] = flux_density_ext_self_rms
    decomp_results['flux_density_extended_model'] = flux_density_extended_model
    decomp_results['flux_density_compact'] = flux_density_compact
    decomp_results['flux_density_model'] = (flux_density_compact +
                                            flux_density_extended_model)
    decomp_results['flux_density_res'] = flux_res



    # print('r_half_light (old vs new) = {:0.2f} vs {:0.2f}'.format(p.r_half_light, p_copy.r_half_light))
    ax5.annotate(r"$S_\nu^{\rm core-comp}/S_\nu^{\rm total}=$"+'{:0.2f}'.format(flux_density_compact/flux_data),
                (0.33, 0.32), xycoords='figure fraction', fontsize=18)
    ax5.annotate(r"$S_\nu^{\rm ext}/S_\nu^{\rm total}\ \ \ =$"+'{:0.2f}'.format(flux_density_ext/flux_data),
                (0.33, 0.29), xycoords='figure fraction', fontsize=18)
    ax5.annotate(r"$S_\nu^{\rm ext \ model}/S_\nu^{\rm total}\ \ \ =$"+'{:0.2f}'.format(flux_density_extended_model/flux_data),
                (0.33, 0.26), xycoords='figure fraction', fontsize=18)
    ax5.annotate(r"$S_\nu^{\rm res}/S_\nu^{\rm total}\ \ \ =$"+'{:0.2f}'.format(flux_res/flux_data),
                (0.33, 0.23), xycoords='figure fraction', fontsize=18)
    
    # plt.tight_layout()
    plt.savefig(
        imagename.replace('.fits', '_ext_vs_comp'+special_name+'.jpg'),
        dpi=300,
        bbox_inches='tight')

    # save_data = True
    # if save_data == True:
    if obs_type != 'radio':
        exteded_file_name = imagename.replace('.fits', '') + \
                            special_name + '_extended.fits'
        pf.writeto(exteded_file_name,extended,overwrite=True)
        copy_header(imagename,exteded_file_name)
        compact_file_name = imagename.replace('.fits', '') + \
                            special_name + '_conv_compact.fits'
        pf.writeto(compact_file_name,compact,overwrite=True)
        copy_header(imagename,compact_file_name)
        decomp_results['compact_model_image'] = compact_file_name
        decomp_results['extended_model_image'] = exteded_file_name
    

    return(decomp_results)



def plot_interferometric_decomposition(imagename0, imagename,
                                       modelname, residualname,
                                       crop=False, box_size=512,
                                       max_percent_lowlevel=99.0,
                                       max_percent_highlevel=99.9999,
                                       NAME=None, EXT='.pdf',
                                       vmin0=None,
                                       run_phase = '1st',
                                       vmin_factor=3,vmax_factor=0.1,
                                       SPECIAL_NAME='', show_figure=True):
    """

    """
    fig = plt.figure(figsize=(16, 16))
    try:
        g = pf.getdata(imagename)
        I1 = pf.getdata(imagename0)
        if len(np.shape(g) == 4):
            g = g[0][0]
            I1 = I1[0][0]
        m = pf.getdata(modelname)
        r = pf.getdata(residualname)
    except:
        I1 = load_fits_data(imagename0)
        g = load_fits_data(imagename)
        m = load_fits_data(modelname)
        r = load_fits_data(residualname)

    dx1 = g.shape[0]/2
    dx = g.shape[0]/2
    if crop == True:
        xin, xen, yin, yen = do_cutout(imagename, box_size=box_size,
                                       center=None, return_='box')
        #         I1 = I1[int(2*xin):int(xen/2),int(2*yin):int(yen/2)]
        I1 = I1[int(xin + box_size / 1.25):int(xen - box_size / 1.25),
             int(yin + box_size / 1.25):int(yen - box_size / 1.25)]
        dx1 = I1.shape[0]/2
        # g = g[xin:xen,yin:yen]
        # m = m[xin:xen,yin:yen]
        # r = r[xin:xen,yin:yen]

    if mad_std(I1) == 0:
        std0 = I1.std()
    else:
        std0 = mad_std(I1)

    if mad_std(g) == 0:
        std = g.std()
    else:
        std = mad_std(g)

    if mad_std(r) == 0:
        std_r = r.std()
    else:
        std_r = mad_std(r)

    if mad_std(m) == 0:
        std_m = m.std()
    else:
        std_m = mad_std(m)

    #     print(I1)
    vmin0 = 3 * std  # 0.5*g.min()#
    if vmin0 is not None:
        vmin0 = vmin0
    else:
        vmin0 = 3 * std

    vmax0 = 1.0 * g.max()
    vmin = vmin_factor * std  # 0.5*g.min()#
    vmax = vmax_factor * g.max()
    vmin_r = 0.5 * r.min()  # 1*std_r
    vmax_r = 1.0 * r.max()
    vmin_m = 1 * mad_std(m)  # vmin#0.01*std_m#0.5*m.min()#
    vmax_m = m.max()  # vmax#0.5*m.max()

    # levels_I1 = np.geomspace(2*I1.max(), 1.5 * np.std(I1), 6)
    levels_I1 = np.geomspace(2 * I1.max(), vmin0, 6)
    levels_g = np.geomspace(2*g.max(), 3 * std, 6)
    levels_m = np.geomspace(2*m.max(), 20 * std_m, 6)
    levels_r = np.geomspace(2*r.max(), 3 * std_r, 6)
    levels_neg = np.asarray([-3]) * std
    script_R = "\u211B"
    script_M = "\u2133"
    if run_phase == '1st':
        title_labels = [r"$I_1^{\rm mask}[\sigma_{\mathrm{opt}}]$",
                        r"$I_2$",
                        r""+script_M+r"$_{1,2} = I_{1}^{\rm mask}[\sigma_{\mathrm{"
                                     r"opt}}] * "
                                     r"\theta_2$",
                        r""+script_R+r"$_{1,2} = I_2 -$"+script_M+r"$_{1,2}$"
                        ]

    if run_phase == '2nd':
        # title_labels = [r""+script_R+r"$_{1,2}$",
        #                 r"$I_3$",
        #                 r"$I_{1}^{\rm mask} * \theta_3 + $"+script_R+r"$_{1,2} "
        #                                                              r"* \theta_3$",
        #                 r""+script_R+r"$_{T}$"
        #                 ]
        title_labels = [r""+script_R+r"$_{1,2}$",
                        r"$I_3$",
                        r""+script_M+"$_{1,3} + $"+script_M+r"$_{2,3}$",
                        r""+script_R+r"$_{T}$"
                        ]

    if run_phase == 'compact':
        title_labels = [r""+script_R+r"$_{1,2}$",
                        r"$I_3$",
                        r"$I_{1}^{\rm mask} * \theta_3$",
                        r"$I_3 - I_{1}^{\rm mask} * \theta_3$"
                        ]

    # colors = [(0, 0, 0), (1, 1, 1)]
    # cmap_name = 'black_white'
    # import matplotlib.colors as mcolors
    # cm = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors,
    #                                                N=len(levels_g))
    cm = 'gray'
    #     norm = simple_norm(g,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    norm = visualization.simple_norm(g, stretch='linear',
                                     max_percent=max_percent_lowlevel)
    norm0 = simple_norm(abs(I1), min_cut=0.5 * np.std(I1), max_cut=vmax,
                        stretch='sqrt')  # , max_percent=max_percent_highlevel)
    norm2 = simple_norm(abs(g), min_cut=vmin, max_cut=vmax,
                        stretch='asinh',asinh_a=0.02)  # , max_percent=max_percent_highlevel)
    CM = 'magma_r'
    ax = fig.add_subplot(1, 4, 1)

    #     im = ax.imshow(I1, cmap='gray_r',norm=norm,alpha=0.2)

    #     im_plot = ax.imshow(g, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax, vmin=vmin)#norm=norm
    im_plot = ax.imshow(I1, cmap='magma_r', extent=[-dx1,dx1,-dx1,dx1],
                        origin='lower', alpha=1.0, norm=norm0)

    ax.set_title(title_labels[0])

    cmap_magma_r = plt.cm.get_cmap('magma_r')
    # contour_palette = ['#000000', '#444444', '#888888', '#DDDDDD']
    # contour_palette = ['#000000', '#222222', '#444444', '#666666', '#888888',
    #                    '#AAAAAA', '#CCCCCC', '#EEEEEE', '#FFFFFF']
    contour_palette = ['#000000', '#444444', '#666666', '#EEEEEE', '#EEEEEE',
                       '#FFFFFF']


    ax.contour(I1, levels=levels_I1[::-1], colors=contour_palette,
               extent=[-dx1, dx1, -dx1, dx1],
               linewidths=1.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)

    cell_size = get_cell_size(imagename0)

    xticks = np.linspace(-dx1, dx1, 5)
    xticklabels = np.linspace(-dx1*cell_size, +dx1*cell_size, 5)
    xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
    ax.set_yticks(xticks,xticklabels)
    ax.set_xticks(xticks,xticklabels)
    ax.set_xlabel(r'Offset [arcsec]')
    ax = add_beam_to_image(imagename0, ax=ax, dx=dx1,
                           cell_size=cell_size)
    # ax.set_yticks([])
    # ax.set_yticklabels([])

    #     cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.0,0.38,0.02,0.23]))#,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')
    ax = fig.add_subplot(1, 4, 2)
    #     im = ax.imshow(g, cmap='gray_r',norm=norm,alpha=0.2)

    #     im_plot = ax.imshow(g, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax, vmin=vmin)#norm=norm
    im_plot = ax.imshow(g, cmap='magma_r',extent=[-dx,dx,-dx,dx],
                        origin='lower', alpha=1.0, norm=norm2)

    ax.set_title(title_labels[1])


    ax.contour(g, levels=levels_g[::-1], colors=contour_palette,
               extent=[-dx,dx,-dx,dx],
               linewidths=1.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)

    # cb = plt.colorbar(mappable=plt.gca().images[0],
    #                   cax=fig.add_axes([-0.0, 0.40, 0.02,0.19]))  # ,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')

    cb = plt.colorbar(mappable=plt.gca().images[0],
                      cax=fig.add_axes([0.07, 0.40, 0.02,0.19]),
                      orientation='vertical',shrink=1, aspect='auto',
                      pad=1, fraction=1.0,
                      drawedges=False, ticklocation='left')
    cb.formatter = CustomFormatter(factor=1000, useMathText=True)
    cb.update_ticks()
#     print('++++++++++++++++++++++')
#     print(plt.gca().images[0])
    cb.set_label(r'Flux Density [mJy/beam]', labelpad=1)
    cb.ax.xaxis.set_tick_params(pad=1)
    cb.ax.tick_params(labelsize=12)
    cb.outline.set_linewidth(1)
    # cb.dividers.set_color('none')

    cell_size = get_cell_size(imagename)
    # ax = add_beam_to_image(imagename, ax=ax, dx=dx,
    #                        cell_size=cell_size)
    xticks = np.linspace(-dx, dx, 5)
    xticklabels = np.linspace(-dx*cell_size, +dx*cell_size, 5)
    xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
    ax.set_yticks(xticks,xticklabels)
    ax.set_xticks(xticks,xticklabels)
    # ax.set_yticks([])
    # ax.set_yticklabels([])

    ax = plt.subplot(1, 4, 3)

    #     im_plot = ax.imshow(m, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax_m, vmin=vmin_m)#norm=norm
    im_plot = ax.imshow(m, cmap='magma_r',extent=[-dx,dx,-dx,dx],
                        origin='lower', alpha=1.0, norm=norm2)
    ax.set_title(title_labels[2])
    ax.contour(m, levels=levels_g[::-1],
               colors=contour_palette,
               extent=[-dx,dx,-dx,dx],
               linewidths=1.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)
    #     cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.08,0.3,0.02,0.4]))#,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')

    cell_size = get_cell_size(modelname)
    xticks = np.linspace(-dx, dx, 5)
    xticklabels = np.linspace(-dx*cell_size, +dx*cell_size, 5)
    xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
    ax.set_yticks(xticks,xticklabels)
    ax.set_xticks(xticks,xticklabels)
    # ax = add_beam_to_image(modelname, ax=ax, dx=dx,
    #                        cell_size=cell_size)
    # ax.set_yticks([])
    # ax.set_yticklabels([])



    ax = plt.subplot(1, 4, 4)
    norm_re = simple_norm(r, min_cut=vmin, max_cut=vmax, stretch='sqrt')  # , max_percent=max_percent_highlevel)
    #     norm = simple_norm(r,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    #     ax.imshow(r,origin='lower',cmap='magma_r',alpha=1.0,vmax=vmax_r, vmin=vmin)#norm=norm
    ax.imshow(r, origin='lower',extent=[-dx,dx,-dx,dx],
              cmap='magma_r', alpha=1.0, norm=norm2)
    #     ax.imshow(r, cmap='magma_r',norm=norm,alpha=0.3,origin='lower')

    ax.contour(r, levels=levels_r[::-1],
               extent=[-dx,dx,-dx,dx],
               colors=contour_palette,
               linewidths=1.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)
    ax.contour(r, levels=levels_neg[::-1],extent=[-dx,dx,-dx,dx],
               colors='k', linewidths=1.0,
               alpha=1.0)

    ax = add_beam_to_image(imagename, ax=ax, dx=dx,
                           cell_size=cell_size)

    cell_size = get_cell_size(residualname)
    xticks = np.linspace(-dx, dx, 5)
    xticklabels = np.linspace(-dx*cell_size, +dx*cell_size, 5)
    xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
    ax.set_yticks(xticks,xticklabels)
    ax.set_xticks(xticks,xticklabels)
    # ax.set_yticks([])
    # ax.set_yticklabels([])


    ax.set_title(title_labels[3])
    #     cb1=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([0.91,0.40,0.02,0.19]))
    #     cb1.set_label(r'Flux [Jy/beam]',labelpad=1)
    #     cb1.ax.xaxis.set_tick_params(pad=1)
    #     cb1.ax.tick_params(labelsize=12)
    #     cb1.outline.set_linewidth(1)
    if NAME is not None:
        plt.savefig(NAME.replace('.fits', '') + SPECIAL_NAME + EXT, dpi=300,
                    bbox_inches='tight')
        plt.savefig(NAME.replace('.fits', '') + SPECIAL_NAME + '.jpg', dpi=300,
                    bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()

def fast_plot(imagename0, imagename, modelname, residualname, crop=False, box_size=512,
              max_percent_lowlevel=99.0, max_percent_highlevel=99.9999,
              NAME=None, EXT='.pdf', SPECIAL_NAME='', show_figure=True):

    """
    Fast plotting of image <> model <> residual images.
    TO BE REMOVED

    CHECK >> PLOT_INT_DEC
    """
    fig = plt.figure(figsize=(16, 16))
    try:
        g = pf.getdata(imagename)
        I1 = pf.getdata(imagename0)
        if len(np.shape(g) == 4):
            g = g[0][0]
            I1 = I1[0][0]
        m = pf.getdata(modelname)
        r = pf.getdata(residualname)
    except:
        I1 = load_fits_data(imagename0)
        g = load_fits_data(imagename)
        m = load_fits_data(modelname)
        r = load_fits_data(residualname)

    if crop == True:
        xin, xen, yin, yen = do_cutout(imagename, box_size=box_size, center=None, return_='box')
        #         I1 = I1[int(2*xin):int(xen/2),int(2*yin):int(yen/2)]
        I1 = I1[int(xin + box_size / 1.25):int(xen - box_size / 1.25),
             int(yin + box_size / 1.25):int(yen - box_size / 1.25)]
        # g = g[xin:xen,yin:yen]
        # m = m[xin:xen,yin:yen]
        # r = r[xin:xen,yin:yen]

    if mad_std(I1) == 0:
        std0 = I1.std()
    else:
        std0 = mad_std(I1)

    if mad_std(g) == 0:
        std = g.std()
    else:
        std = mad_std(g)

    if mad_std(r) == 0:
        std_r = r.std()
    else:
        std_r = mad_std(r)

    if mad_std(m) == 0:
        std_m = m.std()
    else:
        std_m = mad_std(m)

    #     print(I1)
    vmin0 = 3 * std  # 0.5*g.min()#
    vmax0 = 1.0 * g.max()
    vmin = 0.5 * std  # 0.5*g.min()#
    vmax = 1.0 * g.max()
    vmin_r = 0.5 * r.min()  # 1*std_r
    vmax_r = 1.0 * r.max()
    vmin_m = 1 * mad_std(m)  # vmin#0.01*std_m#0.5*m.min()#
    vmax_m = m.max()  # vmax#0.5*m.max()

    levels_I1 = np.geomspace(I1.max(), 1.5 * np.std(I1), 7)
    levels_g = np.geomspace(g.max(), 3 * std, 7)
    levels_m = np.geomspace(m.max(), 20 * std_m, 7)
    levels_r = np.geomspace(r.max(), 3 * std_r, 7)

    #     norm = simple_norm(g,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    norm = visualization.simple_norm(g, stretch='linear', max_percent=max_percent_lowlevel)
    norm0 = simple_norm(abs(I1), min_cut=0.5 * np.std(I1), max_cut=vmax,
                        stretch='sqrt')  # , max_percent=max_percent_highlevel)
    norm2 = simple_norm(abs(g), min_cut=vmin, max_cut=vmax, stretch='sqrt')  # , max_percent=max_percent_highlevel)
    CM = 'magma_r'
    ax = fig.add_subplot(1, 4, 1)

    #     im = ax.imshow(I1, cmap='gray_r',norm=norm,alpha=0.2)

    #     im_plot = ax.imshow(g, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax, vmin=vmin)#norm=norm
    im_plot = ax.imshow(I1, cmap='magma_r', origin='lower', alpha=1.0, norm=norm0)

    ax.set_title(r'$I_1^{\rm mask}$')

    ax.contour(I1, levels=levels_I1[::-1], colors='#009E73', linewidths=0.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)
    #     cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.0,0.38,0.02,0.23]))#,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')

    ax = fig.add_subplot(1, 4, 2)
    #     im = ax.imshow(g, cmap='gray_r',norm=norm,alpha=0.2)

    #     im_plot = ax.imshow(g, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax, vmin=vmin)#norm=norm
    im_plot = ax.imshow(g, cmap='magma_r', origin='lower', alpha=1.0, norm=norm2)

    ax.set_title(r'$I_2$')

    ax.contour(g, levels=levels_g[::-1], colors='#009E73', linewidths=0.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)
    cb = plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.0, 0.40, 0.02,
                                                                      0.19]))  # ,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')
    cb.set_label(r'Flux [Jy/beam]', labelpad=1)
    cb.ax.xaxis.set_tick_params(pad=1)
    cb.ax.tick_params(labelsize=12)
    cb.outline.set_linewidth(1)
    # cb.dividers.set_color('none')

    ax = plt.subplot(1, 4, 3)

    #     im_plot = ax.imshow(m, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax_m, vmin=vmin_m)#norm=norm
    im_plot = ax.imshow(m, cmap='magma_r', origin='lower', alpha=1.0, norm=norm2)
    ax.set_title(r'$I_{1}^{\rm mask} * \theta_2$')
    ax.contour(m, levels=levels_g[::-1], colors='#009E73', linewidths=0.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)
    #     cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.08,0.3,0.02,0.4]))#,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')

    ax = plt.subplot(1, 4, 4)
    norm_re = simple_norm(r, min_cut=vmin, max_cut=vmax, stretch='sqrt')  # , max_percent=max_percent_highlevel)
    #     norm = simple_norm(r,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    #     ax.imshow(r,origin='lower',cmap='magma_r',alpha=1.0,vmax=vmax_r, vmin=vmin)#norm=norm
    ax.imshow(r, origin='lower', cmap='magma_r', alpha=1.0, norm=norm2)
    #     ax.imshow(r, cmap='magma_r',norm=norm,alpha=0.3,origin='lower')

    ax.contour(r, levels=levels_r[::-1], colors='#009E73', linewidths=0.2, alpha=1.0)  # cmap='Reds', linewidths=0.75)
    ax.set_yticks([])
    ax.set_title(r'$R_{12}$')
    #     cb1=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([0.91,0.40,0.02,0.19]))
    #     cb1.set_label(r'Flux [Jy/beam]',labelpad=1)
    #     cb1.ax.xaxis.set_tick_params(pad=1)
    #     cb1.ax.tick_params(labelsize=12)
    #     cb1.outline.set_linewidth(1)
    if NAME is not None:
        plt.savefig(NAME.replace('.fits', '') + SPECIAL_NAME + EXT, dpi=300,
                    bbox_inches='tight')
        plt.savefig(NAME.replace('.fits', '') + SPECIAL_NAME + '.jpg', dpi=300,
                    bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()


def plot_data_model_res(imagename, modelname, residualname, reference_image,
                        crop=False,box_size=512, NAME=None, CM='magma_r',
                        vmin_factor=3.0,vmax_factor=0.1,
                        max_percent_lowlevel=99.0, max_percent_highlevel=99.9999,
                        obs_type = 'radio',
                        ext='.pdf', show_figure=True):
    """
    Plots fitting results: image <> model <> residual images.

    Parameters
    ----------
    imagename : str
        Path to the image.
    modelname : str
        Path to the model.
    residualname : str
        Path to the residual.
    reference_image : str
        Path to the reference image.
    crop : bool, optional
        Crop the image to the box_size. The default is False.
    box_size : int, optional
        Size of the box to crop the image. The default is 512.
    NAME : str, optional
        Name of the output file. The default is None.
    CM : str, optional
        Colormap. The default is 'magma_r'.
    vmin_factor : float, optional
        Factor to multiply the standard deviation of the image to set the
        minimum value of the colormap. The default is 3.0.
    vmax_factor : float, optional
        Factor to multiply the maximum value of the image to set the
        maximum value of the colormap. The default is 0.1.
    max_percent_lowlevel : float, optional
        Maximum percentile to set the low level of the colormap. The default is 99.0.
    max_percent_highlevel : float, optional
        Maximum percentile to set the high level of the colormap. The default is 99.9999.
    ext : str, optional
        Extension of the output file. The default is '.pdf'.
    show_figure : bool, optional
        Show the figure. The default is True.
    """
    fig = plt.figure(figsize=(12, 12))
    try:
        try:
            g = pf.getdata(imagename)
            # if len(np.shape(g)==4):
            #     g = g[0][0]
            m = pf.getdata(modelname)
            r = pf.getdata(residualname)
        except:
            g = imagename
            m = modelname
            r = residualname
    except:
        g = load_fits_data(imagename)
        m = load_fits_data(modelname)
        r = load_fits_data(residualname)

    if crop == True:
        xin, xen, yin, yen = do_cutout(reference_image, box_size=box_size,
                                       center=None, return_='box')
        #         I1 = I1[int(2*xin):int(xen/2),int(2*yin):int(yen/2)]
        g = g[xin:xen, yin:yen]
        m = m[xin:xen, yin:yen]
        r = r[xin:xen, yin:yen]

    if mad_std(g) == 0:
        std = g.std()
    else:
        std = mad_std(g)

    if residualname is not None:
        if mad_std(r) == 0:
            std_r = r.std()
        else:
            std_r = mad_std(r)
    else:
        std_r = std

    if mad_std(m) == 0:
        std_m = m.std()
    else:
        std_m = mad_std(m)

    dx = g.shape[0]/2
    try:
        cell_size = get_cell_size(reference_image)
        axis_units_label = r'Offset [arcsec]'
    except:
        print('No cell or pixel size information in the image wcs/header. '
              'Setting cell/pixel size = 1.')
        cell_size = 1
        axis_units_label = r'Offset [px]'

    #     print(I1)
    vmin = vmin_factor * std  # 0.5*g.min()#
    vmax = vmax_factor * g.max()
    vmin_r = vmin  # 1.0*r.min()#1*std_r
    vmax_r = vmax #1.0 * r.max()
    vmin_m = vmin  # 1*mad_std(m)#vmin#0.01*std_m#0.5*m.min()#
    vmax_m = vmax  # 0.5*m.max()#vmax#0.5*m.max()

    levels_g = np.geomspace(3*g.max(), 3 * std, 6)
    levels_m = np.geomspace(3*m.max(), 10 * std_m, 6)
    levels_r = np.geomspace(3*r.max(), 3 * std_r, 6)

    #     norm = simple_norm(g,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    norm = visualization.simple_norm(g, stretch='linear',
                                     max_percent=max_percent_lowlevel)
    norm2 = simple_norm(abs(g), min_cut=vmin, max_cut=vmax,
                        stretch='asinh',asinh_a=0.05)  # , max_percent=max_percent_highlevel)

    ax = fig.add_subplot(1, 3, 1)

    #     im = ax.imshow(g, cmap='gray_r',norm=norm,alpha=0.2)

    #     im_plot = ax.imshow(g, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax, vmin=vmin)#norm=norm
    im_plot = ax.imshow(g, cmap=CM,
                        origin='lower', extent=[-dx,dx,-dx,dx],
                        alpha=1.0, norm=norm2)

    ax = add_beam_to_image(imagename=reference_image, ax=ax,
                           dx=dx,cell_size=cell_size)

    ax.set_title(r'Data')

    contour_palette = ['#000000', '#444444', '#666666', '#EEEEEE',
                       '#EEEEEE', '#FFFFFF']

    ax.contour(g, levels=levels_g[::-1], colors=contour_palette,
               linewidths=1.2,extent=[-dx, dx, -dx, dx],
               alpha=1.0)  # cmap='Reds', linewidths=0.75)
    cb1 = plt.colorbar(mappable=plt.gca().images[0],
                       cax=fig.add_axes([0.91, 0.40, 0.02, 0.19]))
    if obs_type == 'radio':
        cb1.formatter = CustomFormatter(factor=int(1000/vmax_factor), 
                                        useMathText=True)
        cb1.update_ticks()
        cb1.set_label(r'Flux Density [mJy/beam]', labelpad=1)
        cb1.ax.xaxis.set_tick_params(pad=1)
        cb1.ax.tick_params(labelsize=12)
        cb1.outline.set_linewidth(1)
    else:
        cb1.formatter = CustomFormatter(factor=int(1/vmax_factor), 
                                        useMathText=True)
        cb1.update_ticks()
        cb1.set_label(r'Pixel Intensity]', labelpad=1)
        cb1.ax.xaxis.set_tick_params(pad=1)
        cb1.ax.tick_params(labelsize=12)
        cb1.outline.set_linewidth(1)
        
    """
    # No need for this additional colorbar.
    cb = plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes(
        [-0.0, 0.40, 0.02,
         0.19]))  # ,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')
    cb.formatter = CustomFormatter(factor=1000, useMathText=True)
    cb.update_ticks()
    cb.set_label(r'Flux [mJy/beam]', labelpad=1)
    cb.ax.xaxis.set_tick_params(pad=1)
    cb.ax.tick_params(labelsize=12)
    cb.outline.set_linewidth(1)
    # cb.dividers.set_color('none')
    """
    xticks = np.linspace(-dx, dx, 4)
    xticklabels = np.linspace(-dx*cell_size, +dx*cell_size, 4)
    xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
    ax.set_yticks(xticks,xticklabels)
    ax.set_xticks(xticks,xticklabels)
    ax.set_xlabel(axis_units_label)
    # ax.set_yticks([])
    # ax.set_yticklabels([])


    ax = plt.subplot(1, 3, 2)

    #     im_plot = ax.imshow(m, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax_m, vmin=vmin_m)#norm=norm
    norm_mod = simple_norm(m, min_cut=vmin, max_cut=vmax,
                           stretch='asinh',asinh_a=0.05)  # , max_percent=max_percent_highlevel)

    im_plot = ax.imshow(m, cmap=CM,
                        origin='lower',extent=[-dx,dx,-dx,dx],
                        alpha=1.0,
                        norm=norm_mod)
    ax.set_title(r'Model')
    ax.contour(m, levels=levels_g[::-1],
               extent=[-dx, dx, -dx, dx],
               colors=contour_palette, linewidths=1.2,
               alpha=1.0)  # cmap='Reds', linewidths=0.75)
    #     cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.08,0.3,0.02,0.4]))#,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')
    ax.set_yticks(xticks,xticklabels)
    ax.set_xticks(xticks,xticklabels)
    ax.set_xlabel(axis_units_label)
    # ax.set_yticks([])
    ax.set_yticklabels([])

    ax = plt.subplot(1, 3, 3)
    norm_re = simple_norm(r, min_cut=vmin, max_cut=vmax,
                          stretch='asinh',asinh_a=0.05)  # , max_percent=max_percent_highlevel)
    #     norm = simple_norm(r,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    #     ax.imshow(r,origin='lower',cmap='magma_r',alpha=1.0,vmax=vmax_r, vmin=vmin)#norm=norm
    ax.imshow(r, origin='lower',extent=[-dx, dx, -dx, dx],
              cmap=CM, alpha=1.0, norm=norm_re)
    #     ax.imshow(r, cmap='magma_r',norm=norm,alpha=0.3,origin='lower')

    ax.contour(r, levels=levels_g[::-1],
               extent=[-dx, dx, -dx, dx],
               colors=contour_palette, linewidths=1.2,
               alpha=1.0)  # cmap='Reds', linewidths=0.75)
    levels_neg = np.asarray([-6 * std])

    ax.contour(r, levels=levels_neg[::-1],
               extent=[-dx, dx, -dx, dx],
               colors='k', linewidths=1.0,
               alpha=1.0)
    ax.set_yticks(xticks, xticklabels)
    ax.set_xticks(xticks, xticklabels)
    ax.set_xlabel(axis_units_label)
    # ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_title(r'Residual')
    # cb1.dividers.set_color('none')
    if NAME != None:
        plt.savefig(NAME + ext, dpi=300, bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()


def plot_image_model_res(imagename, modelname, residualname, reference_image, crop=False,
               box_size=512, NAME=None, CM='magma_r',
               vmin_factor=3.0,vmax_factor=0.1,
               max_percent_lowlevel=99.0, max_percent_highlevel=99.9999,
               ext='.pdf', show_figure=True):
    """
    Fast plotting of image <> model <> residual images.

    """
    fig = plt.figure(figsize=(12, 12))
    try:
        try:
            g = pf.getdata(imagename)
            # if len(np.shape(g)==4):
            #     g = g[0][0]
            m = pf.getdata(modelname)
            r = pf.getdata(residualname)
        except:
            g = imagename
            m = modelname
            r = residualname
    except:
        g = load_fits_data(imagename)
        m = load_fits_data(modelname)
        r = load_fits_data(residualname)

    if crop == True:
        xin, xen, yin, yen = do_cutout(reference_image, box_size=box_size,
                                       center=None, return_='box')
        #         I1 = I1[int(2*xin):int(xen/2),int(2*yin):int(yen/2)]
        g = g[xin:xen, yin:yen]
        m = m[xin:xen, yin:yen]
        r = r[xin:xen, yin:yen]

    if mad_std(g) == 0:
        std = g.std()
    else:
        std = mad_std(g)

    if mad_std(r) == 0:
        std_r = r.std()
    else:
        std_r = mad_std(r)

    if mad_std(m) == 0:
        std_m = m.std()
    else:
        std_m = mad_std(m)

    #     print(I1)
    vmin = vmin_factor * std  # 0.5*g.min()#
    vmax = vmax_factor * g.max()
    vmin_r = vmin  # 1.0*r.min()#1*std_r
    vmax_r = vmax #1.0 * r.max()
    vmin_m = vmin  # 1*mad_std(m)#vmin#0.01*std_m#0.5*m.min()#
    vmax_m = vmax  # 0.5*m.max()#vmax#0.5*m.max()

    levels_g = np.geomspace(g.max(), 3 * std, 7)
    levels_m = np.geomspace(m.max(), 10 * std_m, 7)
    levels_r = np.geomspace(r.max(), 3 * std_r, 7)

    #     norm = simple_norm(g,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    norm = visualization.simple_norm(g, stretch='linear',
                                     max_percent=max_percent_lowlevel)
    norm2 = simple_norm(abs(g), min_cut=vmin, max_cut=vmax,
                        stretch='asinh',asinh_a=0.05)  # , max_percent=max_percent_highlevel)

    ax = fig.add_subplot(2, 3, 1)

    #     im = ax.imshow(g, cmap='gray_r',norm=norm,alpha=0.2)

    #     im_plot = ax.imshow(g, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax, vmin=vmin)#norm=norm
    im_plot = ax.imshow(g, cmap=CM, origin='lower', alpha=1.0, norm=norm2)

    ax.set_title(r'Image')

    ax.contour(g, levels=levels_g[::-1], colors='#009E73', linewidths=0.2,
               alpha=1.0)  # cmap='Reds', linewidths=0.75)
    cb = plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes(
        [-0.0, 0.40, 0.02,
         0.19]))  # ,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')
    cb.set_label(r'Flux [Jy/beam]', labelpad=1)
    cb.ax.xaxis.set_tick_params(pad=1)
    cb.ax.tick_params(labelsize=12)
    cb.outline.set_linewidth(1)
    # cb.dividers.set_color('none')

    ax = plt.subplot(2, 3, 2)

    #     im_plot = ax.imshow(m, cmap='magma_r',origin='lower',alpha=1.0,vmax=vmax_m, vmin=vmin_m)#norm=norm
    norm_mod = simple_norm(m, min_cut=vmin, max_cut=vmax,
                           stretch='asinh',asinh_a=0.05)  # , max_percent=max_percent_highlevel)

    im_plot = ax.imshow(m, cmap=CM, origin='lower', alpha=1.0,
                        norm=norm_mod)
    ax.set_title(r'Model')
    ax.contour(m, levels=levels_g[::-1], colors='#009E73', linewidths=0.2,
               alpha=1.0)  # cmap='Reds', linewidths=0.75)
    #     cb=plt.colorbar(mappable=plt.gca().images[0], cax=fig.add_axes([-0.08,0.3,0.02,0.4]))#,format=ticker.FuncFormatter(fmt))#cax=fig.add_axes([0.01,0.7,0.5,0.05]))#, orientation='horizontal')

    ax = plt.subplot(2, 3, 3)
    norm_re = simple_norm(r, min_cut=vmin, max_cut=vmax,
                          stretch='asinh',asinh_a=0.05)  # , max_percent=max_percent_highlevel)
    #     norm = simple_norm(r,stretch='asinh',asinh_a=0.01)#,vmin=vmin,vmax=vmax)
    #     ax.imshow(r,origin='lower',cmap='magma_r',alpha=1.0,vmax=vmax_r, vmin=vmin)#norm=norm
    ax.imshow(r, origin='lower', cmap=CM, alpha=1.0, norm=norm_re)
    #     ax.imshow(r, cmap='magma_r',norm=norm,alpha=0.3,origin='lower')

    ax.contour(r, levels=levels_r[::-1], colors='#009E73', linewidths=0.2,
               alpha=1.0)  # cmap='Reds', linewidths=0.75)
    ax.set_yticks([])
    ax.set_title(r'Residual')
    cb1 = plt.colorbar(mappable=plt.gca().images[0],
                       cax=fig.add_axes([0.91, 0.40, 0.02, 0.19]))
    cb1.set_label(r'Flux [Jy/beam]', labelpad=1)
    cb1.ax.xaxis.set_tick_params(pad=1)
    cb1.ax.tick_params(labelsize=12)
    cb1.outline.set_linewidth(1)
    return(ax,plt,fig)
    # cb1.dividers.set_color('none')
    # if NAME != None:
    #     plt.savefig(NAME + ext, dpi=300, bbox_inches='tight')
    #     if show_figure == True:
    #         plt.show()
    #     else:
    #         plt.close()


def eimshow(imagename, crop=False, box_size=128, center=None, with_wcs=True,
            vmax=None,fig=None,
            dx_shift = 0,
            dy_shift = 0,
            vmax_factor=0.5, neg_levels=np.asarray([-3]), CM='magma_r',
            cmap_cont='magma_r',
            rms=None, plot_title=None, apply_mask=False,
            add_contours=True, extent=None, projection='offset', add_beam=False,
            vmin_factor=3, plot_colorbar=False, cbar_orientation=None,pad=-0.2,
            figsize=(5, 5), aspect=None,n_pts_labels = 5,cbar_n_points=6,num_contours=6,
            show_axis='on',flux_units='mJy',add_frequency=False,freq_label=None,
            source_distance=None, scalebar_length=250 * u.pc,
            ax=None, return_extent=False,no_show=False,
            save_name=None, special_name='',
            text_annotation=None, text_position=(0.3, 0.8),
            verbose = 0):
    """
    Customised imshow function for plotting images. It includes the option to
    automatically add contour levels, custom projections, colorbars, scalebars, etc.

    Parameters
    ----------
    imagename : str, np.ndarray
        Path to the image.
    crop : bool, optional
        Crop the image to the box_size. The default is False.
    box_size : int, optional
        Size of the box to crop the image. The default is 128.
    center : tuple, optional
        Center of the image. The default is None.
    with_wcs : bool, optional
        Use the wcs information of the image. The default is True.
    vmax : float, optional
        Maximum value of the colormap. The default is None.
    vmax_factor : float, optional
        Factor to multiply the maximum value of the image to set the
        maximum value of the colormap. The default is 0.5.
    neg_levels : np.ndarray, optional
        Negative levels for the contours. The default is np.asarray([-3]).
    CM : str, optional
        Colormap. The default is 'magma_r'.
    cmap_cont : str, optional
        Colormap for the contours. The default is 'terrain'.
    rms : float, optional
        RMS of the image. The default is None.
    plot_title : str, optional
        Title of the plot. The default is None.
    apply_mask : bool, optional
        Apply a mask to the image. The default is False.
    add_contours : bool, optional
        Add contours to the image. The default is True.
    extent : list, optional
        Extent of the image. The default is None.
    projection : str, optional
        Projection of the image. The default is 'offset'.
    add_beam : bool, optional
        Add the beam to the image. The default is False.
    vmin_factor : int, optional
        Factor to multiply the standard deviation of the image to set the
        minimum value of the colormap. The default is 3.
    plot_colorbar : bool, optional
        Plot the colorbar. The default is True.
    figsize : tuple, optional
        Size of the figure. The default is (5, 5).
    aspect : str, optional
        Aspect of the image. The default is None.
    show_axis : str, optional
        Show the axis. The default is 'on'.
    flux_units : str, optional
        Units of the flux. The default is 'Jy'.
    source_distance : float, optional
        Distance to the source. The default is None.
    scalebar_length : float, optional
        Length of the scalebar. The default is 250 * u.pc.
    ax : matplotlib.pyplot.axis, optional
        Axis of the plot. The default is None.
    save_name : str, optional
        Name of the output file. The default is None.
    special_name : str, optional
        Special name for the output file. The default is ''.
    """
    try:
        import cmasher as cmr
        # print('Imported cmasher for density maps.'
        #       'If you would like to use, examples:'
        #       'CM = cmr.ember,'
        #       'CM = cmr.flamingo,'
        #       'CM = cmr.gothic'
        #       'CM = cmr.lavender')
        """
        ... lilac,rainforest,sepia,sunburst,torch.
        Diverging: copper,emergency,fusion,infinity,pride'
        """
    except:
        print('Error importing cmasher. If you want '
              'to use its colormaps, install it. '
              'Then you can use for example:'
              'CM = cmr.flamingo')
    if ax is None:
        fig = plt.figure(figsize=figsize)
        # if isinstance(box_size, int):
        #     fig = plt.figure(figsize=figsize)
        # else:
        #     scale_fig_x = box_size[0]/box_size[1]
        #     fig = plt.figure(figsize=(figsize[0]*scale_fig_x,figsize[1]))
    else:
        if fig is None:
            fig = plt.figure(figsize=figsize)
        else:
            pass
    if isinstance(imagename, str) == True:

        if with_wcs == True:
            hdu = pf.open(imagename)
            ww = WCS(hdu[0].header, naxis=2)
            try:
                if len(np.shape(hdu[0].data) == 2):
                    g = hdu[0].data[0][0]
                else:
                    g = hdu[0].data
            except:
                g = load_fits_data(imagename)
        if with_wcs == False:
            g = load_fits_data(imagename)

        if crop == True:
            yin, yen, xin, xen = do_cutout(imagename, box_size=box_size,
                                           center=center, return_='box')
            g = g[xin:xen, yin:yen]
            # crop = False

        if apply_mask == True:
            _, mask_d = mask_dilation(imagename, cell_size=None,
                                      sigma=6, rms=None,
                                      dilation_size=None,
                                      iterations=2, dilation_type='disk',
                                      PLOT=False, show_figure=False)
            print('Masking emission....')
            g = g * mask_d[xin:xen, yin:yen]
        
        g = np.nan_to_num(g,nan=0)
    
    else:
        g = np.nan_to_num(imagename,nan=0)
        mask_d = 1
        # print('3', g)

        if crop == True:
            xin, xen, yin, yen = do_cutout(imagename, box_size=box_size,
                                           center=center, return_='box')
            g = g[xin:xen, yin:yen]
            if apply_mask == True:
                print('Masking emission....')
                g = g * mask_d[xin:xen, yin:yen]

    if rms is not None:
        std = rms
    else:
        if mad_std(g) == 0:
            """
            About std:
                mad_std is much more robust than np.std.
                But:
                    if mad_std is applied to a masked image, with zero
                    values outside the emission region, mad_std(image) is zero!
                    So, in that case, np.std is a good option.
            """
            # print('5', g)
            std = np.nanstd(g)
        else:
            std = mad_std(g)

    if isinstance(imagename, str) == True:
        try:
            cell_size = get_cell_size(imagename)
            axis_units_label = r'Offset [arcsec]'
        except:
            print(
                'No cell or pixel size information in the image wcs/header. '
                'Setting cell/pixel size = 1.')
            cell_size = 1
            axis_units_label = r'Offset [px]'
    else:
        cell_size = 1
        axis_units_label = r'Offset [px]'
    
    dx = g.shape[1] / 2
    dy = g.shape[0] / 2
    if ax is None:

        if (projection == 'celestial') and (with_wcs == True) and (isinstance(
                imagename, str) == True):
            ax = fig.add_subplot(projection=ww.celestial)
            ax.set_xlabel('RA', fontsize=14)
            ax.set_ylabel('DEC', fontsize=14)
            ax.grid()

        if isinstance(imagename, str) == False:
            projection = 'px'

        if projection == 'offset':
            ax = fig.add_subplot()
            # dx = g.shape[0] / 2
            axis_units_label = r'Offset [arcsec]'
            ax.set_xlabel(axis_units_label, fontsize=14)

        dx = g.shape[1] / 2
        dy = g.shape[0] / 2
        if projection == 'px':
            ax = fig.add_subplot()
            cell_size = 1
            # dx = g.shape[0] / 2
            ax.set_xlabel('x pix')
            ax.set_ylabel('y pix')
            axis_units_label = r'Offset [px]'
            ax.set_xlabel(axis_units_label, fontsize=14)
            ax.set_ylabel(axis_units_label, fontsize=14)

    else:
        projection = 'offset'
        if projection == 'offset':
            # ax = fig.add_subplot()
            # dx = g.shape[0] / 2
            axis_units_label = r'Offset [arcsec]'
            ax.set_xlabel(axis_units_label, fontsize=14)

    
    xticks = np.linspace(-dx+dx_shift/2, dx+dx_shift/2, n_pts_labels)
    yticks = np.linspace(-dy+dy_shift/2, dy+dy_shift/2, n_pts_labels)
    xticklabels = np.linspace(-(dx-dx_shift/2) * cell_size, +(dx+dx_shift/2) * cell_size, n_pts_labels)
    yticklabels = np.linspace(-(dy-dy_shift/2) * cell_size, +(dy+dy_shift/2) * cell_size, n_pts_labels)
    # if dx < 10:
    #     xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
    # else:
    #     xticklabels = ['{:.0f}'.format(xtick) for xtick in xticklabels]
    # if dy < 10:
    #     yticklabels = ['{:.2f}'.format(ytick) for ytick in yticklabels]
    # else:
    #     yticklabels = ['{:.0f}'.format(ytick) for ytick in yticklabels]

    if (projection =='offset') or (projection == 'celestial'):
        xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
        yticklabels = ['{:.2f}'.format(ytick) for ytick in yticklabels]
    else:
        xticklabels = ['{:.0f}'.format(xtick) for xtick in xticklabels]
        yticklabels = ['{:.0f}'.format(ytick) for ytick in yticklabels]

    ax.set_yticks(yticks, yticklabels)
    ax.set_xticks(xticks, xticklabels)
    ax.set_aspect('equal')

    ax.tick_params(axis='y', which='both', labelsize=16, color='black',
                   pad=5)
    ax.tick_params(axis='x', which='both', labelsize=16, color='black',
                   pad=5)
    # ax2_x.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    # ax3_y.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    if projection != 'celestial':
        ax.grid(which='both', axis='both', color='gray', linewidth=0.6,
                alpha=0.7)
        ax.grid(which='both', axis='both', color='gray', linewidth=0.6,
                alpha=0.7)
    else:
        ax.grid()

    ax.axis(show_axis)
    # ax.axis(show_axis)

    vmin = vmin_factor * std
    if extent is None:
        extent = [-dx+dx_shift, dx+dx_shift, -dy+dy_shift, dy+dy_shift]
    #     print(g)

    if vmax is not None:
        vmax = vmax
    else:
        if vmax_factor is not None:
            vmax = vmax_factor * g.max()
        else:
            vmax = 0.95 * g.max()


    norm0 = simple_norm(g, stretch='linear', max_percent=99.0)
    norm = simple_norm(g, stretch='sqrt', asinh_a=0.02, min_cut=vmin,
                       max_cut=vmax)
    
    if no_show:
        # # plot the first normalization (low level, transparent)
        im_plot = ax.imshow(g, origin='lower',aspect=aspect,
                            cmap='gray', norm=norm0, alpha=0.0, extent=extent)

        # cm = copy.copy(plt.cm.get_cmap(CM))
        # cm.set_under((0, 0, 0, 0))
        im_plot = ax.imshow((g), cmap=CM, origin='lower', alpha=0.0, extent=extent,
                            norm=norm,
                            aspect=aspect)  # ,vmax=vmax, vmin=vmin)#norm=norm

    else:
        # # plot the first normalization (low level, transparent)
        im_plot = ax.imshow(g, origin='lower',aspect=aspect,
                            cmap='gray', norm=norm0, alpha=0.5, extent=extent)

        # cm = copy.copy(plt.cm.get_cmap(CM))
        # cm.set_under((0, 0, 0, 0))
        im_plot = ax.imshow((g), cmap=CM, origin='lower', alpha=1.0, extent=extent,
                            norm=norm,
                            aspect=aspect)  # ,vmax=vmax, vmin=vmin)#norm=norm
        
    


    if add_contours:
        try:
            levels_g = np.geomspace(2.0 * g.max(), 5 * std, num_contours)
            levels_low = np.asarray([4 * std, 3 * std])
            levels_black = np.geomspace(vmin_factor * std + 0.00001, 2.5 * g.max(), num_contours)
            levels_neg = neg_levels * std
            levels_white = np.geomspace(g.max(), 0.1 * g.max(), num_contours)

            contour_palette_ = ['#000000', '#444444', '#666666', '#EEEEEE',
                               '#EEEEEE', '#FFFFFF']
            if '_r' in CM:
                contour_palette = contour_palette_
            else:
                contour_palette = contour_palette_[::-1]

            if no_show:
                contour = ax.contour(g, levels=levels_g[::-1],
                                    #  colors=contour_palette[::-1], #for dark cmap
                                    colors=contour_palette, #for light cmap
                                    # aspect=aspect,
                                    linewidths=1.2, extent=extent,
                                    alpha=0.0)

                contour = ax.contour(g, levels=levels_low[::-1],
                                    colors='brown',
                                    # aspect=aspect,
                                    # linestyles=['dashed', 'dashdot'],
                                    linewidths=1.0, extent=extent,
                                    alpha=0.0)

            else:
                contour = ax.contour(g, levels=levels_g[::-1],
                                    #  colors=contour_palette[::-1], #for dark cmap
                                    colors=contour_palette, #for light cmap
                                    # aspect=aspect,
                                    linewidths=1.2, extent=extent,
                                    alpha=1.0)

                contour = ax.contour(g, levels=levels_low[::-1],
                                    colors='brown',
                                    # aspect=aspect,
                                    # linestyles=['dashed', 'dashdot'],
                                    linewidths=1.0, extent=extent,
                                    alpha=1.0)
            # ax.clabel(contour, inline=1, fontsize=10)
        except:
            pass
        try:
            if no_show:
                ax.contour(g, levels=levels_neg[::-1], colors='k',
                        linewidths=1.0, extent=extent,
                        alpha=0.0)

            else:
                ax.contour(g, levels=levels_neg[::-1], colors='k',
                        linewidths=1.0, extent=extent,
                        alpha=1.0)
        except:
            pass
        
    if add_frequency:
        try:
            if freq_label is None:
                if isinstance(imagename, str) == True:
                    frequency = f"{(getfreqs([imagename])[0]/1e9):.1f} GHz"
            else:
                frequency = freq_label
            ax.annotate(frequency,
                        xy=(0.70, 0.82), xycoords='axes fraction',
                        fontsize=12,
                        bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                        color='red')
        except:
            pass
        
    
    if plot_colorbar:
        
        divider = make_axes_locatable(ax)
        
        try:
            # cb = plt.colorbar(mappable=plt.gca().images[0],
            #                   cax=fig.add_axes([0.90, 0.15, 0.05, 0.70]))

            # cb = plt.colorbar(im_plot, 
            #                   ax=ax,
            #                   cax=fig.add_axes([0.90, 0.15, 0.05, 0.70]))
            # cbar_ax = fig.add_axes([0.15, 0.85, 0.7, 0.05])  # [left, bottom, width, height]
            if cbar_orientation == 'horizontal':
                cax = divider.append_axes("top", size="7%", pad=0.05)
                cb = fig.colorbar(im_plot, 
                                #   pad=pad,
                                  cax=cax,fraction=0.046,
                                  orientation='horizontal')
            else:
                # cb = plt.colorbar(im_plot, 
                #                 ax=ax,
                #                 cax=fig.add_axes([0.90, 0.15, 0.05, 0.70])
                #                 )
                cax = divider.append_axes("right", size="7%", pad=0.05)
                cb = plt.colorbar(im_plot, 
                                ax=ax,fraction=0.046,
                                cax=cax
                                )
            # cb = plt.colorbar(im_plot, 
            #                   ax=ax,
            #                   orientation='horizontal', fraction=0.046, pad=0.04)

            if flux_units == 'Jy':
                cb.set_label(r"Flux Density [Jy/Beam]", labelpad=10, fontsize=16)
                cb.formatter = CustomFormatter(factor=1, useMathText=True)
                cb.update_ticks()
            if flux_units == 'mJy':
                cb.set_label(r"Flux Density [mJy/Beam]", labelpad=10, fontsize=16)
                cb.formatter = CustomFormatter(factor=1000, useMathText=True)
                cb.update_ticks()
            if flux_units == 'any':
                cb.set_label(r"Pixel Intensity", labelpad=10, fontsize=16)
                cb.formatter = CustomFormatter(factor=1, useMathText=True)
                cb.update_ticks()
            if flux_units == 'sfr':
                cb.set_label(r"${\rm M}_{\rm \odot} \ {\rm yr^{-1}} \ {\rm beam^{-1}}$", labelpad=10, fontsize=16)
                cb.formatter = CustomFormatter(factor=1, useMathText=True)
                cb.update_ticks()

            levels_colorbar2 = np.geomspace(1.0 * vmax, 3 * std,
                                            cbar_n_points)
            cb.set_ticks(levels_colorbar2)

            cb.ax.yaxis.set_tick_params(labelleft=True, labelright=False,
                                        tick1On=False, tick2On=False)
            cb.ax.yaxis.tick_right()
            # cb.set_ticks(levels_colorbar2)
            # cb.set_label(r'Flux [mJy/beam]', labelpad=10, fontsize=16)
            #         cbar.ax.xaxis.set_tick_params(pad=0.1,labelsize=10)
            cb.ax.tick_params(labelsize=16)
            cb.outline.set_linewidth(1)
            # cbar.dividers.set_color(None)

            # Make sure the color bar has ticks and labels at the top, since the bar is on the top as well.
            cb.ax.xaxis.set_ticks_position('top')
            cb.ax.xaxis.set_label_position('top')
        except:
            pass

    if plot_title is not None:
        ax.set_title(plot_title)

    if add_beam == True:

        if isinstance(imagename, str) == True:
            try:
                from matplotlib.patches import Ellipse
                imhd = imhead(imagename)
                a = imhd['restoringbeam']['major']['value']
                b = imhd['restoringbeam']['minor']['value']
                pa = imhd['restoringbeam']['positionangle']['value']
                if projection == 'px':
                    el = Ellipse((-(dx-dx_shift) * 0.85, 
                                  -(dy-dy_shift) * 0.85), 
                                 b, a, angle=pa,
                                 facecolor='black', alpha=1.0)
                else:
                    el = Ellipse((-(dx-dx_shift) * 0.85, 
                                  -(dy-dy_shift) * 0.85), 
                                 b / cell_size,
                                 a / cell_size,
                                 angle=pa, facecolor='black', alpha=1.0)

                ax.add_artist(el, )

                Oa = '{:.2f}'.format(a)
                Ob = '{:.2f}'.format(b)

                blabel_pos_y, blabel_pos_x = g.shape
                blabel_pos_x = (blabel_pos_x + dx+dx_shift)# * (dx/dy)
                blabel_pos_y = (blabel_pos_y + dy+dy_shift) 

                #         ax.annotate(r'$' + Oa +'\\times'+Ob+'$',
                #                     xy=(blabel_pos_x* 0.77, blabel_pos_y * 0.58), xycoords='data',
                #                     fontsize=15,bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                #                     color='red')
                ax.annotate(r"$" + Oa + "''\\times" + Ob + "''$",
                            xy=(0.63, 0.08), xycoords='axes fraction',
                            fontsize=12,
                            bbox=dict(boxstyle='round', facecolor='white',
                                      alpha=0.9),
                            color='red')

                el.set_clip_box(ax.bbox)
            except:
                print('Error adding beam.')

    if source_distance is not None:
        # try:
        ww.wcs.radesys = 'icrs'
        radesys = ww.wcs.radesys
        # distance = source_distance * u.Mpc
        distance = angular_distance_cosmo(source_distance)  # * u.Mpc
        #         scalebar_length = scalebar_length
        scalebar_loc = (0.99, 0.99)  # y, x
        left_side = coordinates.SkyCoord(
            *ww.celestial.wcs_pix2world(
                g.shape[0],
                g.shape[1],
                0) * u.deg,
            frame=radesys.lower())

        length = (scalebar_length / distance).to(u.arcsec,
                                                 u.dimensionless_angles())

        scale_bar_length_pixels = length.value / cell_size
        scale_bar_position = (-(dx-1*dx_shift) * 0.50, -(dy-1*dy_shift) * 0.9)

        ax.annotate('',
                    xy=(scale_bar_position[0] + scale_bar_length_pixels,
                        scale_bar_position[1]),
                    # xy=(0.1, 0.1), ##
                    xytext=scale_bar_position, arrowprops=dict(arrowstyle='-',
                                                               color='black', #navy
                                                               lw=5))

        ax.text(scale_bar_position[0] + scale_bar_length_pixels / 2,
                scale_bar_position[1] + scale_bar_length_pixels / 15,
                f'{scalebar_length}', fontsize=20,
                color='black', ha='center',weight='bold',
                va='bottom')
        # except:
        #     print('Error adding scalebar.')
    if text_annotation:
        ax.annotate(text_annotation, xy=text_position, xycoords='figure fraction',
                    fontsize=18, color='black', fontweight='bold')

    if save_name != None:
    #         if not os.path.exists(save_name+special_name+'.jpg'):
        plt.savefig(save_name + special_name + '.jpg', dpi=600,
                    bbox_inches='tight')
        # plt.savefig(save_name + special_name + '.pdf', dpi=600,
        #             bbox_inches='tight')
    # if no_show:
    #     plt.clf()
    #     plt.close()
    if return_extent:
        return ax,extent
    else:
        return ax

def plot_alpha_map(alphaimage,alphaimage_error,radio_map,frequencies,
                   mask_good_alpha = None,
                   levels_g = None,
                   vmin_factor = 3,neg_levels=np.asarray([-3]),
                   vmin=None,vmax=None,rms=None,figsize=(6, 6),
                   extent = None,crop=False,centre=None,
                   plot_title='',label_colorbar='',
                   cmap='magma_r',n_contours = 8,
                   projection='offset',cell_size = None,
                   box_size=None,plot_alpha_error_points=False,
                   save_name=None):

    from matplotlib.colors import ListedColormap
    
    fig = plt.figure(figsize=figsize)
    

    # fig, ax = plt.subplots(figsize=figsize)

    if mask_good_alpha is None:
        mask_good_alpha = np.ones_like(alphaimage)

    if isinstance(radio_map, str) == True:
        g = load_fits_data(radio_map)
    else:
        g = radio_map

    if crop:
        if centre is None:
            centre = (int(alphaimage.shape[0]/2),int(alphaimage.shape[1]/2))
        else:
            centre = centre
        
        
        mask_good_alpha = do_cutout_2D(mask_good_alpha, 
                                       box_size=box_size,
                                       center=centre, 
                                       return_='data')
            
        _alphaimage = do_cutout_2D(alphaimage, 
                                   box_size=box_size, 
                                   center=centre, 
                                   return_='data')

        _alphaimage_error = do_cutout_2D(alphaimage_error, 
                                         box_size=box_size, 
                                          center=(centre[0],centre[1]), 
                                          return_='data')
        
        
        
        _g = do_cutout_2D(g, 
                          box_size=box_size, 
                          center=centre, 
                          return_='data')
        
        

    else:
        _alphaimage = alphaimage.copy()
        _alphaimage_error = alphaimage_error.copy()
        _g = g


    # _alphaimage[~mask_good_alpha] = np.nan
    # _alphaimage_error[~mask_good_alpha] = np.nan
    # _g[~mask_good_alpha] = np.nan

    # Access the 'magma' colormap
    base_cmap = plt.colormaps[cmap]

    # Create a new colormap with 3 colors from the base colormap
    new_colors = base_cmap(np.linspace(0, 1, 10))
    custom_cmap = ListedColormap(new_colors)
    # custom_cmap = 'magma_r'


    if vmin is None:
        vmin = np.nanmin(_alphaimage) - 1*abs(np.nanmedian(_alphaimage))
    else:
        vmin = vmin
    if vmax is None:
        vmax = np.nanmax(_alphaimage) + 1*abs(np.nanmedian(_alphaimage))
    else:
        vmax = vmax
    
    if projection == 'offset':
        ax,extent = eimshow(radio_map,return_extent = True,figsize=figsize,
                        add_contours=False,CM='magma',no_show=True,
                            crop=crop,center=centre,box_size=box_size)
        if cell_size is None:
            if isinstance(radio_map, str):
                cell_size = get_cell_size(radio_map)
            else:
                raise ValueError(f"Please, if selecting offset projection mode," 
                                 f"provide the cell_size of the image.")
    else:
        ax = fig.add_subplot(111)  # 111 means 1x1 grid, first subplot
        cell_size = 1
    
    img = plt.imshow(_alphaimage, origin='lower', 
                     vmin=vmin, 
                     extent=extent,
    #                  vmax=0.0,
                     vmax=vmax, 
                     cmap=custom_cmap,
    #                  cmap='plasma'
                    )
    if rms is None:
        std = mad_std(_g)
    else:
        std = rms 
    
    # print(extent)
    if levels_g is None:
        contour_palette = ['#000000', '#444444', '#666666', '#EEEEEE',
                        '#EEEEEE', '#FFFFFF']
    else:
        contour_palette = 'darkgray' #darkkhaki, tan, thistle, olive

    
    if levels_g is None:
        levels_g = np.geomspace(3.0 * _g.max(), 5 * std, n_contours)
    else:
        levels_g = levels_g
    

    # levels_low = np.asarray([4 * std, 3 * std])
    # levels_black = np.geomspace(vmin_factor * std + 0.00001, 2.5 * _g.max(), 6)
    # levels_neg = neg_levels * std
    # levels_white = np.geomspace(_g.max(), 0.1 * _g.max(), 6)



    contour = ax.contour(_g, levels=levels_g[::-1],
                         colors=contour_palette,
                          origin='lower', 
                        #   cmap='Greys',
                         # aspect=aspect,
                         linewidths=1.2, extent=extent,
                         alpha=1.0)

    if plot_alpha_error_points:
        # Downsample the data for scatter plot
        downscale_factor = 10  # Adjust this factor as needed
        x_downsampled, y_downsampled = np.meshgrid(
            np.arange(0, _alphaimage.shape[1], downscale_factor),
            np.arange(0, _alphaimage.shape[0], downscale_factor)
        )
        
        marker_sizes = 50 * abs(_alphaimage_error[y_downsampled, x_downsampled]/_alphaimage[y_downsampled, x_downsampled])
        min_error_pos = np.where(marker_sizes==np.nanmin(marker_sizes))
        max_error_pos = np.where(marker_sizes==np.nanmax(marker_sizes))
        # print(x_downsampled*cell_size,y_downsampled*cell_size)
        sc = ax.scatter(x_downsampled*cell_size, y_downsampled*cell_size, s=marker_sizes, 
                     c='cyan', alpha=0.9)
        _ = ax.scatter(x_downsampled[min_error_pos]*cell_size, 
                        y_downsampled[min_error_pos]*cell_size, 
                        s=marker_sizes[min_error_pos], 
                         c='black', alpha=0.8,
                        label=fr'Min err $\alpha$ = {(marker_sizes[min_error_pos][0]/50):.2f}')
        _ = ax.scatter(x_downsampled[max_error_pos]*cell_size, 
                        y_downsampled[max_error_pos]*cell_size, 
                        s=marker_sizes[max_error_pos], 
                        c='black', alpha=0.8,
        #                 marker='x',
                        label=fr'Max err $\alpha$ = {(marker_sizes[max_error_pos][0]/50):.2f}')
        plt.legend()
    if label_colorbar == '':
        label_colorbar = fr"Spectral Index $\alpha$ [{(frequencies[0]/1e9):.1f} ~ {(frequencies[-1]/1e9):.1f} GHz]"
    else:
        _label_colorbar = fr"[{(frequencies[0]/1e9):.1f} ~ {(frequencies[-1]/1e9):.1f} GHz]"
        label_colorbar = label_colorbar + _label_colorbar

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%",
                              pad=0.05)  # Adjust 'pad' for space and 'size' for width
    cbar = fig.colorbar(img, cax=cax)
    cbar.set_label(label=label_colorbar)
    # cbar.formatter = CustomFormatter(factor=1000, useMathText=True)
    # cbar.update_ticks()
    
    
    # cbar = plt.colorbar(img, 
    #                     label=label_colorbar)
    ax.grid(which='both', axis='both', color='gray', linewidth=0.6,
                alpha=0.7)
    
    ax.set_title(f"{plot_title} ({len(frequencies)} images)")
    # plt.savefig('alpha_map_and_errors_MCG12_S_C_X_new_fit_in_linear.pdf',dpi=300,bbox_inches='tight')
    # plt.show()
    if save_name != None:
        plt.savefig(save_name, dpi=300,
                    bbox_inches='tight')
        
    pass
    

def add_beam_to_image(imagename, ax, dx, cell_size):
    if isinstance(imagename, str) == True:
        try:
            from matplotlib.patches import Ellipse

            imhd = imhead(imagename)
            a = imhd['restoringbeam']['major']['value']
            b = imhd['restoringbeam']['minor']['value']
            pa = imhd['restoringbeam']['positionangle']['value']
            el = Ellipse((-dx * 0.85, -dx * 0.85), b / cell_size, a / cell_size,
                         angle=pa, facecolor='black', alpha=1.0)

            ax.add_artist(el, )

            Oa = '{:.2f}'.format(a)
            Ob = '{:.2f}'.format(b)

            blabel_pos_x, blabel_pos_y = 2*dx, 2*dx
            blabel_pos_x = blabel_pos_x + dx
            blabel_pos_y = blabel_pos_y + dx

            #         ax.annotate(r'$' + Oa +'\\times'+Ob+'$',
            #                     xy=(blabel_pos_x* 0.77, blabel_pos_y * 0.58), xycoords='data',
            #                     fontsize=15,bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
            #                     color='red')
            ax.annotate(r'$' + Oa + '\\times' + Ob + '$',
                        xy=(0.60, 0.06), xycoords='axes fraction',
                        fontsize=15, bbox=dict(boxstyle='round', facecolor='white',
                                               alpha=0.9),
                        color='red')

            el.set_clip_box(ax.bbox)
        except:
            print('Error adding beam.')
    return ax


def plot_image(image, residual_name=None, box_size=200, box_size_inset=60,
               center=None, rms=None, add_inset='auto', add_beam=True,
               vmin_factor=3, max_percent_lowlevel=99.0,
               levels_neg_factors=np.asarray([-3]),
               max_percent_highlevel=99.9999, vmax=None,
               do_cut=False, CM='magma_r', cbar_axes=[0.03, 0.11, 0.05, 0.77],
               # cbar_axes=[0.9, 0.15, 0.04, 0.7],
               source_distance=None, save_name=None, special_name='',
               scalebar_length=1.0 * u.kpc,
               show_axis='on', plot_color_bar=True, figsize=(5, 5),
               cmap_cont='terrain',
               projection='offset'):
    import cmasher as cmr
    import astropy.io.fits as fits
    from astropy import coordinates
    import matplotlib.ticker as mticker
    from matplotlib.ticker import FuncFormatter, FormatStrFormatter

    hdu = fits.open(image)
    ww = WCS(hdu[0].header)
    #     ww.wcs.ctype = [ 'XOFFSET' , 'YOFFSET' ]
    imhd = imhead(image)
    if do_cut == True:
        if center is None:
            st = imstat(image)
            # print('  >> Center --> ', st['maxpos'])
            yin, yen, xin, xen = st['maxpos'][0] - box_size, st['maxpos'][
                0] + box_size, st['maxpos'][1] - box_size, \
                                 st['maxpos'][1] + box_size
        else:
            yin, yen, xin, xen = center[0] - box_size, center[0] + box_size, \
                                 center[1] - box_size, center[1] + box_size
    else:
        yin, yen, xin, xen = 0, -1, 0, -1

    #     cmr is a package for additional density plot maps, if you would like to use
    #        >> https://cmasher.readthedocs.io/user/sequential/rainforest.html#rainforest
    #     CM = 'cmr.horizon'
    #     cm = cmr.neon

    fontsize = 12
    tick_fontsize = 12
    fig = plt.figure(figsize=figsize)
    scalling = 1

    # ax = fig.add_subplot(projection=ww.celestial[cutout,cutout])

    pixel_scale = (ww.pixel_scale_matrix[1, 1] * 3600)

    # if projection == 'celestial':
    ax = fig.add_subplot(projection=ww.celestial)

    #     xoffset_in = (xen_cut-xin_cut) * pixel_scale / 2
    #     yoffset_in = (yen_cut-yin_cut) * pixel_scale / 2

    extent = [xin, xen, yin, yen]
    # improve image visualization with normalizations
    norm = visualization.simple_norm(
        hdu[0].data.squeeze()[xin:xen, yin:yen] * scalling, stretch='linear',
        max_percent=max_percent_lowlevel)
    # plot the first normalization (low level, transparent)
    im = ax.imshow(hdu[0].data.squeeze()[xin:xen, yin:yen], origin='lower',
                   cmap=CM, norm=norm, alpha=0.2, extent=extent)

    cm = copy.copy(plt.cm.get_cmap(CM))
    cm.set_under((0, 0, 0, 0))

    data_range = hdu[0].data.squeeze()[xin:xen, yin:yen] * scalling

    # set min and max levels
    if rms != None:  # in case you want to set the min level manually.
        vmin = vmin_factor * rms
        std = rms
    else:
        if mad_std(data_range) == 0:
            std = data_range.std()
        else:
            std = mad_std(data_range)
        vmin = vmin_factor * std
    if vmax is not None:
        vmax = vmax
    else:
        vmax = 1.0 * data_range.max()

    # set second normalization
    norm2 = simple_norm(abs(hdu[0].data.squeeze())[xin:xen, yin:yen] * scalling,
                        min_cut=vmin,
                        max_cut=vmax, stretch='asinh',
                        asinh_a=0.05)  # , max_percent=max_percent_highlevel)
    norm2.vmin = vmin

    # this is what is actually shown on the final plot, better contrast.
    im = ax.imshow(hdu[0].data.squeeze()[xin:xen, yin:yen] * scalling,
                   origin='lower',
                   norm=norm2, cmap=cm, aspect='auto', extent=extent)

    levels_colorbar = np.geomspace(2.0 * vmax, 5 * std,
                                   8)  # draw contours only until 5xstd level.
    levels_neg = np.asarray(levels_neg_factors * std)
    levels_low = np.asarray([4 * std, 3 * std])
    #     levels_neg = np.asarray([])
    #     levels_low = np.asarray([])

    #     levels_colorbar = np.append(levels_neg,levels_pos)

    levels_colorbar2 = np.geomspace(1.0 * vmax, 3 * std,
                                    5)  # draw contours only until 5xstd level.
    ax.contour(hdu[0].data.squeeze()[xin:xen, yin:yen], extent=extent,
               levels=levels_colorbar[::-1], cmap=cmap_cont,
               linewidths=1.25)  # ,alpha=0.6)

    ax.contour(hdu[0].data.squeeze()[xin:xen, yin:yen], levels=levels_low[::-1],
               colors='brown', linestyles=['dashed', 'dashdot'], extent=extent,
               linewidths=1.0, alpha=0.5)
    ax.contour(hdu[0].data.squeeze()[xin:xen, yin:yen], levels=levels_neg[::-1],
               colors='gray', extent=extent,
               linewidths=1.5, alpha=1.0)

    ww.wcs.radesys = 'icrs'
    radesys = ww.wcs.radesys
    if source_distance is not None:
        # distance = source_distance * u.Mpc
        distance = angular_distance_cosmo(source_distance)  # * u.Mpc
        img = hdu[0].data.squeeze()[xin:xen, yin:yen]
        #         scalebar_length = scalebar_length
        scalebar_loc = (0.82, 0.1)  # y, x
        left_side = coordinates.SkyCoord(
            *ww.celestial[xin:xen, yin:yen].wcs_pix2world(
                scalebar_loc[1] * img.shape[1],
                scalebar_loc[0] * img.shape[0],
                0) * u.deg,
            frame=radesys.lower())

        length = (scalebar_length / distance).to(u.arcsec,
                                                 u.dimensionless_angles())
        make_scalebar(ax, left_side, length, color='red', linestyle='-',
                      label=f'{scalebar_length:0.1f}',
                      text_offset=0.1 * u.arcsec, fontsize=24)

    #     _ = ax.set_xlabel(f"Right Ascension {radesys}")
    #     _ = ax.set_ylabel(f"Declination {radesys}")
    _ = ax.set_xlabel(f"Right Ascension")
    _ = ax.set_ylabel(f"Declination")

    if projection == 'celestial':
        ra = ax.coords['ra']
        ra.set_major_formatter('hh:mm:ss.s')
        dec = ax.coords['dec']
        #         ra.set_axislabel(f"RA ({radesys})", fontsize=fontsize)
        #         dec.set_axislabel(f"Dec ({radesys})", fontsize=fontsize, minpad=0.0)
        ra.set_axislabel(f"RA", fontsize=fontsize)
        dec.set_axislabel(f"Dec", fontsize=fontsize, minpad=0.0)

        ra.ticklabels.set_fontsize(tick_fontsize)
        ra.set_ticklabel(exclude_overlapping=True)
        dec.ticklabels.set_fontsize(tick_fontsize)
        dec.set_ticklabel(exclude_overlapping=True)
        ax.tick_params(axis="y", direction="in", pad=-25)
        ax.tick_params(axis="x", direction="in", pad=-25)
        ax.axis(show_axis)
    if projection == 'offset':
        """
        This is a workaround to set the axes in terms of offsets [arcsec] 
        relative to the image center. I found this to work because I must use 
        the wcs coordinates (projection=celestial) if the scale bar is ploted. 
        So this uses the coordinates, but then the axes are removed and on top of 
        that the relative offsets are added.   
        However, this assumes that the source is at the center of the image. 
        """
        ax.axis('off')
        xoffset = (hdu[0].data.shape[0] - xin * 2) * pixel_scale / 2
        yoffset = (hdu[0].data.shape[1] - yin * 2) * pixel_scale / 2
        ax2_x = ax.twinx()
        ax3_y = ax2_x.twiny()
        ax2_x.set_ylabel('Offset [arcsec]', fontsize=14)
        ax3_y.set_xlabel('Offset [arcsec]', fontsize=14)
        ax2_x.yaxis.set_ticks(np.linspace(-xoffset, xoffset, 7))
        ax3_y.xaxis.set_ticks(np.linspace(-yoffset, yoffset, 7))

        ax2_x.tick_params(axis='y', which='both', labelsize=16, color='black',
                          pad=-30)
        ax3_y.tick_params(axis='x', which='both', labelsize=16, color='black',
                          pad=-25)
        ax2_x.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax3_y.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        ax2_x.grid(which='both', axis='both', color='gray', linewidth=0.6,
                   alpha=0.7)
        ax3_y.grid(which='both', axis='both', color='gray', linewidth=0.6,
                   alpha=0.7)
        ax2_x.axis(show_axis)
        ax3_y.axis(show_axis)

    freq = '{:.2f}'.format(imhd['refval'][2] / 1e9)
    label_pos_x, label_pos_y = hdu[0].data.squeeze()[xin:xen, yin:yen].shape
    label_pos_x = label_pos_x + xin
    label_pos_y = label_pos_y + yin
    #     ax.annotate(r'' + freq + 'GHz', xy=(0.35, 0.05),  xycoords='figure fraction', fontsize=14,
    #                 color='red')
    ax.annotate(r'' + freq + 'GHz',
                xy=(0.55, 0.82), xycoords='axes fraction',
                fontsize=18,
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                color='red')

    if plot_color_bar == True:
        def format_func(value, tick_number, scale_density='mJy'):
            # Use the custom formatter for the colorbar ticks
            mantissa = value * 1000
            return r"${:.1f}$".format(mantissa)

        cax = fig.add_axes(cbar_axes)

        cbar = fig.colorbar(im, cax=cax, orientation='vertical',
                            shrink=1, aspect='auto', pad=10, fraction=1.0,
                            drawedges=False, ticklocation='left')
        #         cbar.ax.yaxis.set_major_formatter(FuncFormatter(format_func))
        #         cbar.formatter = mticker.ScalarFormatter(useMathText=True)
        cbar.formatter = CustomFormatter(factor=1000, useMathText=True)
        cbar.update_ticks()

        #         cbar.update_ticks()
        # cbar_axes = [0.12, 0.95, 0.78, 0.06]
        # cax = fig.add_axes(cbar_axes)
        #
        # cbar = fig.colorbar(im, cax=cax, orientation='horizontal', format='%.0e',
        #                     shrink=1.0, aspect='auto', pad=0.1, fraction=1.0, drawedges=False
        #                     )
        cbar.ax.yaxis.set_tick_params(labelleft=True, labelright=False,
                                      tick1On=False, tick2On=False)
        cbar.ax.yaxis.tick_left()
        cbar.set_ticks(levels_colorbar2)
        cbar.set_label(r'Flux [mJy/beam]', labelpad=10, fontsize=16)
        #         cbar.ax.xaxis.set_tick_params(pad=0.1,labelsize=10)
        cbar.ax.tick_params(labelsize=16)
        cbar.outline.set_linewidth(1)
        # cbar.dividers.set_color(None)

        # Make sure the color bar has ticks and labels at the top, since the bar is on the top as well.
        cbar.ax.xaxis.set_ticks_position('top')
        cbar.ax.xaxis.set_label_position('top')

    imhd = imhead(image)
    if add_inset == 'auto':
        a = imhd['restoringbeam']['major']['value']
        if a / pixel_scale <= 15:
            add_inset = True
        else:
            add_inset = False

    # add residual inset
    if residual_name is not None:
        residual_data = pf.getdata(residual_name)
        axins_re = inset_axes(ax, width="40%", height="40%", loc='lower right',
                              bbox_to_anchor=(0.00, -0.10, 1.0, 1.0),
                              bbox_transform=ax.transAxes)
        xoffset = hdu[0].data.shape[0] * pixel_scale / 2
        yoffset = hdu[0].data.shape[1] * pixel_scale / 2
        extent = [-xoffset, xoffset, -yoffset, yoffset]
        vmin_re = np.min(residual_data)
        vmax_re = np.max(residual_data)

        norm_re = visualization.simple_norm(residual_data * scalling,
                                            stretch='linear',
                                            max_percent=max_percent_lowlevel)
        # im = ax.imshow(fh[0].data.squeeze()[cutout,cutout], cmap='gray_r',norm=norm)

        #     if projection is 'offset':
        imre = axins_re.imshow(residual_data, origin='lower',
                               cmap=CM, norm=norm_re, alpha=0.2, extent=extent)

        norm_re2 = simple_norm(residual_data * scalling, min_cut=vmin_re,
                               max_cut=vmax_re, stretch='linear',
                               asinh_a=0.05)  # , max_percent=max_percent_highlevel)

        imre2 = axins_re.imshow(residual_data * scalling, origin='lower',
                                norm=norm_re2, cmap=CM, aspect='auto',
                                extent=extent)
        # axins_re.contour(residual_data, colors='k',extent=extent,
        #               linewidths=1.0, alpha=1.0)
        #         axins_re.imshow(residual_data, cmap=CM,extent= [-xoffset, xoffset, -yoffset, yoffset],origin='lower')
        axins_re.set_ylabel('', fontsize=10)
        axins_re.set_xlabel('', fontsize=10)
        axins_re.xaxis.set_ticks(np.linspace(-xoffset, xoffset, 4))
        axins_re.yaxis.set_ticks([])
        axins_re.tick_params(axis='y', which='both', labelsize=11, color='black')
        axins_re.tick_params(axis='x', which='both', labelsize=11, color='black')
        axins_re.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axins_re.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        #         cbr_re = plt.colorbar(imre2)
        axins_re.axis(show_axis)
    #         axins_re.axis('off')

    if add_inset == True:
        st = imstat(image)
        print('  >> Center --> ', st['maxpos'])
        # sub region of the original image
        yin_cut, yen_cut, xin_cut, xen_cut = st['maxpos'][0] - box_size_inset, \
                                             st['maxpos'][0] + box_size_inset, \
                                             st['maxpos'][1] - box_size_inset, \
                                             st['maxpos'][1] + box_size_inset

        axins = inset_axes(ax, width="40%", height="40%", loc='lower left',
                           bbox_to_anchor=(0.05, -0.05, 1.0, 1.0),
                           bbox_transform=ax.transAxes)

        xoffset_in = (xen_cut - xin_cut) * pixel_scale / 2
        yoffset_in = (yen_cut - yin_cut) * pixel_scale / 2

        extent_inset = [-xoffset_in, xoffset_in, -yoffset_in, yoffset_in]

        Z2 = hdu[0].data.squeeze()[xin_cut:xen_cut, yin_cut:yen_cut]

        vmax_inset = np.max(Z2)
        vmin_inset = 1 * np.std(Z2)

        norm_inset = visualization.simple_norm(Z2, stretch='linear',
                                               max_percent=max_percent_lowlevel)

        norm2_inset = simple_norm(abs(Z2), min_cut=vmin,
                                  max_cut=vmax_inset, stretch='asinh',
                                  asinh_a=0.008)  # , max_percent=max_percent_highlevel)
        norm2_inset.vmin = vmin

        axins.imshow(Z2, cmap=CM, norm=norm_inset, alpha=0.2,
                     extent=extent_inset,
                     aspect='auto', origin='lower')
        axins.imshow(Z2, norm=norm2_inset, extent=extent_inset,
                     cmap=cm,
                     origin="lower", alpha=1.0, aspect='auto')

        axins.set_ylabel('', fontsize=10)
        axins.set_xlabel('', fontsize=10)
        axins.xaxis.set_ticks(np.linspace(-xoffset_in, xoffset_in, 4))
        axins.yaxis.set_ticks(np.linspace(-yoffset_in, yoffset_in, 4))
        axins.tick_params(axis='y', which='both', labelsize=14, color='black')
        axins.tick_params(axis='x', which='both', labelsize=14, color='black')
        axins.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
        axins.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))

        axins.grid(True, alpha=0.9)
        levels_inset = np.geomspace(3.0 * Z2.max(), 5 * std,
                                    6)  # draw contours only until 5xstd level.
        levels_inset_neg = np.asarray([-3 * std])
        levels_inset_low = np.asarray([4 * std, 3 * std])

        csi = axins.contour(Z2, levels=levels_colorbar[::-1], cmap=cmap_cont,
                            extent=extent_inset,
                            linewidths=1.25, alpha=1.0)
        #         axins.clabel(csi, inline=False, fontsize=8, manual=False, zorder=99)
        axins.contour(Z2, levels=levels_inset_low[::-1], colors='brown',
                      extent=extent_inset, linestyles=['dashed', 'dashdot'],
                      linewidths=1.0, alpha=1.0)
        axins.contour(Z2, levels=levels_inset_neg[::-1], colors='gray',
                      extent=extent_inset,
                      linewidths=1.5, alpha=1.0)

        # inset zoom box is not working properly, so I am doing it manually.
        import matplotlib.patches as patches
        from matplotlib.patches import FancyArrowPatch
        rect = patches.Rectangle((yin_cut, xin_cut),
                                 xen_cut - xin_cut, yen_cut - yin_cut,
                                 linewidth=2, edgecolor='black',
                                 facecolor='none', alpha=0.5)
        ax.add_patch(rect)
        arrow_start = (yin_cut / 2, xin_cut / 2)
        arrow_end = (yin_cut / 2 + 50, xin_cut / 2 + 50)
        arrow = FancyArrowPatch(arrow_start, arrow_end, arrowstyle='->',
                                mutation_scale=24, linewidth=1.5, color='black')
        ax.add_patch(arrow)

        # manually add lines to connect the rectangle and arrow
        x1, y1 = rect.get_xy()
        x2, y2 = arrow_end
        ax.plot([x1 + 10, x2 + 100], [y1 + 10, y2 + 100], linestyle='--',
                color='black')
        #         ax.indicate_inset_zoom(axins)
        axins.axis(show_axis)

        if add_beam == True:
            from matplotlib.patches import Ellipse
            a_s = imhd['restoringbeam']['major']['value']
            b_s = imhd['restoringbeam']['minor']['value']
            pa_s = imhd['restoringbeam']['positionangle']['value']
            el_s = Ellipse((xoffset_in + 10, yoffset_in + 10), b_s,
                           a_s, angle=pa_s,
                           facecolor='r', alpha=0.5)
            axins.add_artist(el_s)
            el_s.set_clip_box(axins.bbox)

    #         axins.set_xlim(x1, x2)
    #         axins.set_ylim(y1, y2)
    #         axins.set_xticklabels([])
    #         axins.set_yticklabels([])
    #         ax.indicate_inset_zoom(axins, edgecolor="black")

    add_inset2 = False
    if add_inset2 == True:
        #         ax.axis(show_axis)

        st = imstat(image)
        print('  >> Center --> ', st['maxpos'])
        # sub region of the original image
        xin_cut, xen_cut, yin_cut, yen_cut = st['maxpos'][0] - box_size_inset, \
                                             st['maxpos'][0] + box_size_inset, \
                                             st['maxpos'][1] - box_size_inset, \
                                             st['maxpos'][1] + box_size_inset

        print(hdu[0].data.squeeze()[xin_cut:xen_cut, yin_cut:yen_cut].shape)
        Z2 = hdu[0].data.squeeze()[xin_cut:xen_cut, yin_cut:yen_cut]

        vmax_inset = np.max(Z2)
        vmin_inset = 1 * np.std(Z2)

        axins = inset_axes(ax, width="40%", height="40%", loc='lower left',
                           bbox_to_anchor=(0.00, -0.10, 1.0, 1.0),
                           bbox_transform=ax.transAxes)
        # axins = ax.inset_axes([0.00, 0.1, 0.35, 0.4], transform=ax.transData)
        # axins = ax.inset_axes([0.00, 0.1, 0.35, 0.4], transform=ax.transAxes)
        x1, x2, y1, y2 = xin_cut, xen_cut, yin_cut, yen_cut
        print(xen_cut - xin_cut)
        extent = np.asarray([xin_cut,
                             xen_cut,
                             yin_cut,
                             yen_cut
                             ]) * pixel_scale / 2

        xoffset_in = (xen_cut - xin_cut)  # * pixel_scale / 2
        yoffset_in = (yen_cut - yin_cut)  # * pixel_scale / 2
        #         xoffset_in = (xin_cut - xen_cut)
        #         yoffset_in = (yin_cut - yen_cut)
        # extent_arcsec = [-xoffset_in, xoffset_in, -yoffset_in, yoffset_in]
        extent_arcsec = extent
        #         extent_arcsec= [xin_cut* pixel_scale, xen_cut* pixel_scale,
        #                         yin_cut* pixel_scale, yen_cut* pixel_scale]
        norm_inset = visualization.simple_norm(Z2, stretch='linear',
                                               max_percent=max_percent_lowlevel)

        norm2_inset = simple_norm(abs(Z2), min_cut=vmin,
                                  max_cut=vmax_inset, stretch='asinh',
                                  asinh_a=0.05)  # , max_percent=max_percent_highlevel)
        norm2_inset.vmin = vmin

        axins.imshow(Z2, cmap=CM, norm=norm_inset, alpha=0.2,
                     extent=extent_arcsec,
                     aspect='auto', origin='lower')
        axins.imshow(Z2, norm=norm2_inset, extent=extent_arcsec,
                     cmap=cm,
                     origin="lower", alpha=1.0, aspect='auto')

        axins.grid(True, alpha=0.3)
        levels_inset = np.geomspace(3.0 * Z2.max(), 5 * std,
                                    8)  # draw contours only until 5xstd level.
        levels_inset_neg = np.asarray([-3 * std])

        csi = axins.contour(Z2, levels=levels_colorbar[::-1], colors='grey',
                            extent=extent_arcsec,
                            linewidths=1.0, alpha=1.0)
        axins.clabel(csi, inline=False, fontsize=8, manual=False, zorder=99)
        axins.contour(Z2, levels=levels_inset_neg[::-1], colors='k',
                      extent=extent_arcsec,
                      linewidths=1.0, alpha=1.0)

        if add_beam == True:
            from matplotlib.patches import Ellipse
            a_s = imhd['restoringbeam']['major']['value']
            b_s = imhd['restoringbeam']['minor']['value']
            pa_s = imhd['restoringbeam']['positionangle']['value']
            el_s = Ellipse((xin_cut + 10, yin_cut + 10), b_s / pixel_scale,
                           a_s / pixel_scale, angle=pa_s,
                           facecolor='r', alpha=0.5)
            axins.add_artist(el_s)
            el_s.set_clip_box(axins.bbox)

        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticklabels([])
        axins.set_yticklabels([])

        ax.indicate_inset_zoom(axins, edgecolor="black")

    if add_beam == True:
        from matplotlib.patches import Ellipse
        a = imhd['restoringbeam']['major']['value']
        b = imhd['restoringbeam']['minor']['value']
        pa = imhd['restoringbeam']['positionangle']['value']
        el = Ellipse((hdu[0].data.shape[0] - 60 - xin, 60 + xin),
                     b / pixel_scale, a / pixel_scale, angle=pa,
                     facecolor='black', alpha=1.0)
        ax.add_artist(el, )

        Oa = '{:.2f}'.format(a)
        Ob = '{:.2f}'.format(b)

        blabel_pos_x, blabel_pos_y = hdu[0].data.squeeze()[xin:xen,
                                     yin:yen].shape
        blabel_pos_x = blabel_pos_x + xin
        blabel_pos_y = blabel_pos_y + yin

        #         ax.annotate(r'$' + Oa +'\\times'+Ob+'$',
        #                     xy=(blabel_pos_x* 0.77, blabel_pos_y * 0.58), xycoords='data',
        #                     fontsize=15,bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
        #                     color='red')
        ax.annotate(r'$' + Oa + '\\times' + Ob + '$',
                    xy=(0.55, 0.12), xycoords='axes fraction',
                    fontsize=15,
                    bbox=dict(boxstyle='round', facecolor='white', alpha=0.9),
                    color='red')

        el.set_clip_box(ax.bbox)

    #     plt.tight_layout()
    if save_name != None:
        #         if not os.path.exists(save_name+special_name+'.jpg'):
        plt.savefig(save_name + special_name + '.jpg', dpi=300,
                    bbox_inches='tight')
        plt.savefig(save_name + special_name + '.pdf', dpi=600,
                    bbox_inches='tight')
    #         else:
    #             print('Skiping save')
    #     plt.show()
    return (plt, ax)


def plot_compare_conv_deconv(conv_image, deconv_image):
    # Create figure and subplots
    fig, ax = plt.subplots(figsize=(6, 3))
    img_size_x = conv_image.shape[0]
    img_size_y = conv_image.shape[1]

    # Plot image on left side
    im1 = ax.imshow(deconv_image, cmap='viridis', extent=[-500, 0, -500, 500], aspect='auto')
    ax.contour(deconv_image, extent=[-500, 0, -500, 500],
               cmap='Greys')  # ,levels=np.geomspace(mad_std(Z1),Z1.max(),3),extent=[-5, 0, -5, 5])

    # Plot image on right side
    im2 = ax.imshow(conv_image, cmap='plasma', extent=[0, 500, -500, 500], aspect='auto')
    ax.contour(conv_image, extent=[0, 500, -500, 500],
               cmap='Greys')  # ,levels=np.geomspace(mad_std(Z1),Z1.max(),3),extent=[-5, 0, -5, 5])
    # # Draw diagonal line
    # line = plt.Line2D([0, 5], [-5, 0], color='white', linewidth=2, linestyle='--')
    # ax.add_line(line)

    # Set axis limits and labels
    ax.set_xlim([-500, 500])
    ax.set_ylim([-500, 500])
    ax.set_xlabel('X')
    ax.set_ylabel('Y')

    # # Add colorbars
    # cbar1 = plt.colorbar(im1, ax=ax, shrink=0.7)
    # cbar1.set_label('Z1')
    # cbar2 = plt.colorbar(im2, ax=ax, shrink=0.7)
    # cbar2.set_label('Z2')

    plt.show()





def add_ellipse(ax, x0, y0, d_r, q, PA, label=None, show_center=True,
                label_offset=5, **kwargs):
    """
    Add an elliptical aperture to a plot with a distance label.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Axes object to which the ellipse should be added.
    x0 : float
        The x-coordinate of the center of the ellipse.
    y0 : float
        The y-coordinate of the center of the ellipse.
    d_r : float
        The radial distance from the center of the ellipse to its edge.
    q : float
        The axis ratio of the ellipse (b/a).
    PA : float
        The position angle of the major axis of the ellipse, in degrees.
    label : str or None, optional
        The label to use for the distance. If None, no label will be added.
    label_offset : float, optional
        The distance between the label and the ellipse edge, in pixels.
    **kwargs : optional
        Additional arguments that are passed to the Ellipse constructor.
    """
    a = d_r / (1 - q ** 2) ** 0.5
    # a = d_r
    b = a * q
    theta = np.deg2rad(PA)
    ellipse = Ellipse(xy=(x0, y0), width=a * 2, height=b * 2, angle=PA,
                      linestyle='-.', **kwargs)
    ax.add_artist(ellipse)
    if show_center:
        ax.plot(x0, y0, '+', color='white', ms=10)

    # Add distance label
    if label is not None:
        dx = label_offset * np.cos(theta)
        dy = label_offset * np.sin(theta)
        label_x, label_y = x0 + (a + label_offset) * np.cos(theta), y0 + (
                    b + label_offset) * np.sin(theta)
        print(dx, dy)
        ax.annotate(label, xy=(label_x, label_y),
                    xytext=(label_x + dx, label_y + dy),
                    #                     ha='center', va='center',
                    fontsize=14
                    #                     rotation=PA,
                    #                     textcoords='offset pixels'
                    )


def eimshow_with_cgrid(image, pix_to_pc, rms=None, dx=400, dy=400, do_cut=False, data_2D=None,
                       vmin_factor=3.0, vmax_factor=0.1, add_contours=True,
                       neg_levels=np.asarray([-3]),
                       centre=None, apply_mask=False, sigma_mask=6, CM='magma_r',
                       cmap_cont='terrain',
                       dilation_size=None, iterations=2,
                       r_d=200, circles=np.asarray([1, 3, 9])):
    fig = plt.figure(figsize=(4, 4))
    fig.subplots_adjust(wspace=0, hspace=0)
    if data_2D is not None:
        image_data = data_2D
    else:
        image_data = load_fits_data(image)
    cell_size = get_cell_size(image)

    if apply_mask:
        _, mask_dilated = mask_dilation(image, cell_size=cell_size,
                                        sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations, rms=rms,
                                        PLOT=False)
        image_data = image_data * mask_dilated

    pos_x, pos_y = nd.maximum_position(image_data)

    ax1 = fig.add_subplot(1, 1, 1)
    #     rms= mad_std(load_fits_data(residuallist_comb[-1]))
    if do_cut == True:
        image_data_plot = image_data[int(pos_x - dx):int(pos_x + dx),
                          int(pos_y - dy):int(pos_y + dy)]
        centre_plot = nd.maximum_position(image_data_plot)
    else:
        centre_plot = pos_x, pos_y
        image_data_plot = image_data
    # im1 = ax1.imshow(deconv_image_cut, cmap='magma', aspect='auto')
    ax1 = eimshow(image_data_plot, rms=rms, vmin_factor=vmin_factor, ax=ax1, CM=CM,
                  extent=[-dx, dx, -dx, dx],
                  neg_levels=neg_levels, cmap_cont=cmap_cont,
                  add_contours=add_contours, vmax_factor=vmax_factor)

    #     if centre:
    #     ax1.plot(centre_plot[0]-(pos_x-dx),centre_plot[1]-(pos_x-dx), marker=r'x', color='black',ms=10)
    ax1.plot(centre_plot[0] * cell_size * 0, centre_plot[1] * cell_size * 0, marker=r'x',
             color='white', ms=10)
    #     ax1.set_title('Deconvolved',size=14)
    ax1 = add_circle_grids(ax=ax1, image=image_data_plot, pix_to_pc=pix_to_pc,
                           #                            center=centre_plot[::-1],
                           center=(centre_plot[1] * cell_size, centre_plot[0] * cell_size),
                           circles=circles, r_d=r_d,  # in pc
                           add_labels=False)

    xticks = np.linspace(-dx, dx, 5)
    xticklabels = np.linspace(-dx * cell_size, +dx * cell_size, 5)
    xticklabels = ['{:.2f}'.format(xtick) for xtick in xticklabels]
    ax1.set_yticks(xticks, [])
    ax1.set_xticks(xticks, xticklabels)

    ax1.set_xlabel('offset [arcsec]')

    ax1.grid(False)
    ax1.axis('on')
    return (ax1, image_data_plot)


def add_circle_grids(ax, image, pix_to_pc, r_d=200, add_labels=True,
                     extent=None, center=None,
                     circles=np.asarray([1, 3, 6, 12])):
    # Set the center of the image
    size_max = int(np.sqrt((0.5 * image.shape[0]) ** 2.0 +
                           (0.5 * image.shape[1]) ** 2.0))
    if center is None:
        center = (image.shape[0] // 2, image.shape[1] // 2)

    # Set the radial distance between circles

    r_d_pix = r_d / pix_to_pc
    # Ni = int(size_max/r_d)
    # Create a figure and axis object
    # fig, ax = plt.subplots()
    #     ax = eimshow(load_fits_data(imagelist_vla[-1]),vmin_factor=3,rms=rms)

    # Plot the image
    # ax.imshow(image)

    # Add the circular dashed lines
    prev_text = None
    angle = -90

    for i in range(len(circles)):
        radius = circles[i] * r_d_pix
        circle = Circle(center, circles[i] * r_d_pix, color='limegreen', lw=2,
                        linestyle='-.',
                        fill=False)
        ax.add_patch(circle)
        label = f'{circles[i] * r_d} pc'
        if add_labels:
            x = center[0] + radius * np.cos(np.deg2rad(angle)) + 50
            y = center[1] + radius * np.sin(np.deg2rad(angle)) + 50
            text = Text(x, y, label, ha='center', va='center', color='green')
            ax.add_artist(text)
            prev_text = text
        angle = angle + 90
    # Show the plot
    ax.axis(False)
    return (ax)


def add_contours(ax, image, levels, scale, colors='white', linewidths=1.0):
    """
    Add contours to a plot with distance labels.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The Axes object to which the contours should be added.
    image : numpy.ndarray
        The image data.
    levels : list or array-like
        The contour levels to plot.
    scale : float
        The scale of the image, in kpc/pixel.
    colors : str or list or array-like, optional
        The color(s) to use for the contours.
    linewidths : float or list or array-like, optional
        The width(s) of the contour lines.
    """
    # Plot the contours
    cs = ax.contour(image, levels=levels, colors=colors, linewidths=linewidths)

    # Add distance labels
    for i, level in enumerate(levels):
        if len(cs.collections) <= i:
            continue
        c = cs.collections[i]
        label = f"{level * scale:.1f} kpc"
        try:
            x, y = c.get_paths()[0].vertices.mean(axis=0)
        except IndexError:
            continue
        ax.text(x, y, label, color='red', ha='center', va='center')

    # Add a colorbar
    ax.figure.colorbar(cs, ax=ax)


def plot_fluxes(x_data,y_data,yerr,
                add_fit_legend=True,
                title_text=None,
                log_plot=True,
                labels = [],
                ecolors = [],
                symbols = [],
                figsize=(6,3),
                fig_save_name=None
               ):

    fig, ax1 = plt.subplots(1, 1, sharex=True, figsize=figsize,
                                   # gridspec_kw={'height_ratios': [1, 1]}
                                 )
    for i in range(y_data.shape[0]):
    # fig, ax = plt.subplots()
    # x_resample = np.linspace(np.min(x_data)*0.7, np.max(y_data)*1.3, 500)
        
        if labels != []:
            label = labels[i]
        else:
            label = None
            
        if ecolors != []:
            ecolor = ecolors[i]
        else:
            ecolor='black'

        if symbols != []:
            fmt = symbols[i]
        else:
            fmt = 'd'
            
        ax1.errorbar(x_data, y_data[i], yerr=yerr[i], markeredgewidth=3.0,
                     fmt=fmt,label=label, color=ecolor, ecolor=ecolor,alpha=0.5)


    if add_fit_legend == True:
        ax1.legend(loc=(0.05, 0.05),frameon=True,prop={'size': 11})
    if np.nanmax(y_data) < 5:
        ax1.set_ylim(0.05*np.nanmin(y_data),1.2*np.nanmax(y_data))
    else:
        ax1.set_ylim(0.05*np.nanmin(y_data),1.2*np.nanmax(y_data))
    
    if title_text is not None:
        ax1.set_title(title_text)

    if log_plot == True:
        ax1.semilogx()
        ax1.semilogy()

    plt.subplots_adjust(hspace=0.05)
    # ax1.legend(loc=(0.7, 0.55),prop={'size': 14},frameon=False)
    ax1.legend(
        # loc=(0.7, 0.80),
               prop={'size': 14},ncol=2,
               frameon=True)
    ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
    ax1.set_xlabel(r'$\nu$ [GHz]')
    if fig_save_name is not None:
        plt.savefig(f"{fig_save_name}.pdf", bbox_inches='tight', pad_inches = 0.01)
    # plt.show()
    return(fig,ax1)
