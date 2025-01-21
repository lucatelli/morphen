def write_data_for_SED_fitter(freq, flux, flux_err, theta1, theta2, flag, filename):
    with open(filename, 'w') as file:
        # Writing the header
        file.write("# Space for commenting the references for the fluxes\n")
        file.write("#--------------------------------------------------------------------------\n")
        file.write("#Freq     Flux       Flux_err   theta1    theta2   Flag\n")
        file.write("# GHz     mJy        mJy        arcsec    arcsec   0=data,1=UL,2=fake\n")
        file.write("#--------------------------------------------------------------------------\n")
        
        # Writing the data
        for f, _flux, _flux_err, t1, t2, fl in zip(freq, flux, flux_err, theta1, theta2, flag):
            file.write(f"{f:.5f}   {_flux:.5f}   {_flux_err:.5f}     {t1:.2f}      {t2:.2f}      {fl}\n")

def write_data_for_SED_fitter_v2(freq, flux, flux_err, flux_err_3, theta1, theta2, flag, filename):
    with open(filename, 'w') as file:
        # Writing the header
        file.write("# Space for commenting the references for the fluxes\n")
        file.write("#--------------------------------------------------------------------------\n")
        file.write("#Freq     Flux        Flux_err   Flux_err_3  theta1   theta2   Flag\n")
        file.write("# GHz     mJy         mJy        mJy         arcsec   arcsec   0=data,1=UL,2=fake\n")
        file.write("#--------------------------------------------------------------------------\n")
        
        # Writing the data
        for f, _flux, _flux_err, _flux_err_3, t1, t2, fl in zip(freq, flux, flux_err, flux_err_3, theta1, theta2, flag):
            file.write(f"{f:.5f}   {_flux:.5f}    {_flux_err:.5f}    {_flux_err_3:.5f}     {t1:.2f}     {t2:.2f}     {fl}\n")

def prepare_data(root_paths,prefix_images):
    MFS_images = []
    MFS_residuals = []
    imagelist = []
    residuallist = []
    freqlist = []
    for i in range(len(root_paths)):
        root_path_i = root_paths[i]
        prefix_images_i = prefix_images[i]
        imagelist_i = glob.glob(root_path_i+prefix_images_i)
        imagelist_i, residuallist_i, freqs_i = sort_list_by_freq(imagelist_i)
        MFS_images.append(imagelist_i[0].replace('0000-image','MFS-image'))
        MFS_residuals.append(residuallist_i[0].replace('0000-residual','MFS-residual'))

        imagelist.extend(imagelist_i)
        residuallist.extend(residuallist_i)
        freqlist.extend(freqs_i)
        
    freqlist_MFS = getfreqs(MFS_images)


    imagelist_beam,residuallist_beam = \
        sort_list_by_beam_size(imagelist=MFS_images,
                                     return_df = False)
#     ax=eimshow(imagelist_beam[0],crop=True,add_beam=True)
    return(np.asarray(MFS_images),
           np.asarray(MFS_residuals),
           np.asarray(imagelist),
           np.asarray(residuallist),
           freqlist_MFS,
           np.asarray(freqlist))


def compute_image_stats_wsclean(path,
                        image_list,
                        image_statistics,
                        prefix='',
                        show_figure=True):
    """
    This function will compute statistics of a cleaned image from a wsclean run
    (provided an image prefix). It will also store
    associated model and residual images.

    Parameters
    ----------
    path : str
        Path to the image files.
    image_list : list
        List to store  the image names of each self-cal step.
    image_statistics : dict
        Dictionary to store the statistics of the images at a given self-cal step.
        It can be an existing dictionary.
    prefix : str
        Prefix of the image files.

    """
    file_list = glob.glob(f"{path}*{prefix}*MFS-image.fits")
    file_list.sort(key=os.path.getmtime, reverse=False)
    try:
        image_list[prefix] = file_list[-1]
    except:
        image_list[prefix] = file_list
    image_list[prefix+'_residual'] = image_list[prefix].replace(
        'MFS-image.fits', 'MFS-residual.fits')
    image_list[prefix+'_model'] = image_list[prefix].replace(
        'MFS-image.fits', 'MFS-model.fits')

    sigma = 6.0
    

    level_stats = level_statistics(image_list[prefix],sigma=sigma)
    image_stats = get_image_statistics(imagename=image_list[prefix],
                                             residual_name = image_list[prefix+'_residual'],
                                             dic_data=level_stats,
                                             sigma_mask=sigma)
    img_props = compute_image_properties(image_list[prefix],
                                               image_list[prefix+'_residual'],
                                               results = image_stats,
                                               sigma_mask = sigma,
                                               show_figure=show_figure)[-1]



    image_statistics[prefix] = img_props
    return(image_statistics,image_list)

def create_wsclean_mask(imagename,rms_mask,sigma_mask,mask_grow_iterations,PLOT=False):
    original_mask,dilated_mask = mask_dilation(imagename,
                                                     PLOT=PLOT,
                                                     rms=rms_mask,
                                                     sigma=sigma_mask,
                                                     iterations=mask_grow_iterations)

    mask = dilated_mask
    mask_wslclean = mask * 1.0  # mask in wsclean is inverted??? But how this works?
    mask_name = imagename.replace('.fits', '') + '_mask.fits'
    pf.writeto(mask_name, mask_wslclean, overwrite=True)
    
    #copy the fits header from the original image into the mask.
    copy_header(imagename,
                      mask_name,
                      mask_name)
    
    return(mask_name,dilated_mask)


def run_wsclean(g_name, imsize='2048', imsizey='2048',cell='0.06asec',
                opt_args = "''",
                robust=0.5,base_name=None,
                savemodel=False,shift=None,
                nsigma_automask='4.0', nsigma_autothreshold='2.0',
                datacolumn='DATA',mask=None,
                niter=1000,quiet=True,
                with_multiscale=False, scales='None',
                uvtaper=[],nc = 4,
                image_list={},image_statistics={},
                calculate_subband_fluxes=True):


    g_vis = g_name + '.ms'
    if imsizey is None:
        imsizey = imsize
    if base_name is None:
        base_name  = 'final_image'
    else:
        base_name = base_name

    imaging_script_path = os.path.dirname(os.path.abspath(__file__))
    print(' >> Imaging script path:',imaging_script_path)
    
    # os.system("export OPENBLAS_NUM_THREADS=1 && python /mnt/ext10TB/GitHub/morphen/selfcal/imaging_with_wsclean.py --f " +
    os.system("export OPENBLAS_NUM_THREADS=1 && python "+imaging_script_path+"/imaging_with_wsclean.py --f " +
              g_name + " --sx "
              + str(imsize) + " --sy " + str(imsizey) + " --niter "
              + str(niter) + " --data " + datacolumn + " --cellsize " + cell
              + ' --nsigma_automask ' + nsigma_automask + ' --mask '+str(mask)
              + ' --nsigma_autothreshold ' + nsigma_autothreshold
              + ' --scales ' + scales + ' --nc ' + str(nc)
              +' --opt_args '+ opt_args
              +' --quiet '+ str(quiet) + ' --with_multiscale '+str(with_multiscale)
              + ' --shift ' + str(shift)
              + " --r " + str(robust) + " --t "+str(uvtaper)
              + " --update_model " + str(savemodel) + " --save_basename " + base_name)

    if calculate_subband_fluxes:
        image_statistics,image_list = compute_image_stats_wsclean(path=os.path.dirname(g_name) + '/',
                                                        image_list=image_list,
                                                        image_statistics=image_statistics,
                                                        prefix=base_name)

        print('Using residual map as RMS estimator:',image_list[base_name+'_residual'])
        rms_mask = mad_std(load_fits_data(image_list[base_name+'_residual']))
        mask_name,dilated_mask = create_wsclean_mask(image_list[base_name],
                                             rms_mask=rms_mask,
                                             sigma_mask=6.0,
                                             mask_grow_iterations=2,
                                             PLOT=True)


    #     eimshow(image_list[base_name],
    #                   rms=mad_std(load_fits_data(image_list[base_name+'_residual'])),
    #                   # crop=True,box_size=300,
    #                   crop=True,box_size=int(4*image_statistics[base_name]['C95radii']),
    #                   add_beam=True)


        try:
            centre = nd.maximum_position(load_fits_data(image_list[base_name]))[::-1]
            # centre = (int(image_statistics[prefix]['y0']),int(image_statistics[prefix]['x0']))
            fig = plt.figure(figsize=(16, 8))
            ax0 = fig.add_subplot(1, 2, 2)
            ax0 = eimshow(imagename = image_list[base_name],
                                center=centre,
                                projection = 'offset',plot_colorbar=True,
                                vmax_factor=0.8,
                                rms=mad_std(load_fits_data(image_list[base_name + '_residual'])),
                                # crop=True,box_size=300,
                                figsize=(8, 8), ax=ax0,fig=fig,
                                crop=True, box_size=int(4 * image_statistics[base_name]['C95radii']),
                                # save_name=image_list[prefix].replace('.fits', '_map'),
                                add_beam=True)
            ax0.set_title(f'Radio Map')
            ax1 = fig.add_subplot(1, 2, 1)
            ax1 = eimshow(imagename = image_list[base_name + '_residual'],
                                center=centre,
                                projection='offset',
                                vmin_factor=-3.0, vmax_factor=0.99,
                                add_contours=False,
                                figsize=(8, 8), ax=ax1,fig=fig,
                                crop=True, box_size=int(4 * image_statistics[base_name]['C95radii']),
                                save_name=image_list[base_name].replace('.fits', '_map'),
                                plot_title = f'Residual Map', plot_colorbar=False,
                                add_beam=False)

        except:
            print('--==>> Error on plotting radio map.')
            pass


        sub_band_images = glob.glob(image_list[base_name].replace('-MFS-image.fits','')+'-????-image.fits')
        sub_band_residuals = []
        for i in range(len(sub_band_images)):
            sub_band_residuals.append(sub_band_images[i].replace('-image.fits','-residual.fits'))



        _FLUXES = []
        _FLUXES_err = []
    #     mask_MFS = load_fits_data(mask_name)
        for i in range(len(sub_band_images)):
            print('++>> Computing flux density on sub-band image ',sub_band_images[i])
            print('++>> Associated sub-band  residual image is ',sub_band_residuals[i])

            img_props = compute_image_properties(sub_band_images[i],
                                                    sub_band_residuals[i],
                                                    sigma_mask = 6.0,
                                                    last_level=1.0,
                                                    mask=dilated_mask,
                                                    verbose=1,
                                                    show_figure=True)[-1]

            # flux_density,flux_density_err = img_props['total_flux_levels'], img_props['flux_error_res_3']
            flux_density,flux_density_err = img_props['total_flux_mask'], img_props['flux_error_res_3']

            print('Flux density = ', flux_density)
            _FLUXES.append(flux_density)
            _FLUXES_err.append(flux_density_err)
        FLUXES = np.asarray(_FLUXES)
        FLUXES_err = np.asarray(_FLUXES_err)
        freqlist = getfreqs(sub_band_images)

        try:
            do_mcmc_fit = True
            mini,result_1,param_dict = \
                do_fit_spec_RC_linear(freqlist,FLUXES*1000,FLUXES_err*1000,
                                    plot_errors_shade = True,
                                    do_mcmc_fit=do_mcmc_fit,
                                    basename_save=image_list[base_name],
                                    verbose=2)
        except:
            do_mcmc_fit = False
            mini,result_1,param_dict = \
                do_fit_spec_RC_linear(freqlist,FLUXES*1000,FLUXES_err*1000,
                                    plot_errors_shade = True,
                                    do_mcmc_fit=do_mcmc_fit,
                                    basename_save=image_list[base_name],
                                    verbose=2)
        print('---------------------------------------------')
        print('------------ RESTORING BEAM SIZE ------------')
        Omaj,Omin, _, _, _ = beam_shape(image_list[base_name])
        print(f"                {Omaj:.3f} arcec            ")
        print('---------------------------------------------')
    else:
        image_list = {}
        image_statistics = {}
        Omaj = None
    
    return(image_list,image_statistics,Omaj)


# import builtins
# import inspect

# # # Override the built-in print function
# # original_print = builtins.print

# # def custom_print(*args, **kwargs):
# #     frame = inspect.currentframe().f_back
# #     info = inspect.getframeinfo(frame)
# #     file_name = info.filename
# #     line_number = info.lineno
# #     original_print(f'[{file_name}:{line_number}]', *args, **kwargs)

# # builtins.print = custom_print


# # Save the original print function
# original_print = builtins.print

# def custom_print(*args, **kwargs):
#     frame = inspect.currentframe().f_back
#     info = inspect.getframeinfo(frame)
#     file_name = info.filename
#     line_number = info.lineno
#     # Call the original print function with the additional file and line info
#     original_print(f'[{file_name}:{line_number}]', *args, **kwargs)

# # Override the built-in print function
# builtins.print = custom_print


def beam_area(Omaj, Omin, cellsize):
    '''
    Computes the estimated projected beam area (theroetical),
    given the semi-major and minor axis
    and the cell size used during cleaning.
    Return the beam area in pixels.
    '''
    BArea = ((np.pi * Omaj * Omin) / (4 * np.log(2))) / (cellsize ** 2.0)
    return (BArea)


def beam_area2(image, cellsize=None):
    '''
    Computes the estimated projected beam area (theroetical),
    given the semi-major and minor axis
    and the cell size used during cleaning.
    Return the beam area in pixels.
    '''
    if cellsize is None:
        try:
            cellsize = get_cell_size(image)
        except:
            print('Unable to read cell size from image header. '
                  'Please, provide the cell size of the image!')
            pass
    imhd = imhead(image)
    Omaj = imhd['restoringbeam']['major']['value']
    Omin = imhd['restoringbeam']['minor']['value']
    BArea = ((np.pi * Omaj * Omin) / (4 * np.log(2))) / (cellsize ** 2.0)
    return (BArea)

def getfreqs(fitslist):
    freqs = []
    for fitsfile in fitslist:
        hdu = fits.open(fitsfile)
        hdr = hdu[0].header
        freq = hdr['CRVAL3']
        freqs.append(freq)
    _freqs = np.array(freqs)
    return _freqs

def beam_shape(image):
    '''
    Return the beam shape (bmin,bmaj,pa) from given image.
    It uses CASA's function `imhead`.

    '''
    import numpy as np
    from astropy import units as u
    cell_size = get_cell_size(image)
    imhd = imhead(image)
    Omaj = imhd['restoringbeam']['major']['value']
    Omin = imhd['restoringbeam']['minor']['value']
    PA = imhd['restoringbeam']['positionangle']['value']
    freq = imhd['refval'][2] / 1e9
    """
    bmaj,bmin,PA,freq = beam_shape(crop_image)
    """
    bmaj = Omaj*u.arcsec
    bmin = Omin*u.arcsec
    freq_ = freq * u.GHz

    fwhm_to_sigma = 1./(8*np.log(2))**0.5
    BAarcsec = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)
    # # BA
    # equiv = u.brightness_temperature(freq_)
    # (0.0520*u.Jy/BA).to(u.K, equivalencies=equiv)
    return (Omaj,Omin,PA,freq,BAarcsec)

def sort_list_by_beam_size(imagelist, residuallist=None,return_df=False):
    """
    Sort a list of images files by beam size.

    If no residual list is provided, it will assume that the residual files are in the same
    directory as the image files, and that the residual files have the same name prefix as the
    image files.


    Parameters
    ----------
    imagelist : list
        List of image files.
    residuallist : list, optional
        List of associated residual files.
    return_df : bool, optional
        Return a pandas dataframe with the beam sizes.

    Returns
    -------
    imagelist_sort : list
        Sorted list of image files.
    residuallist_sort : list
        Sorted list of residual files.
    df_beam_sizes_sorted : pandas dataframe
        Dataframe with the beam sizes.
    """
    beam_sizes_list = []
    for i in tqdm(range(len(imagelist))):
        beam_sizes = {}
        aO, bO, _, _, _ = beam_shape(imagelist[i])
        beam_size_arcsec = np.sqrt(aO*bO)
        # beam_size_px, _, _ = get_beam_size_px(imagelist[i])
        beam_sizes['imagename'] = imagelist[i]
        if residuallist is not None:
            beam_sizes['residualname'] = residuallist[i]
        else:
            beam_sizes['residualname'] = \
                imagelist[i].replace('/MFS_images/','/MFS_residuals/')\
                            .replace('-image','-residual')
        beam_sizes['id'] = i
        beam_sizes['B_size_arcsec'] = beam_size_arcsec
        beam_sizes_list.append(beam_sizes)

    df_beam_sizes = pd.DataFrame(beam_sizes_list)
    df_beam_sizes_sorted = df_beam_sizes.sort_values('B_size_arcsec')
    imagelist_sort = np.asarray(df_beam_sizes_sorted['imagename'])
    residuallist_sort = np.asarray(df_beam_sizes_sorted['residualname'])
    # i = 0
    # for image in imagelist_sort:
    #     print(i,'>>',os.path.basename(image))
    #     i=i+1
    if return_df==True:
        return(imagelist_sort,residuallist_sort,df_beam_sizes_sorted)
    if return_df==False:
        return (imagelist_sort, residuallist_sort)


def sort_list_by_freq(imagelist, return_df=False):
    """
    Sort a list of images files by frequency.

    If no residual list is provided, it will assume that the residual files are in the same
    directory as the image files, and that the residual files have the same name prefix as the
    image files.

    Parameters
    ----------
    imagelist : list
        List of image files.
    residuallist : list, optional
        List of associated residual files.
    return_df : bool, optional
        Return a pandas dataframe with the beam sizes.

    Returns
    -------
    imagelist_sort : list
        Sorted list of image files.
    residuallist_sort : list
        Sorted list of residual files.
    df_beam_sizes_sorted : pandas dataframe
        Dataframe with the beam sizes.
    """
    freqs = getfreqs(imagelist)

    # Get indices that would sort the frequency array
    sorted_indices = np.argsort(freqs)
    # Sort the frequency array and image list using the sorted indices
    sorted_freq_array = freqs[sorted_indices]
    sorted_image_list = [imagelist[i] for i in sorted_indices]

    residuallist = []
    for i in range(len(sorted_image_list)):
        residuallist_i = \
            sorted_image_list[i].replace('/MFS_images/', '/MFS_residuals/') \
                .replace('-image.', '-residual.').replace('-image-pb', '-residual-pb')
        residuallist.append(residuallist_i)

    imagelist_sort = np.asarray(sorted_image_list)
    residuallist_sort = np.asarray(residuallist)
    i = 0
    for image in imagelist_sort:
        print(i, '>>', os.path.basename(image))
        i = i + 1
    return (imagelist_sort, residuallist_sort, sorted_freq_array)


def get_beam_size_px(imagename):
    aO,bO,_,_,_ = beam_shape(imagename)
    cs = get_cell_size(imagename)
    aO_px = aO/cs
    bO_px = bO/cs
    beam_size_px = np.sqrt(aO_px * bO_px)
    return(beam_size_px,aO_px,bO_px)

def beam_physical_area(imagename,z):
    '''
    Return the beam shape (bmin,bmaj,pa) given an image.
    '''
    import numpy as np
    from astropy import units as u
    cell_size = get_cell_size(imagename)
    imhd = imhead(imagename)
    Omaj = imhd['restoringbeam']['major']['value']
    Omin = imhd['restoringbeam']['minor']['value']
    """
    bmaj,bmin,PA,freq = beam_shape(crop_image)
    """
    pc_scale = arcsec_to_pc(z=z,cell_size=cell_size)
    bmaj = Omaj*pc_scale
    bmin = Omin*pc_scale

    fwhm_to_sigma = 1./(8*np.log(2))**0.5
    BAarcsec = 2.*np.pi*(bmaj*bmin*fwhm_to_sigma**2)

    return (pc_scale,bmaj,bmin,BAarcsec)


def get_phase_centre(vis):
    from astropy.coordinates import SkyCoord
    import astropy.units as u

    msmd = casatools.msmetadata()
    ms = casatools.ms()
    tb = casatools.table()

    msmd.open(vis)
    ra_radians = msmd.phasecenter()['m0']['value']
    dec_radians = msmd.phasecenter()['m1']['value']
    msmd.close()
    # Convert to SkyCoord object
    coord = SkyCoord(ra=ra_radians * u.radian, dec=dec_radians * u.radian, frame='icrs')

    # Format the output using 'hmsdms'
    formatted_coord = coord.to_string('hmsdms')
    formatted_ra, formatted_dec = formatted_coord.split()

    formatted_ra_hms = formatted_ra.replace('h', ':').replace('m', ':').replace('s', '')
    formatted_dec_dms = formatted_dec.replace('d', '.').replace('m', '.').replace('s', '')

    formatted_output = "J2000 {} {}".format(formatted_ra_hms, formatted_dec_dms)

    print(formatted_output)
    return (formatted_output)

def tcreate_beam_psf(imname, cellsize=None,size=(128,128),app_name='',
                     aspect=None):
    """
    From an interferometric image, reconstruct the restoring beam as a PSF
    gaussian image, with the same size as the original image.

    Parameters
    ----------
    imname : str
        Image name.
    cellsize : float
        Cell size of the image in arcsec.
    size : tuple
        Size of the PSF image.
    app_name : str
        Name to append to the PSF image.
    aspect : str
        Aspect ratio of the PSF image. If 'equal', the PSF will be circular.
        If None, the PSF will be elliptical.
    Returns
    -------
    psf_name : str
        PSF image name.

    """
    if cellsize is None:
        try:
            cellsize = get_cell_size(imname)
        except:
            print('Please, provide a cellsize for the image.')
            return (ValueError)

    msmd = casatools.msmetadata()
    ms = casatools.ms()
    tb = casatools.table()
    cl = casatools.componentlist()
    ia = IA()
    qa = casatools.quanta()
    # ia.open(image)

    imst = imstat(imname)
    imhd = imhead(imname)
    tb.close()
    direction = "J2000 10h00m00.0s -30d00m00.0s"
    cl.done()
    freq = str(imhd['refval'][2] / 1e9) + 'GHz'
    if (aspect=='circular') or (aspect=='equal'):
        print('WARNING: Using circular Gaussian for Gaussian beam convolution.')
        minoraxis = str(imhd['restoringbeam']['minor']['value']) + str(
            imhd['restoringbeam']['major']['unit'])
        majoraxis = str(imhd['restoringbeam']['minor']['value']) + str(
            imhd['restoringbeam']['major']['unit'])

    if (aspect == 'elliptical') or (aspect is None):
    # else:
        print('INFO: Using Elliptical Gaussian for Gaussian beam convolution.')
        minoraxis = str(imhd['restoringbeam']['minor']['value']) + str(
            imhd['restoringbeam']['minor']['unit'])
        majoraxis = str(imhd['restoringbeam']['major']['value']) + str(
            imhd['restoringbeam']['major']['unit'])
    # print(f"++==>>  PSF major/minor axis = ', {majoraxis} X {minoraxis}")

    pa = str(imhd['restoringbeam']['positionangle']['value']) + str(
        imhd['restoringbeam']['positionangle']['unit'])
    cl.addcomponent(dir=direction, flux=1.0, fluxunit='Jy', freq=freq,
                    shape="Gaussian", majoraxis=majoraxis, minoraxis=minoraxis,
                    positionangle=pa)
    ia.fromshape(imname.replace('.fits', '_beampsf'+app_name+'.im'),
                 [size[0],size[1],1,1], overwrite=True)
    cs = ia.coordsys()
    cs.setunits(['rad', 'rad', '', 'Hz'])
    cell_rad = qa.convert(qa.quantity(str(cellsize) + "arcsec"), "rad")['value']
    cs.setincrement([-cell_rad, cell_rad], 'direction')
    cs.setreferencevalue([qa.convert("10h", 'rad')['value'],
                          qa.convert("-30deg", 'rad')['value']],
                         type="direction")
    cs.setreferencevalue(freq, 'spectral')
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.modify(cl.torecord(), subtract=False)
    exportfits(
        imagename=imname.replace('.fits', '_beampsf'+app_name+'.im'),
        fitsimage=imname.replace('.fits', '_beampsf'+app_name+'.fits'),
        overwrite=True)
    psf_name = imname.replace('.fits', '_beampsf'+app_name+'.fits')
    cl.close()
    return(psf_name)


def create_box_around_peak(imagename, fractions=None):
    """
    Create a box with 25% (or specified fraction) of the image
    around the peak.
    """
    ihl = imhead(imagename, mode='list')
    ih = imhead(imagename)
    st = imstat(imagename)

    M = ihl['shape'][0]
    N = ihl['shape'][1]
    if fractions is None:
        frac_X = int(0.25 * M)
        frac_Y = int(0.25 * N)
    else:
        frac_X = int(fractions[0] * M)
        frac_Y = int(fractions[1] * N)
    # slice_pos_X = 0.15 * M
    # slice_pos_Y = 0.85 * N
    slice_pos_X = st['maxpos'][0]
    slice_pos_Y = st['maxpos'][1]

    box_edge = np.asarray([slice_pos_X - frac_X,
                           slice_pos_Y - frac_Y,
                           slice_pos_X + frac_X,
                           slice_pos_Y + frac_Y]).astype(int)

    box_edge_str = str(box_edge[0]) + ',' + str(box_edge[1]) + ',' + \
                   str(box_edge[2]) + ',' + str(box_edge[3])

    return (box_edge_str)


def create_box(imagename, fracX=0.15, fracY=0.15):
    """
    Create a box with 20% of the image
    at an edge (upper left) of the image.
    To be used with casa tasks.
    """
    ihl = imhead(imagename, mode='list')
    ih = imhead(imagename)

    M = ihl['shape'][0]
    N = ihl['shape'][1]
    frac_X = int(fracX * M)
    frac_Y = int(fracY * N)
    slice_pos_X = int(0.02 * M)
    slice_pos_Y = int(0.98 * N)

    box_edge = np.asarray([slice_pos_X,
                           slice_pos_Y - frac_Y,
                           slice_pos_X + frac_X,
                           slice_pos_Y]).astype(int)

    box_edge_str = str(box_edge[0]) + ',' + str(box_edge[1]) + ',' + \
                   str(box_edge[2]) + ',' + str(box_edge[3])

    return (box_edge_str, ih)


def create_box_np(imagename, fracX=0.15, fracY=0.15):
    """
    Create a box with 20% of the image
    at an edge (upper left) of the image.
    To be used on a np.array.
    """
    ihl = imhead(imagename, mode='list')
    ih = imhead(imagename)

    M = ihl['shape'][0]
    N = ihl['shape'][1]

    frac_X = int(fracX * M)
    frac_Y = int(fracY * N)
    slice_pos_X = int(0.02 * M)
    slice_pos_Y = int(0.98 * N)

    cut = [[slice_pos_X, slice_pos_X + frac_X],
           [slice_pos_Y - frac_Y, slice_pos_Y]
           ]

    return (cut)



def convolve_2D_smooth(imagename, imagename2=None,
                       mode='same', add_prefix='_convolve2D'):
    """
    CASA's way of convolving images.
    Very slow, but accurate.

    This function will be removed in the future, as all convolution operations
    will be migrated into pure python tasks.

    Parameters
    ----------
        imagename (str): The name of the image to convolve.
        imagename2 (str): The name of the image to use as the restoring beam.
        mode (str): The mode of convolution. Can be 'same' or 'transfer'.
        add_prefix (str): The prefix to add to the output image name.

    """
    if mode == 'same':
        imhd = imhead(imagename)
        imsmooth(imagename=imagename,
                 outfile=imagename.replace('.fits', '_convolved2D.fits'), overwrite=True,
                 major=imhd['restoringbeam']['major'],
                 minor=imhd['restoringbeam']['minor'],
                 pa=imhd['restoringbeam']['positionangle'])
        return (imagename.replace('.fits', '_convolved2D.fits'))

    if mode == 'transfer' and imagename2 != None:
        '''
        Use restoring beam from image1 to convolve image2.
        '''
        imhd = imhead(imagename2)
        outfile = imagename2.replace('.fits', add_prefix + '.fits')
        imsmooth(imagename=imagename,
                 outfile=outfile, overwrite=True,
                 major=imhd['restoringbeam']['major'],
                 minor=imhd['restoringbeam']['minor'],
                 pa=imhd['restoringbeam']['positionangle'])
        return (outfile)




"""
 ____  _               _          
|  _ \| |__  _   _ ___(_) ___ ___ 
| |_) | '_ \| | | / __| |/ __/ __|
|  __/| | | | |_| \__ \ | (__\__ \
|_|   |_| |_|\__, |___/_|\___|___/
             |___/  
             
#Physics
"""
def T_B(theta_maj, theta_min, freq, I):
    # https://science.nrao.edu/facilities/vla/proposing/TBconv
    # https://www.cv.nrao.edu/~sransom/web/Ch2.html
    # theta_maj/min are gaussian beam maj and min axes in units of arcsecs
    # freq is the central frequency in units of GHz
    # I is the total flux measured by the beam - peak flux
    # divifing I/theta_maj & min converts I(mJy/beam) to S/solid angle
    brightness_temp = 1000*1222 * I / ((freq ** 2) * (theta_maj * theta_min))

    return brightness_temp/1e5


def Tb_source(Snu,freq,theta1,theta2,z,
              Snu_err=None,theta1_err=None,theta2_err=None):
    """
    Compute the brightness temperature, provided the deconvolved model parameters
    having semi-major and semi-minor FWHM's theta1 and theta2.
    
    Parameters
    ----------
    Snu : float
        The flux density at the frequency nu.
    freq : float
        The frequency at which the flux density is measured.
    theta1 : float
        The deconvolved FWHM axis of the source.
    theta2 : float
        The deconvolved FWHM axis of the source.
    z : float
        The redshift of the source.
    
    Returns
    -------
    Tb : float
        The brightness temperature of the source divided by 1e5.
    """
    const = 1.8e9 * (1+z)*1000
    Tb = ((const * Snu)/(freq*freq*(theta1*1000)*(theta2*1000)))/1e5
    if (Snu_err is not None) and (theta1_err is not None) and (theta2_err is not None):
        Tb_err = Tb * np.sqrt((Snu_err/Snu)**2.0 
                              + (theta1_err/theta1)**2.0 
                              + (theta2_err/theta2)**2.0)
        return(Tb,Tb_err)
    else:
        return(Tb)

def Olim(Omaj,SNR):
    Olim_ = Omaj * np.sqrt((4*np.log(2)/np.pi) * np.log(SNR/(SNR-1)))
    return(Olim_)

def deconv_R50(R50conv,theta12):
    return(np.sqrt(4*R50conv**2.0 - theta12**2.0 )/2)

def phi_source(OH,SpH,OL,SpL):
    """
    ## Source Sizes
    If the circular Gaussian source is imaged with two different resolutions
    $\theta_H$ and $\theta_L$, the ratio of the image peak brightnesses is
    \begin{equation}
        \frac{S_p^{H}}{S_p^{L}} =
        \left(
        1 + \frac{\phi^2}{\theta_L^2}
        \right)
        \left(
        1 + \frac{\phi^2}{\theta_H^2}
        \right)^{-1}
    \end{equation}
    This equation can be solved for the source size
    \begin{align}
        \phi =
        \left[
            \frac{\theta_L^2 \theta_H^2 (S_p^L - S_p^H)}{\theta_L^2 S_p^H - \theta_H^2 S_p^L}
        \right]^{1/2}
    \end{align}

    """
    nume = ((OL**2) * (OH**2)) * (SpL - SpH)
    deno = OL*OL*SpH - OH*OH*SpL
    phi_size = np.sqrt(nume/deno)
    return(phi_size)

def get_size_params(df,pix_to_pc):
    Bb = (df['bmin_pc']/pix_to_pc)*df['cell_size']
    Ba = (df['bmaj_pc']/pix_to_pc)*df['cell_size']
    max_im = df['max']
    Bsize = np.sqrt(Bb*Ba)
    return(Bsize,max_im)

def LT(In,Rn,n):
    """
    Total luminosity of a Sersic Function

    """
    bn = 2*n - 1.0/3.0 + 4/(405*n)
    num = 2*np.pi*In*Rn*Rn*n*np.exp(bn)*scipy.special.gamma(2*n)
    den = bn**(2*n)
    return(num/den)


def D_Cfit(z):
    h = 0.687
    H0 = 100 * h  # km/(s * Mpc)
    c = 299792.458  # km/s
    # comoving distance Dc
    D_H0 = c / H0  # Mpc
    a = 1.0 / (1 + z)
    n1 = a / (1 - a)
    n2 = (1 - a) / (0.785 + a)
    n3 = (1 - a) / ((0.312 + a) ** 2.0)
    D_C = D_H0 / (n1 + 0.2278 + 0.2070 * n2 - 0.0158 * n3)
    DC_MPC = D_C
    return (DC_MPC)


def compute_Lnu(flux, z, alpha):
    dist_conversion_factor = 3.08567758128 * 1e24  # m #3.08567758128*(10**24) #cm
    # lum_conversion_factor =  10**(-26) # W/(m^2 Hz Jy)
    lum_conversion_factor = 1e-23
    # D_L = luminosity_distance_cosmo(z=z)
    D_C = comoving_distance_cosmo(z=z) * dist_conversion_factor
    # D_Lnu = D_C * ((1 + z) ** ((alpha + 1) / 2)) * dist_conversion_factor
    D_L = D_C * (1+z)
    # flux is in Jy
    L_nu = (4.0 * np.pi * (D_L ** 2.0)/((1+z)**(alpha+1))) * flux * lum_conversion_factor
    # print(f"Comoving distance DC= {D_C/dist_conversion_factor} Mpc")
    # print(f"Luminosity distance DL= {D_L/dist_conversion_factor} Mpc")
    # print(f"Spetral Luminosity Lnu = {L_nu}")
    # Ld = 4.0 * np.pi * D_L**2.0 * flux
    # L_nu_error = 0.0
    return (L_nu)

def compute_Lnu_v2(flux, z, alpha,freq,nu0=1.0):
    dist_conversion_factor = 3.08567758128 * 1e24  # m #3.08567758128*(10**24) #cm
    # lum_conversion_factor =  10**(-26) # W/(m^2 Hz Jy)
    lum_conversion_factor = 1e-23
    # D_L = luminosity_distance_cosmo(z=z)
    D_C = comoving_distance_cosmo(z=z) * dist_conversion_factor
    # D_Lnu = D_C * ((1 + z) ** ((alpha + 1) / 2)) * dist_conversion_factor
    D_L = D_C * (1+z)
    # flux is in Jy
    L_nu = (4.0 * np.pi * (D_L ** 2.0)/((1+z)**(alpha+1))) * ((freq/nu0)**alpha) *flux * lum_conversion_factor
    # print(f"Comoving distance DC= {D_C/dist_conversion_factor} Mpc")
    # print(f"Luminosity distance DL= {D_L/dist_conversion_factor} Mpc")
    # print(f"Spetral Luminosity Lnu = {L_nu}")
    # Ld = 4.0 * np.pi * D_L**2.0 * flux
    # L_nu_error = 0.0
    return (L_nu)


def tabatabaei(nu, al, z):
    thermal_frac = 1. / (1. + 13. * ((nu * (1 + z)) ** (0.1 + al)))
    return (thermal_frac)


# Calculate a luminosity, star formation rate, and uncertainties given a flux density
# Using equation relating synchrotron emission to star formation rate given in Murphy et. al (2011)
# Also using Condon & Matthews (2018) to calculate spectral luminosity distance
def calc_params(flux, flux_error, redshift, redshift_error, freqs, alphas=-0.85,
                allowthermal=False):
    # Defining symbols (sympy)
    z = Symbol('z')  # redshift
    zunc = Symbol('zunc')  # redshift uncertainty
    nu = Symbol('nu')  # frequency
    nuunc = Symbol('nuunc')  # frequency uncertainty
    al = Symbol('al')  # non-thermal spectral index alpha
    alunc = Symbol('alunc')  # alpha uncertainty
    f = Symbol('f')  # flux
    f_unc = Symbol('f_unc')  # flux uncertainty
    Ho = Symbol('Ho')  # hubble constant
    Ho_unc = Symbol('Ho_unc')  # hubble uncertainty

    # define speed of light and Hubble distance
    c = 299792.458  # km/s
    Dho = c / Ho

    # Define symbolic formluas for desired quantities
    # convenience definition from Murphy et al. paper
    a = 1 / (1 + z)
    # Comoving distance formula
    # 3E24 factor converts between cm and Mpc
    Dc = (Dho / (a / (1 - a) + 0.2278 + 0.2070 * (1 - a) / (
                0.785 + a) - 0.0158 * (1 - a) / (
                             (0.312 + a) ** 2))) * 3.08567758128 * (
                     10 ** 24)  # cm
    # luminosity distance formula
    Dl = (1 + z) * Dc  # cm
    # spectral luminosity distance
    Dl_nu = Dl * ((1 + z) ** (-(al + 1) / 2))
    # inverse square law to get luminosity
    Lumform = (4 * np.pi * f * (10 ** -23) * Dl_nu ** 2)  # ergs/s
    # SFR formula in solar masses/yr (Murphy et. al)
    kroupa_to_salpeter = 1.5
    # SFRform = kroupa_to_salpeter*(6.64e-29*(nu**(-al))*Lumform)
    if allowthermal:
        L_NT = Lumform / (1 + 1 / 10. * ((1 + z) * nu) ** (-0.1 - al))
        SFRform = 6.64e-29 * (nu ** (-al)) * L_NT
    else:
        SFRform = 6.64e-29 * (nu ** (-al)) * Lumform
    # luminosity uncertainty formula - simple error propagation
    Lum_stat_unc = ((diff(Lumform, f) * f_unc) ** 2) ** 0.5
    Lum_syst_unc = ((diff(Lumform, z) * zunc) ** 2 + (
                diff(Lumform, Ho) * Ho_unc) ** 2) ** 0.5
    # SFR uncertainty formula
    SFR_stat_uncertainty = ((diff(SFRform, f) * f_unc) ** 2) ** .5
    SFR_syst_uncertainty = ((diff(SFRform, z) * zunc) ** 2 + (
                diff(SFRform, Ho) * Ho_unc) ** 2 + (
                                        diff(SFRform, al) * alunc) ** 2) ** .5

    # Define constants
    Hubble = 67.8  # km/s/Mpc
    Hubble_unc = 2
    # freqs = 1.51976491105  # GHz
    freqsigs = 0
    alphasig = 0.05

    output = []

    # substitute in values into symbolic expressions
    # SFR values
    SF = SFRform.subs({nu: freqs, al: alphas, f: flux, z: redshift, Ho: Hubble})
    SF_stat = SFR_stat_uncertainty.subs(
        {nu: freqs, al: alphas, f: flux, z: redshift,
         f_unc: flux_error, Ho: Hubble})
    SF_syst = SFR_syst_uncertainty.subs(
        {nu: freqs, al: alphas, f: flux, z: redshift, zunc: redshift_error,
         alunc: alphasig, f_unc: flux_error, Ho: Hubble, Ho_unc: Hubble_unc})
    # luminosity values
    Lum = Lumform.subs({f: flux, z: redshift, al: alphas, Ho: Hubble})
    Lum_stat = Lum_stat_unc.subs(
        {f: flux, z: redshift, f_unc: flux_error, Ho: Hubble, al: alphas})
    Lum_syst = Lum_syst_unc.subs(
        {f: flux, z: redshift, f_unc: flux_error, zunc: redshift_error,
         Ho: Hubble, Ho_unc: Hubble_unc})

    output.append(Lum)
    output.append(Lum_stat)
    output.append(SF)
    output.append(SF_stat)

    return output


# def compute_SFR_general(flux, frequency, z, alpha, alpha_NT=-0.85, flux_error=None,
#                    calibration_kind='Murphy12', return_with_error=False):
#     '''
#         To do:
#             [ ] - Implement error estimates
#             [ ] - Check the thermal contribution
#             [ ] - Use spectral index from the image
#             [ ] - Explore different Te's (electron temperature)
#     '''

#     if calibration_kind == 'Murphy11':
#         Lnu_NT, Lnu_NT_error = compute_Lnu(flux, z,
#                                            alpha)  # 0.0014270422727500343
#         SFR = 6.64 * (1e-29) * ((frequency) ** (-alpha_NT)) * Lnu_NT
#         if flux_error is not None:
#             Lnu_NT_error, Lnu_NT_error2 = compute_Lnu(flux_error, z, alpha)
#             SFR_error = 6.64 * (1e-29) * ((frequency) ** (-alpha_NT)) * Lnu_NT_error
#         else:
#             SFR_error = 0.0

#     if calibration_kind == 'Tabatabaei2017':
#         '''
#         There is something wrong for this kind!
#         '''
#         Lnu_NT, Lnu_NT_error = compute_Lnu(flux, z,
#                                            alpha)  # 0.0014270422727500343
#         SFR = 1.11 * 1e-37 * 1e9 * frequency * Lnu_NT
#         if flux_error is not None:
#             Lnu_NT_error, Lnu_NT_error2 = compute_Lnu(flux_error, z, alpha)
#             SFR_error = 1.11 * 1e-37 * 1e9 * frequency * Lnu_NT_error
#         else:
#             SFR_error = 0.0

#     if calibration_kind == 'Murphy12':
#         Te = 1e4
#         Lnu_NT, Lnu_NT_error = compute_Lnu(flux, z, alpha)
#         SFR = 1e-27 * ( \
#                     (2.18 * ((Te / (1e4)) ** 0.45) * (frequency ** (-0.1)) + \
#                      15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
#             ) * Lnu_NT

#         if flux_error is not None:
#             Lnu_NT_error, Lnu_NT_error2 = compute_Lnu(flux_error, z, alpha)
#             SFR_error = 1e-27 * ( \
#                         (2.18 * ((Te / (1e4)) ** 0.45) * (frequency ** (-0.1)) + \
#                          15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
#                 ) * Lnu_NT_error
#         else:
#             SFR_error = 0.0
#     # print('SFR =', SFR, '+/-', SFR_error, 'Mo/yr')
#     if return_with_error:
#         return (SFR, SFR_error)
#     else:
#         return (SFR)

def compute_SFR_general(flux, frequency, z, alpha, alpha_NT=-0.85, 
                        nu0 = 1.0,
                        flux_error=None,
                        calibration_kind='Murphy12'):
    '''
        To do:
            [ ] - Implement error estimates
            [ ] - Check the thermal contribution
            [ ] - Use spectral index from the image
            [ ] - Explore different Te's (electron temperature)
    '''

    Lnu,Lnu_error,SFR,SFR_error = None, None, None, None
    if calibration_kind == 'Murphy11':
        Lnu_NT = compute_Lnu(flux, z, alpha)  # 0.0014270422727500343
        SFR = 6.64 * (1e-29) * ((frequency) ** (-alpha_NT)) * Lnu_NT
        if flux_error is not None:
            Lnu_NT_error = compute_Lnu(flux_error, z, alpha)
            SFR_error = 6.64 * (1e-29) * ((frequency) ** (-alpha_NT)) * Lnu_NT_error
        else:
            SFR_error = 0.0

    # if calibration_kind == 'Tabatabaei2017':
    #     '''
    #     There is something wrong for this kind!
    #     '''
    #     Lnu_NT = compute_Lnu(flux, z, alpha)  # 0.0014270422727500343
    #     SFR = 1.11 * 1e-37 * 1e9 * frequency * Lnu_NT
    #     if flux_error is not None:
    #         Lnu_NT_error = compute_Lnu(flux_error, z, alpha)
    #         SFR_error = 1.11 * 1e-37 * 1e9 * frequency * Lnu_NT_error
    #     else:
    #         SFR_error = 0.0

    # if calibration_kind == 'Murphy12':
    #     Te = 1e4
    #     Lnu_NT = compute_Lnu(flux, z, alpha_NT)
    #     SFR = 1e-27 * ( \
    #                 (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
    #                  15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
    #         ) * Lnu_NT

    #     if flux_error is not None:
    #         Lnu_NT_error = compute_Lnu(flux_error, z, alpha_NT)
    #         SFR_error = 1e-27 * ( \
    #                     (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
    #                      15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
    #             ) * Lnu_NT_error
    #     else:
    #         SFR_error = 0.0
    
    if calibration_kind == 'Murphy12':
        Te = 1e4
        Lnu = compute_Lnu(flux, z, alpha_NT)
        SFR = 1*1e-27 * ( \
                    (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                     15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
            ) * Lnu

        if flux_error is not None:
            Lnu_error = compute_Lnu(flux_error, z, alpha_NT)
            SFR_error = 1*1e-27 * ( \
                        (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                         15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
                ) * Lnu_error
        else:
            SFR_error = 0.0

    if calibration_kind == 'Tabatabaei2017':
        Te = 1e4
        Lnu = compute_Lnu(flux, z, alpha_NT)
        SFR = 1*1e-27 * ( \
                    (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                     15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
            ) * Lnu

        if flux_error is not None:
            Lnu_error = compute_Lnu(flux_error, z, alpha_NT)
            SFR_error = 1*1e-27 * ( \
                        (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                         15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
                ) * Lnu_error
        else:
            SFR_error = 0.0
    
    if calibration_kind == 'Murphy12_NT':
        Lnu = compute_Lnu(flux, z, alpha_NT)
        SFR = 6.64*(1e-29) * ((frequency/nu0)**(-alpha_NT)) * Lnu
        
        if flux_error is not None:
            Lnu_error = compute_Lnu(flux_error, z, alpha_NT)
            SFR_error = 6.64*1e-29 * ((frequency/nu0)**(-alpha_NT)) * Lnu_error
        else:
            SFR_error = 0.0
            
    if calibration_kind == 'Murphy12_T':
        Te = 1e4
        Lnu = compute_Lnu(flux, z, alpha)
        SFR = 4.59*1e-28 * ((frequency/nu0)**(-alpha)) * ((Te / (1e4)) ** (-0.45)) * Lnu
        
        if flux_error is not None:
            Lnu_error = compute_Lnu(flux_error, z, alpha)
            SFR_error = 4.59*1e-28 * ((frequency/nu0)**(-alpha)) * ((Te / (1e4)) ** (-0.45)) * Lnu_error
        else:
            SFR_error = 0.0
    # print(f'nu =', frequency)
    # print(f'alpha =', alpha)
    # print(f'SFR =', SFR, '+/-', SFR_error, 'Mo/yr')
    print(f'Lnu =', Lnu, '+/-', Lnu_error, 'erg/s/Hz')
    # print('SFR =', SFR, '+/-', SFR_error, 'Mo/yr')
    return (SFR, SFR_error)



def compute_SFR_general_v2(flux, frequency, z, alpha, alpha_NT=-0.85, 
                        nu0 = 1.0,
                        flux_error=None,
                        calibration_kind='Murphy12'):
    '''
        To do:
            [ ] - Implement error estimates
            [ ] - Check the thermal contribution
            [ ] - Use spectral index from the image
            [ ] - Explore different Te's (electron temperature)
    '''
    Lnu,Lnu_error,SFR,SFR_error = None, None, None, None
    
    if calibration_kind == 'Murphy11':
        Lnu_NT = compute_Lnu(flux, z, alpha)  # 0.0014270422727500343
        SFR = 6.64 * (1e-29) * ((frequency) ** (-alpha_NT)) * Lnu_NT
        if flux_error is not None:
            Lnu_NT_error = compute_Lnu(flux_error, z, alpha)
            SFR_error = 6.64 * (1e-29) * ((frequency) ** (-alpha_NT)) * Lnu_NT_error
        else:
            SFR_error = 0.0

    # if calibration_kind == 'Tabatabaei2017':
    #     '''
    #     There is something wrong for this kind!
    #     '''
    #     Lnu_NT = compute_Lnu(flux, z, alpha)  # 0.0014270422727500343
    #     SFR = 1.11 * 1e-37 * 1e9 * frequency * Lnu_NT
    #     if flux_error is not None:
    #         Lnu_NT_error = compute_Lnu(flux_error, z, alpha)
    #         SFR_error = 1.11 * 1e-37 * 1e9 * frequency * Lnu_NT_error
    #     else:
    #         SFR_error = 0.0

    # if calibration_kind == 'Murphy12':
    #     Te = 1e4
    #     Lnu_NT = compute_Lnu(flux, z, alpha_NT)
    #     SFR = 1e-27 * ( \
    #                 (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
    #                  15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
    #         ) * Lnu_NT

    #     if flux_error is not None:
    #         Lnu_NT_error = compute_Lnu(flux_error, z, alpha_NT)
    #         SFR_error = 1e-27 * ( \
    #                     (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
    #                      15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
    #             ) * Lnu_NT_error
    #     else:
    #         SFR_error = 0.0
    
    if calibration_kind == 'Murphy12':
        Te = 1e4
        # Lnu_tot = compute_Lnu_v2(flux, z, alpha_NT,frequency,nu0)
        Lnu = compute_Lnu(flux, z, alpha_NT)
        SFR = 1*1e-27 * ( \
                    (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                     15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
            ) * Lnu

        if flux_error is not None:
            # Lnu_tot_error = compute_Lnu_v2(flux_error, z, alpha_NT,frequency,nu0)
            Lnu_error = compute_Lnu(flux, z, alpha_NT)
            SFR_error = 1*1e-27 * ( \
                        (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                         15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
                ) * Lnu_error
        else:
            SFR_error = 0.0

    if calibration_kind == 'Tabatabaei2017':
        Te = 1e4
        # Lnu_tot = compute_Lnu_v2(flux, z, alpha_NT,frequency,nu0)
        Lnu = compute_Lnu(flux, z, alpha_NT)
        SFR = 1*1e-27 * ( \
                    (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                     15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
            ) * Lnu

        if flux_error is not None:
            # Lnu_tot_error = compute_Lnu_v2(flux_error, z, alpha_NT,frequency,nu0)
            Lnu_error = compute_Lnu(flux, z, alpha_NT)
            SFR_error = 1*1e-27 * ( \
                        (2.18 * ((Te / (1e4)) ** (0.45)) * (frequency ** (-0.1)) + \
                         15.1 * (frequency ** (alpha_NT))) ** (-1.00) \
                ) * Lnu_error
        else:
            SFR_error = 0.0
    
    if calibration_kind == 'Murphy12_NT':
        Lnu = compute_Lnu_v2(flux, z, alpha_NT,frequency,nu0)
        SFR = 6.64*(1e-29) * ((frequency)**(-alpha_NT)) * Lnu
        
        if flux_error is not None:
            Lnu_error = compute_Lnu_v2(flux_error, z, alpha_NT,frequency,nu0)
            SFR_error = 6.64*1e-29 * ((frequency)**(-alpha_NT)) * Lnu_error
        else:
            SFR_error = 0.0
            
    if calibration_kind == 'Murphy12_T':
        Te = 1e4
        Lnu = compute_Lnu_v2(flux, z, alpha,frequency,nu0)
        SFR = 4.59*1e-28 * ((frequency)**(-alpha)) * ((Te / (1e4)) ** (-0.45)) * Lnu
        
        if flux_error is not None:
            Lnu_error = compute_Lnu_v2(flux_error, z, alpha,frequency,nu0)
            SFR_error = 4.59*1e-28 * ((frequency)**(-alpha)) * ((Te / (1e4)) ** (-0.45)) * Lnu_error
        else:
            SFR_error = 0.0
    print(f'nu =', frequency)
    print(f'alpha =', alpha)
    print(f'SFR =', SFR, '+/-', SFR_error, 'Mo/yr')
    print(f'Lnu =', Lnu, '+/-', Lnu_error, 'erg/s/Hz')
    return (SFR, SFR_error)


def interferometric_decomposition(image1, image2, image3=None,
                                  residual1=None, residual2=None, residual3=None,
                                  iterations=2,
                                  dilation_size=None, std_factor=1.0,sub_bkg=False,
                                  ref_mask=None,sigma=6):
    """
    Peform an interferometric image decomposition using the e-MERLIN and JVLA.

    Parameters
    ----------
    image1 : str
        Path to an e-MERLIN image.
    image2 : str
        Path to a combined image.
    image3 : str, optional
        Path to a JVLA image.
    residual1 : str, optional
        Path to the residual image of the e-MERLIN image.
    residual2 : str, optional
        Path to the residual image of the combined image.
    residual3 : str, optional
        Path to the residual image of the JVLA image.
    iterations : int, optional
        Number of iterations to perform the sigma masking.
    dilation_size : int, optional
        Size of the dilation kernel.
    std_factor : float, optional
        Factor to multiply the standard deviation of the residual image.
    ref_mask : array, optional
        Reference JVLA mask to be used in the masking process.

    """
    def fit_sigma(image1_data, image2_data):
        '''
        Functiont to optmize sigma masking.
        This minimizes the negative values when subtracting e-MERLIN and VLA images.
        '''

        def opt_sigma_old(params):
            """
            This the old implementation of the minimisation function.

            Note: This function works very well, and it is exact. However, it is very slow.
            """
            sigma_level = params['sigma_level']
            # print('Current sigma is:', sigma_level)
            gomask, gmask = mask_dilation(g, cell_size=None,
                                          iterations=iterations,
                                          sigma=sigma_level,
                                          dilation_size=dilation_size_i,
                                          PLOT=False)

            I1mask = gmask * g  # - convolve_fft(gmask,kernel)
            I1mask_name = image_cut_i.replace('.fits', '_I1mask.fits')
            pf.writeto(I1mask_name,
                       I1mask, overwrite=True)

            copy_header(image_cut_i, I1mask_name, I1mask_name)

            M12 = convolve_2D_smooth(I1mask_name,
                                     mode='transfer',
                                     imagename2=image_cut_j,
                                     add_prefix='_M12')
            R12_ = g_next - load_fits_data(M12)  # + np.mean(g[g>3*mad_std(g)])
            #             residual_mask = R12_ - std_factor*std_level
            residual_mask = R12_ - std_factor * np.std(R12_)
            return (residual_mask)

        def opt_sigma(params):
            '''
            Improved version, do not use casa, so do not need to write and read files.
            Using convolve_fft seems to give **almost** the same results as CASA's
            function smooth.

            ** There is a small difference in the relative amplitude when convolving image2 with
            the beam of image1. The factor is exactly the same as the ratio of the beam areas:
            beam_area2(image2) / beam_area2(image1)
            '''
            sigma_level = params['sigma_level']
            # std_factor = params['std_factor']
            # print(f"Current opt values = sigma={sigma_level.value}, std_factor={std_factor.value}")
            print(f"Current opt values = sigma={sigma_level.value}")
            _, dilated_mask = mask_dilation(image1_data, cell_size=None,
                                                 iterations=2,
                                          sigma=sigma_level,
                                          dilation_size=dilation_size_i,
                                          PLOT=False)
            masked_image = dilated_mask * image1_data * ref_mask
            M12 = _fftconvolve_jax(masked_image, PSF_DATA_j)
            # M12 = scipy.signal.fftconvolve(masked_image, PSF_DATA, mode='same')
            mask_M12 = (M12 > 1e-6)
            _R12 = (image2_data - M12 + offset_2) * ref_mask
            # R12 = _R12 + abs(jnp.nanmedian(_R12[_R12<0]))
            # R12 = (image2_data - (M12 + 0*std_factor * bkg_2))*ref_mask
            # residual_mask = 1000*(R12_ - std_factor * std_level) * ref_mask  # *mask_M12#avoid
            return (np.array(_R12).copy())

        bounds_i, bounds_f = 6, 100
        sigma_i, sigma_f = bounds_i, bounds_f
        std_fac_bounds_i, std_fac_bounds_f = 0.99, 1.01
        std_fac_i, std_fac_f = std_fac_bounds_i, std_fac_bounds_f
        pars0 = 6.0, 1.0
        sigma_0, std_fac_0 = pars0
        fit_params = lmfit.Parameters()
        fit_params.add("sigma_level", value=sigma_0, min=sigma_i, max=sigma_f)
        # fit_params.add("std_factor", value=std_fac_0, min=std_fac_i, max=std_fac_f)
        #         fit_params.add("offset", value=1.0, min=0.5, max=2)
        solver_method = 'nelder'
        #         mini = lmfit.Minimizer(opt_sigma,fit_params,max_nfev=5000,nan_policy='omit',reduce_fcn='neglogcauchy')
        mini = lmfit.Minimizer(opt_sigma, fit_params, max_nfev=50000,
                               nan_policy='omit', reduce_fcn='neglogcauchy')
        result_1 = mini.minimize(method=solver_method,tol=1e-12,
                                 options={'xatol': 1e-12,
                                          'fatol': 1e-8,
                                          'adaptive': True})

        result = mini.minimize(method=solver_method, params=result_1.params,
                               tol=1e-8,
                               options={'xatol': 1e-12,
                                        'fatol': 1e-12,
                                        'adaptive': True})
        # result_1 = mini.minimize(method='least_squares',
        #                          tr_solver="exact",
        #                          tr_options={'regularize': True},
        #                          # x_scale='jac', loss="cauchy",
        #                          ftol=1e-14, xtol=1e-14, gtol=1e-14,
        #                          f_scale=1,
        #                          max_nfev=200000,
        #                          verbose=2)
        #
        #
        # result = mini.minimize(method='least_squares',
        #                        params=result_1.params,
        #                        tr_solver="exact",
        #                        tr_options={'regularize': True},
        #                        # x_scale='jac', loss="cauchy",
        #                        ftol=1e-14, xtol=1e-14, gtol=1e-14,
        #                        f_scale=1,
        #                        max_nfev=200000,
        #                        verbose=2)

        # result_1 = mini.minimize(method='differential_evolution',
        #                          # tr_solver="exact",
        #                          # tr_options={'regularize': False},
        #                          # # x_scale='jac', loss="cauchy",
        #                          # ftol=1e-14, xtol=1e-14, gtol=1e-14,
        #                          # f_scale=1,
        #                          # max_nfev=200000,
        #                          verbose=2)
        #
        #
        # result = mini.minimize(method='differential_evolution',
        #                        params=result_1.params,
        #                        # tr_solver="exact",
        #                        # tr_options={'regularize': False},
        #                        # # x_scale='jac', loss="cauchy",
        #                        # ftol=1e-14, xtol=1e-14, gtol=1e-14,
        #                        # f_scale=1,
        #                        # max_nfev=200000,
        #                        verbose=2)
        #         result = mini.minimize(method=solver_method,params=result_1.params,
        #                                  max_nfev=30000, #x_scale='jac',  # f_scale=0.5,
        #                                  tr_solver="exact",
        #                                  tr_options={'regularize': True},
        #                                  ftol=1e-14, xtol=1e-14, gtol=1e-14, verbose=2)

        parFit = result.params['sigma_level'].value #, result.params['std_factor'].value
        Err_parFit = result.params['sigma_level'].stderr #, result.params['std_factor'].stderr
        # resFit = result.residual
        # chisqr = result.chisqr
        # redchi = result.redchi
        # pars = parFit
        # sigma_opt_fit, std_factor_opt_fit = pars
        return (parFit, Err_parFit, result)

    # set image names, they must have the same cell_size,image size AND MUST BE ALIGNED!!!
    image_cut_i = image1  # highest resolution image
    image_cut_j = image2  # intermediate resolution image

    # read data files
    image1_data = load_fits_data(image_cut_i)
    image2_data = load_fits_data(image_cut_j)

    if sub_bkg:
        if residual2 is not None:
            bkg_2 = load_fits_data(residual2)
            offset_2 = 0.5 * bkg_2
        else:
            bkg_1 = sep_background(image_cut_i, apply_mask=False,
                                   show_map=False,use_beam_fraction=True).back()
            bkg_2 = sep_background(image_cut_j, apply_mask=False,
                                   show_map=False,use_beam_fraction=True).back()
            offset_2 = 0.5 * bkg_2

    else:
        if residual2 is not None:
            bkg_2 = load_fits_data(residual2)
            offset_2 = 0.5 * bkg_2
        else:
            bkg_1 = mad_std(image1_data)
            bkg_2 = mad_std(image2_data)
            offset_2 = 0.5*np.std(image2_data)


    psf_name_j = tcreate_beam_psf(image_cut_j,
                                  size=image1_data.shape,
                                  aspect=None)

    correction_factor_ij = beam_area2(image_cut_j) / beam_area2(image_cut_i)
    PSF_DATA_j = load_fits_data(psf_name_j)*correction_factor_ij  # this will result similar results as

    omaj_i, omin_i, _, _, _ = beam_shape(image_cut_i)
    dilation_size_i = int(np.sqrt(omaj_i * omin_i) / (2 * get_cell_size(image_cut_i)))
    print('Usig mask dilution size of (for image i)=', dilation_size_i)

    results = {}
    results_short = {}

    results['imagename_i'] = os.path.basename(image_cut_i)
    results['imagename_j'] = os.path.basename(image_cut_j)

    parFit, Err_parFit, result_mini = fit_sigma(image1_data, image2_data)
    print('Optmized Sigma=', parFit, '+/-', Err_parFit)
    # sigma_level, std_factor = parFit
    sigma_level = parFit

    #Dilate the mask referent to image1 (high-res image).
    I01mask, I1mask = mask_dilation(image1_data, cell_size=get_cell_size(image_cut_i),
                                    iterations=iterations,
                                    sigma=sigma_level,
                                    dilation_size=dilation_size_i, PLOT=True)

    I1mask_data = I1mask * image1_data * ref_mask # - convolve_fft(gmask,kernel)
    # gg = gmask - gaussian_filter(gmask,150)
    I1mask_name = image_cut_i.replace('.fits', '_I1mask.fits')
    pf.writeto(I1mask_name, I1mask_data, overwrite=True)
    copy_header(image_cut_i, I1mask_name, I1mask_name)

    pf.writeto(image_cut_i.replace('.fits', '_mask_bool.fits'), I1mask * 1.0,
               overwrite=True)
    copy_header(image_cut_i, image_cut_i.replace('.fits', '_mask_bool.fits'),
                image_cut_i.replace('.fits', '_mask_bool.fits'))


    M12_data = scipy.signal.fftconvolve(I1mask_data, PSF_DATA_j, mode='same') #+offset_2
    M12 = image_cut_j.replace('.fits', '_M12opt.fits')
    pf.writeto(M12, M12_data, overwrite=True)
    copy_header(image_cut_j, M12, M12)

    # subtract image2 from theta2*image1_mask and save images
    R12_ = load_fits_data(image_cut_j) - M12_data # + 1*
    # std_factor*std_level# std_factor*std_level# + np.mean(g[g>3*mad_std(g)])
    R12 = image_cut_j.replace('.fits', '_R12opt.fits')
    pf.writeto(R12, R12_, overwrite=True)
    copy_header(image_cut_j, R12, R12)
    ###################################
    ###################################
    ###################################
    std_1 = mad_std(load_fits_data(residual1))
    std_2 = mad_std(load_fits_data(residual2))
    _, I1mask = mask_dilation(image1, PLOT=False, dilation_type='disk', rms=std_1,
                                     sigma=result_mini.params['sigma_level'].value,
                                     iterations=2, dilation_size=None)
    _, mask_I1 = mask_dilation(image1, PLOT=False, dilation_type='disk', rms=std_1,
                                     sigma=6.0,
                                     iterations=2, dilation_size=None)
    _, mask_I2 = mask_dilation(image2, PLOT=False, dilation_type='disk', rms=std_2,
                                     sigma=6.0,
                                     iterations=2, dilation_size=None)
    _, mask_M12 = mask_dilation(M12, PLOT=False, dilation_type='disk', rms=std_2,
                                      sigma=6, iterations=2, dilation_size=None)
    _, mask_R12 = mask_dilation(R12, PLOT=False, dilation_type='disk', rms=std_2,
                                      sigma=6, iterations=2, dilation_size=None)

    results['I1_name'] = image1
    results['I2_name'] = image2
    results['R12_name'] = R12
    results['M12_name'] = M12

    print(f"  ++==>> Computing image statistics on Image1...")
    I1_props = compute_image_properties(img=image1,
                                        residual=residual1,
                                        rms=std_1,
                                        apply_mask=False,
                                        show_figure=False,
                                        mask=mask_I1,
                                        last_level=2.0)[-1]

    print(f"  ++==>> Computing image statistics on Image-Mask1...")
    I1mask_props = compute_image_properties(img=I1mask_name,
                                            residual=residual1,
                                            rms=std_1,
                                            apply_mask=False,
                                            show_figure=False,
                                            mask=I1mask,
                                            last_level=2.0)[-1]

    print(f"  ++==>> Computing image statistics on Image2...")
    I2_props = compute_image_properties(img=image2,
                                        residual=residual2,
                                        rms=std_2,
                                        apply_mask=False,
                                        show_figure=False,
                                        mask=mask_I2,
                                        last_level=2.0)[-1]

    print(f"  ++==>> Computing image statistics on R12...")
    R12_props = compute_image_properties(img=R12,
                                         residual=residual2,
                                         rms=std_2,
                                         apply_mask=False,
                                         show_figure=False,
                                         mask=mask_R12,
                                         last_level=2.0)[-1]
    print(f"  ++==>> Computing image statistics on M12...")
    M12_props = compute_image_properties(img=M12,
                                         residual=residual2,
                                         rms=std_2,
                                         show_figure=False,
                                         apply_mask=False,
                                         mask=mask_M12,
                                         last_level=2.0)[-1]

    """
    Store measured properties:
        - total flux density on: I1, I1mask, R12, M12 and I2
        - peak brightness on: I1, I1mask, R12, M12 and I2
        - sizes on: I1, I1mask, R12, M12 and I2
    """

    results['S_I1'] = I1_props['total_flux_mask']
    results['S_I1mask'] = I1mask_props['total_flux_mask']
    results['S_R12'] = R12_props['total_flux_mask']
    results['S_M12'] = M12_props['total_flux_mask']
    results['S_I2'] = I2_props['total_flux_mask']

    results['Speak_I1'] = I1_props['peak_of_flux']
    results['Speak_I1mask'] = I1mask_props['peak_of_flux']
    results['Speak_R12'] = R12_props['peak_of_flux']
    results['Speak_M12'] = M12_props['peak_of_flux']
    results['Speak_I2'] = I2_props['peak_of_flux']

    results['C50radii_I1'] = I1_props['C50radii']
    results['C50radii_I1mask'] = I1mask_props['C50radii']
    results['C50radii_R12'] = R12_props['C50radii']
    results['C50radii_M12'] = M12_props['C50radii']
    results['C50radii_I2'] = I2_props['C50radii']

    results['C95radii_I1'] = I1_props['C95radii']
    results['C95radii_I1mask'] = I1mask_props['C95radii']
    results['C95radii_R12'] = R12_props['C95radii']
    results['C95radii_M12'] = M12_props['C95radii']
    results['C95radii_I2'] = I2_props['C95radii']


    results['diff_SI2_SI1'] = results['S_I2'] - results['S_I1']
    results['ratio_SI1_SI2'] = results['S_I1'] / results['S_I2']
    results['diff_SpeakI2_SpeakI1'] = results['Speak_I2'] - results['Speak_I1']
    results['ratio_SpeakI1_SpeakI2'] = results['Speak_I1'] / results['Speak_I2']

    # omask_i, mask_i = mask_dilation(image_cut_i,
    #                                 cell_size=get_cell_size(image_cut_i),
    #                                 sigma=sigma, iterations=2,
    #                                 dilation_size=dilation_size_i, PLOT=True)

    omaj_j, omin_j, _, _, _ = beam_shape(image_cut_j)
    dilation_size_j = int(
        np.sqrt(omaj_j * omin_j) / (2 * get_cell_size(image_cut_j)))
    print('Usig mask dilution size of (for image j)=', dilation_size_j)

    omask_j, mask_j = mask_dilation(image_cut_j,
                                    cell_size=get_cell_size(image_cut_j),
                                    iterations=2, sigma=sigma,
                                    dilation_size=dilation_size_j, PLOT=False)

    plot_interferometric_decomposition(I1mask_name, image_cut_j,
                                       M12, R12,
                                       vmax_factor=2.0, vmin_factor=0.5,
                                       run_phase='1st',
                                       crop=True, NAME=image_cut_i,
                                       SPECIAL_NAME='_I1_2_M12_R12')

    print('####################################################')
    print('----------------- **** REPORT **** -----------------')
    print('####################################################')
    print(f"Flux Density I1 (high-res)              = {results['S_I1'] * 1000:.2f} mJy")
    print(f"Flux Density I1mask (high-res compact)  = {results['S_I1mask'] * 1000:.2f} mJy")
    print(f"Flux Density M12 (mid-res compact conv) = {results['S_M12'] * 1000:.2f} mJy")
    print(f"Flux Density I2 (mid-res)               = {results['S_I2'] * 1000:.2f} mJy")
    print(f"Flux Density R12 (mid-res extended)     = {results['S_R12'] * 1000:.2f} mJy")
    print('------------------ ************** ------------------')
    print(f"Diff Flux Density I2-I1                 = {results['diff_SI2_SI1'] * 1000:.2f} mJy")
    print(f"Ratio Flux Density I1/I2                = {results['ratio_SI1_SI2']:.2f}")
    print(f"Diff Peak SpeakI2 - SpeakI1             = "
          f"{results['diff_SpeakI2_SpeakI1']*1000:.2f} mJy/beam")
    print(f"Ratio Peak SpeakI1 / SpeakI2            = "
          f"{results['ratio_SpeakI1_SpeakI2']:.2f}")
    print('####################################################')

    if image3 is not None:

        def image_sum(image_cut_k, M13, M23,offset_3):
            def opt_res(params):
                # print(f"Params: a={params['a'].value},"
                #       f"b={params['b'].value},"
                #       f"c={params['c'].value}")
                optimisation = (params['a'] * M23_data + params['b'] * M13_data + params['c'] *
                                offset_3)
                I3ext = (image3_data - optimisation)
                #     I3ext = image_data - (params['a']*M23_data-params['c']*zero_off)- (params['b']*M13_data-params['d']*zero_off)
                # return (I3ext+abs(np.nanmin(I3ext)))
                return (I3ext)

            def opt_sub_3(image_data, M23_data, M13_data, bkg_3):
                bounds_i = 0.1, 0.1, -1.0  # ,-3.0
                bounds_f = 10.01, 10.01, 1.0  # ,+3.0
                pars0 = 1.0, 1.0, 0.5  # ,0.0
                ai, bi, ci = bounds_i
                af, bf, cf = bounds_f
                a0, b0, c0 = pars0
                fit_params = lmfit.Parameters()
                fit_params.add("a", value=a0, min=ai, max=af)
                fit_params.add("b", value=b0, min=bi, max=bf)
                fit_params.add("c", value=c0, min=ci, max=cf)
                #     fit_params.add("d", value=d0, min=di, max=df)
                solver_method = 'least_squares'
                mini = lmfit.Minimizer(opt_res, fit_params, max_nfev=50000,
                                       nan_policy='omit', reduce_fcn='neglogcauchy')

                result_1 = mini.minimize(method='nelder', tol=1e-8,
                                         options={'xatol': 1e-8,
                                                  'fatol': 1e-8,
                                                  'adaptive': True})

                result = mini.minimize(method='nelder', params=result_1.params,
                                       tol=1e-8,
                                       options={'xatol': 1e-8,
                                                'fatol': 1e-8,
                                                'adaptive': True})

                #             result_1 = mini.minimize(method=solver_method)
                # result_1 = mini.minimize(method=solver_method,
                #                          tr_solver="exact",
                #                          tr_options={'regularize': True},
                #                          x_scale='jac', loss="cauchy",
                #                          ftol=1e-12, xtol=1e-12, gtol=1e-12,
                #                          f_scale=1.0,
                #                          max_nfev=5000, verbose=2)
                # result = mini.minimize(method=solver_method,
                #                        params=result_1.params,
                #                        tr_solver="exact",
                #                        tr_options={'regularize': True},
                #                        x_scale='jac', loss="cauchy",
                #                        ftol=1e-12, xtol=1e-12, gtol=1e-12,
                #                        f_scale=1.0,
                #                        max_nfev=5000, verbose=2)

                parFit = result.params['a'].value, result.params['b'].value, \
                    result.params['c'].value  # ,result.params['d'].value
                Err_parFit = result.params['a'].stderr, result.params['b'].stderr, \
                    result.params['c'].stderr  # ,result.params['d'].stderr
                # resFit = result.residual
                # chisqr = result.chisqr
                # redchi = result.redchi
                # pars = parFit
                # res_opt = pars
                return (result, parFit, Err_parFit)

            image3_data = mask_k * load_fits_data(image_cut_k)
            M13_data = load_fits_data(M13)
            M23_data = load_fits_data(M23)

            result_mini_I3, parFit, Err_parFit = \
                opt_sub_3(image3_data, M23_data, M13_data, offset_3)
            a, b, c = parFit
            a_err, b_err, c_err = Err_parFit
            print(f"Linear Image Combination Optmisation: "
                  f"        (a,b,c)={parFit}+/-{Err_parFit}")

            model_total = a * load_fits_data(M23) + b * load_fits_data(M13) + 1 * c * offset_3
            M123 = image_cut_j.replace('.fits', '_M123_opt.fits')
            pf.writeto(M123, model_total, overwrite=True)
            copy_header(image_cut_k, M123, M123)

            model_comp = 0 * a * load_fits_data(M23) + b * load_fits_data(M13) + 0 * c * offset_3 / 2
            Mcomp = image_cut_j.replace('.fits', '_M13_opt.fits')
            pf.writeto(Mcomp, model_comp, overwrite=True)
            copy_header(image_cut_k, Mcomp, Mcomp)

            model_extended = a * load_fits_data(M23)# + c * bkg_3 / 2
            Mext = image_cut_j.replace('.fits', '_M23_opt.fits')
            pf.writeto(Mext, model_extended, overwrite=True)
            copy_header(image_cut_k, Mext, Mext)

            I3re = load_fits_data(image_cut_k) - (a * load_fits_data(M23) + b * load_fits_data(M13) + 0 * c * offset_3)
            I3_RT = image_cut_j.replace('.fits', '_RT.fits')
            I3ext = load_fits_data(image_cut_k) - (0 * a * load_fits_data(M23) + b * load_fits_data(M13) +  0 * c * offset_3 / 2)
            I3ext_name = image_cut_j.replace('.fits', '_residual_extended.fits')
            I3comp = load_fits_data(image_cut_k) - (a * load_fits_data(M23) + 0 * b * load_fits_data(M13) + 0 * c * offset_3 / 2)
            I3comp_name = image_cut_j.replace('.fits', '_residual_comp.fits')

            pf.writeto(I3_RT, I3re, overwrite=True)
            pf.writeto(I3ext_name, I3ext, overwrite=True)
            pf.writeto(I3comp_name, I3comp, overwrite=True)
            copy_header(image_cut_k, I3_RT, I3_RT)
            copy_header(image_cut_k, I3ext_name, I3ext_name)
            copy_header(image_cut_k, I3comp_name, I3comp_name)

            return (M123, Mcomp, Mext, I3_RT, I3ext_name, I3comp_name, result_mini_I3, parFit, Err_parFit)


        print('Running decomposition for Image3...')
        image_cut_k = image3  # lowest resolution image
        results['imagename_k'] = os.path.basename(image_cut_k)

        omaj_k, omin_k, _, _, _ = beam_shape(image_cut_k)
        dilation_size_k = int(
            np.sqrt(omaj_k * omin_k) / (2 * get_cell_size(image_cut_k)))
        print('Usig mask dilution size of (for image k)=', dilation_size_k)

        omask_k, mask_k = mask_dilation(image_cut_k,
                                        cell_size=get_cell_size(image_cut_k),
                                        sigma=6, iterations=2,
                                        dilation_size=dilation_size_k,
                                        PLOT=False)

        image3_data = load_fits_data(image_cut_k)
        bak = beam_area2(image_cut_k, cellsize=get_cell_size(image_cut_k))
        #         image_cut_k = imagelist[k]#cut_image2(imagelist[k],size=(1024,1024))

        psf_name_k = tcreate_beam_psf(image_cut_j,
                                      size=image1_data.shape,
                                      aspect=None)

        correction_factor_jk = beam_area2(image_cut_k) / beam_area2(image_cut_j)
        correction_factor_ik = beam_area2(image_cut_k) / beam_area2(image_cut_i)
        PSF_DATA_k = load_fits_data(psf_name_k)  # this will result similar results as

        if sub_bkg:
            bkg_3 = sep_background(image_cut_k, apply_mask=False,
                                   show_map=False, use_beam_fraction=True).back()
            offset_3 = 0.5 * bkg_3

        else:
            if residual3 is not None:
                bkg_3 = load_fits_data(residual3)
                offset_3 = 0.5 * bkg_3
            else:
                bkg_3 = mad_std(image3_data)
                offset_3 = 0.5 * bkg_3


        #read the result of R12, but without the offest component!
        R12_data = load_fits_data(R12)# - offset_2
        R12_mask = mask_j  # (R12_data>1*mad_std(R12_data))
        # gg = gmask - convolve_fft(gmask,kernel)
        R12_data_mask = R12_mask * R12_data  # - gaussian_filter(gmask,10)

        #Generate new file for a masked region of R12.
        R12mask = image_cut_j.replace('.fits', '_R12_mask.fits')
        pf.writeto(R12mask, R12_data_mask, overwrite=True)
        copy_header(image_cut_j, R12mask, R12mask)

        # Convolve R12mask with the beam of image3.
        M23_data = scipy.signal.fftconvolve(R12_data_mask,
                                            PSF_DATA_k*correction_factor_jk,
                                            mode='same')
        M23 = image_cut_j.replace('.fits', '_M23.fits')
        pf.writeto(M23, M23_data, overwrite=True)
        copy_header(image_cut_k, M23, M23)

        # # Convolving I1mask with the beam of image3.
        # M13_data = scipy.signal.fftconvolve(I1mask_data,
        #                                     PSF_DATA_k*correction_factor_ik,
        #                                     mode='same')
        # M13 = image_cut_k.replace('.fits', '_M13.fits')
        # pf.writeto(M13, M13_data, overwrite=True)
        # copy_header(image_cut_k, M13, M13)

        # Convolving I1mask with the beam of image3.
        M13_data = scipy.signal.fftconvolve(M12_data,
                                            PSF_DATA_k*correction_factor_jk,
                                            mode='same')
        M13 = image_cut_k.replace('.fits', '_M13.fits')
        pf.writeto(M13, M13_data, overwrite=True)
        copy_header(image_cut_k, M13, M13)

        M123_opt, M13_opt, M23_opt, I3_RT, \
            I3ext_name, I3comp_name, result_mini_I3, parFit, Err_parFit = \
            image_sum(image_cut_k, M13, M23,offset_3)

        a, b, c = parFit
        a_err, b_err, c_err = Err_parFit

        results['I3_name'] = image3
        results['M23_name'] = M23_opt
        results['M13_name'] = M13_opt
        results['M123_name'] = M123_opt
        results['I3ext_name'] = I3ext_name
        results['I3_RT_name'] = I3_RT

        results['a'] = a
        results['b'] = b
        results['c'] = c
        results['a_err'] = a_err
        results['b_err'] = b_err
        results['c_err'] = c_err

        std_3 = mad_std(load_fits_data(residual3))

        _, mask_M13 = mask_dilation(M13_opt, PLOT=False, dilation_type='disk', rms=std_3,
                                    sigma=6.0, iterations=2, dilation_size=None)
        _, mask_I3 = mask_dilation(image3, PLOT=False, dilation_type='disk', rms=std_3,
                                   sigma=6.0,
                                   iterations=2, dilation_size=None)
        _, mask_I3ext = mask_dilation(I3ext_name, PLOT=False, dilation_type='disk', rms=std_3,
                                   sigma=6.0,
                                   iterations=2, dilation_size=None)

        print(f"  ++==>> Computing image statistics on Image3...")
        I3_props = compute_image_properties(img=image3,
                                            residual=residual3,
                                            rms=std_3,
                                            apply_mask=False,
                                            show_figure=False,
                                            mask=mask_I3,
                                            last_level=2.0)[-1]
        print(f"  ++==>> Computing image statistics on extended emission...")
        I3ext_props = compute_image_properties(img=I3ext_name,
                                               residual=residual3,
                                               rms=std_3,
                                               apply_mask=False,
                                               show_figure=False,
                                               mask=mask_I3,
                                               last_level=2.0)[-1]

        print(f"  ++==>> Computing image statistics on M23...")
        M23_props = compute_image_properties(img=M23_opt,
                                               residual=residual3,
                                               rms=std_3,
                                               apply_mask=False,
                                               show_figure=False,
                                               mask=mask_I3,
                                               last_level=2.0)[-1]
        print(f"  ++==>> Computing image statistics on M13...")
        M13_props = compute_image_properties(img=M13_opt,
                                             residual=residual3,
                                             rms=std_3,
                                             apply_mask=False,
                                             show_figure=False,
                                             mask=mask_M13,
                                             last_level=2.0)[-1]
        print(f"  ++==>> Computing image statistics on residual RT...")
        I3_RT_props = compute_image_properties(img=I3_RT,
                                             residual=residual3,
                                             rms=std_3,
                                             apply_mask=False,
                                             show_figure=True,
                                             mask=mask_I3,
                                             last_level=2.0)[-1]

        results['S_I3'] = I3_props['total_flux_mask']
        results['S_I3ext'] = I3ext_props['total_flux_mask']
        results['S_M23'] = M23_props['total_flux_mask']
        results['S_M13'] = M13_props['total_flux_mask']
        results['S_I3RT'] = I3_RT_props['total_flux_mask']

        results['Speak_I3'] = I3_props['peak_of_flux']
        results['Speak_I3ext'] = I3ext_props['peak_of_flux']
        results['Speak_M23'] = M23_props['peak_of_flux']
        results['Speak_M13'] = M13_props['peak_of_flux']
        results['Speak_I3RT'] = I3_RT_props['peak_of_flux']

        results['C50radii_I3'] = I3_props['C50radii']
        results['C50radii_I3ext'] = I3ext_props['C50radii']
        results['C50radii_M23'] = M23_props['C50radii']
        results['C50radii_M13'] = M13_props['C50radii']
        results['C50radii_I3RT'] = I3_RT_props['C50radii']

        results['C95radii_I3'] = I3_props['C95radii']
        results['C95radii_I3ext'] = I3ext_props['C95radii']
        results['C95radii_M23'] = M23_props['C95radii']
        results['C95radii_M13'] = M13_props['C95radii']
        results['C95radii_I3RT'] = I3_RT_props['C95radii']

        results['diff_SI3_SI1'] = results['S_I3'] - results['S_I1']
        results['ratio_SI1_SI3'] = results['S_I1'] / results['S_I3']
        results['diff_SpeakI3_SpeakI1'] = results['Speak_I3'] - results['Speak_I1']
        results['ratio_SpeakI1_SpeakI3'] = results['Speak_I1'] / results['Speak_I3']

        print('####################################################')
        print('----------------- **** REPORT **** -----------------')
        print('####################################################')
        print(f"Flux Density I3 (low-res)             = {results['S_I3'] * 1000:.2f} mJy")
        print(f"Flux Density Comp M13 (low-res)       = {results['S_M13'] * 1000:.2f} mJy")
        print(f"Flux Density Ext M23 (low-res)        = {results['S_M23'] * 1000:.2f} mJy")
        print(f"Flux Density Extended (I3 low-res)    = {results['S_I3ext'] * 1000:.2f} mJy")
        print(f"Flux Density Residual (RT low-res)    = {results['S_I3RT'] * 1000:.2f} mJy")
        print('------------------ ************** ------------------')
        print(f"Diff Flux Density I3-I1                 = {results['diff_SI3_SI1'] * 1000:.2f} mJy")
        print(f"Ratio Flux Density I1/I3                = {results['ratio_SI1_SI3']:.2f}")
        print(f"Diff Peak SpeakI3 - SpeakI1             = "
              f"{results['diff_SpeakI3_SpeakI1'] * 1000:.2f} mJy/beam")
        print(f"Ratio Peak SpeakI1 / SpeakI3            = "
              f"{results['ratio_SpeakI1_SpeakI3']:.2f}")
        print('####################################################')

        #
        # results['S_I3_23_res_full_py'] = np.sum(I3_residual_23) / bak
        # results['S_I3_23_res_mask_py'] = np.sum(mask_k * I3_residual_23) / bak
        # results['a'] = a
        # results['b'] = b
        # results['c'] = c
        #
        df = pd.DataFrame.from_dict(results, orient='index').T
        df.to_csv(image3.replace('.fits', '_interf_decomposition.csv'),
                  header=True,
                  index=False)

        plot_interferometric_decomposition(R12mask, image_cut_k,
                                           M123_opt,  # M23,
                                           I3_RT,
                                           vmin0=mad_std(load_fits_data(image_cut_j)),
                                           vmax_factor=0.5, vmin_factor=2.0,
                                           run_phase='2nd',
                                           crop=False, box_size=512,
                                           NAME=image_cut_j,
                                           SPECIAL_NAME='_R12_I2_M123_I3re')

        plot_interferometric_decomposition(R12mask, image_cut_k, M13_opt,
                                           # M23,
                                           I3ext_name,
                                           vmax_factor=0.5, vmin_factor=2.0,
                                           crop=False, box_size=512,
                                           run_phase='compact',
                                           NAME=image_cut_j,
                                           SPECIAL_NAME='_R12_I3_Mcomp_I3ext')

    if image3 is None:
        return (result_mini, results, results_short, I1mask_data, I1mask_name, R12, M12)
    else:
        return (result_mini, result_mini_I3, results, I1mask, I1mask_name, R12, M12, M123_opt,
                M13_opt, M23_opt, I3_RT, I3ext_name)

image_decomposition = deprecated("image_decomposition",
                              "interferometric_decomposition")(interferometric_decomposition)

def perform_interferometric_decomposition(imagelist_em, imagelist_comb,
                                          imagelist_vla, residuallist_vla,
                                          idx_em, idx_vla, idxs_em, z,
                                          std_factor=0.5, ref_mask=None,
                                          sigma=6):

    """
    Perform interferometric decomposition of the images in the list.

    Parameters
    ----------
    imagelist_em : list
        List of e-merlin images.
    imagelist_comb : list
        List of combined images.
    imagelist_vla : list
        List of jvla images.
    residuallist_vla : list
        List of jvla residual images.
    idx_em : int
        Index of the e-merlin image.
    idx_vla : int
        Index of the jvla image.
    idxs_em : list
        List of indices of the e-merlin images.
    z : float
        Redshift of the source.
    std_factor : float; default = 0.5
        Factor to multiply the standard deviation of the image.
        Do not modify this parameter, unless you know what you are doing. See
        Fig. 5 in the paper for more details.
    ref_mask : array
        Reference mask.
    sigma : float; default = 6
        Sigma level for the mask dilation.

    """
    int_results = []
    I1_MASK_LIST = []
    R12_MASK_LIST = []
    R12_LIST = []
    M12_LIST = []
    M123_LIST = []
    M23_Mext_opt_LIST = []
    I3RE_RT_LIST = []
    M13_Mcomp_opt_LIST = []
    I3RE_LIST = []
    I3EXT_LIST = []
    IMAGELIST_ERROR = []
    imagelist_em_short = imagelist_em[idxs_em]
    rms = mad_std(load_fits_data(residuallist_vla[idx_vla]))
    EXTENDED_PROPS = []
    EXTENDED_2nd_PROPS = []
    COMP_EM_PROPS = []
    COMP_PROPS = []
    COMP_VLA_PROPS = []
    for l in range(0, len(imagelist_em_short)):
        for j in tqdm(range(0, len(imagelist_comb))):
            try:
                #                 i = idx_em # e-merlin image
                k = idx_vla  # almost pure jvla image, but needs to have the same cellsize as of the e-merlin one
                result_mini, results, results_short, Imask, I1mask_name, I1mask_2, R12, R12conv, M12, \
                    M123_opt, Mcomp_opt, Mext_opt, I3re_name, I3ext_name, I3comp_name, I3_residual_23_name = image_decomposition(
                    #                 image1 = imagelist_comb[2],
                    image1=imagelist_em_short[l],
                    image2=imagelist_comb[j],
                    image3=imagelist_vla[k],
                    std_factor=std_factor, dilation_size=None, iterations=2,
                    ref_mask=ref_mask, sigma=sigma
                )
                int_results.append(results)
                I1_MASK_LIST.append(I1mask_name)
                R12_MASK_LIST.append(R12conv)
                M12_LIST.append(M12)
                R12_LIST.append(R12)
                M123_LIST.append(M123_opt)
                M23_Mext_opt_LIST.append(Mext_opt)
                M13_Mcomp_opt_LIST.append(Mcomp_opt)
                I3RE_RT_LIST.append(I3re_name)
                I3EXT_LIST.append(I3ext_name)

                rms_i = mad_std(load_fits_data(imagelist_em_short[l]))
                rms_j = mad_std(load_fits_data(R12))

                _, mask_extended = mask_dilation(imagelist_comb[j],
                                                 # imagelist_vla[k],
                                                 sigma=sigma, iterations=2,
                                                 dilation_size=None, PLOT=True)
                #                 _, mask_extended_2nd = mask_dilation(I3ext_name,#imagelist_vla[k],
                #                              sigma=6, iterations=2,
                #                                 dilation_size=None, PLOT=True)

                _, mask_compact = mask_dilation(M12, rms=rms,
                                                # we have to give a rms for model-based images
                                                sigma=sigma, iterations=2,
                                                dilation_size=None, PLOT=True)
                results_extended = run_analysis_list([R12], residuallist_vla[k],
                                                     imagelist_vla[k],
                                                     z,
                                                     mask_extended, rms_j, sigma=sigma)
                results_extended_2nd = run_analysis_list([I3ext_name],
                                                         residuallist_vla[k],
                                                         imagelist_vla[k],
                                                         z,
                                                         ref_mask, rms, sigma=sigma)

                results_compact = run_analysis_list([M12], residuallist_vla[k],
                                                    imagelist_vla[k],
                                                    z,
                                                    mask_compact, rms_j, sigma=sigma)
                results_compact_vla = run_analysis_list([Mcomp_opt],
                                                        residuallist_vla[k],
                                                        imagelist_vla[k],
                                                        z,
                                                        ref_mask, rms, sigma=sigma)

                results_compact_EM = run_analysis_list([I1mask_name],
                                                       residuallist_vla[k],
                                                       imagelist_vla[k],
                                                       z,
                                                       mask_compact, rms_i, sigma=sigma)
                EXTENDED_PROPS.append(results_extended)
                EXTENDED_2nd_PROPS.append(results_extended_2nd)
                COMP_EM_PROPS.append(results_compact_EM)
                COMP_PROPS.append(results_compact)
                COMP_VLA_PROPS.append(results_compact_vla)

            except:
                print('Error on minimising image:', imagelist_comb[j])
                IMAGELIST_ERROR.append(imagelist_comb[j])
    return (
    int_results, EXTENDED_PROPS, EXTENDED_2nd_PROPS, COMP_PROPS, COMP_EM_PROPS,
    COMP_VLA_PROPS,I1_MASK_LIST,
    R12_MASK_LIST, M12_LIST, R12_LIST, M123_LIST, M23_Mext_opt_LIST,
    M13_Mcomp_opt_LIST, I3RE_RT_LIST, I3EXT_LIST)





"""
#Utils
"""

def get_vis_amp_uvwave(vis,spwid,avg_in_time=True,wantedpol='RR,LL'):
    msmd.open(vis)
    scans = msmd.scansforintent('*TARGET*').tolist()
    msmd.close()

    ms.open(vis)
    ms.selectinit(datadescid=spwid)
    ms.select({'scan_number': scans})
    ms.selectpolarization(wantedpol=wantedpol)
    # ms.selectchannel(nchan=1)
    mydata = ms.getdata(['time', 'amplitude', 'axis_info',
                         'u', 'v',
                         'flag'], ifraxis=True)
    ms.close()

    freq_axis = mydata['axis_info']['freq_axis']['chan_freq']
    mydata['amplitude'][mydata['flag']] = np.nan
    antsel = np.ones_like(mydata['axis_info']['ifr_axis']['ifr_shortname'], dtype='bool')

    amp_avg_time = np.nanmean(mydata['amplitude'].T, axis=0)
    amp_avg_chan = np.nanmean(mydata['amplitude'].T, axis=2)
    amp_avg_chan_corr = np.nanmean(amp_avg_chan, axis=2)
    amp_avg_chan_corr_time = np.nanmean(amp_avg_chan_corr, axis=0)

    # freq_axis
    lightspeed = 299792458.0  # speed of light in m/s
    wavelen = (lightspeed / np.mean(freq_axis, axis=0)) * 1e3

    uu = mydata['u'].T / wavelen
    vv = mydata['v'].T / wavelen
    uvwave = np.sqrt(vv ** 2.0 + uu ** 2.0)

    if avg_in_time==True:
        uvwave_final = np.nanmean(uvwave,axis=0)
        amp_avg_final = amp_avg_chan_corr_time
    else:
        uvwave_final = uvwave
        amp_avg_final = amp_avg_chan_corr
    return(uvwave_final[antsel],amp_avg_final[antsel])




def read_spws(vis):
    tb.open(vis + "/DATA_DESCRIPTION")
    SPECTRAL_WINDOW_ID = tb.getcol("SPECTRAL_WINDOW_ID")
    tb.close()
    # print(SPECTRAL_WINDOW_ID)
    return (SPECTRAL_WINDOW_ID)

def read_spws_v2(vis):
    spwtab = tb.open(vis+'::SPECTRAL_WINDOW')
    tb.close()
    return(spwtab)

def read_colsnames(vis):
    tb.open(vis)
    colnames = tb.colnames()
    tb.close()
    # print(colnames)
    return (colnames)


def query_ms(vis, spwid, query_args=['UVW', 'FLAG']):
    ms.open(vis)
    # select the key
    ms.selectinit(datadescid=spwid)
    query_results = ms.getdata(query_args)
    ms.selectinit(reset=True)
    ms.close()
    return (query_results)

def get_uvwave(vis, spw_id):
    from astropy.constants import c
    ms.open(vis)
    ms.selectinit(spw_id)
    d = ms.getdata(["uvw"])
    ms.done()
    msmd.open(vis)
    chan_freq = msmd.chanfreqs(spw_id)
    msmd.done()
    nchan = len(chan_freq)
    # print(nchan)
    # print(chan_freq)

    # d["uvw"] is an array of float64 with shape [3, nvis]
    u, v, w = d["uvw"]  # unpack into len nvis vectors
    broadcast = np.ones((nchan, 1))
    uu = u * broadcast
    vv = v * broadcast
    ww = w * broadcast

    wavelengths = c.value / chan_freq[:, np.newaxis]  # m
    uw = 1e-3 * uu / wavelengths  # [klambda]
    vw = 1e-3 * vv / wavelengths  # [klambda]
    return (uw, vw, wavelengths)

def get_blinfo():
    return dict([((x, y), np.where((ant1 == x) & (ant2 == y))[0]) 
                 for x in ant_uniq for y in ant_uniq if y > x])


def get_uvwave_tab(vis,index_freq=None):
    # spw info
    def get_blinfo():
        return dict([((x, y), np.where((ant1 == x) & (ant2 == y))[0]) 
                     for x in ant_uniq for y in ant_uniq if y > x])

    lightspeed = 299792458.0 # speed of light in m/s
    
    # spwtab = tb.open(vis+'::SPECTRAL_WINDOW')
    # chan_freq = tb.getcol('CHAN_FREQ').flatten()
    # tb.close()
    ms.open(vis)
    mydata = ms.getdata(['axis_info'])
    chan_freq = mydata['axis_info']['freq_axis']['chan_freq'].flatten()
    ms.close()

    #vis and antenna info
    tb.open(vis)
    uvw = tb.getcol('UVW')
    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    flag = tb.getcol('FLAG')
    flag_row = tb.getcol('FLAG_ROW')
    tb.close()
    tb.open(vis+'::ANTENNA')
    station = tb.getcol('STATION')
    tb.close()
    
    ant_uniq = np.unique(np.hstack((ant1, ant2)))
    bldict = get_blinfo()
    if index_freq is not None:
        wavelen = (lightspeed/chan_freq[index_freq])*1e3 # in Klambda
    else:
        wavelen = (lightspeed/chan_freq.mean())*1e3 # in Klambda
    # flag_final = np.logical_or(flag, flag_row[:,np.newaxis,np.newaxis])
    flag_final = np.logical_or(flag, flag_row[np.newaxis,np.newaxis,:])
    # flag_final_VLA = np.logical_or(flag_VLA, flag_row_VLA[:,np.newaxis,np.newaxis])
    # wavelen = (lightspeed/chan_freq.mean())*1e3
    uwave = uvw[0]/wavelen
    vwave = uvw[1]/wavelen
    uvwave = np.sqrt(uwave**2.0 + vwave**2.0)
    return(uvw,uvwave,wavelen)

    
    
def plot_uv_wavelength(vis, color='black', label=None,
                       chunk_size=int(4), downsample_factor=int(50),
                      fig=None, ax=None, alpha=1.0, figsize=(10,10)):
    """
    Plot UV wavelength (uwave vs vwave) for radio interferometric data.
    
    Parameters:
    -----------
    vis : str
        Path to measurement set (MS) file
    color : str, optional
        Color of the plot points
    chunk_size : int, optional
        Size of chunks for frequency averaging
    downsample_factor : int, optional
        Factor by which to downsample the data points
    fig : matplotlib.figure.Figure, optional
        Existing figure to plot on
    ax : matplotlib.axes.Axes, optional
        Existing axes to plot on
    alpha : float, optional
        Transparency of plot points
    figsize : tuple, optional
        Figure size if creating new figure
        
    Returns:
    --------
    fig, ax : tuple
        Matplotlib figure and axes objects
    """
    LIGHT_SPEED = 299792458.0  # speed of light in m/s
    
    def get_baseline_info(ant1, ant2, ant_uniq):
        """Create dictionary mapping baseline pairs to data indices."""
        return {(x, y): np.where((ant1 == x) & (ant2 == y))[0]
                for x in ant_uniq for y in ant_uniq if y > x}
    
    msmd = casatools.msmetadata()
    ms = casatools.ms()
    tb = casatools.table()
    print(f"++==> Reading info of visibility: {vis}")
    ms.open(vis)
    mydata = ms.getdata(['axis_info'])
    chan_freq = mydata['axis_info']['freq_axis']['chan_freq'].flatten()
    ms.close()
    
    #vis and antenna info
    tb.open(vis)
    uvw = tb.getcol('UVW')
    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    flag = tb.getcol('FLAG')
    flag_row = tb.getcol('FLAG_ROW')
    tb.close()
    tb.open(vis+'::ANTENNA')
    station = tb.getcol('STATION')
    tb.close()
    
    # Get unique antennas and create baseline dictionary
    ant_uniq = np.unique(np.hstack((ant1, ant2)))
    baseline_dict = get_baseline_info(ant1, ant2, ant_uniq)
    
    # Create station labels dictionary
    labels = {bl: f"{station[bl[0]]}-{station[bl[1]]}" 
             for bl in baseline_dict.keys()}
    
    # Setup plotting
    if fig is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    
    # Prepare frequency chunks for averaging
    n_chunks = len(chan_freq) // chunk_size
    freq_chunks = chan_freq.reshape((n_chunks, chunk_size))
    chunk_freqs = freq_chunks.mean(axis=1)
    wavelengths = (LIGHT_SPEED / chunk_freqs.reshape(-1, 1)) * 1e3  # in Klambda

    total_baselines = len(list(baseline_dict.keys()))
    print(f"      Processing baselines for uv plot.")
    print(f"      Total number of baselines = {total_baselines}")
    # Plot each baseline
    kk = 0 
    for baseline in tqdm(list(baseline_dict.keys())):
        # Get baseline indices
        bl_indices = baseline_dict[baseline]
        
        # Extract u,v coordinates
        u = uvw[0, bl_indices]
        v = uvw[1, bl_indices]
        
        # Create conjugate points (including a NaN separator)
        # u_points = np.tile(np.hstack([u, np.nan, -u]), (n_chunks, 1)) / wavelengths
        # v_points = np.tile(np.hstack([v, np.nan, -v]), (n_chunks, 1)) / wavelengths
        uwave_points = np.tile(np.hstack([u, -u]), (n_chunks, 1)) / wavelengths
        vwave_points = np.tile(np.hstack([v, -v]), (n_chunks, 1)) / wavelengths
        
        if kk < 1:
            # n_points = 2*len(chan_freq)*uvw.shape[-1]*total_baselines
            n_points = 2*vwave_points.shape[-1]*total_baselines
            print(f"      Total number of uv-lambda points = {n_points/1e6} million")
            print(f"      Total number of downsampled uv-lambda points to plot = {n_points/(chunk_size*downsample_factor*1e6)} million.")
        # Plot with downsampling
        ax.plot(uwave_points[:, ::downsample_factor], 
                vwave_points[:, ::downsample_factor],
                '.', markersize=0.1, color=color, alpha=alpha)
        kk = kk+1
    
    # if label is not None:
        """
        Plot the last set of points to create a label.
        """
        # ax.plot(uwave_points[0, 0], 
        #         vwave_points[0, 0],
        #         '.', markersize=0.1, color=color, 
        #         label = label,
        #         alpha=alpha)
        # ax.legend(loc='upper right')
        # On the last iteration, add a custom legend entry with a larger marker size
        
        # legend_marker = Line2D([0], [0], marker='.', color='w', label=label, 
        #                     markerfacecolor=color, markersize=5, alpha=1.0)
        # ax.add_artist(ax.legend(handles=[legend_marker], loc='upper right', framealpha=0.5))
        
    # ax.set_xlabel('U (kλ)')
    # ax.set_ylabel('V (kλ)')
    # ax.set_title('UV Coverage')
    # ax.grid(True, alpha=0.3)
    # ax.set_aspect('equal')
    
    return fig, ax
    
def plot_uwave_vwave(vis,color=None,fig=None,ax=None,alpha=1.0,
                     chunk_size = int(8), #channel average
                     downsample_factor = 200,
                     figsize=(8,8),
                     uv_ranges=[],
                     plot_axes = True,
                     save_figure=False,
                     show_figure=False,
                     title_text = None,
                     figure_save_name='uvwave_coverage.jpg'
                    ):
    lightspeed = 299792458.0 # speed of light in m/s
    def get_blinfo():
        return dict([((x, y), np.where((ant1 == x) & (ant2 == y))[0]) 
                    for x in ant_uniq for y in ant_uniq if y > x])
    
#     spwtab = tb.open(vis+'::SPECTRAL_WINDOW')
#     chan_freq = tb.getcol('CHAN_FREQ').flatten()
#     tb.close()
    msmd = casatools.msmetadata()
    ms = casatools.ms()
    tb = casatools.table()
    print(f"++==> Reading info of visibility: {vis}")
    
    msmd.open(vis)
    bandwidth = msmd.bandwidths()
    nspw = len(bandwidth)
    chan_freqs_all = np.empty(nspw, dtype=object)
    spws_freq = np.zeros(nspw)

    for nch in range(nspw):
        chan_freqs_all[nch] = msmd.chanfreqs(nch)
        spws_freq[nch] = np.mean(chan_freqs_all[nch])
    msmd.done()
    chan_freq = np.concatenate(chan_freqs_all)
    
    #vis and antenna info
    tb.open(vis)
    uvw = tb.getcol('UVW')
    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    flag = tb.getcol('FLAG')
    flag_row = tb.getcol('FLAG_ROW')
    tb.close()
    tb.open(vis+'::ANTENNA')
    station = tb.getcol('STATION')
    tb.close()

    ant_uniq = np.unique(np.hstack((ant1, ant2)))
    bldict = get_blinfo()
    wavelen = (lightspeed/chan_freq.mean())*1e6 # in Klambda
    # flag_final = np.logical_or(flag, flag_row[:,np.newaxis,np.newaxis])
    flag_final = np.logical_or(flag, flag_row[np.newaxis,np.newaxis,:])
    # flag_final_VLA = np.logical_or(flag_VLA, flag_row_VLA[:,np.newaxis,np.newaxis])


    labels = {}
    for bl in bldict.keys():
        labels[bl] = station[bl[0]]+'-'+station[bl[1]]

    # plt.ion() 
    if fig is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    
    def chunk_average(arr, chunk_size):
        # Reshape the array to have chunks of size `chunk_size`
        reshaped = arr.reshape(arr.shape[0], -1, chunk_size)
        # Average over each chunk along the last axis
        return reshaped.mean(axis=2)
    
    def downsample(arr, factor):
        return arr[:, ::factor]
    
    chunk_size_plot = 1000
    
    
    
    reshaped_size = int(len(chan_freq)/chunk_size)
    avg_chan = chan_freq.reshape((reshaped_size, chunk_size))
    chunk_chan_mean = avg_chan.mean(axis=1)
    chunk_chan_mean_broad = np.repeat(chunk_chan_mean, chunk_size)
    chunk_chan_mean_column = chunk_chan_mean.reshape(-1, 1)
    UUU = []
    VVV = []
    kk = 0
    # colors = [
    #     '#1F77B4', '#FF7F0E', '#2CA02C', '#D62728',
    #     '#9467BD', '#8C564B', '#E377C2', '#7F7F7F',
    #     '#BCBD22', '#17BECF']
    colors = [
        '#E69F00', '#56B4E9', '#009E73', '#F0E442',
        '#0072B2', '#D55E00', '#CC79A7', '#000000']
    u_points_max = []
    v_points_max = []
    print(f"      Processing baselines for uv plot.")
    print(f"      Total number of baselines = {len(list(bldict.keys()))}")
    kk = 0 
    wavelen = (lightspeed/chunk_chan_mean_column)*1e3 # in Klambda
    for bl in tqdm(list(bldict.keys())[::1]):
        # print(len(list(bldict.keys())))
        
        u = uvw[0,bldict[bl]]
        v = uvw[1,bldict[bl]]
        uuu = np.tile(np.hstack([u, -u]), (reshaped_size, 1))/wavelen
        # print(uuu.shape[0] *uuu.shape[1])
        vvv = np.tile(np.hstack([v, -v]), (reshaped_size, 1))/wavelen
        # Reduce the datasets
        nspw = uuu.shape[0]
        # colors = plt.cm.viridis(np.linspace(0, 1, nspw))
        
        UUU_downsampled = downsample(uuu, downsample_factor)
        VVV_downsampled = downsample(vvv, downsample_factor)
        # print(UUU_downsampled)
        if color is None:
            color = colors[kk % len(colors)]
        # ax.plot(UUU_downsampled[i], VVV_downsampled[i], 'o', markersize=1.0, color=color, alpha=0.7)
        # if kk < 1:
        #     n_points = 2*len(chan_freq)*uvw.shape[-1]*len(list(bldict.keys()))
        #     print(f"      Total number of uv points = {n_points/1e9} billion")
        #     print(f"      Total number of downsampled uv points to plot = {n_points/(chunk_size*downsample_factor*1e9)} billion.")
            
        ax.plot(UUU_downsampled.T, VVV_downsampled.T, '.', markersize=0.5, color=color, alpha=alpha)
        
        max_u_point = np.nanmax(UUU_downsampled)
        max_v_point = np.nanmax(VVV_downsampled)
        u_points_max.append(max_u_point)
        v_points_max.append(max_v_point)
        kk = kk + 1
        # for i in range(nspw):
        #     color = colors[i % len(colors)]
        #     ax.plot(UUU_downsampled[i], VVV_downsampled[i], 'o', markersize=1.0, color=color, alpha=0.7)
#         UUU_reduced = chunk_average(uuu, chunk_size_plot)
#         VVV_reduced = chunk_average(vvv, chunk_size_plot)
#         ax.plot(UUU_reduced.T, VVV_reduced.T, '.', markersize=0.5, color=color, alpha=0.7)
#         # ax.clear()  # Clear previous plot
#         ax.plot(uuu, vvv, '.',markersize=0.1,color=color,alpha=alpha)#, label=labels[bl])
#     #     fig.canvas.draw()
#     #     fig.canvas.flush_events()
#     # plt.ioff()
#     # plt.show()
    
    # xy_limits = [np.nanmin([np.nanmin(VVV_downsampled), np.nanmin(UUU_downsampled)]),
    #              np.nanmax([np.nanmax(VVV_downsampled), np.nanmax(UUU_downsampled)])]
    xy_limits = [-np.nanmax([np.nanmax(v_points_max), np.nanmax(u_points_max)]),
                 np.nanmax([np.nanmax(v_points_max), np.nanmax(u_points_max)])]
    
    if plot_axes == True:
        ax.set_xlabel('$u$ [k$\lambda$]')
        ax.set_ylabel('$v$ [k$\lambda$]')
        # ax.grid()
        if title_text is not None:
            ax.set_title(title_text)
    uv_ranges.append(xy_limits)
    # ax.set_ylim(xy_limits[0]*1.1,xy_limits[1]*1.1)
    # ax.set_xlim(xy_limits[0]*1.1,xy_limits[1]*1.1)
    ax.set_ylim(np.nanmin(uv_ranges)*1.1,np.nanmax(uv_ranges)*1.1)
    ax.set_xlim(np.nanmin(uv_ranges)*1.1,np.nanmax(uv_ranges)*1.1)
    ax.axis('equal')
    if save_figure:
        plt.savefig(figure_save_name,bbox_inches='tight',dpi=300)
        plt.clf()
        plt.close()
    if show_figure:
        plt.show()
        plt.clf()
        plt.close()
    # return(fig,ax,uuu,UUU_downsampled,VVV_downsampled,uvw)
    return(fig,ax,uv_ranges)


def plot_uwave_vwave_v3(vis,color=None,fig=None,ax=None,alpha=1.0,
                     chunk_size = int(8), #channel average
                     downsample_factor = 200,
                     chunk_size_plot = 100,
                     figsize=(8,8),
                     uv_ranges=[],
                     plot_axes = True,
                     save_figure=False,
                     show_figure=False,
                     title_text = None,
                     figure_save_name='uvwave_coverage.jpg'
                    ):
    lightspeed = 299792458.0 # speed of light in m/s
    def get_blinfo():
        return dict([((x, y), np.where((ant1 == x) & (ant2 == y))[0]) 
                    for x in ant_uniq for y in ant_uniq if y > x])
    

    msmd = casatools.msmetadata()
    ms = casatools.ms()
    tb = casatools.table()
    print(f"++==> Reading info of visibility: {vis}")
    
    msmd.open(vis)
    bandwidth = msmd.bandwidths()
    nspw = len(bandwidth)
    chan_freqs_all = np.empty(nspw, dtype=object)
    spws_freq = np.zeros(nspw)

    for nch in range(nspw):
        chan_freqs_all[nch] = msmd.chanfreqs(nch)
        spws_freq[nch] = np.mean(chan_freqs_all[nch])
    msmd.done()
    chan_freq = np.concatenate(chan_freqs_all)
    
    #vis and antenna info
    tb.open(vis)
    uvw = tb.getcol('UVW')
    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    flag = tb.getcol('FLAG')
    flag_row = tb.getcol('FLAG_ROW')
    tb.close()
    tb.open(vis+'::ANTENNA')
    station = tb.getcol('STATION')
    tb.close()

    ant_uniq = np.unique(np.hstack((ant1, ant2)))
    bldict = get_blinfo()
    wavelen = (lightspeed/chan_freq.mean())*1e6 # in Klambda
    # flag_final = np.logical_or(flag, flag_row[:,np.newaxis,np.newaxis])
    flag_final = np.logical_or(flag, flag_row[np.newaxis,np.newaxis,:])
    # flag_final_VLA = np.logical_or(flag_VLA, flag_row_VLA[:,np.newaxis,np.newaxis])


    labels = {}
    for bl in bldict.keys():
        labels[bl] = station[bl[0]]+'-'+station[bl[1]]

    # plt.ion() 
    if fig is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    
    
    def chunk_average(arr, chunk_size):
        # Reshape the array to have chunks of size `chunk_size`
        reshaped = arr[:, :arr.shape[1] - arr.shape[1] % chunk_size].reshape(arr.shape[0], -1, chunk_size)
        # Average over each chunk along the last axis
        return reshaped.mean(axis=2)


    
    def downsample(arr, factor):
        return arr[:, ::factor]
    
    
    
    
    reshaped_size = len(chan_freq) // chunk_size_plot
    avg_chan = chan_freq[:reshaped_size * chunk_size_plot].reshape(reshaped_size, chunk_size_plot)
    chunk_chan_mean = avg_chan.mean(axis=1)
    chunk_chan_mean_broad = np.repeat(chunk_chan_mean, chunk_size_plot)
    chunk_chan_mean_column = chunk_chan_mean.reshape(-1, 1)
    UUU = []
    VVV = []
    kk = 0

    colors = [
        '#E69F00', '#56B4E9', '#009E73', '#F0E442',
        '#0072B2', '#D55E00', '#CC79A7', '#000000']
    u_points_max = []
    v_points_max = []
    print(f"      Processing baselines for uv plot.")
    print(f"      Total number of baselines = {len(list(bldict.keys()))}")
    kk = 0 
    for bl in tqdm(list(bldict.keys())[::1]):
        # print(len(list(bldict.keys())))
        wavelen = (lightspeed / chunk_chan_mean_column) * 1e3  # in Klambda
        u = uvw[0, bldict[bl]]
        v = uvw[1, bldict[bl]]
        uuu = np.tile(np.hstack([u, np.nan, -u]), (reshaped_size, 1)) / wavelen
        # print(uuu.shape[0] *uuu.shape[1])
        vvv = np.tile(np.hstack([v, np.nan, -v]), (reshaped_size, 1)) / wavelen
        # Reduce the datasets
        nspw = uuu.shape[0]
        
        UUU_downsampled = chunk_average(uuu, downsample_factor)
        VVV_downsampled = chunk_average(vvv, downsample_factor)
        # print(UUU_downsampled)
        
        if color is None:
            color = colors[kk % len(colors)]
            
        if UUU_downsampled.size > 0 and UUU_downsampled.size > 0:
            # Plot the averaged points if they are not empty
            ax.plot(UUU_downsampled.T, VVV_downsampled.T, '.', markersize=0.5, color=color, alpha=alpha)
            max_u_point = np.nanmax(UUU_downsampled)
            max_v_point = np.nanmax(VVV_downsampled)
            u_points_max.append(max_u_point)
            v_points_max.append(max_v_point)
            
        kk = kk + 1
        
    xy_limits = [-np.nanmax([np.nanmax(v_points_max), np.nanmax(u_points_max)]),
                 np.nanmax([np.nanmax(v_points_max), np.nanmax(u_points_max)])]
    
    if plot_axes == True:
        ax.set_xlabel('$u$ [k$\lambda$]')
        ax.set_ylabel('$v$ [k$\lambda$]')
        # ax.grid()
        if title_text is not None:
            ax.set_title(title_text)
    uv_ranges.append(xy_limits)

    ax.set_ylim(np.nanmin(uv_ranges)*1.1,np.nanmax(uv_ranges)*1.1)
    ax.set_xlim(np.nanmin(uv_ranges)*1.1,np.nanmax(uv_ranges)*1.1)
    ax.axis('equal')
    if save_figure:
        plt.savefig(figure_save_name,bbox_inches='tight',dpi=300)
        plt.clf()
        plt.close()
    if show_figure:
        plt.show()
        plt.clf()
        plt.close()
    return(fig,ax,uv_ranges)



def plot_uwave_vwave_v4(vis, color=None, fig=None, ax=None, alpha=1.0,
                        chunk_size=int(8),  # channel average
                        downsample_factor=200,
                        chunk_size_plot=100,
                        figsize=(8, 8),
                        uv_ranges=[],
                        plot_axes=True,
                        save_figure=False,
                        show_figure=False,
                        title_text=None,
                        figure_save_name='uvwave_coverage.jpg'):
    
    # Constants
    lightspeed = 299792458.0  # speed of light in m/s

    def get_blinfo():
        # Dictionary mapping antenna pairs to the indices of their baselines
        return dict([((x, y), np.where((ant1 == x) & (ant2 == y))[0])
                    for x in ant_uniq for y in ant_uniq if y > x])
    
    # Reading metadata and visibility information
    msmd = casatools.msmetadata()
    tb = casatools.table()
    print(f"++==> Reading info of visibility: {vis}")
    
    # Get spectral information
    msmd.open(vis)
    nspw = msmd.nspw()
    chan_freqs_all = np.array([msmd.chanfreqs(spw) for spw in range(nspw)], dtype=object)
    msmd.done()
    chan_freq = np.concatenate(chan_freqs_all)
    
    # Get visibility and antenna data
    tb.open(vis)
    uvw = tb.getcol('UVW')
    ant1 = tb.getcol('ANTENNA1')
    ant2 = tb.getcol('ANTENNA2')
    flag = tb.getcol('FLAG')
    flag_row = tb.getcol('FLAG_ROW')
    tb.close()
    
    tb.open(vis + '::ANTENNA')
    station = tb.getcol('STATION')
    tb.close()

    # Get unique antennas and baselines
    ant_uniq = np.unique(np.hstack((ant1, ant2)))
    bldict = get_blinfo()

    # Handle flags
    flag_final = np.logical_or(flag, flag_row[np.newaxis,np.newaxis,:])

    # Wavelength calculation for each channel
    wavelengths = lightspeed / chan_freq

    # Average the channels for better performance
    reshaped_size = len(chan_freq) // chunk_size_plot
    avg_chan = chan_freq[:reshaped_size * chunk_size_plot].reshape(reshaped_size, chunk_size_plot)
    chunk_chan_mean = avg_chan.mean(axis=1)
    wavelengths_mean = lightspeed / chunk_chan_mean

    # Initialize plot if not provided
    if fig is None:
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111)
    
    colors = ['#E69F00', '#56B4E9', '#009E73', '#F0E442', '#0072B2', '#D55E00', '#CC79A7', '#000000']
    
    u_points_max = []
    v_points_max = []
    
    print(f"      Processing baselines for uv plot.")
    print(f"      Total number of baselines = {len(list(bldict.keys()))}")
    
    # Loop through the baselines
    for kk, bl in enumerate(tqdm(list(bldict.keys()))):
        u = uvw[0, bldict[bl]]
        v = uvw[1, bldict[bl]]
        
        # Apply flagging
        # flagged_u = np.ma.masked_array(u, flag_final[bldict[bl]])
        # flagged_v = np.ma.masked_array(v, flag_final[bldict[bl]])

        # # Skip this baseline if there's no unflagged data
        # if flagged_u.size == 0 or flagged_v.size == 0:
        #     continue
        
        # # Average the u and v points over chunks
        # valid_size = flagged_u.size // chunk_size_plot * chunk_size_plot
        # if valid_size == 0:
        #     continue
        
        # Reshape and average
        u_chunked = u.reshape(-1, chunk_size_plot).mean(axis=1)
        v_chunked = v.reshape(-1, chunk_size_plot).mean(axis=1)

        # Convert u and v to uwave and vwave using averaged wavelengths
        uwave = u_chunked / wavelengths_mean
        vwave = v_chunked / wavelengths_mean

        # Downsample the data for efficient plotting
        uwave_downsampled = uwave[::downsample_factor]
        vwave_downsampled = vwave[::downsample_factor]

        # Choose color for the baseline
        plot_color = color if color is not None else colors[kk % len(colors)]
        
        # Plot the downsampled data if non-empty
        if uwave_downsampled.size > 0 and vwave_downsampled.size > 0:
            ax.plot(uwave_downsampled, vwave_downsampled, '.', markersize=0.5, color=plot_color, alpha=alpha)
            u_points_max.append(np.nanmax(np.abs(uwave_downsampled)))
            v_points_max.append(np.nanmax(np.abs(vwave_downsampled)))

    # Set plot limits
    if u_points_max and v_points_max:
        max_limit = np.nanmax([np.nanmax(u_points_max), np.nanmax(v_points_max)])
        xy_limits = [-max_limit, max_limit]
        ax.set_xlim(xy_limits)
        ax.set_ylim(xy_limits)

    # Set plot labels and title
    if plot_axes:
        ax.set_xlabel('$u$ [kλ]')
        ax.set_ylabel('$v$ [kλ]')
        if title_text is not None:
            ax.set_title(title_text)
    ax.axis('equal')
    if save_figure:
        plt.savefig(figure_save_name, bbox_inches='tight', dpi=300)
        plt.clf()
        plt.close()
    
    if show_figure:
        plt.show()
        plt.clf()
        plt.close()
    
    return fig, ax, uv_ranges




# def plot_uwave_vwave_v2(vis,color='black',fig=None,ax=None,alpha=1.0):
#     lightspeed = 299792458.0 # speed of light in m/s
#     def get_blinfo():
#         return dict([((x, y), np.where((ant1 == x) & (ant2 == y))[0]) 
#                     for x in ant_uniq for y in ant_uniq if y > x])
#     spwtab = tb.open(vis+'::SPECTRAL_WINDOW')
#     chan_freq = tb.getcol('CHAN_FREQ').flatten()
#     tb.close()

#     #vis and antenna info
#     tb.open(vis)
#     uvw = tb.getcol('UVW')
#     ant1 = tb.getcol('ANTENNA1')
#     ant2 = tb.getcol('ANTENNA2')
#     flag = tb.getcol('FLAG')
#     flag_row = tb.getcol('FLAG_ROW')
#     tb.close()
#     tb.open(vis+'::ANTENNA')
#     station = tb.getcol('STATION')
#     tb.close()

#     ant_uniq = np.unique(np.hstack((ant1, ant2)))
#     bldict = get_blinfo()
#     wavelen = (lightspeed/chan_freq.mean())*1e6 # in Klambda
#     # flag_final = np.logical_or(flag, flag_row[:,np.newaxis,np.newaxis])
#     flag_final = np.logical_or(flag, flag_row[np.newaxis,np.newaxis,:])
#     # flag_final_VLA = np.logical_or(flag_VLA, flag_row_VLA[:,np.newaxis,np.newaxis])


#     labels = {}
#     for bl in bldict.keys():
#         labels[bl] = station[bl[0]]+'-'+station[bl[1]]

#     # plt.ion() 
#     if fig is None:
#         fig = plt.figure()
#         ax = fig.add_subplot(111)
    
#     def chunk_average(arr, chunk_size):
#         # Reshape the array to have chunks of size `chunk_size`
#         reshaped = arr.reshape(arr.shape[0], -1, chunk_size)
#         # Average over each chunk along the last axis
#         return reshaped.mean(axis=2)
    
#     # Function to downsample an array along its columns (data points)
#     def downsample(arr, factor):
#         return arr[:, ::factor]

#     # Function to average an array's rows (channels)
#     def average_channels(arr, factor):
#         # Reshape the array to have rows of size `factor`
#         reshaped = arr.reshape(-1, factor, arr.shape[1])
#         # Average over each set of rows
#         return reshaped.mean(axis=1)
    
#     # Downsampling and averaging factors
#     downsample_factor = 1000  # For data points
#     downsample_channel_factor = 64  # For channels
    
#     chunk_size = int(12) #channel average
#     reshaped_size = int(len(chan_freq)/chunk_size)
#     avg_chan = chan_freq.reshape((reshaped_size, chunk_size))
#     chunk_chan_mean = avg_chan.mean(axis=1)
#     chunk_chan_mean_broad = np.repeat(chunk_chan_mean, chunk_size)
#     chunk_chan_mean_column = chunk_chan_mean.reshape(-1, 1)
#     UUU = []
#     VVV = []
#     for bl in tqdm(list(bldict.keys())[::1]):
#         wavelen = (lightspeed/chunk_chan_mean_column)*1e3 # in Klambda
#         u = uvw[0,bldict[bl]]
#         v = uvw[1,bldict[bl]]
#         uuu = np.tile(np.hstack([u, np.nan, -u]), (reshaped_size, 1))/wavelen
#         # print(uuu.shape[0] *uuu.shape[1])
#         vvv = np.tile(np.hstack([v, np.nan, -v]), (reshaped_size, 1))/wavelen
#         # Reduce the datasets
        
#         # Step 1: Downsample data points
#         UUU_downsampled = downsample(uuu, downsample_factor)
#         VVV_downsampled = downsample(vvv, downsample_factor)

#         # Step 2: Average channels
#         UUU_reduced = average_channels(UUU_downsampled, downsample_channel_factor)
#         VVV_reduced = average_channels(VVV_downsampled, downsample_channel_factor)

#         ax.plot(UUU_reduced.T, VVV_reduced.T, '.', markersize=0.5, color=color, alpha=0.7)
# #         UUU_reduced = chunk_average(uuu, chunk_size_plot)
# #         VVV_reduced = chunk_average(vvv, chunk_size_plot)
# #         ax.plot(UUU_reduced.T, VVV_reduced.T, '.', markersize=0.5, color=color, alpha=0.7)
# #         # ax.clear()  # Clear previous plot
# #         ax.plot(uuu, vvv, '.',markersize=0.1,color=color,alpha=alpha)#, label=labels[bl])
# #     #     fig.canvas.draw()
# #     #     fig.canvas.flush_events()
# #     # plt.ioff()
# #     # plt.show()
#     return(fig,ax)
