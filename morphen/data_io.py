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

def load_fits_data(image):
    '''
        Name origin:
        ctn > casa to numpy
        Function that read fits files, using casa IA.open or astropy.io.fits.
        Note: For some reason, IA.open returns a rotated mirroed array, so we need
        to undo it by a rotation.
        '''
    if isinstance(image, str) == True:
        try:
            ia = IA()
            ia.open(image)
            try:
                numpy_array = ia.getchunk()[:, :, 0, 0]
            except:
                numpy_array = ia.getchunk()[:, :]
            ia.close()
            # casa gives a mirroed and 90-degree rotated image :(
            data_image = np.rot90(numpy_array)[::-1, ::]
            return (data_image)
        except:
            try:
                _data_image = pf.getdata(image)
                if len(_data_image.shape) == 4:
                    data_image = _data_image[0][0]
                else:
                    data_image = _data_image
                return (data_image)
            except:
                print('Input image is not a fits file.')
                return(ValueError)
    else:
        raise ValueError("Input image is not a string.")

ctn = deprecated("ctn","load_fits_data")(load_fits_data)



def cut_image(img, center=None, size=(1024, 1024),
              cutout_filename=None, special_name=''):
    """
    Cut images keeping updated header/wcs.
    This function is a helper to cut both image and its associated residual
    (from casa or wsclean).
    It saves both images with a cutout prefix.
    If the centre is not given, the peak position will be selected.
    If the size is not defined, cut a standard size of (1024 x 1024).
    It updates the wcs of the croped image.

    To do: crop a fration of the image.
    """
    if center is None:
        """
        Better to include momments instead of peak.
        """
        imst = imstat(img)
        position = (imst['maxpos'][0], imst['maxpos'][1])
    else:
        position = center

    hdu = pf.open(img)[0]
    wcs = WCS(hdu.header, naxis=2)

    cutout = Cutout2D(hdu.data[0][0], position=position, size=size, wcs=wcs)
    hdu.data = cutout.data
    hdu.header.update(cutout.wcs.to_header())
    if cutout_filename is None:
        cutout_filename = img.replace('-image.fits', '-image_cutout' +
                                      special_name + '.fits')
    hdu.writeto(cutout_filename, overwrite=True)

    # do the same for the residual image
    hdu2 = pf.open(img)[0]
    wcs2 = WCS(hdu2.header, naxis=2)

    hdu_res = pf.open(img.replace('-image.fits', '-residual.fits'))[0]
    # plt.imshow(hdu_res.data[0][0])
    wcs_res = WCS(hdu_res.header, naxis=2)
    cutout_res = Cutout2D(hdu_res.data[0][0],
                          position=position, size=size, wcs=wcs2)
    hdu2.data = cutout_res.data
    hdu2.header.update(cutout_res.wcs.to_header())
    # if cutout_filename is None:
    cutout_filename_res = img.replace('-image.fits', '-residual_cutout' +
                                      special_name + '.fits')
    hdu2.writeto(cutout_filename_res, overwrite=True)
    return(cutout_filename,cutout_filename_res)


def do_cutout(image, box_size=(200,200), center=None, return_='data'):
    """
    Fast cutout of a numpy array.

    Returs: numpy data array or a box for that cutout, if asked.
    """
    if isinstance(box_size, int):
        box_size = (box_size,box_size)

    if center is None:
        if isinstance(image, str) == True:
            # imhd = imhead(image)
            st = imstat(image)
            print('  >> Center --> ', st['maxpos'])
            xin, xen, yin, yen = (st['maxpos'][0] - box_size[0],
                                  st['maxpos'][0] + box_size[0],
                                  st['maxpos'][1] - box_size[1],
                                  st['maxpos'][1] + box_size[1])
            data_cutout = load_fits_data(image)[xin:xen, yin:yen]

        else:
            try:
                max_x, max_y = np.where(load_fits_data(image) == load_fits_data(image).max())
                xin = max_x[0] - box_size[0]
                xen = max_x[0] + box_size[0]
                yin = max_y[0] - box_size[1]
                yen = max_y[0] + box_size[1]
            except:
                max_x, max_y = np.where(image == image.max())
                xin = max_x[0] - box_size[0]
                xen = max_x[0] + box_size[0]
                yin = max_y[0] - box_size[1]
                yen = max_y[0] + box_size[1]

            data_cutout = image[xin:xen, yin:yen]


    else:
        xin, xen, yin, yen = (center[0] - box_size[0], center[0] + box_size[0],
                              center[1] - box_size[1], center[1] + box_size[1])
        if isinstance(image, str) == True:
            data_cutout = load_fits_data(image)[xin:xen, yin:yen]
        else:
            data_cutout = image[xin:xen, yin:yen]

    if return_ == 'data':
        return (data_cutout)
    if return_ == 'box':
        box = xin, xen, yin, yen  # [xin:xen,yin:yen]
        return (box)


# def do_cutout_2D(image_data, box_size=300, center=None, return_='data'):
#     """
#     Fast cutout of a numpy array.
#
#     Returs: numpy data array or a box for that cutout, if asked.
#     """
#
#     if center is None:
#         x0, y0= nd.maximum_position(image_data)
#         print('  >> Center --> ', x0, y0)
#         if x0-box_size>1:
#             xin, xen, yin, yen = x0 - box_size, x0 + box_size, \
#                                  y0 - box_size, y0 + box_size
#         else:
#             print('Box size is larger than image!')
#             return ValueError
#     else:
#         xin, xen, yin, yen = center[0] - box_size, center[0] + box_size, \
#             center[1] - box_size, center[1] + box_size
#     if return_ == 'data':
#         data_cutout = image_data[xin:xen, yin:yen]
#         return (data_cutout)
#     if return_ == 'box':
#         box = xin, xen, yin, yen  # [xin:xen,yin:yen]
#         return(box)

def do_cutout_2D(image_data, box_size=(300,300), 
                 center=None, centre_mode='max',
                 return_='data'):
    """
    Fast cutout of a numpy array.

    Returs: numpy data array or a box for that cutout, if asked.
    """
    if isinstance(box_size, int):
        box_size = (box_size,box_size)

    if center is None:
        if centre_mode == 'max':
            x0, y0= nd.maximum_position(image_data)
        if centre_mode == 'image_centre':
            x0, y0 = image_data.shape[0]//2, image_data.shape[1]//2
            
        # print('  >> Center --> ', x0, y0)
        if x0-box_size[0]>1:
            xin, xen, yin, yen = x0 - box_size[0], x0 + box_size[0], \
                                 y0 - box_size[1], y0 + box_size[1]
        else:
            print('Box size is larger than image!')
            return ValueError
    else:
        xin, xen, yin, yen = center[0] - box_size[0], center[0] + box_size[0], \
            center[1] - box_size[1], center[1] + box_size[1]
    if return_ == 'data':
        data_cutout = image_data[xin:xen, yin:yen]
        return (data_cutout)
    if return_ == 'box':
        box = xin, xen, yin, yen  # [xin:xen,yin:yen]
        return(box)




def copy_header(image, image_to_copy, file_to_save=None):
    """
    For image files with no wcs, copy the header from a similar/equal image to
    the wanted file.
    Note: This is intended to be used to copy headers from images to their
    associated models and residuals.
    Note: Residual CASA images do not have Jy/Beam units, so this function
        can be used to copy the header/wcs information to the wanted file
        in order to compute the total flux in residual maps after the
        header has been copied.
    """
    from astropy.io import fits
    if file_to_save is None:
        file_to_save = image_to_copy.replace('.fits', 'header.fits')
    # Open the source image and get its header
    with fits.open(image) as hdu1:
        header = hdu1[0].header
        # Open the target image and update its header
        with fits.open(image_to_copy, mode='update') as hdu2:
            hdu2[0].header.update(header)
            hdu2.flush()
            hdu2.close()
    pass

def get_frequency(imagename):
    """
    Get the frequency of a radio observation from the wcs of a fits image.
    """
    from astropy.io import fits
    from astropy.wcs import WCS

    # Open the FITS file
    with fits.open(imagename) as hdulist:
        header = hdulist[0].header

    # Extract WCS information
    wcs_info = WCS(header)

    for i in range(1, wcs_info.naxis + 1):
        if 'FREQ' in header.get(f'CTYPE{i}', ''):
            freq_ref = header.get(f'CRVAL{i}')
            frequency = freq_ref/1e9
    return(frequency)


def pad_psf(imagename,psfnasme):
    import numpy as np
    psf_data = pf.getdata(psfnasme)
    image_data = load_fits_data(imagename)
    # Assume that 'psf' and 'image' are the original psf and image arrays, respectively
    psf_image_size = psf_data.shape[0] # assuming the psf is square
    image_size = image_data.shape[0] # assuming the image is square
    padding_size = (image_size - psf_image_size + 1) // 2

    # Create a new array of zeros with the desired padded size
    padded_psf = np.zeros((image_size, image_size))
    start_idx = image_size // 2 - psf_image_size // 2
    end_idx = start_idx + psf_image_size

    # Copy the original psf into the center of the new array
    padded_psf[start_idx:end_idx, start_idx:end_idx] = psf_data

    # Copy the original psf into the center of the new array
#     padded_psf[padding_size:-padding_size, padding_size:-padding_size] = psf_data
    pf.writeto(imagename.replace('.fits','_psf.fits'),padded_psf,overwrite=True)
    return(imagename.replace('.fits','_psf.fits'))



def format_coords(dec_raw):
    deg, rest = dec_raw.split('.', 1)
    min, rest = rest.split('.', 1)
    sec = rest
    dec_formatted = f"{deg} {min} {sec}"
    return (dec_formatted)


def conver_str_coords(ra, dec):
    from astropy.coordinates import Angle, SkyCoord
    import astropy.units as u

    #     # Example hour-angle coordinates
    # #     ha_str = '13h15m34.9461s'
    # #     dec_str = '+62d07m28.6912s'

    #     # Convert the hour-angle and declination to angles
    #     ha = Angle(ha_str,unit='hourangle')
    #     print(ha)
    #     dec = Angle(dec_str,unit='deg')

    #     # Create a SkyCoord object with the coordinates and convert to ICRS frame
    #     coords = SkyCoord(ha, dec, unit=(u.hourangle, u.deg), frame='icrs')

    #     # Get the RA and Dec in degrees
    #     ra_deg = coords.ra.deg
    #     dec_deg = coords.dec.deg
    coor = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    ra_deg = coor.ra.degree
    dec_deg = coor.dec.degree
    # print(ra_deg, dec_deg)
    return (ra_deg, dec_deg)



def calculate_centroid_and_imagesize(coordinates, cellsize):
    """
    Calculate the centroid and imagesize for interferometric imaging given a list of celestial coordinates.
    
    
    # Example: Fake coordinates
    coordinates = [
        '20:37:31.075 +25:33:04.896',
        '20:37:28.675 +25:34:14.596',
        '20:37:34.475 +25:32:54.096',
        '20:37:30.875 +25:31:44.896',
    ]
    cellsize = 0.2  # Arcseconds

    # Calculate centroid and imagesize
    result = calculate_centroid_and_imagesize(coordinates, cellsize)

    # Extract phase centre and image dimensions
    phase_centre = SkyCoord(result["phase_centre"], unit=(u.hourangle, u.deg))
    ra_imagesize, dec_imagesize = result["imagesize"]

    # Visualize
    sky_coords = SkyCoord(coordinates, unit=(u.hourangle, u.deg))
    plt.figure(figsize=(8, 6))
    plt.scatter(sky_coords.ra.deg, sky_coords.dec.deg, color="blue", label="Source Positions")
    plt.scatter(phase_centre.ra.deg, phase_centre.dec.deg, color="red", marker="x", s=100, label="Centroid (Phase Centre)")

    # Draw rectangle
    margin_factor = 1.1
    ra_min = min(sky_coords.ra.deg) - (margin_factor - 1) * (max(sky_coords.ra.deg) - min(sky_coords.ra.deg))
    ra_max = max(sky_coords.ra.deg) + (margin_factor - 1) * (max(sky_coords.ra.deg) - min(sky_coords.ra.deg))
    dec_min = min(sky_coords.dec.deg) - (margin_factor - 1) * (max(sky_coords.dec.deg) - min(sky_coords.dec.deg))
    dec_max = max(sky_coords.dec.deg) + (margin_factor - 1) * (max(sky_coords.dec.deg) - min(sky_coords.dec.deg))

    plt.plot(
        [ra_min, ra_min, ra_max, ra_max, ra_min],
        [dec_min, dec_max, dec_max, dec_min, dec_min],
        color="orange", label="Imaging Region", linestyle="--"
    )

    plt.xlabel("RA (deg)")
    plt.ylabel("Dec (deg)")
    plt.legend()
    plt.title("Calculated Imaging Region and Centroid")
    plt.grid()
    plt.show()

    
    Parameters:
    coordinates (list): List of celestial coordinates as strings (e.g., '20:37:31.075 +25.33.04.896').
    cellsize (float): Pixel size in arcseconds.
    
    
    Returns:
    dict: Phase centre and imagesize with margins.
    """
    # Parse coordinates
    sky_coords = SkyCoord(coordinates, unit=(u.hourangle, u.deg))
    ra_values = sky_coords.ra.deg
    dec_values = sky_coords.dec.deg

    # Compute centroid
    ra_centroid = np.mean(ra_values)
    dec_centroid = np.mean(dec_values)
    phase_centre = SkyCoord(ra_centroid, dec_centroid, unit=u.deg)

    # Compute the angular extent
    ra_extent = (max(ra_values) - min(ra_values)) * np.cos(np.radians(dec_centroid))
    dec_extent = max(dec_values) - min(dec_values)

    # Add margins
    margin_factor = 1.1
    ra_extent_margined = ra_extent * margin_factor
    dec_extent_margined = dec_extent * margin_factor

    # Convert angular size to pixel size
    ra_imagesize = int(np.ceil(ra_extent_margined * 3600 / cellsize))
    dec_imagesize = int(np.ceil(dec_extent_margined * 3600 / cellsize))

    # Return results
    return {
        "phase_centre": phase_centre.to_string('hmsdms'),
        "phase_centre_formatted": phase_centre.to_string(style='hmsdms', sep=':', precision=3),
        "imagesize": (ra_imagesize, dec_imagesize),
    }



def find_offsets(reference_coords, target_coords):
    # Compute the pixel offsets between target and reference coordinates
    offset_x = reference_coords[0] - target_coords[0]
    offset_y = reference_coords[1] - target_coords[1]
    return offset_x, offset_y


def cutout_2D_radec_backup(imagename, residualname=None, ra_f=None, dec_f=None, cutout_size=1024,
                    special_name='',correct_shift=False,ref_cutout_image=None):
    from astropy.io import fits
    import os
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord
    # load image data and header
    if ra_f is None:
        imst = imstat(imagename)
        print('maxpos = ', imst['maxpos'])
        print('maxposf = ', imst['maxposf'])
        coords = imst['maxposf'].split(',')
        ra = coords[0]
        dec = format_coords(coords[1])
        # print(ra, dec)
        ra_f, dec_f = conver_str_coords(ra, dec)
        print('ra_f,dec_f = ', ra_f,dec_f)

    with fits.open(imagename) as hdul:
        image_data, header = hdul[0].data, hdul[0].header

        # create a WCS object from the header
        wcs = WCS(header, naxis=2)
        # wcs.wcs.radesys = 'icrs'

        # set the center and size of the cutout
        #     ra_f,dec_f = conver_str_coords(ra,dec)
        center_ra = ra_f  # center RA in degrees
        center_dec = dec_f  # center Dec in degrees

        # center = SkyCoord(ra=center_ra, dec=center_dec, unit='deg',from)
        center = SkyCoord(ra=center_ra * u.degree, dec=center_dec * u.degree,
                          frame='icrs')
        # print(center)
        # create a Cutout2D object
        cutout = Cutout2D(image_data[0][0], center, cutout_size, wcs=wcs)
        # apply shift
        if correct_shift == True:
            if ref_cutout_image is not None:
                ref_image_cutout_data = load_fits_data(ref_cutout_image)
                x_ref, y_ref = nd.maximum_position(ref_image_cutout_data)[::-1]
                reference_source_coords = (x_ref, y_ref)
                offset_x, offset_y = find_offsets(reference_source_coords,
                                                  nd.maximum_position(cutout.data)[::-1])
                print(f" !!!! Offsets of peak position are: {offset_x, offset_y}.")
                aligned_target_image = shift(cutout.data, (offset_y, offset_x),
                                             mode='constant')
                new_hdul = fits.HDUList(
                    [fits.PrimaryHDU(header=hdul[0].header, data=aligned_target_image)])
            else:
                #the code must stop
                print('No reference image was provided. '
                      'No shift correction will be applied.')
        else:
            new_hdul = fits.HDUList(
                [fits.PrimaryHDU(header=hdul[0].header, data=cutout.data)])

        new_hdul[0].header.update(cutout.wcs.to_header())
        savename_img = os.path.dirname(imagename) + '/' + os.path.basename(imagename).replace(
            '.fits', '.cutout' + special_name + '.fits')
        new_hdul.writeto(savename_img, overwrite=True)

    if residualname is not None:
        with fits.open(residualname) as hdul:
            image_data, header = hdul[0].data, hdul[0].header

            # create a WCS object from the header
            wcs = WCS(header, naxis=2)
            # wcs.wcs.radesys = 'icrs'

            # set the center and size of the cutout
            #     ra_f,dec_f = conver_str_coords(ra,dec)
            center_ra = ra_f  # center RA in degrees
            center_dec = dec_f  # center Dec in degrees

            # center = SkyCoord(ra=center_ra, dec=center_dec, unit='deg',from)
            center = SkyCoord(ra=center_ra * u.degree, dec=center_dec * u.degree,
                              frame='icrs')

            # create a Cutout2D object
            cutout = Cutout2D(image_data[0][0], center, cutout_size, wcs=wcs)
            if correct_shift == True:
                if ref_cutout_image is not None:
                    # ref_image_cutout_data = load_fits_data(ref_cutout_image)
                    # x_ref, y_ref = nd.maximum_position(ref_image_cutout_data)[::-1]
                    # reference_source_coords = (x_ref, y_ref)
                    # offset_x, offset_y = find_offsets(reference_source_coords,
                    #                                   nd.maximum_position(cutout.data)[::-1])
                    # print(f" !!!! Offsets of peak position are: {offset_x, offset_y}.")
                    aligned_target_image = shift(cutout.data, (offset_y, offset_x),
                                                 mode='constant')
                    new_hdul = fits.HDUList(
                        [fits.PrimaryHDU(header=hdul[0].header, data=aligned_target_image)])
                else:
                    # the code must stop
                    print('No reference image was provided. '
                          'No shift correction will be applied.')
            else:
                new_hdul = fits.HDUList(
                    [fits.PrimaryHDU(header=hdul[0].header, data=cutout.data)])

            # new_hdul = fits.HDUList(
            #     [fits.PrimaryHDU(header=hdul[0].header, data=cutout.data)])
            new_hdul[0].header.update(cutout.wcs.to_header())
            savename_res = os.path.dirname(residualname) + '/' + os.path.basename(
                residualname).replace('.fits', '.cutout' + special_name + '.fits')
            new_hdul.writeto(savename_res, overwrite=True)
    return (ra_f, dec_f,savename_img,wcs)

def cutout_2D_radec(imagename, residualname=None, ra_f=None, dec_f=None, cutout_size=1024,
                    special_name='', correct_shift=False, ref_cutout_image=None, pixel_coords=None):
    from astropy.io import fits
    import os
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord
    
    # load image data and header
    with fits.open(imagename) as hdul:
        image_data, header = hdul[0].data, hdul[0].header
        wcs = WCS(header, naxis=2)

        if pixel_coords is not None:
            # Convert pixel coordinates to RA/Dec
            ra_f, dec_f = wcs.pixel_to_world(*pixel_coords).ra.degree, wcs.pixel_to_world(*pixel_coords).dec.degree
            print(f"Converted pixel coordinates {pixel_coords} to RA/Dec: ({ra_f}, {dec_f})")
        
        elif ra_f is None or dec_f is None:
            imst = imstat(imagename)
            print('maxpos = ', imst['maxpos'])
            print('maxposf = ', imst['maxposf'])
            coords = imst['maxposf'].split(',')
            ra = coords[0]
            dec = format_coords(coords[1])
            ra_f, dec_f = conver_str_coords(ra, dec)
            print('ra_f, dec_f = ', ra_f, dec_f)
        
        # set the center and size of the cutout
        center = SkyCoord(ra=ra_f * u.degree, dec=dec_f * u.degree, frame='icrs')
        # print('Centre SkyCoord = ',center)
        # create a Cutout2D object
        cutout = Cutout2D(image_data[0][0], center, cutout_size, wcs=wcs)
        
        # apply shift
        if correct_shift:
            if ref_cutout_image is not None:
                ref_image_cutout_data = load_fits_data(ref_cutout_image)
                x_ref, y_ref = nd.maximum_position(ref_image_cutout_data)[::-1]
                reference_source_coords = (x_ref, y_ref)
                offset_x, offset_y = find_offsets(reference_source_coords,
                                                  nd.maximum_position(cutout.data)[::-1])
                print(f" !!!! Offsets of peak position are: {offset_x, offset_y}.")
                aligned_target_image = shift(cutout.data, (offset_y, offset_x),
                                             mode='constant')
                new_hdul = fits.HDUList(
                    [fits.PrimaryHDU(header=hdul[0].header, data=aligned_target_image)])
            else:
                # the code must stop
                print('No reference image was provided. '
                      'No shift correction will be applied.')
        else:
            new_hdul = fits.HDUList(
                [fits.PrimaryHDU(header=hdul[0].header, data=cutout.data)])

        new_hdul[0].header.update(cutout.wcs.to_header())
        savename_img = os.path.dirname(imagename) + '/' + os.path.basename(imagename).replace(
            '.fits', '.cutout.' + special_name + '.fits')
        new_hdul.writeto(savename_img, overwrite=True)

    if residualname is not None:
        with fits.open(residualname) as hdul:
            image_data, header = hdul[0].data, hdul[0].header
            wcs = WCS(header, naxis=2)

            center = SkyCoord(ra=ra_f * u.degree, dec=dec_f * u.degree, frame='icrs')

            cutout = Cutout2D(image_data[0][0], center, cutout_size, wcs=wcs)
            if correct_shift:
                if ref_cutout_image is not None:
                    aligned_target_image = shift(cutout.data, (offset_y, offset_x),
                                                 mode='constant')
                    new_hdul = fits.HDUList(
                        [fits.PrimaryHDU(header=hdul[0].header, data=aligned_target_image)])
                else:
                    print('No reference image was provided. '
                          'No shift correction will be applied.')
            else:
                new_hdul = fits.HDUList(
                    [fits.PrimaryHDU(header=hdul[0].header, data=cutout.data)])

            new_hdul[0].header.update(cutout.wcs.to_header())
            savename_res = os.path.dirname(residualname) + '/' + os.path.basename(
                residualname).replace('.fits', '.cutout.' + special_name + '.fits')
            new_hdul.writeto(savename_res, overwrite=True)
    
    return ra_f, dec_f, savename_img



def cutout_2D_radec_v2(imagename, residualname=None, ra_f=None, dec_f=None, cutout_size=1024,
                       cut_model=False,special_name='', correct_shift=False, 
                       ref_cutout_image=None, pixel_coords=None,
                       custom_save_path=None, custom_save_name=None):
    from astropy.io import fits
    import os
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord
    
    if os.path.exists(custom_save_path) is False:
        os.makedirs(custom_save_path)

    prefix_image_add = f"-{imagename.split('-image')[0].split('-')[-1]}"
    
    with fits.open(imagename) as hdul:
        image_data, header = hdul[0].data, hdul[0].header
        wcs = WCS(header, naxis=2)

        if pixel_coords is not None:
            ra_f, dec_f = wcs.pixel_to_world(*pixel_coords).ra.degree, wcs.pixel_to_world(*pixel_coords).dec.degree
            print(f"Converted pixel coordinates {pixel_coords} to RA/Dec: ({ra_f}, {dec_f})")
        elif ra_f is None or dec_f is None:
            imst = imstat(imagename)
            print('maxpos = ', imst['maxpos'])
            print('maxposf = ', imst['maxposf'])
            coords = imst['maxposf'].split(',')
            ra = coords[0]
            dec = format_coords(coords[1])
            ra_f, dec_f = conver_str_coords(ra, dec)
            print('ra_f, dec_f = ', ra_f, dec_f)
            
        center = SkyCoord(ra=ra_f * u.degree, dec=dec_f * u.degree, frame='icrs')
        # print('Centre SkyCoord = ', center)
        cutout = Cutout2D(image_data[0][0] if len(image_data.shape) == 4 else image_data, center, cutout_size, wcs=wcs)
        
        if len(image_data.shape) == 4:
            new_data = np.zeros((image_data.shape[0], image_data.shape[1], cutout_size, cutout_size))
            new_data[0,0] = cutout.data
        else:
            new_data = cutout.data

        if correct_shift:
            if ref_cutout_image is not None:
                ref_image_cutout_data = load_fits_data(ref_cutout_image)
                x_ref, y_ref = nd.maximum_position(ref_image_cutout_data)[::-1]
                reference_source_coords = (x_ref, y_ref)
                offset_x, offset_y = find_offsets(reference_source_coords,
                                                nd.maximum_position(cutout.data)[::-1])
                print(f" !!!! Offsets of peak position are: {offset_x, offset_y}.")
                if len(image_data.shape) == 4:
                    new_data[0,0] = shift(new_data[0,0], (offset_y, offset_x), mode='constant')
                else:
                    new_data = shift(new_data, (offset_y, offset_x), mode='constant')
            else:
                raise ValueError(f"No reference image was provided."
                                 f"No shift correction will be applied.")

        header.update(cutout.wcs.to_header())
        
        if len(image_data.shape) == 4:
            header['NAXIS3'] = image_data.shape[1]
            header['NAXIS4'] = image_data.shape[0]
            header['NAXIS1'] = cutout_size
            header['NAXIS2'] = cutout_size

        hdu = fits.PrimaryHDU(data=new_data, header=header)
        if custom_save_path is not None and custom_save_name is not None:
            savename_img = custom_save_path + '/' + custom_save_name+prefix_image_add+'-image.fits'
        else:
            savename_img = os.path.dirname(imagename) + '/' + os.path.basename(imagename).replace('.fits', '.cutout.' + special_name + '.fits')
        hdu.writeto(savename_img, overwrite=True)

    if residualname is not None:
        with fits.open(residualname) as hdul:
            image_data, header = hdul[0].data, hdul[0].header
            wcs = WCS(header, naxis=2)
            
            cutout = Cutout2D(image_data[0][0] if len(image_data.shape) == 4 else image_data, center, cutout_size, wcs=wcs)
            
            if len(image_data.shape) == 4:
                new_data = np.zeros((image_data.shape[0], image_data.shape[1], cutout_size, cutout_size))
                new_data[0,0] = cutout.data
            else:
                new_data = cutout.data

            if correct_shift and ref_cutout_image is not None:
                if len(image_data.shape) == 4:
                    new_data[0,0] = shift(new_data[0,0], (offset_y, offset_x), mode='constant')
                else:
                    new_data = shift(new_data, (offset_y, offset_x), mode='constant')
            header.update(cutout.wcs.to_header())
            if len(image_data.shape) == 4:
                header['NAXIS3'] = image_data.shape[1]
                header['NAXIS4'] = image_data.shape[0]
                header['NAXIS1'] = cutout_size
                header['NAXIS2'] = cutout_size

            hdu = fits.PrimaryHDU(data=new_data, header=header)
            
            if custom_save_path is not None and custom_save_name is not None:
                savename_res = custom_save_path + '/' + custom_save_name+prefix_image_add+'-residual.fits'
            else:
                savename_res = os.path.dirname(residualname) + '/' + os.path.basename(residualname).replace('.fits', '.cutout.' + special_name + '.fits')
            
            hdu.writeto(savename_res, overwrite=True)
    if cut_model:
        model_name = imagename.replace('-image','-model')
        with fits.open(model_name) as hdul:
            image_data, header = hdul[0].data, hdul[0].header
            wcs = WCS(header, naxis=2)
            
            cutout = Cutout2D(image_data[0][0] if len(image_data.shape) == 4 else image_data, center, cutout_size, wcs=wcs)
            
            if len(image_data.shape) == 4:
                new_data = np.zeros((image_data.shape[0], image_data.shape[1], cutout_size, cutout_size))
                new_data[0,0] = cutout.data
            else:
                new_data = cutout.data

            if correct_shift and ref_cutout_image is not None:
                if len(image_data.shape) == 4:
                    new_data[0,0] = shift(new_data[0,0], (offset_y, offset_x), mode='constant')
                else:
                    new_data = shift(new_data, (offset_y, offset_x), mode='constant')
            header.update(cutout.wcs.to_header())
            if len(image_data.shape) == 4:
                header['NAXIS3'] = image_data.shape[1]
                header['NAXIS4'] = image_data.shape[0]
                header['NAXIS1'] = cutout_size
                header['NAXIS2'] = cutout_size

            hdu = fits.PrimaryHDU(data=new_data, header=header)
            
            if custom_save_path is not None and custom_save_name is not None:
                savename_model = custom_save_path + '/' + custom_save_name+prefix_image_add+'-model.fits'
            else:
                savename_model = os.path.dirname(model_name) + '/' + os.path.basename(model_name).replace('.fits', '.cutout.' + special_name + '.fits')
            
            hdu.writeto(savename_model, overwrite=True)
   
    return ra_f, dec_f, savename_img

        
def calculate_pixel_distance(image_file, x0, y0):
    """
    Calculate the relative distance of a pixel from the center of an image
    and return its celestial coordinates.

    Parameters:
    -----------
    image_file : str
        Path to the FITS file.
    x0, y0 : float
        Pixel coordinates of the target position.

    Returns:
    --------
    relative_distance_arcsec : float
        Relative distance of the pixel from the center of the image in arcseconds.
    pixel_coords : str
        Celestial coordinates of the input pixel in the format "RA Dec".
    """
    with fits.open(image_file) as hdul:
        wcs = WCS(hdul[0].header)
        data = hdul[0].data
    celestial_wcs = wcs.celestial

    ny, nx = data.shape[-2:]
    x_center, y_center = nx / 2, ny / 2
    center_coords = celestial_wcs.pixel_to_world(x_center, y_center)
    target_coords = celestial_wcs.pixel_to_world(x0, y0)
    
    relative_distance_arcsec = center_coords.separation(target_coords).arcsecond
    
    pixel_coords = target_coords.to_string('hmsdms')
    
    return relative_distance_arcsec, pixel_coords    
            

"""
 ____              _
/ ___|  __ ___   _(_)_ __   __ _
\___ \ / _` \ \ / / | '_ \ / _` |
 ___) | (_| |\ V /| | | | | (_| |
|____/ \__,_| \_/ |_|_| |_|\__, |
                           |___/
#Saving
"""


def save_results_csv(result_mini, save_name, ext='.csv', save_corr=True,
                     save_params=True):
    values = result_mini.params.valuesdict()
    if save_corr:
        try:
            covariance = result_mini.covar
            covar_df = pd.DataFrame(covariance, index=values.keys(),
                                    columns=values.keys())
            covar_df.to_csv(save_name + '_mini_corr' + ext, index_label='parameter')
        except:
            print('Error saving covariance matrix. Skiping...')

    if save_params:
        try:
            stderr = [result_mini.params[name].stderr for name in values.keys()]
            df = pd.DataFrame({'value': list(values.values()), 'stderr': stderr},
                              index=values.keys())
            df.to_csv(save_name + '_mini_params' + ext, index_label='parameter')
        except:
            print('Errors not present in mini, saving only parameters.')
            df = pd.DataFrame(result_mini.params.valuesdict(), index=['value'])
            df.T.to_csv(save_name + '_mini_params' + ext,
                        index_label='parameter')