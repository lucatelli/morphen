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

def do_cutout(image, box_size=(200, 200), center=None, return_='data'):
    """
    Perform a cutout of a 2D astronomical image array.
    
    COORDINATE CONVENTION:
    - Input coordinates follow (x, y) convention where x=column, y=row
    - Array indexing follows numpy convention: array[y, x] or array[row, column]
    - Image shape is (height, width) corresponding to (n_rows, n_columns)
    
    Parameters
    ----------
    image : str or numpy.ndarray
        Either a file path to a FITS image or a 2D numpy array
    box_size : int or tuple of int, default (200, 200)
        Size of the cutout box. If int, creates a square box.
        Values represent half-width in each dimension.
    center : tuple of int, optional
        Center coordinates (x, y) for the cutout, where x is the column index
        and y is the row index. If None, uses the position of the maximum value in the image.
    return_ : {'data', 'box'}, default 'data'
        What to return: 'data' returns the cutout array,
        'box' returns the bounding box coordinates (xmin, xmax, ymin, ymax)
    
    Returns
    -------
    numpy.ndarray or tuple
        Either the cutout data array or bounding box coordinates
        
    Raises
    ------
    ValueError
        If inputs are invalid or cutout extends beyond image boundaries
    TypeError
        If image type is not supported
    """
    
    # Input validation and normalization
    if isinstance(box_size, int):
        box_size = (box_size, box_size)
    elif not (isinstance(box_size, (tuple, list)) and len(box_size) == 2):
        raise ValueError("box_size must be an integer or a tuple/list of two integers")
    
    box_size = tuple(int(s) for s in box_size)
    if any(s <= 0 for s in box_size):
        raise ValueError("box_size values must be positive")
    
    if return_ not in ['data', 'box']:
        raise ValueError("return_ must be 'data' or 'box'")
    
    # Load data if image is a file path
    if isinstance(image, str):
        try:
            data = load_fits_data(image)
        except Exception as e:
            raise ValueError(f"Failed to load image from path '{image}': {e}")
    elif isinstance(image, np.ndarray):
        data = image
    else:
        raise TypeError("image must be a string (file path) or numpy array")
    
    # Validate data is 2D
    if data.ndim != 2:
        raise ValueError(f"Image data must be 2D, got {data.ndim}D array")
    
    height, width = data.shape
    
    # Determine center coordinates
    if center is None:
        # Find position of maximum value (replacing CASA imstat functionality)
        max_indices = np.unravel_index(np.nanargmax(data), data.shape)
        center_y, center_x = max_indices
        print(f'  >> Center (auto-detected at max value) --> ({center_x}, {center_y})')
    else:
        if not (isinstance(center, (tuple, list)) and len(center) == 2):
            raise ValueError("center must be a tuple/list of two numbers")
        center_x, center_y = int(center[0]), int(center[1])
    
    # Calculate bounding box coordinates
    half_width_x, half_width_y = box_size
    xmin = center_x - half_width_x
    xmax = center_x + half_width_x
    ymin = center_y - half_width_y
    ymax = center_y + half_width_y
    
    # Bounds checking and adjustment
    original_bounds = (xmin, xmax, ymin, ymax)
    
    # Clamp coordinates to image boundaries
    xmin = max(0, xmin)
    xmax = min(width, xmax)
    ymin = max(0, ymin)
    ymax = min(height, ymax)
    
    # Check if any adjustment was needed
    if (xmin, xmax, ymin, ymax) != original_bounds:
        print(f"Warning: Cutout box adjusted to fit within image boundaries")
        print(f"  Original: x=[{original_bounds[0]}:{original_bounds[1]}], "
              f"y=[{original_bounds[2]}:{original_bounds[3]}]")
        print(f"  Adjusted: x=[{xmin}:{xmax}], y=[{ymin}:{ymax}]")
    
    # Validate final box size
    if xmin >= xmax or ymin >= ymax:
        raise ValueError(f"Invalid cutout box: center ({center_x}, {center_y}) "
                        f"with box_size {box_size} results in empty or invalid region")
    
    # Return requested output
    if return_ == 'box':
        return (xmin, xmax, ymin, ymax)
    else:
        # Extract cutout (note: numpy indexing is [y, x])
        cutout = data[ymin:ymax, xmin:xmax]
        
        if cutout.size == 0:
            raise ValueError("Cutout resulted in empty array")
            
        return cutout


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


def copy_header(source_image, target_image, file_to_save=None):
    """
    Copy header information from a source FITS file to a target FITS file.
    Specially designed to work with HST and other multi-extension FITS files.
    
    This function preserves the exact structure of the source file while
    replacing the data with that from the target file.
    
    Parameters
    ----------
    source_image : str
        Path to the source FITS file (from which headers will be copied).
    target_image : str
        Path to the target FITS file (to which headers will be copied).
    file_to_save : str, optional
        If provided, save the result to this path instead of modifying target_image.
        Default is None (modify target_image directly).
        
    Returns
    -------
    str
        Path to the updated FITS file
    """
    from astropy.io import fits
    import numpy as np
    
    # Determine output file
    output_file = file_to_save if file_to_save is not None else target_image
    
    # Read in source and target data
    with fits.open(source_image) as source_hdul, fits.open(target_image) as target_hdul:
        # Create a new HDUList for the output
        new_hdus = fits.HDUList()
        
        # Get the target data
        target_data = None
        if len(target_hdul) > 0 and target_hdul[0].data is not None:
            target_data = target_hdul[0].data
        elif len(target_hdul) > 1:
            # If primary HDU has no data, try the first extension
            target_data = target_hdul[1].data
        
        # Copy each HDU from source, but replace data with target data
        for i, hdu in enumerate(source_hdul):
            # Create a copy of the source HDU
            new_hdu = hdu.copy()
            
            # Replace the data if this is the HDU that should contain it
            if i == 0 and source_hdul[0].data is not None:
                # Primary HDU has data in source, so put target data here
                new_hdu.data = target_data
            elif i == 1 and (source_hdul[0].data is None or len(source_hdul) == 1):
                # First extension should get the data if primary doesn't have it
                new_hdu.data = target_data
            
            # Add to our new HDUList
            new_hdus.append(new_hdu)
    
    # Write to output file
    new_hdus.writeto(output_file, overwrite=True)
    
    return output_file

def get_cell_size_old(imagename):
    """
    Get the cell size/pixel size in arcsec from an image header wcs.
    """
    hdu = pf.open(imagename)
    ww = WCS(hdu[0].header)
    pixel_scale = (ww.pixel_scale_matrix[1,1]*3600)
    cell_size =  pixel_scale.copy()
    return(cell_size)


def get_cell_size(imagename, ext=None):
    """
    Get the pixel scale (cell size) from a FITS header in arcseconds.
    
    Parameters
    ----------
    imagename : str
        Path to the FITS file.
    ext : int, optional
        Extension number to use. If None, will find the best extension.
        
    Returns
    -------
    float
        Cell size in arcseconds.
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    
    with fits.open(imagename) as hdul:
        if ext is None:
            # Find the extension with valid WCS
            for i, hdu in enumerate(hdul):
                try:
                    if hdu.header and ('CD1_1' in hdu.header or 'CDELT1' in hdu.header):
                        ext = i
                        break
                except Exception:
                    pass
            
            # If still None, use primary header
            if ext is None:
                ext = 0
        
        header = hdul[ext].header
        
        # Try different keywords for pixel scale
        if 'CD1_1' in header:
            cell_size = abs(header['CD1_1']) * 3600.0  # deg to arcsec
        elif 'CDELT1' in header:
            cell_size = abs(header['CDELT1']) * 3600.0  # deg to arcsec
        else:
            # Try to get from WCS
            try:
                wcs = WCS(header, naxis=2)
                if wcs.has_celestial:
                    cell_size = wcs.proj_plane_pixel_scales()[0].to('arcsec').value
                else:
                    cell_size = 1.0
            except Exception:
                cell_size = 1.0
                
    return cell_size

def find_sci_extension(hdul):
    """
    Find the science extension with valid data in a FITS file.
    
    Parameters
    ----------
    hdul : HDUList
        The HDUList from a FITS file.
        
    Returns
    -------
    tuple
        (extension_index, data, wcs)
    """
    from astropy.wcs import WCS
    
    # First, try to find SCI extension which is standard for HST
    for i, hdu in enumerate(hdul):
        if hasattr(hdu, 'name') and hdu.name == 'SCI' and hdu.data is not None:
            try:
                wcs = WCS(hdu.header, naxis=2)
                if wcs.has_celestial:
                    return i, hdu.data, wcs
            except Exception:
                pass
    
    # If no SCI extension found, try any extension with data and valid WCS
    for i, hdu in enumerate(hdul):
        if hdu.data is not None:
            try:
                wcs = WCS(hdu.header, naxis=2)
                if wcs.has_celestial:
                    return i, hdu.data, wcs
            except Exception:
                pass
    
    # Fall back to primary HDU
    try:
        wcs = WCS(hdul[0].header, naxis=2)
        return 0, hdul[0].data, wcs
    except Exception:
        return 0, hdul[0].data, None



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


def peak_image_alignment(reference_image, target_image, 
                         mask=None,
                         apply_filter=True):
    """
    Align images using phase correlation with robust preprocessing.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    mask_threshold : float, optional
        Threshold in units of standard deviation for creating masks
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    """
    import numpy as np
    from scipy import fftpack
    from scipy import ndimage
    
    def preprocess_image(image):
        # Normalize image
        image = image - np.nanmean(image)
        image = image / np.nanstd(image)
        
        # # Create mask if threshold provided
        # if mask_threshold is not None:
        #     mask = image > mask_threshold * np.nanstd(image)
        #     # Slightly expand mask to include structure edges
        #     mask = ndimage.binary_dilation(mask, iterations=2)
        #     image = image * mask
        if mask is not None:
            image = image * mask
            
        
        # Apply gentle Gaussian smoothing to reduce noise while preserving structure
        image = ndimage.gaussian_filter(image, sigma=1.0)
        
        return image
    
    # Preprocess both images
    ref_processed = preprocess_image(reference_image)
    target_processed = preprocess_image(target_image)
    
    # Apply windowing to reduce edge effects
    window = np.outer(np.hanning(ref_processed.shape[0]), 
                     np.hanning(ref_processed.shape[1]))
    ref_processed *= window
    target_processed *= window
    
    # Compute FFTs
    F1 = fftpack.fft2(ref_processed)
    F2 = fftpack.fft2(target_processed)
    
    # Compute cross-power spectrum with optional high-freq dampening
    eps = 1e-10  # Small number to prevent division by zero
    cross_power = (F1 * F2.conjugate()) / (np.abs(F1 * F2.conjugate()) + eps)
    
    if apply_filter:    
        # Apply frequency domain filtering to focus on significant structure
        freq_filter = ndimage.gaussian_filter(np.abs(F1 * F2.conjugate()), sigma=2)
        freq_filter = freq_filter / np.max(freq_filter)
        cross_power *= freq_filter
    
    
    # Inverse FFT
    cc = np.real(fftpack.ifft2(cross_power))
    
    # Find peak with sub-pixel precision using center of mass of peak region
    max_loc = np.unravel_index(np.argmax(cc), cc.shape)
    peak_region = ndimage.binary_dilation(cc == np.max(cc), iterations=2)
    refined_y, refined_x = ndimage.center_of_mass(cc * peak_region)
    
    # Convert to shifts
    shift_y = refined_y if refined_y < ref_processed.shape[0]//2 \
        else refined_y - ref_processed.shape[0]
    shift_x = refined_x if refined_x < ref_processed.shape[1]//2 \
        else refined_x - ref_processed.shape[1]
    
    return shift_y, shift_x


def apply_alignment(reference_image, target_image, return_shifted=False):
    """
    Align target image to reference image and optionally return the shifted image.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    return_shifted : bool
        If True, returns the aligned image
        
    Returns:
    --------
    tuple or tuple, array
        (shift_y, shift_x) or ((shift_y, shift_x), aligned_image)
    """
    from scipy.ndimage import shift
    
    # Get optimal shifts
    shift_y, shift_x = robust_image_alignment(reference_image, target_image)
    
    if return_shifted:
        # Apply shift with spline interpolation for better accuracy
        aligned_image = shift(target_image, (shift_y, shift_x), 
                            mode='constant', order=3)
        return (shift_y, shift_x), aligned_image
    
    return shift_y, shift_x





def structural_image_alignment(reference_image, target_image,
                               mask=None,
                               apply_filter=True
                              ):
    """
    Align images using phase correlation based on overall structure,
    without any dependence on peak positions.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    """
    import numpy as np
    from scipy import fftpack
    from scipy import ndimage
    
    def preprocess_image(image):
        # Normalize image
        image = image - np.nanmean(image)
        image = image / np.nanstd(image)
        
        # Apply gentle Gaussian smoothing to reduce noise while preserving structure
        image = ndimage.gaussian_filter(image, sigma=1.0)
        if mask is not None:
            image = image * mask
        
        return image
    
    # Preprocess both images
    ref_processed = preprocess_image(reference_image)
    target_processed = preprocess_image(target_image)

    if apply_filter:
        # Apply windowing to reduce edge effects
        window = np.outer(np.hanning(ref_processed.shape[0]), 
                         np.hanning(ref_processed.shape[1]))
        ref_processed *= window
        target_processed *= window
    
    # Compute FFTs
    F1 = fftpack.fft2(ref_processed)
    F2 = fftpack.fft2(target_processed)
    
    # Compute normalized cross-power spectrum
    eps = 1e-10
    cross_power = (F1 * F2.conjugate()) / (np.abs(F1 * F2.conjugate()) + eps)
    
    # Inverse FFT to get correlation surface
    cc = np.real(fftpack.ifft2(cross_power))
    
    # Find the shift from the entire correlation surface
    shift = np.unravel_index(np.argmax(cc), cc.shape)
    
    # Convert to relative shifts
    shift_y = shift[0] if shift[0] < ref_processed.shape[0]//2 \
        else shift[0] - ref_processed.shape[0]
    shift_x = shift[1] if shift[1] < ref_processed.shape[1]//2 \
        else shift[1] - ref_processed.shape[1]
    
    return shift_y, shift_x


def align_by_phase_lmfit_enhanced(reference_image, target_image, search_window=20):
    """
    Enhanced version of phase correlation alignment using lmfit.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    search_window : int
        Maximum pixel shift to consider in each direction
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    lmfit.ModelResult
        The full fit result object with additional statistics
    """
    import numpy as np
    from scipy import fftpack
    from scipy.ndimage import shift, gaussian_filter
    from lmfit import Parameters, minimize, Minimizer
    
    # Preprocess images
    def preprocess(image):
        # Normalize
        norm = (image - np.mean(image)) / np.std(image)
        
        # Apply Hanning window
        window = np.outer(np.hanning(image.shape[0]), 
                         np.hanning(image.shape[1]))
        return norm * window
    
    ref_processed = preprocess(reference_image)
    target_processed = preprocess(target_image)
    # Pre-compute FFT of reference
    F1 = fftpack.fft2(ref_processed)
    # F2 = fftpack.fft2(target_processed)
    
    # Create frequency weight matrix (emphasize mid frequencies)
    fy = fftpack.fftfreq(F1.shape[0])[:, np.newaxis]
    fx = fftpack.fftfreq(F1.shape[1])[np.newaxis, :]
    freq_weight = np.exp(-(fx**2 + fy**2) / 0.1)  # Adjust 0.1 as needed
    
    def objective(params):
        """
        Enhanced objective function using weighted phase correlation.
        """
        shift_y = params['shift_y'].value
        shift_x = params['shift_x'].value
        
        # Shift and preprocess target
        shifted = shift(target_processed, (shift_y, shift_x), 
                       mode='constant', cval=0, order=3)
        # target_processed = preprocess(shifted)
        
        # Compute FFT of shifted target
        # F2 = fftpack.fft2(target_processed)
        F2 = fftpack.fft2(shifted)
        # Compute weighted cross-power spectrum
        # cross_power = F1 * F2.conjugate() * freq_weight
        cross_power = F1 * F2.conjugate()
        
        # Normalize
        eps = 1e-10
        normalized_cross_power = cross_power / (np.abs(cross_power) + eps)
        
        # Get correlation and apply Gaussian smoothing for stability
        correlation = np.abs(fftpack.ifft2(normalized_cross_power))
        # correlation = gaussian_filter(correlation, sigma=1.0)
        
        # Compute metric
        return -np.max(correlation)
    
    # Set up parameters
    params = Parameters()
    params.add('shift_y', value=0, min=-search_window, max=search_window)
    params.add('shift_x', value=0, min=-search_window, max=search_window)

    mini = Minimizer(objective, params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')
    
    # Perform the minimization
    # result = mini.minimize(method='nelder')
    result = mini.minimize(method='least_squares', loss="cauchy",
                           tr_solver="exact",verbose=0)
    
    # # Try different optimization methods
    # methods = ['nelder', 'powell', 'cobyla']
    # best_result = None
    # best_metric = float('inf')
    
    # for method in methods:
    #     try:
    #         result = minimize(objective, params, method=method)
    #         if result.success and result.residual < best_metric:
    #             best_result = result
    #             best_metric = result.residual
    #     except:
    #         continue
    
    # if best_result is None:
    #     raise ValueError("None of the optimization methods succeeded")
    
    # shift_y = best_result.params['shift_y'].value
    # shift_x = best_result.params['shift_x'].value
    shift_y = result.params['shift_x'].value
    shift_x = result.params['shift_y'].value
    
    return shift_y, shift_x


def align_by_lmfit(reference_image, target_image, search_window=10):
    """
    Align images using lmfit to minimize the squared difference between them.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    search_window : int
        Maximum pixel shift to consider in each direction
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    lmfit.ModelResult
        The full fit result object with additional statistics
    """
    import numpy as np
    from scipy.ndimage import shift
    from lmfit import Parameters, minimize, Minimizer
    
    # Normalize images
    ref_norm = (reference_image - np.mean(reference_image)) / np.std(reference_image)
    target_norm = (target_image - np.mean(target_image)) / np.std(target_image)
    
    def objective(params):
        """
        Objective function to minimize.
        Returns array of residuals (will be squared internally by lmfit).
        """
        shift_y = params['shift_y'].value
        shift_x = params['shift_x'].value
        
        # Shift the target image
        shifted = shift(target_norm, (shift_y, shift_x), mode='constant', cval=0)
        
        # Return flattened residuals
        return (ref_norm - shifted).ravel()
    
    # Set up parameters with bounds
    params = Parameters()
    params.add('shift_y', value=0, min=-search_window, max=search_window)
    params.add('shift_x', value=0, min=-search_window, max=search_window)

    mini = Minimizer(objective, params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')
    
    # Perform the minimization
    result = mini.minimize(method='least_squares', loss="cauchy",
                           tr_solver="exact",verbose=0)
    # result = minimize(objective, params, method='leastsq')
    
    # Get the optimal shifts
    shift_y = result.params['shift_y'].value
    shift_x = result.params['shift_x'].value
    
    return shift_y, shift_x



def align_by_ncc_lmfit(reference_image, target_image, search_window=25):
    """
    Align images using lmfit to maximize normalized cross correlation.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    search_window : int
        Maximum pixel shift to consider in each direction
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    lmfit.ModelResult
        The full fit result object with additional statistics
    """
    import numpy as np
    from scipy.ndimage import shift
    from lmfit import Parameters, minimize, Minimizer
    
    # Normalize images
    ref_norm = (reference_image - np.mean(reference_image)) / np.std(reference_image)
    target_norm = (target_image - np.mean(target_image)) / np.std(target_image)
    
    def objective(params):
        """
        Objective function to minimize (negative NCC).
        Returns scalar value.
        """
        shift_y = params['shift_y'].value
        shift_x = params['shift_x'].value
        
        # Shift the target image
        shifted = shift(target_norm, (shift_y, shift_x), mode='constant', cval=0)
        
        # Compute normalized cross correlation
        numerator = np.sum(ref_norm * shifted)
        denominator = np.sqrt(np.sum(ref_norm**2) * np.sum(shifted**2))
        ncc = numerator / denominator if denominator != 0 else 0
        
        # Return negative NCC (since we're minimizing)
        return -ncc
    
    # Set up parameters with bounds
    params = Parameters()
    params.add('shift_y', value=0, min=-search_window, max=search_window)
    params.add('shift_x', value=0, min=-search_window, max=search_window)
    
    mini = Minimizer(objective, params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')
    
    # Perform the minimization
    result = mini.minimize(method='least_squares', loss="cauchy",
                           tr_solver="exact",verbose=0)
    
    # Get the optimal shifts
    shift_y = result.params['shift_y'].value
    shift_x = result.params['shift_x'].value
    
    return shift_y, shift_x


def align_by_phase_lmfit(reference_image, target_image, search_window=20):
    """
    Align images using lmfit to optimize phase correlation in Fourier space.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    search_window : int
        Maximum pixel shift to consider in each direction
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    lmfit.ModelResult
        The full fit result object with additional statistics
    """
    import numpy as np
    from scipy import fftpack
    from scipy.ndimage import shift
    from lmfit import Parameters, minimize, Minimizer
    
    # Normalize images
    ref_norm = (reference_image - np.mean(reference_image)) / np.std(reference_image)
    target_norm = (target_image - np.mean(target_image)) / np.std(target_image)
    
    # Compute FFT of reference image once
    F1 = fftpack.fft2(ref_norm)
    
    def objective(params):
        """
        Objective function to minimize negative phase correlation.
        """
        shift_y = params['shift_y'].value
        shift_x = params['shift_x'].value
        
        # Shift target image
        shifted = shift(target_norm, (shift_y, shift_x), mode='constant', cval=0)
        
        # Compute FFT of shifted target
        F2 = fftpack.fft2(shifted)
        
        # Compute cross-power spectrum
        cross_power = F1 * F2.conjugate()
        
        # Normalize to get only phase information
        eps = 1e-10  # Small number to prevent division by zero
        normalized_cross_power = cross_power / (np.abs(cross_power) + eps)
        
        # Inverse FFT to get correlation
        correlation = np.abs(fftpack.ifft2(normalized_cross_power))
        
        # Compute metric to minimize
        # We use negative correlation since we're minimizing
        return -np.max(correlation)
    
    # Set up parameters with bounds
    params = Parameters()
    params.add('shift_y', value=0, min=-search_window, max=search_window)
    params.add('shift_x', value=0, min=-search_window, max=search_window)
    
    # Perform the minimization
    # Using 'nelder' method as it works well with this type of objective function
    # result = minimize(objective, params, method='nelder')
    mini = Minimizer(objective, params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')
    
    # Perform the minimization
    result = mini.minimize(method='least_squares', loss="cauchy",
                           tr_solver="exact",verbose=0)
    
    # Get the optimal shifts
    shift_y = result.params['shift_x'].value
    shift_x = result.params['shift_y'].value
    
    return shift_y, shift_x

def align_by_minimization(reference_image, target_image, search_window=10):
    """
    Align images by minimizing the squared difference between them.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    search_window : int
        Maximum pixel shift to consider in each direction
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    """
    import numpy as np
    from scipy.ndimage import shift
    from scipy.optimize import minimize
    
    # Normalize images
    ref_norm = (reference_image - np.mean(reference_image)) / np.std(reference_image)
    target_norm = (target_image - np.mean(target_image)) / np.std(target_image)
    
    def compute_ssd(params):
        """
        Compute sum of squared differences for given shift parameters.
        """
        shift_y, shift_x = params
        # Shift the target image
        shifted = shift(target_norm, (shift_y, shift_x), mode='constant', cval=0)
        # Compute squared difference
        diff = (ref_norm - shifted) ** 2
        return np.sum(diff)
    
    # Initial guess (can be [0,0] or from a coarse estimation)
    initial_guess = [0, 0]
    
    # Set bounds for the optimization
    bounds = [(-search_window, search_window), (-search_window, search_window)]
    
    # Minimize the SSD
    result = minimize(compute_ssd, 
                     initial_guess,
                     bounds=bounds,
                     method='L-BFGS-B')
    
    return result.x[0], result.x[1]


def align_by_ncc_minimization(reference_image, target_image, search_window=10):
    """
    Align images by maximizing normalized cross correlation.
    
    Parameters:
    -----------
    reference_image : 2D numpy array
        The reference image
    target_image : 2D numpy array
        The image to be aligned
    search_window : int
        Maximum pixel shift to consider in each direction
        
    Returns:
    --------
    tuple : (shift_y, shift_x)
        The optimal shift to align target with reference
    """
    import numpy as np
    from scipy.ndimage import shift
    from scipy.optimize import minimize
    
    # Normalize images
    ref_norm = (reference_image - np.mean(reference_image)) / np.std(reference_image)
    target_norm = (target_image - np.mean(target_image)) / np.std(target_image)
    
    def compute_ncc(params):
        """
        Compute negative normalized cross correlation (negative because we minimize)
        """
        shift_y, shift_x = params
        # Shift the target image
        shifted = shift(target_norm, (shift_y, shift_x), mode='constant', cval=0)
        
        # Compute normalized cross correlation
        numerator = np.sum(ref_norm * shifted)
        denominator = np.sqrt(np.sum(ref_norm**2) * np.sum(shifted**2))
        ncc = numerator / denominator if denominator != 0 else 0
        
        # Return negative because minimize looks for minimum
        return -ncc
    
    # Initial guess
    initial_guess = [0, 0]
    
    # Set bounds for the optimization
    bounds = [(-search_window, search_window), (-search_window, search_window)]
    
    # Minimize negative NCC (equivalent to maximizing NCC)
    result = minimize(compute_ncc, 
                     initial_guess,
                     bounds=bounds,
                     method='L-BFGS-B')
    
    return result.x[0], result.x[1]


def cutout_2D_radec(imagename, residualname=None, ra_f=None, dec_f=None, cutout_size=1024,
                    special_name='', correct_shift=False, ref_cutout_image=None, pixel_coords=None,
                    mask=None, apply_filter=True,shift_correction_mode='peak'):
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
        if mask is not None:
            # mask = do_cutout_2D(mask, box_size=cutout_size, 
            #                     center=center, return_='data')
            mask = Cutout2D(mask.astype(int), center, cutout_size, wcs=wcs)
            mask=(mask.data).astype(bool)
        
        # apply shift
        if correct_shift:
            if ref_cutout_image is not None:
                ref_image_cutout_data = load_fits_data(ref_cutout_image)
                # x_ref, y_ref = nd.maximum_position(ref_image_cutout_data)[::-1]
                # reference_source_coords = (x_ref, y_ref)
                # offset_x, offset_y = find_offsets(reference_source_coords,
                #                                   nd.maximum_position(cutout.data)[::-1])
                if shift_correction_mode == 'peak':
                    # print(" ++==>> Applying peak-based image alignment.")
                    offset_y, offset_x = \
                        peak_image_alignment(ref_image_cutout_data,
                                            cutout.data,
                                            mask=mask,
                                            apply_filter=apply_filter
                        )
                elif shift_correction_mode == 'structural':
                    # print(" ++==>> Applying structural-based image alignment.")
                    offset_y, offset_x = \
                        structural_image_alignment(ref_image_cutout_data,
                                            cutout.data,
                                            mask=mask,
                        )
                elif shift_correction_mode == 'image_diff':
                    # print(" ++==>> Applying image difference-based image alignment.")
                    raise ValueError("Image difference-based alignment not implemented yet.")

                else:
                    raise ValueError("Invalid shift correction mode. Choose 'peak' or 'structural'.")
                
                print(f"        > Offset of image position is: ({int(offset_x)}, {int(offset_y)}).")
                
                aligned_target_image = shift(cutout.data, (int(offset_y), int(offset_x)),
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
                    aligned_target_image = shift(cutout.data, (int(offset_y), int(offset_x)),
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
                       mask=None, apply_filter=True,shift_correction_mode='peak',
                       custom_save_path=None, custom_save_name=None):
    from astropy.io import fits
    import os
    from astropy.wcs import WCS
    from astropy.nddata import Cutout2D
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord
    
    if custom_save_path is not None:
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
        if mask is not None:
            # mask = do_cutout_2D(mask, box_size=cutout_size, 
            #                     center=center, return_='data')
            mask = Cutout2D(mask.astype(int), center, cutout_size, wcs=wcs)
            mask=(mask.data).astype(bool)
        
        if len(image_data.shape) == 4:
            new_data = np.zeros((image_data.shape[0], image_data.shape[1], cutout_size, cutout_size))
            new_data[0,0] = cutout.data
        else:
            new_data = cutout.data

        if correct_shift:
            if ref_cutout_image is not None:
                ref_image_cutout_data = load_fits_data(ref_cutout_image)
                # x_ref, y_ref = nd.maximum_position(ref_image_cutout_data)[::-1]
                # reference_source_coords = (x_ref, y_ref)
                # offset_x, offset_y = find_offsets(reference_source_coords,
                #                                 nd.maximum_position(cutout.data)[::-1])
                # print(f" !!!! Offsets of peak position are: {offset_x, offset_y}.")
                if shift_correction_mode == 'peak':
                    # print(" ++==>> Applying peak-based image alignment.")
                    offset_y, offset_x = \
                        peak_image_alignment(ref_image_cutout_data,
                                            cutout.data,
                                            mask=mask,
                                            apply_filter=apply_filter
                        )
                elif shift_correction_mode == 'structural':
                    # print(" ++==>> Applying structural-based image alignment.")
                    offset_y, offset_x = \
                        structural_image_alignment(ref_image_cutout_data,
                                            cutout.data,
                                            mask=mask,
                        )
                elif shift_correction_mode == 'image_diff':
                    # print(" ++==>> Applying image difference-based image alignment.")
                    raise ValueError("Image difference-based alignment not implemented yet.")

                else:
                    raise ValueError("Invalid shift correction mode. Choose 'peak' or 'structural'.")
                
                print(f"        > Offset of image position is: ({int(offset_x)}, {int(offset_y)}).")
                if len(image_data.shape) == 4:
                    new_data[0,0] = shift(new_data[0,0], (int(offset_y), int(offset_x)), mode='constant')
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