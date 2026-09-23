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

def load_fits_data(image, ext=None, plane=None):
    '''
    Read a FITS image into a 2D numpy array. Pure astropy -- no CASA.

    Name origin:
    ctn > casa to numpy. The function used to open images through
    `casatools.image` and undo the result with `np.rot90(a)[::-1, :]`, described
    at the time as undoing a "mirrored and 90-degree rotated" array. That
    expression is in fact exactly `a.T`, and all it did was convert CASA's
    `getchunk()` [x, y] axis order into the [y, x] order that astropy already
    returns. Verified on radio 4D, HSC, Legacy Survey, PSF, HST MEF and JWST MEF
    files: both routes give bit-identical arrays. So astropy needs no correction
    and the CASA dependency bought nothing.

    The returned array is indexed [row, column] = [y, x], as `astropy.io.fits`
    gives it, which is the convention every consumer in this codebase assumes.

    Parameters
    ----------
    image : str or numpy.ndarray
        Path to a FITS file. An array is returned unchanged, so call sites may
        pass either without guarding on the type.
    ext : int or str, optional
        Extension to read, by index or by EXTNAME. Default: the extension named
        'SCI' if one holds data (the HST/JWST convention), else the first
        extension that holds data.
    plane : int, optional
        Which plane to take when the data has a real (non-degenerate) extra
        axis, e.g. an RGB composite or a spectral cube. Default: plane 0, with a
        warning naming the shape. Degenerate axes are unwrapped silently.

    Returns
    -------
    numpy.ndarray
        A 2D in-memory array (not a memmap).

    Raises
    ------
    TypeError
        `image` is neither a string nor an array.
    FileNotFoundError
        No such path, or it is a CASA image directory rather than a FITS file.
    ValueError
        The file exists but cannot be read as FITS, or the requested extension
        or plane does not exist.
    '''
    if isinstance(image, np.ndarray):
        return image
    if not isinstance(image, str):
        raise TypeError(f"load_fits_data: expected a file path or a numpy "
                        f"array, got {type(image).__name__}.")

    if not os.path.exists(image):
        raise FileNotFoundError(f"load_fits_data: no such file {image!r}.")
    if os.path.isdir(image):
        # The old CASA branch could open these; this one deliberately cannot.
        # Nothing else in the pipeline can either -- get_cell_size,
        # compute_image_properties and the plotting all go through astropy --
        # so a CASA image has to be converted before it is of any use here.
        raise FileNotFoundError(
            f"load_fits_data: {image!r} is a directory, not a FITS file. "
            f"CASA images are no longer read directly; convert it first with "
            f"exportfits(imagename={image!r}, fitsimage={image!r}+'.fits').")

    def _has_data(hdu):
        """True if this HDU holds an image, from the header alone.

        Deliberately header-only: touching `hdu.data` to find out would both
        force a read of every extension scanned and raise outright on scaled
        integer data (see the memmap note below).
        """
        return int(hdu.header.get('NAXIS', 0) or 0) >= 2

    def _select(hdul):
        if ext is not None:
            try:
                hdu = hdul[ext]
            except (KeyError, IndexError) as exc:
                raise ValueError(f"load_fits_data: no extension {ext!r} in "
                                 f"{os.path.basename(image)}.") from exc
            if not _has_data(hdu):
                raise ValueError(f"load_fits_data: extension {ext!r} of "
                                 f"{os.path.basename(image)} holds no data.")
            return hdu
        # HST and JWST products put the image in an extension named SCI and
        # leave the primary HDU empty, alongside WHT/ERR/DQ/CON planes of
        # the same shape. Name it explicitly rather than trusting position.
        hdu = next((h for h in hdul
                    if getattr(h, 'name', '') == 'SCI' and _has_data(h)), None)
        if hdu is None:
            hdu = next((h for h in hdul if _has_data(h)), None)
        if hdu is None:
            raise ValueError(f"load_fits_data: no extension of "
                             f"{os.path.basename(image)} holds any data.")
        return hdu

    try:
        hdul = pf.open(image, memmap=True)
    except Exception as exc:
        raise ValueError(f"load_fits_data: could not read {image!r} as "
                         f"FITS ({exc}).") from exc

    # try/finally rather than `with hdul:` -- the block may REBIND hdul below,
    # and a `with` would close the original handle while leaking the new one.
    try:
        hdu = _select(hdul)
        # Scaled data (BZERO/BSCALE/BLANK, as written for unsigned integers and
        # by some segmentation maps) cannot be memory-mapped -- astropy refuses
        # with "Cannot load a memory-mapped image ... Set memmap=False" -- so
        # reopen without it. Decided from the header, before any read, because
        # the check itself is what would raise.
        if any(k in hdu.header for k in ('BZERO', 'BSCALE', 'BLANK')):
            hdul.close()
            try:
                hdul = pf.open(image, memmap=False)
            except Exception as exc:
                raise ValueError(f"load_fits_data: could not read {image!r} as "
                                 f"FITS ({exc}).") from exc
            hdu = _select(hdul)

        data = hdu.data
        if data is None:
            raise ValueError(f"load_fits_data: the selected extension of "
                             f"{os.path.basename(image)} holds no data.")
        original_shape = data.shape
        if data.ndim < 2:
            raise ValueError(f"load_fits_data: {os.path.basename(image)} holds "
                             f"a {data.ndim}D array {original_shape}, not an "
                             f"image.")

        # Peel leading axes until the array is 2D. Radio images carry degenerate
        # Stokes/frequency axes -- (1, 1, ny, nx) and (1, ny, nx) -- which are
        # pure packaging and go silently. A leading axis longer than 1 is real
        # data being discarded, so it is announced. Only leading axes are
        # touched, so a legitimately narrow (ny, 1) image stays 2D, which
        # np.squeeze would not respect.
        while data.ndim > 2:
            if data.shape[0] == 1:
                data = data[0]
            else:
                idx = 0 if plane is None else int(plane)
                if not -data.shape[0] <= idx < data.shape[0]:
                    raise ValueError(
                        f"load_fits_data: plane {idx} is out of range for "
                        f"{os.path.basename(image)} with shape "
                        f"{original_shape}.")
                if plane is None:
                    print(f"!! WARNING !! load_fits_data: "
                          f"{os.path.basename(image)} has shape "
                          f"{original_shape}; taking plane {idx} of "
                          f"{data.shape[0]}. Pass plane=<i> to choose.")
                data = data[idx]

        # FITS stores data big-endian, so astropy hands back a byte-swapped
        # dtype ('>f4'). JAX refuses any array whose dtype is not in native
        # byte order -- jit raises "Error interpreting argument ... as an
        # abstract array" -- and the CASA branch this replaced never exposed
        # that, because casacore always returned a native float64 copy.
        # Normalise here, once, instead of at each of the ~195 call sites.
        #
        # Floating-point data additionally comes back as float64, which is the
        # precision the pipeline was written against. Integer data keeps its
        # own type in native order: casacore forced those to float64 too, but
        # corrupted the values doing it (checked on a uint32 segmentation map),
        # so there is nothing there worth reproducing.
        data = np.asarray(data)
        if data.dtype.kind == 'f':
            out_dtype = np.float64
        elif not data.dtype.isnative:
            out_dtype = data.dtype.newbyteorder('=')
        else:
            out_dtype = data.dtype
        # np.array, not np.asarray: this must materialise the memmap before the
        # file closes even when no conversion is needed.
        return np.array(data, dtype=out_dtype)
    finally:
        hdul.close()

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

# def do_cutout(image, box_size=(200,200), center=None, return_='data'):
#     """
#     Fast cutout of a numpy array.

#     Returs: numpy data array or a box for that cutout, if asked.
#     """
#     if isinstance(box_size, int):
#         box_size = (box_size,box_size)

#     if center is None:
#         if isinstance(image, str) == True:
#             # imhd = imhead(image)
#             st = imstat(image)
#             print('  >> Center --> ', st['maxpos'])
#             xin, xen, yin, yen = (st['maxpos'][0] - box_size[0],
#                                   st['maxpos'][0] + box_size[0],
#                                   st['maxpos'][1] - box_size[1],
#                                   st['maxpos'][1] + box_size[1])
#             data_cutout = load_fits_data(image)[xin:xen, yin:yen]

#         else:
#             try:
#                 max_y, max_x = np.where(load_fits_data(image) == load_fits_data(image).max())
#                 xin = max_x[0] - box_size[0]
#                 xen = max_x[0] + box_size[0]
#                 yin = max_y[0] - box_size[1]
#                 yen = max_y[0] + box_size[1]
#             except:
#                 max_y, max_x = np.where(image == image.max())
#                 xin = max_x[0] - box_size[0]
#                 xen = max_x[0] + box_size[0]
#                 yin = max_y[0] - box_size[1]
#                 yen = max_y[0] + box_size[1]

#             data_cutout = image[xin:xen, yin:yen]


#     else:
#         xin, xen, yin, yen = (center[0] - box_size[0], center[0] + box_size[0],
#                               center[1] - box_size[1], center[1] + box_size[1])
#         if isinstance(image, str) == True:
#             data_cutout = load_fits_data(image)[xin:xen, yin:yen]
#         else:
#             data_cutout = image[xin:xen, yin:yen]

#     if return_ == 'data':
#         return (data_cutout)
#     if return_ == 'box':
#         box = xin, xen, yin, yen  # [xin:xen,yin:yen]
#         return (box)


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




def copy_header_old(image, image_to_copy, file_to_save=None):
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


# def copy_header(source_image, target_image, file_to_save=None, exclude_keywords=None):
#     """
#     Copy header information from a source FITS file to a target FITS file,
#     with automatic handling for complex files like HST that have multiple headers.
    
#     Parameters
#     ----------
#     source_image : str
#         Path to the source FITS file (from which headers will be copied).
#     target_image : str
#         Path to the target FITS file (to which headers will be copied).
#     file_to_save : str, optional
#         If provided, save the result to this path instead of modifying target_image.
#         Default is None (modify target_image directly).
#     exclude_keywords : list, optional
#         List of header keywords to exclude from copying.
#         Default is None (copy all keywords).
        
#     Returns
#     -------
#     str
#         Path to the updated FITS file
    
#     Notes
#     -----
#     - This function automatically handles HST files by copying both primary and 
#       extension headers to maintain astrometric information.
#     - For HST files, both the primary header and the first extension header are copied.
#     - For standard FITS files, only the primary header is copied.
#     """
#     from astropy.io import fits
#     import os
#     import warnings
    
#     # Set up default parameters
#     if exclude_keywords is None:
#         exclude_keywords = ['HISTORY', 'COMMENT']  # Often these should be excluded
    
#     # Create a new file if file_to_save is specified
#     output_file = file_to_save if file_to_save is not None else target_image
    
#     # If we're creating a new file, copy the target first
#     if file_to_save is not None:
#         with fits.open(target_image) as hdul:
#             hdul.writeto(file_to_save, overwrite=True)
    
#     # Open both files
#     with fits.open(source_image) as source_hdul:
#         with fits.open(output_file, mode='update') as target_hdul:
            
#             # Check if this is likely an HST file by looking for typical HST keywords
#             is_hst = False
#             if len(source_hdul) > 1:
#                 primary_header = source_hdul[0].header
#                 if ('TELESCOP' in primary_header and primary_header['TELESCOP'] == 'HST') or \
#                    ('INSTRUME' in primary_header and primary_header['INSTRUME'] in 
#                     ['ACS', 'WFC3', 'STIS', 'NICMOS', 'WFPC2', 'COS', 'FOS', 'GHRS']):
#                     is_hst = True
            
#             # Handle HST files specially
#             if is_hst:
#                 print(f"Detected HST file. Copying both primary and extension headers.")
                
#                 # First, ensure the target has the right structure (at least 2 HDUs)
#                 if len(target_hdul) < 2:
#                     warnings.warn("Target file doesn't have enough HDUs for HST data. "
#                                  "Attempting to match source structure.")
                    
#                     # Create a new fits file with the same structure as the source
#                     temp_file = output_file + ".temp"
#                     with fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU()]) as new_hdul:
#                         # Copy data from target
#                         if len(target_hdul) > 0:
#                             new_hdul[0].data = target_hdul[0].data
#                         new_hdul.writeto(temp_file, overwrite=True)
                    
#                     # Close current handles and reopen with new structure
#                     target_hdul.close()
#                     os.replace(temp_file, output_file)
#                     target_hdul = fits.open(output_file, mode='update')
                
#                 # Copy primary header (HDU 0)
#                 for key, value in source_hdul[0].header.items():
#                     if key not in exclude_keywords:
#                         target_hdul[0].header[key] = value
                
#                 # Copy first extension header (HDU 1) - contains vital WCS info for HST
#                 if len(source_hdul) > 1 and len(target_hdul) > 1:
#                     for key, value in source_hdul[1].header.items():
#                         if key not in exclude_keywords:
#                             target_hdul[1].header[key] = value
#                 else:
#                     warnings.warn("Could not copy extension header - missing in source or target")
            
#             # Standard case - just copy the primary header
#             else:
#                 for key, value in source_hdul[0].header.items():
#                     if key not in exclude_keywords:
#                         target_hdul[0].header[key] = value
            
#             # Save changes
#             target_hdul.flush()
    
#     return output_file

# def copy_header(source_image, target_image, file_to_save=None, 
#                 hdu_mapping=None, exclude_keywords=None, overwrite_mode='replace'):
#     """
#     Copy header information from a source FITS file to a target FITS file,
#     with support for multiple HDUs and customizable behavior.
    
#     Parameters
#     ----------
#     source_image : str
#         Path to the source FITS file (from which headers will be copied).
#     target_image : str
#         Path to the target FITS file (to which headers will be copied).
#     file_to_save : str, optional
#         If provided, save the result to this path instead of modifying target_image.
#         Default is None (modify target_image directly).
#     hdu_mapping : dict, optional
#         Dictionary mapping source HDU indices to target HDU indices.
#         For example, {0: 0, 1: 1} copies the 0th HDU header from source to 0th HDU
#         in target, and 1st HDU header from source to 1st HDU in target.
#         Default is None (copy all HDUs that exist in both files with matching indices).
#     exclude_keywords : list, optional
#         List of header keywords to exclude from copying.
#         Default is None (copy all keywords).
#     overwrite_mode : str, optional
#         How to handle existing header keywords in target.
#         Options:
#             'update' - Update target with source, keeping target keywords not in source
#             'replace' - Replace target header completely with source header
#             'preserve' - Only add keywords from source that don't exist in target
#         Default is 'update'.
        
#     Returns
#     -------
#     str
#         Path to the updated FITS file
    
#     Notes
#     -----
#     - This function handles complex FITS files with multiple headers/HDUs.
#     - Common usage is for copying WCS information from data to model/residual images.
#     - For HST and similar files with multiple headers, use the hdu_mapping parameter
#       to specify which headers should be copied where.
#     """
#     from astropy.io import fits
#     import os
#     import warnings
    
#     # Set up default parameters
#     if exclude_keywords is None:
#         exclude_keywords = []
    
#     # Create a new file if file_to_save is specified
#     if file_to_save is not None:
#         # Copy target_image to file_to_save first
#         with fits.open(target_image) as hdul:
#             hdul.writeto(file_to_save, overwrite=True)
#         output_file = file_to_save
#     else:
#         output_file = target_image
    
#     # Open the source image
#     with fits.open(source_image) as source_hdul:
#         source_hdu_count = len(source_hdul)
        
#         # Open the target image for updating
#         with fits.open(output_file, mode='update') as target_hdul:
#             target_hdu_count = len(target_hdul)
            
#             # Determine which HDUs to copy
#             if hdu_mapping is None:
#                 # Default: copy matching HDUs
#                 hdu_mapping = {i: i for i in range(min(source_hdu_count, target_hdu_count))}
            
#             # Process each HDU according to the mapping
#             for source_idx, target_idx in hdu_mapping.items():
#                 # Check if indices are valid
#                 if source_idx >= source_hdu_count:
#                     warnings.warn(f"Source HDU index {source_idx} does not exist. Skipping.")
#                     continue
#                 if target_idx >= target_hdu_count:
#                     warnings.warn(f"Target HDU index {target_idx} does not exist. Skipping.")
#                     continue
                
#                 source_header = source_hdul[source_idx].header
#                 target_header = target_hdul[target_idx].header
                
#                 # Handle different overwrite modes
#                 if overwrite_mode == 'replace':
#                     # Create a new header object to avoid reference issues
#                     new_header = fits.Header()
#                     for key, value in source_header.items():
#                         if key not in exclude_keywords:
#                             new_header[key] = value
#                     target_hdul[target_idx].header = new_header
                
#                 elif overwrite_mode == 'update':
#                     # Update target with source, preserving target keys not in source
#                     for key, value in source_header.items():
#                         if key not in exclude_keywords:
#                             target_header[key] = value
                
#                 elif overwrite_mode == 'preserve':
#                     # Only add keys that don't exist in target
#                     for key, value in source_header.items():
#                         if key not in exclude_keywords and key not in target_header:
#                             target_header[key] = value
            
#             # Save changes
#             target_hdul.flush()
    
#     return output_file

# def copy_header(source_image, target_image, file_to_save=None, exclude_keywords=None):
#     """
#     Copy header information from a source FITS file to a target FITS file,
#     preserving the multi-extension structure essential for HST and other complex FITS files.
    
#     Parameters
#     ----------
#     source_image : str
#         Path to the source FITS file (from which headers will be copied).
#     target_image : str
#         Path to the target FITS file (to which headers will be copied).
#     file_to_save : str, optional
#         If provided, save the result to this path instead of modifying target_image.
#         Default is None (modify target_image directly).
#     exclude_keywords : list, optional
#         List of header keywords to exclude from copying.
#         Default is None (excludes only HISTORY and COMMENT by default).
        
#     Returns
#     -------
#     str
#         Path to the updated FITS file
#     """
#     from astropy.io import fits
#     import os
#     import warnings
    
#     # Set default exclusions
#     if exclude_keywords is None:
#         exclude_keywords = ['HISTORY', 'COMMENT']
    
#     # Determine output file
#     output_file = file_to_save if file_to_save is not None else target_image
    
#     # Step 1: Open both files to analyze their structure
#     with fits.open(source_image) as source_hdul:
#         with fits.open(target_image) as target_hdul:
            
#             # Step 2: Create a new HDUList that will preserve the target data with source headers
#             new_hdul = fits.HDUList()
            
#             # Step 3: For each HDU in the source, create a corresponding HDU in the new list
#             for i, source_hdu in enumerate(source_hdul):
#                 # If we have a corresponding HDU in target, use its data
#                 if i < len(target_hdul):
#                     # Create appropriate HDU type with target data but empty header
#                     if isinstance(source_hdu, fits.PrimaryHDU):
#                         new_hdu = fits.PrimaryHDU(data=target_hdul[i].data)
#                     elif isinstance(source_hdu, fits.ImageHDU):
#                         new_hdu = fits.ImageHDU(data=target_hdul[i].data)
#                     elif isinstance(source_hdu, fits.BinTableHDU):
#                         new_hdu = fits.BinTableHDU(data=target_hdul[i].data)
#                     elif isinstance(source_hdu, fits.TableHDU):
#                         new_hdu = fits.TableHDU(data=target_hdul[i].data)
#                     else:
#                         # For other HDU types, create a generic image HDU
#                         new_hdu = fits.ImageHDU(data=target_hdul[i].data)
                
#                 # If target doesn't have this HDU, create empty one with same type as source
#                 else:
#                     if isinstance(source_hdu, fits.PrimaryHDU):
#                         new_hdu = fits.PrimaryHDU()
#                     elif isinstance(source_hdu, fits.ImageHDU):
#                         new_hdu = fits.ImageHDU()
#                     elif isinstance(source_hdu, fits.BinTableHDU):
#                         new_hdu = fits.BinTableHDU()
#                     elif isinstance(source_hdu, fits.TableHDU):
#                         new_hdu = fits.TableHDU()
#                     else:
#                         new_hdu = fits.ImageHDU()
                
#                 # Step 4: Copy header from source to the new HDU
#                 for key, value in source_hdu.header.items():
#                     if key not in exclude_keywords:
#                         new_hdu.header[key] = value
                
#                 # Add this HDU to our new HDUList
#                 new_hdul.append(new_hdu)
            
#             # Step 5: Add any remaining HDUs from target that don't exist in source
#             for i in range(len(source_hdul), len(target_hdul)):
#                 new_hdul.append(target_hdul[i].copy())
    
#     # Step 6: Write the new HDUList to the output file
#     new_hdul.writeto(output_file, overwrite=True)
    
#     # Optional: Verify the file was created correctly
#     try:
#         with fits.open(output_file) as test_hdul:
#             num_hdus = len(test_hdul)
#             print(f"Successfully created FITS file with {num_hdus} HDUs/extensions")
#     except Exception as e:
#         warnings.warn(f"Created file may have issues: {str(e)}")
    
#     return output_file

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


# def get_cell_size(imagename, ext=None):
#     """
#     Get the pixel scale (cell size) from a FITS header in arcseconds.
    
#     Parameters
#     ----------
#     imagename : str
#         Path to the FITS file.
#     ext : int, optional
#         Extension number to use. If None, will find the best extension.
        
#     Returns
#     -------
#     float
#         Cell size in arcseconds.
#     """
#     from astropy.io import fits
#     from astropy.wcs import WCS
    
#     with fits.open(imagename) as hdul:
#         if ext is None:
#             # Find the extension with valid WCS
#             for i, hdu in enumerate(hdul):
#                 try:
#                     if hdu.header and ('CD1_1' in hdu.header or 'CDELT1' in hdu.header):
#                         ext = i
#                         break
#                 except Exception:
#                     pass
            
#             # If still None, use primary header
#             if ext is None:
#                 ext = 0
        
#         header = hdul[ext].header
        
#         # Try different keywords for pixel scale
#         if 'CD1_1' in header:
#             cell_size = abs(header['CD1_1']) * 3600.0  # deg to arcsec
#         elif 'CDELT1' in header:
#             cell_size = abs(header['CDELT1']) * 3600.0  # deg to arcsec
#         else:
#             # Try to get from WCS
#             try:
#                 wcs = WCS(header, naxis=2)
#                 if wcs.has_celestial:
#                     cell_size = wcs.proj_plane_pixel_scales()[0].to('arcsec').value
#                 else:
#                     cell_size = 1.0
#             except Exception:
#                 cell_size = 1.0
                
#     return cell_size


def get_cell_size(imagename, ext=None):
    """
    Get the pixel scale (cell size) from a FITS header in arcseconds.
    Properly handles CD matrix transformation for JWST and other instruments.
    
    Parameters
    ----------
    imagename : str
        Path to the FITS file.
    ext : int, optional
        Extension number to use. If None, will find the best extension.
        
    Returns
    -------
    float
        Cell size in arcseconds (average of x and y scales if different).
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    import numpy as np
    
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
        
        # Priority 1: Use WCS projection method (most robust for all cases)
        try:
            wcs = WCS(header, naxis=2)
            if wcs.has_celestial:
                # Get pixel scales using proper WCS method
                pixel_scales = wcs.proj_plane_pixel_scales()
                # Average of x and y scales in arcseconds
                cell_size = np.mean([ps.to('arcsec').value for ps in pixel_scales])
                return cell_size
        except Exception as e:
            pass  # Fall back to manual calculation
        
        # Priority 2: Use CD matrix if present (handles rotation properly)
        if 'CD1_1' in header and 'CD2_2' in header:
            cd1_1 = header.get('CD1_1', 0.0)
            cd1_2 = header.get('CD1_2', 0.0)
            cd2_1 = header.get('CD2_1', 0.0)
            cd2_2 = header.get('CD2_2', 0.0)
            
            # Calculate pixel scales accounting for rotation
            # This is the proper way to extract pixel scale from CD matrix
            scale_x = np.sqrt(cd1_1**2 + cd2_1**2)
            scale_y = np.sqrt(cd1_2**2 + cd2_2**2)
            
            # Average of x and y scales, convert to arcseconds
            cell_size = np.mean([scale_x, scale_y]) * 3600.0
            
            return cell_size
        
        # Priority 3: Use CDELT if present and no CD matrix
        # Note: CDELT should NOT be used if CD matrix exists (as in JWST)
        elif 'CDELT1' in header and 'CD1_1' not in header:
            cdelt1 = abs(header['CDELT1'])
            cdelt2 = abs(header.get('CDELT2', cdelt1))
            
            # Average of x and y scales, convert to arcseconds
            cell_size = np.mean([cdelt1, cdelt2]) * 3600.0
            
            return cell_size
        
        # Priority 4: Default fallback
        else:
            print(f"Warning: Could not determine pixel scale for {imagename}")
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


# Keywords that describe how THIS HDU is laid out on disk. They belong to the
# extension they came from and must never be inherited from another one.
_STRUCTURAL_HEADER_KEYS = {
    'SIMPLE', 'XTENSION', 'BITPIX', 'NAXIS', 'EXTEND', 'PCOUNT', 'GCOUNT',
    'EXTNAME', 'EXTVER', 'EXTLEVEL', 'INHERIT', 'BSCALE', 'BZERO', 'BLANK',
    'CHECKSUM', 'DATASUM', 'END',
}


def merge_science_header(hdul, ext):
    """
    Header for a science extension, with the primary header's keywords folded
    in -- everything a standalone copy of that extension needs to stand on its
    own.

    In a multi-extension file the description of the observation is split
    across two headers: HST and JWST put INSTRUME/TELESCOP/FILTER/DATE-OBS in
    the empty PRIMARY and PHOTFNU/PHOTFLAM/PHOTPLAM/BUNIT plus the WCS in SCI.
    Writing a cutout from the SCI header alone therefore produces a file that
    no longer says which instrument or filter made it, which is exactly what
    `read_data.get_jy_conversion_factor` and `get_filter_info` read.

    The science header wins every conflict -- it owns the WCS, the units and
    the data layout -- and structural keywords are never inherited.

    Parameters
    ----------
    hdul : HDUList
    ext : int
        Index of the science extension, e.g. from `find_sci_extension`.

    Returns
    -------
    astropy.io.fits.Header
        A copy; the input is untouched. For ``ext == 0`` this is just
        ``hdul[0].header.copy()``, so single-extension files (radio, Legacy
        Survey, HSC) are unaffected.
    """
    header = hdul[ext].header.copy()
    if ext == 0:
        return header
    for card in hdul[0].header.cards:
        key = card.keyword
        if not key or key in _STRUCTURAL_HEADER_KEYS:
            continue
        if key.startswith('NAXIS'):
            continue
        if key in ('COMMENT', 'HISTORY'):
            # Commentary cards are provenance (drizzle logs, calibration
            # notes). They carry no meaning that can clash, so append rather
            # than test for presence -- `key in header` is meaningless here.
            header.append(card, end=True)
            continue
        if key in header:
            continue
        header.append(card, end=True)
    return header



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
    Deprecated. Kept as a thin wrapper over `estimate_image_shift`.

    The original implementation used *phase* correlation: the cross-power
    spectrum divided by its own modulus. For beam-convolved radio images that is
    the wrong statistic -- the restoring beam is a low-pass filter, so whitening
    the spectrum amplifies exactly the frequencies that hold nothing but noise.
    Measured against a 2.37 px injected offset it returned 2.03 px, and when one
    component of a two-component source was made 10x brighter in the target (an
    inverted-spectrum component, which is routine across a C-to-Ka baseline) it
    was wrong by 11.6 px.

    `apply_filter` is accepted and ignored; there is no frequency-domain filter
    in the replacement.

    Returns
    -------
    (shift_y, shift_x)
        Same sign convention as before: the shift to apply to `target_image`.
    """
    from image_alignment import estimate_image_shift
    result = estimate_image_shift(reference_image, target_image,
                                  method='chi2', mask=mask)
    return result.dy, result.dx




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
    Deprecated. Kept as a thin wrapper over `estimate_image_shift`.

    The original implementation could only ever return ``(0, 0)`` whenever a
    `mask` was supplied. It multiplied *both* images by the *same* mask and then
    took a whitened cross-power spectrum; the mask's autocorrelation dominates
    that surface and pins its maximum at zero shift. It returned exactly
    ``(0.00, 0.00)`` for a synthetic 2.37 px offset and for a real cross-array
    offset alike, so any pipeline that appeared to run with
    ``shift_correction_mode='structural'`` was in fact applying no correction.

    Returns
    -------
    (shift_y, shift_x)
        The shift to apply to `target_image`.
    """
    from image_alignment import estimate_image_shift
    result = estimate_image_shift(reference_image, target_image,
                                  method='chi2', mask=mask)
    return result.dy, result.dx




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


def cutout_2D_radec(imagename, residualname=None, modelname=None, ra_f=None, dec_f=None, cutout_size=1024,
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
            
    if modelname is not None:
        with fits.open(modelname) as hdul:
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
            savename_res = os.path.dirname(modelname) + '/' + os.path.basename(
                modelname).replace('.fits', '.cutout.' + special_name + '.fits')
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

    # Tag the output with the WSClean sub-band/product token (e.g. '-MFS' from
    # '<prefix>-MFS-image.fits'). Guarded because for a name without '-image'
    # the split falls through to the last dash-separated chunk of the *full
    # path*, which then gets spliced into the output filename and makes the
    # write fail. Generic (non-WSClean) names simply get no tag.
    basename = os.path.basename(imagename)
    if '-image' in basename:
        prefix_image_add = f"-{basename.split('-image')[0].split('-')[-1]}"
    else:
        prefix_image_add = ''

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


def t_cutout_2D_radec(imagename, residualname=None, modelname=None,
                      ra_f=None, dec_f=None, pixel_coords=None,
                      cutout_size=1024, mode='trim',
                      correct_shift=False, ref_cutout_image=None,
                      shift_correction_mode='auto', mask=None,
                      max_shift=None, precomputed_shift=None,
                      apply_filter=True,
                      custom_save_path=None, custom_save_name=None,
                      special_name='', overwrite=True,
                      return_paths=False, verbose=1):
    """
    Definitive RA/Dec cutout routine, merging `cutout_2D_radec` and
    `cutout_2D_radec_v2`.

    ``t_`` for *testing*: this is intended to replace both of the above once it
    has been exercised on real data. What it takes from each, and what it fixes:

    From ``cutout_2D_radec`` (v1)
        - An explicit ``modelname=`` argument. v2's ``cut_model=True`` instead
          derived the model path as ``imagename.replace('-image', '-model')``,
          which silently re-cuts the *image* and writes it out as a model
          whenever the name contains no ``-image``, and raises
          ``FileNotFoundError`` whenever no model exists next to the image.

    From ``cutout_2D_radec_v2``
        - 4D (``NAXIS=4``) cubes are written back as 4D with the degenerate
          Stokes/frequency axes preserved. v1 built ``new_data`` for this and
          then never used it, so a 4D input was silently written out as 2D
          with a stale 4D header.
        - ``custom_save_path`` / ``custom_save_name``.

    Fixed here in both
        - **Edge cutouts.** ``Cutout2D`` returns a *trimmed* (smaller) array
          near an image border, but both parents allocated
          ``(cutout_size, cutout_size)`` and assigned into it, raising a shape
          mismatch. The output array now follows the actual cutout shape, and
          ``NAXIS1/2`` follow it too. Pass ``mode='partial'`` to keep the full
          requested box, NaN-padded, instead of trimming.
        - **Shift consistency.** v2 applied an integer shift to the image but a
          sub-pixel shift to the residual and model, leaving them misaligned
          with each other by up to a pixel. All products now get the same
          integer shift.
        - **`correct_shift` without a reference.** v1 left ``new_hdul``
          unassigned and died later with a confusing ``NameError``.
        - **The WSClean filename token.** ``imagename.split('-image')[0]
          .split('-')[-1]`` falls through to the last dash-separated chunk of
          the *full path* for a generic filename, splicing it into the output
          name and making the write fail. It is now basename-scoped and skipped
          entirely when there is no ``-image`` token.

    Parameters
    ----------
    imagename : str
        Image to cut. Required.
    residualname, modelname : str, optional
        Companion products, cut on the identical sky position and box, with the
        identical shift applied.
    ra_f, dec_f : float, optional
        Cutout centre in **degrees**. If both are None and `pixel_coords` is
        None, the image peak is used (via CASA `imstat`).
    pixel_coords : tuple, optional
        (x, y) pixel position, converted to RA/Dec via the image WCS.
    cutout_size : int
        Box side in pixels.
    mode : {'trim', 'partial'}
        Passed to `Cutout2D`. 'trim' (default, the historical behaviour) returns
        a smaller array at image borders; 'partial' keeps the requested size and
        pads with NaN.
    correct_shift : bool
        Align the cutout to `ref_cutout_image` before writing.
    shift_correction_mode : str
        Alignment estimator, passed through to
        `image_alignment.estimate_image_shift`: one of 'chi2', 'mi', 'xcorr',
        'chi2_shift', 'auto' (default) or 'ensemble'. The legacy values 'peak',
        'structural' and 'image_diff' still work and map onto 'chi2', with a
        DeprecationWarning.

        The shift is applied in two parts: the integer pixels by moving the
        cutout window (exact, no interpolation, no zero-padded edge) and the
        sub-pixel remainder by a Fourier phase ramp. The previous version
        rounded the whole shift to an integer, which on sub-pixel offsets --
        the common case for images off a shared grid -- turned a real 0.3 px
        misalignment into an applied 1 px one.
    max_shift : float, optional
        Largest correction the estimator may search for and return, in pixels.
        Defaults to one reference beam (`bmaj / cell`, floored at 3 px), which
        guards against an estimator that has locked onto the wrong component --
        a real failure mode, and the one `component_flip_suspected` reports.
        That guard is wrong when the misalignment is genuinely larger than a
        beam: the estimate is then truncated onto the `max_shift` radius and
        flagged `'clamped'`, so the cutout comes out better aligned but still
        visibly off. Raise it explicitly in that case (the flag, and a
        `dy_err` that collapses to zero because the optimiser is sitting on its
        own bound rather than in a minimum, are the two symptoms).
    apply_filter : bool
        Ignored. Kept so existing call sites do not break; the replacement
        estimators have no frequency-domain filter.
    precomputed_shift : (dy, dx), optional
        Skip estimation and apply this shift instead. Use it to estimate once
        per band on the high-SNR MFS image and reuse the answer for every
        sub-band of that band; sub-bands from a single imaging run share a grid,
        so re-estimating per sub-band only adds scatter.
    custom_save_path, custom_save_name : str, optional
        Output directory and basename. When both are given the products are
        written as ``<name><token>-image.fits`` / ``-residual.fits`` /
        ``-model.fits``; otherwise they land next to their inputs as
        ``<original>.cutout.<special_name>.fits``.
    return_paths : bool
        If True return a dict of every product path plus the centre and WCS,
        instead of the ``(ra_f, dec_f, savename_img)`` triple that both parent
        functions return.

    Returns
    -------
    (ra_f, dec_f, savename_img) or dict
    """
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

    # WSClean sub-band/product token (e.g. '-MFS' from '<pre>-MFS-image.fits'),
    # basename-scoped and skipped for generic names -- see the docstring.
    basename = os.path.basename(imagename)
    if '-image' in basename:
        prefix_image_add = f"-{basename.split('-image')[0].split('-')[-1]}"
    else:
        prefix_image_add = ''

    # --- resolve the cutout centre -----------------------------------------
    # The celestial WCS lives in the science extension, which is HDU 0 only for
    # single-extension files. HST drz and JWST i2d products have an empty
    # PRIMARY (NAXIS=0, no WCS), and `WCS(hdul[0].header, naxis=2)` on one of
    # those silently yields a non-celestial identity WCS whose `pixel_to_world`
    # returns a list of Quantity rather than a SkyCoord.
    with fits.open(imagename) as hdul:
        _sci_ext, _, wcs0 = find_sci_extension(hdul)

    if wcs0 is None or not wcs0.has_celestial:
        raise ValueError(
            f'No celestial WCS found in {os.path.basename(imagename)}; '
            f'cannot place a cutout on the sky.')

    if pixel_coords is not None:
        sky = wcs0.pixel_to_world(*pixel_coords)
        ra_f, dec_f = sky.ra.degree, sky.dec.degree
        if verbose >= 1:
            print(f'Converted pixel coordinates {pixel_coords} to RA/Dec: '
                  f'({ra_f}, {dec_f})')
    elif ra_f is None or dec_f is None:
        imst = imstat(imagename)
        if verbose >= 1:
            print('maxpos = ', imst['maxpos'])
            print('maxposf = ', imst['maxposf'])
        coords = imst['maxposf'].split(',')
        ra_f, dec_f = conver_str_coords(coords[0], format_coords(coords[1]))
        if verbose >= 1:
            print('ra_f, dec_f = ', ra_f, dec_f)

    center = SkyCoord(ra=ra_f * u.degree, dec=dec_f * u.degree, frame='icrs')

    # The shift is computed once, on the image, and reused verbatim for the
    # residual and model so the three stay on a common grid.
    from image_alignment import (estimate_image_shift, apply_pixel_shift,
                                 split_shift)

    # dy/dx are the *rounded* integer parts, kept under their historical names
    # so `return_paths`' 'offset' key means what it always did. The sub-pixel
    # remainder and the full estimate live alongside them.
    offsets = {'dy': 0, 'dx': 0, 'dy_frac': 0.0, 'dx_frac': 0.0,
               'dy_full': 0.0, 'dx_full': 0.0, 'result': None}

    def _cut_one(filename, product):
        """Cut `filename` on the shared centre/box and write it out."""
        with fits.open(filename) as hdul:
            # Same multi-extension handling as the centre resolution above.
            # The cutout is written as a single PrimaryHDU, so the science
            # header is merged with the primary one first -- otherwise an
            # HST/JWST cutout would lose INSTRUME/FILTER and stop reporting
            # its own instrument and units when read back in.
            ext, image_data, wcs = find_sci_extension(hdul)
            if image_data is None or wcs is None or not wcs.has_celestial:
                raise ValueError(
                    f'No science extension with a celestial WCS found in '
                    f'{os.path.basename(filename)}.')
            header = merge_science_header(hdul, ext)

            is_4d = len(image_data.shape) == 4
            plane = image_data[0][0] if is_4d else image_data

            cutout = Cutout2D(plane, center, cutout_size, wcs=wcs, mode=mode,
                              fill_value=np.nan)
            cut_data = cutout.data
            ny, nx = cut_data.shape
            # The nominal cutout defines the output grid for every product,
            # shifted or not -- see the WCS note below.
            out_wcs = cutout.wcs

            local_mask = None
            if mask is not None:
                local_mask = Cutout2D(np.asarray(mask).astype(int), center,
                                      cutout_size, wcs=wcs, mode=mode,
                                      fill_value=0)
                local_mask = (local_mask.data).astype(bool)

            # --- alignment ---------------------------------------------------
            if correct_shift:
                if product == 'image':
                    if precomputed_shift is not None:
                        dy, dx = (float(precomputed_shift[0]),
                                  float(precomputed_shift[1]))
                        if verbose >= 1:
                            print(f'        > Using precomputed shift '
                                  f'({dy:.3f}, {dx:.3f}) px.')
                    else:
                        if ref_cutout_image is None:
                            raise ValueError('correct_shift=True requires '
                                             '`ref_cutout_image` (or '
                                             '`precomputed_shift`).')
                        result = estimate_image_shift(
                            ref_cutout_image, cut_data,
                            method=shift_correction_mode,
                            mask=local_mask, max_shift=max_shift,
                            verbose=verbose)
                        offsets['result'] = result
                        dy, dx = result.dy, result.dx
                        if verbose >= 1:
                            print(f'        > Offset of image position is: '
                                  f'({dx:.3f}, {dy:.3f}) px {result.flags}')
                        if 'clamped' in result.flags and verbose >= 1:
                            # Worth spelling out: the applied shift is a bound,
                            # not a measurement, so the cutout will still be
                            # visibly misaligned and no amount of re-running
                            # changes it.
                            print(f'        > NOTE: truncated at max_shift='
                                  f'{np.hypot(dy, dx):.3f} px; the estimator '
                                  f'wanted ({result.dx_raw:.3f}, '
                                  f'{result.dy_raw:.3f}) px. Pass a larger '
                                  f'max_shift= to apply it in full.')
                    (offsets['dy'], offsets['dx'],
                     offsets['dy_frac'], offsets['dx_frac']) = split_shift(dy, dx)
                    offsets['dy_full'], offsets['dx_full'] = dy, dx

                dy_i, dx_i = offsets['dy'], offsets['dx']
                dy_f, dx_f = offsets['dy_frac'], offsets['dx_frac']

                # The integer part is absorbed by moving the cutout window,
                # which pulls in real neighbouring pixels instead of the zero
                # padding `scipy.ndimage.shift` leaves behind, and interpolates
                # nothing at all. A box centred at (y0 - dy, x0 - dx) holds
                # exactly the same pixels as shifting the nominal box by
                # (+dy, +dx), minus the padding.
                if dy_i or dx_i:
                    x0, y0 = wcs.world_to_pixel(center)
                    recut = Cutout2D(plane, (float(x0) - dx_i, float(y0) - dy_i),
                                     cutout_size, wcs=wcs, mode=mode,
                                     fill_value=np.nan)
                    if recut.data.shape == cut_data.shape:
                        cut_data = recut.data
                    else:
                        # Only near an image border, where mode='trim' clips the
                        # moved box to a different size. Fall back to shifting
                        # the nominal cutout so the output grid stays fixed.
                        cut_data = apply_pixel_shift(cut_data, dy_i, dx_i)
                # ...leaving at most half a pixel, applied as a flux-conserving
                # Fourier phase ramp.
                if dy_f or dx_f:
                    cut_data = apply_pixel_shift(cut_data, dy_f, dx_f)
                ny, nx = cut_data.shape

            # --- rebuild the array with the original dimensionality ---------
            if is_4d:
                new_data = np.zeros((image_data.shape[0], image_data.shape[1],
                                     ny, nx))
                new_data[0, 0] = cut_data
            else:
                new_data = cut_data

            # The output WCS is the *nominal* one -- the grid of the reference
            # cutout -- not the moved window's. Applying an astrometric shift is
            # the assertion that the reference's astrometry is the correct one,
            # so the aligned product must report the reference's sky mapping;
            # otherwise the pixels line up while the headers disagree, which is
            # what the previous version silently produced. The original mapping
            # is preserved in OCRVAL*/OCRPIX* and ASTRSHF*.
            if correct_shift and (offsets['dy_full'] or offsets['dx_full']):
                for key in ('CRVAL1', 'CRVAL2', 'CRPIX1', 'CRPIX2'):
                    if key in header:
                        header['O' + key] = (header[key],
                                             'pre-alignment ' + key)
                header['ASTRSHFY'] = (offsets['dy_full'],
                                      'applied y shift [pix]')
                header['ASTRSHFX'] = (offsets['dx_full'],
                                      'applied x shift [pix]')
                header['ASTRMETH'] = (str(shift_correction_mode),
                                      'alignment estimator')
                _res = offsets.get('result')
                if _res is not None and 'clamped' in _res.flags:
                    # The applied shift is a bound, not a measurement. Say so
                    # in the file, or a later reader sees only a clean number.
                    header['ASTRCLMP'] = (True,
                                          'shift truncated at max_shift')
                header['ASTRREF'] = (os.path.basename(ref_cutout_image)[:60]
                                     if ref_cutout_image else 'precomputed',
                                     'alignment reference')
                header.add_history(
                    f'morphen: aligned by ({offsets["dy_full"]:.3f}, '
                    f'{offsets["dx_full"]:.3f}) px onto the reference grid')
            header.update(out_wcs.to_header())
            # NAXIS follows the *actual* cutout, which differs from
            # `cutout_size` whenever mode='trim' clipped it at an image border.
            header['NAXIS1'] = nx
            header['NAXIS2'] = ny
            if is_4d:
                header['NAXIS3'] = image_data.shape[1]
                header['NAXIS4'] = image_data.shape[0]

            if custom_save_path is not None and custom_save_name is not None:
                savename = os.path.join(
                    custom_save_path,
                    f'{custom_save_name}{prefix_image_add}-{product}.fits')
            else:
                savename = os.path.join(
                    os.path.dirname(filename),
                    os.path.basename(filename).replace(
                        '.fits', '.cutout.' + special_name + '.fits'))

            fits.PrimaryHDU(data=new_data, header=header).writeto(
                savename, overwrite=overwrite)
            return savename, cutout.wcs

    savename_img, cutout_wcs = _cut_one(imagename, 'image')

    savename_res = None
    if residualname is not None:
        savename_res, _ = _cut_one(residualname, 'residual')

    savename_model = None
    if modelname is not None:
        savename_model, _ = _cut_one(modelname, 'model')

    if return_paths:
        return {'ra': ra_f, 'dec': dec_f,
                'image': savename_img,
                'residual': savename_res,
                'model': savename_model,
                'wcs': cutout_wcs,
                'offset': (offsets['dx'], offsets['dy']),
                'shift': (offsets['dy_full'], offsets['dx_full']),
                'shift_result': offsets['result']}

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
            
def get_radec_from_pixel(image_input, pixel_coords=None):
    """
    Convert pixel coordinates to RA/Dec using WCS information from image.
    
    Parameters
    ----------
    image_input : str, astropy.io.fits.Header, or astropy.wcs.WCS
        Either a FITS filename, a FITS header with WCS info, or a WCS object
    pixel_coords : tuple of (x, y), optional
        Pixel coordinates (x, y) in zero-indexed convention.
        If None, uses the image center.
    
    Returns
    -------
    ra_deg : float
        Right Ascension in degrees
    dec_deg : float
        Declination in degrees
    
    Examples
    --------
    >>> # Use image center
    >>> ra, dec = get_radec_from_pixel('image.fits')
    
    >>> # Specific pixel position
    >>> ra, dec = get_radec_from_pixel('image.fits', pixel_coords=(512, 512))
    
    >>> # With WCS object directly
    >>> from astropy.wcs import WCS
    >>> wcs = WCS(header)
    >>> ra, dec = get_radec_from_pixel(wcs, pixel_coords=(100, 200))
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    import numpy as np
    
    # Handle different input types
    if isinstance(image_input, str):
        # Input is a filename
        with fits.open(image_input) as hdul:
            header = hdul[0].header
            data = hdul[0].data
            wcs = WCS(header, naxis=2)
            
            # Get image shape for center calculation
            if len(data.shape) == 4:
                # Radio data format (freq, stokes, y, x)
                img_shape = data.shape[2:]
            elif len(data.shape) == 2:
                # Simple 2D image
                img_shape = data.shape
            else:
                # 3D or other format
                img_shape = data.shape[-2:]
                
    elif isinstance(image_input, WCS):
        # Input is already a WCS object
        wcs = image_input
        # Get shape from WCS if available
        if hasattr(wcs, '_naxis'):
            img_shape = (wcs._naxis[1], wcs._naxis[0])
        else:
            img_shape = None
            
    else:
        # Assume it's a FITS header
        wcs = WCS(image_input, naxis=2)
        # Try to get shape from header
        if 'NAXIS1' in image_input and 'NAXIS2' in image_input:
            img_shape = (image_input['NAXIS2'], image_input['NAXIS1'])
        else:
            img_shape = None
    
    # Determine pixel coordinates to convert
    if pixel_coords is None:
        # Use image center
        if img_shape is None:
            raise ValueError("Cannot determine image center without shape information. "
                           "Please provide pixel_coords explicitly.")
        center_x = img_shape[1] / 2.0
        center_y = img_shape[0] / 2.0
        pixel_coords = (center_x, center_y)
        print(f"Using image center: pixel ({center_x:.1f}, {center_y:.1f})")
    
    # Convert pixel to world coordinates
    sky_coord = wcs.pixel_to_world(*pixel_coords)
    ra_deg = sky_coord.ra.degree
    dec_deg = sky_coord.dec.degree
    
    # print(f"Pixel ({pixel_coords[0]:.1f}, {pixel_coords[1]:.1f}) -> "
    #       f"RA={ra_deg:.6f} deg, Dec={dec_deg:.6f} deg")
    
    return ra_deg, dec_deg


def get_pixel_from_radec(image_input, ra_deg, dec_deg):
    """
    Convert RA/Dec coordinates to pixel position using WCS information.
    
    This is the inverse operation of get_radec_from_pixel.
    
    Parameters
    ----------
    image_input : str, astropy.io.fits.Header, or astropy.wcs.WCS
        Either a FITS filename, a FITS header with WCS info, or a WCS object
    ra_deg : float
        Right Ascension in degrees
    dec_deg : float
        Declination in degrees
    
    Returns
    -------
    x_pixel : float
        X pixel coordinate (zero-indexed)
    y_pixel : float
        Y pixel coordinate (zero-indexed)
    
    Examples
    --------
    >>> x, y = get_pixel_from_radec('image.fits', 150.0, 2.5)
    """
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    
    # Handle different input types
    if isinstance(image_input, str):
        with fits.open(image_input) as hdul:
            header = hdul[0].header
            wcs = WCS(header, naxis=2)
    elif isinstance(image_input, WCS):
        wcs = image_input
    else:
        # Assume it's a FITS header
        wcs = WCS(image_input, naxis=2)
    
    # Create SkyCoord object
    sky_coord = SkyCoord(ra=ra_deg * u.degree, dec=dec_deg * u.degree, frame='icrs')
    
    # Convert to pixel coordinates
    x_pixel, y_pixel = wcs.world_to_pixel(sky_coord)
    
    print(f"RA={ra_deg:.6f} deg, Dec={dec_deg:.6f} deg -> "
          f"Pixel ({x_pixel:.1f}, {y_pixel:.1f})")
    
    return x_pixel, y_pixel




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
                     save_params=True, image_name=None, other_param_dict=None):
    values = result_mini.params.valuesdict()
    
    if save_corr:
        try:
            covariance = result_mini.covar
            covar_df = pd.DataFrame(covariance, index=values.keys(),
                                    columns=values.keys())
            covar_df.to_csv(save_name + '_mini_corr' + ext, index_label='parameter')
        except:
            print('Error saving covariance matrix. Skipping...')

    if save_params:
        try:
            # Create a dictionary for the single row
            row_data = {}
            
            # Add image name as first column
            if image_name is None:
                image_name = save_name  # Use save_name if no image_name provided
            row_data['#imagename'] = os.path.basename(image_name)
            
            # Add each parameter and its error as separate columns
            for param_name in values.keys():
                row_data[param_name] = values[param_name]
                stderr = result_mini.params[param_name].stderr
                row_data[param_name + '_err'] = stderr if stderr is not None else None
            
            row_data['chisq'] = result_mini.chisqr
            row_data['red_chisq'] = result_mini.redchi
            row_data['aic'] = result_mini.aic
            row_data['bic'] = result_mini.bic
            # Add other parameters from the dictionary if provided
            if other_param_dict is not None:
                for key, value in other_param_dict.items():
                    row_data[key] = value
            # Create DataFrame with single row
            df = pd.DataFrame([row_data])
            df.to_csv(save_name + '_mini_params' + ext, index=False)
            
        except Exception as e:
            print(f'Error saving parameters with errors: {e}')
            print('Saving only parameter values without errors.')
            row_data = {}
            row_data['#imagename'] = os.path.basename(image_name) if image_name is not None else os.path.basename(save_name)
            for param_name in values.keys():
                row_data[param_name] = values[param_name]
            df = pd.DataFrame([row_data])
            df.to_csv(save_name + '_mini_params' + ext, index=False)



