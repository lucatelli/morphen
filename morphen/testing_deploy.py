def adaptive_box_size(data_shape, min_boxes=20, max_box_size=256, min_box_size=32):
    """
    Determine appropriate box size based on image dimensions.
    
    Parameters
    ----------
    data_shape : tuple
        Shape of the image data (height, width).
    min_boxes : int
        Minimum number of boxes in each dimension.
    max_box_size : int
        Maximum allowed box size.
    min_box_size : int
        Minimum allowed box size.
        
    Returns
    -------
    box_size : tuple
        Optimal (width, height) for background estimation boxes.
    """
    h, w = data_shape
    
    # Calculate box sizes to achieve at least min_boxes in each dimension
    bh = max(min_box_size, min(max_box_size, h // min_boxes))
    bw = max(min_box_size, min(max_box_size, w // min_boxes))
    
    # Ensure box sizes are even numbers for better grid alignment
    bw = bw + (bw % 2)
    bh = bh + (bh % 2)
    
    return (bw, bh)


def detect_image_type(header=None, data_shape=None, pixel_scale=None):
    """
    Detect if image is likely HST, Legacy Survey, or other type.
    
    Parameters
    ----------
    header : astropy.io.fits.Header, optional
        FITS header with image metadata.
    data_shape : tuple, optional
        Shape of the image data (height, width).
    pixel_scale : float, optional
        Pixel scale in arcseconds, if known.
        
    Returns
    -------
    image_type : str
        'hst', 'legacy', or 'unknown'
    parameters : dict
        Recommended parameters for background estimation.
    """
    parameters = {}
    
    # Check header for telescope information if available
    if header is not None:
        telescope = header.get('TELESCOP', '').lower()
        instrument = header.get('INSTRUME', '').lower()
        
        if 'hst' in telescope or any(ins in instrument for ins in ['acs', 'wfc3', 'wfpc2']):
            # HST images typically need smaller boxes relative to image size
            # but larger filter sizes due to PSF and resolution
            parameters['min_boxes'] = 30
            parameters['filter_fraction'] = 3  # Filter size relative to box size
            parameters['sigma'] = 2.5
            return 'hst', parameters
            
        if any(key in header for key in ['LEGPIPEV', 'TRACTID', 'SURVEY']):
            # Legacy Survey images typically work well with larger boxes
            parameters['min_boxes'] = 15
            parameters['filter_fraction'] = 5
            parameters['sigma'] = 3.0
            return 'legacy', parameters
    
    # Use image dimensions and pixel scale as fallback
    if data_shape is not None:
        h, w = data_shape
        
        # HST images are typically larger
        if h > 2000 and w > 2000:
            parameters['min_boxes'] = 25
            parameters['filter_fraction'] = 3
            parameters['sigma'] = 2.5
            
            # If we also know pixel scale, we can be more confident
            if pixel_scale is not None and pixel_scale < 0.1:  # HST has small pixel scale
                return 'hst', parameters
        else:
            parameters['min_boxes'] = 15
            parameters['filter_fraction'] = 5
            parameters['sigma'] = 3.0
            
            # If we also know pixel scale, we can be more confident
            if pixel_scale is not None and pixel_scale > 0.2:  # Legacy Survey has larger pixel scale
                return 'legacy', parameters
    
    # Default parameters for unknown image types
    parameters['min_boxes'] = 20
    parameters['filter_fraction'] = 4
    parameters['sigma'] = 3.0
    
    return 'unknown', parameters


def improved_background(imagename, method='sep', mask=None, apply_mask=True, 
                       show_map=False, box_size=None, filter_size=None, 
                       adapt_size=True, image_type=None):
    """
    Enhanced background estimation with multiple methods and adaptive sizing.
    
    Parameters
    ----------
    imagename : str
        Path to the image.
    method : str
        Method to use ('sep', 'photutils', 'iterative', 'wavelet').
    mask : array, optional
        Mask to be applied to the image.
    apply_mask : bool
        If True, calculate the mask from the image.
    show_map : bool
        If True, show the background map.
    box_size : int, tuple, or None
        Box size for background estimation. If None, determined automatically.
    filter_size : int, tuple, or None
        Filter size for smoothing the background. If None, determined automatically.
    adapt_size : bool
        If True, automatically adapt box and filter sizes to the image.
    image_type : str, optional
        Force image type ('hst', 'legacy', or None for auto-detection).
        
    Returns
    -------
    background : ndarray
        Estimated background.
    background_rms : ndarray
        Estimated background RMS.
    """
    # Read data
    import fitsio
    from astropy.io import fits
    
    # Try to read the header for image type detection
    header = None
    try:
        with fits.open(imagename) as hdul:
            header = hdul[0].header
    except:
        pass
    
    # Read the data
    data = fitsio.read(imagename)
    if len(data.shape) == 4:
        data = data[0][0]
    
    # Detect image type and get recommended parameters
    if image_type is None and adapt_size:
        detected_type, params = detect_image_type(header, data.shape)
        print(f"Detected image type: {detected_type}")
    else:
        if image_type == 'hst':
            detected_type, params = 'hst', {'min_boxes': 32, 'filter_fraction': 7, 'sigma': 2.0}
        elif image_type == 'legacy':
            detected_type, params = 'legacy', {'min_boxes': 16, 'filter_fraction': 5, 'sigma': 2.0}
        else:
            detected_type, params = 'unknown', {'min_boxes': 16, 'filter_fraction': 5, 'sigma': 2.0}
    
    # Handle NaNs
    data_no_nan = np.copy(data)
    nan_mask = np.isnan(data)
    if nan_mask.any():
        data_no_nan[nan_mask] = 0
    
    # Determine box size if not specified
    if box_size is None and adapt_size:
        box_size = adaptive_box_size(data.shape, min_boxes=params['min_boxes'])
        print(f"Using adaptive box size: {box_size}")
    elif isinstance(box_size, int):
        box_size = (box_size, box_size)
    
    # Determine filter size if not specified
    if filter_size is None and adapt_size:
        if isinstance(box_size, tuple):
            avg_box = sum(box_size) // 2
        else:
            avg_box = box_size
        
        filter_size = max(3, avg_box // params['filter_fraction'])
        if filter_size % 2 == 0:  # Make sure filter size is odd
            filter_size += 1
        print(f"Using filter size: {filter_size}")
    elif isinstance(filter_size, int):
        filter_size = (filter_size, filter_size)
    
    # Create or apply mask
    if mask is None and apply_mask:
        from photutils import detect_sources
        from astropy.stats import sigma_clipped_stats
        from scipy.ndimage import binary_dilation
        
        # Use parameters based on detected image type
        sigma = params['sigma']
        
        # Estimate background using sigma-clipping
        mean, median, std = sigma_clipped_stats(data_no_nan, sigma=sigma)
        threshold = median + (sigma * std)
        
        # Detect sources
        segm = detect_sources(data_no_nan, threshold, npixels=5)
        if segm is not None:
            # Create initial mask
            mask = segm.data > 0
            
            # Dilate mask
            iterations = 2
            if detected_type == 'hst':
                iterations = 3  # More dilation for HST due to PSF
            
            mask = binary_dilation(mask, iterations=iterations)
        else:
            mask = np.zeros(data.shape, dtype=bool)
    
    # Apply selected method
    if method == 'sep':
        import sep
        
        # Unpack box and filter sizes
        if isinstance(box_size, tuple):
            bw, bh = box_size
        else:
            bw = bh = box_size
            
        if isinstance(filter_size, tuple):
            fw, fh = filter_size
        else:
            fw = fh = filter_size
            
        # Ensure data is in proper format for SEP (float32)
        data_sep = np.ascontiguousarray(data_no_nan.astype(np.float32))
        
        # SEP background estimation
        bkg = sep.Background(data_sep, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
        background = bkg.back()
        background_rms = bkg.rms()
        
    elif method == 'photutils':
        from photutils import Background2D, MedianBackground
        from astropy.stats import SigmaClip
        
        sigma_clip = SigmaClip(sigma=params['sigma'])
        
        try:
            bkg = Background2D(data_no_nan, box_size, filter_size=filter_size,
                              sigma_clip=sigma_clip, bkg_estimator=MedianBackground(),
                              mask=mask)
            background = bkg.background
            background_rms = bkg.background_rms
        except Exception as e:
            print(f"Background2D failed: {e}. Falling back to simpler method.")
            from astropy.stats import sigma_clipped_stats
            
            # Fallback to global estimation
            mean, median, std = sigma_clipped_stats(data_no_nan, sigma=params['sigma'], mask=mask)
            background = np.ones_like(data_no_nan) * median
            background_rms = np.ones_like(data_no_nan) * std
        
    elif method == 'iterative':
        def iterative_background(data, mask=None, n_iterations=5, sigma=3.0):
            from astropy.stats import sigma_clipped_stats
            
            # Create working copy
            working_data = data.copy()
            
            # Apply initial mask if provided
            if mask is not None:
                working_mask = mask.copy()
            else:
                working_mask = np.zeros_like(data, dtype=bool)
            
            for i in range(n_iterations):
                # Calculate statistics with sigma clipping
                mean, median, std = sigma_clipped_stats(working_data, sigma=sigma, mask=working_mask)
                
                # Create mask for values significantly above background
                new_mask = working_data > (median + sigma * std)
                working_mask = working_mask | new_mask
                
                # Replace masked values with median for next iteration
                working_data[new_mask] = median
            
            # Final background estimate
            background = np.ones_like(data) * median
            background_rms = np.ones_like(data) * std
            
            return background, background_rms
        
        background, background_rms = iterative_background(data_no_nan, mask=mask, sigma=params['sigma'])
        
    elif method == 'wavelet':
        try:
            import pywt
            
            # Adjust level based on image size
            level = 4
            if data.shape[0] > 2000 or data.shape[1] > 2000:
                level = 5
            elif data.shape[0] < 500 or data.shape[1] < 500:
                level = 3
            
            # Handle image dimensions - must be divisible by 2^level
            pad_h = 0
            pad_w = 0
            if data.shape[0] % (2**level) != 0:
                pad_h = (2**level) - (data.shape[0] % (2**level))
            if data.shape[1] % (2**level) != 0:
                pad_w = (2**level) - (data.shape[1] % (2**level))
            
            if pad_h > 0 or pad_w > 0:
                padded_data = np.pad(data_no_nan, ((0, pad_h), (0, pad_w)), mode='reflect')
            else:
                padded_data = data_no_nan
            
            # Perform wavelet decomposition
            coeffs = pywt.wavedec2(padded_data, 'sym8', level=level)
            
            # Extract approximation coefficients (lowest frequency)
            cA = coeffs[0]
            
            # Reconstruct using only approximation coefficients
            new_coeffs = [cA] + [None] * level
            bg_padded = pywt.waverec2(new_coeffs, 'sym8')
            
            # Crop back to original size
            background = bg_padded[:data.shape[0], :data.shape[1]]
            
            # For RMS, use local standard deviation
            from astropy.stats import sigma_clipped_stats
            
            # Compute standard deviation of residuals
            mean, median, std = sigma_clipped_stats(data_no_nan - background, sigma=params['sigma'], mask=mask)
            background_rms = np.ones_like(background) * std
            
        except ImportError:
            print("PyWavelets not available. Falling back to iterative method.")
            bg_func = iterative_background if 'iterative_background' in locals() else lambda d, m: (np.median(d) * np.ones_like(d), np.std(d) * np.ones_like(d))
            background, background_rms = bg_func(data_no_nan, mask)
    
    # Restore NaN values
    if nan_mask.any():
        background[nan_mask] = np.nan
        background_rms[nan_mask] = np.nan
    
    # Visualization
    if show_map:
        import matplotlib.pyplot as plt
        
        fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        
        # Use arcsinh scaling for visualization
        vmin, vmax = np.nanpercentile(data, [0.5, 99.5])
        
        # Original data
        im1 = ax1.imshow(np.arcsinh(data), origin='lower', cmap='magma_r', 
                        vmin=np.arcsinh(vmin), vmax=np.arcsinh(vmax))
        plt.colorbar(im1, ax=ax1)
        ax1.set_title('Original Data')
        
        # Background
        im2 = ax2.imshow(np.arcsinh(background), origin='lower', cmap='magma_r')
        plt.colorbar(im2, ax=ax2)
        ax2.set_title(f'Background ({method})')
        
        # Background-subtracted
        im3 = ax3.imshow(np.arcsinh(data - background), origin='lower', cmap='magma_r',
                        vmin=np.arcsinh(vmin), vmax=np.arcsinh(vmax))
        plt.colorbar(im3, ax=ax3)
        ax3.set_title('Background Subtracted')
        
        plt.tight_layout()
        plt.show()
        
        # Also show the background mesh
        if method in ['sep', 'photutils']:
            from matplotlib.patches import Rectangle
            
            if isinstance(box_size, tuple):
                bw, bh = box_size
            else:
                bw = bh = box_size
                
            fig, ax = plt.subplots(figsize=(10, 8))
            im = ax.imshow(np.arcsinh(data), origin='lower', cmap='magma_r',
                          vmin=np.arcsinh(vmin), vmax=np.arcsinh(vmax))
            plt.colorbar(im, ax=ax)
            
            # # Draw background mesh grid
            # for i in range(0, data.shape[0], bh):
            #     for j in range(0, data.shape[1], bw):
            #         rect = Rectangle((j, i), bw, bh, fill=False, 
            #                         edgecolor='white', linewidth=0.5)
            #         ax.add_patch(rect)
                    
            ax.set_title(f'Background Estimation Grid (box size: {box_size})')
            plt.tight_layout()
            plt.show()
    
    return background, background_rms

"""
#Data download.
"""

from astroquery.mast import Observations
from astroquery.ipac.ned import Ned
from astropy.coordinates import SkyCoord
from astropy.nddata.utils import Cutout2D
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
import os
from difflib import get_close_matches

def make_hst_cutout(downloaded_files, coord, size_arcsec, output_filename=None, source_name=None, band=None, ref_coordinate=None):
    """
    Create a cutout from downloaded HST FITS files with complete header preservation.
    
    Parameters
    ----------
    downloaded_files : list
        List of paths to downloaded HST FITS files.
    coord : SkyCoord
        Center coordinates for the cutout (used as reference for offset calculation).
    size_arcsec : float
        Size of the cutout in arcseconds.
    output_filename : str, optional
        Name for the output file. If None, will generate based on source_name and band.
    source_name : str, optional
        Name of the source (used for automatic filename generation).
    band : str, optional
        Filter band (used for automatic filename generation).
    ref_coordinate : tuple of floats, optional
        Reference pixel coordinates (x, y) to use as the actual center of the cutout.
        If provided, the cutout will be centered on this pixel position relative to 
        the original image. This is useful for merger systems where the cataloged 
        coordinates point to one galaxy, but you want to center on a different 
        component. The function will compute the offset from coord to ref_coordinate
        and apply it to center the cutout appropriately.
        
    Returns
    -------
    str or None
        Path to the created cutout file, or None if failed.
    """
    # Find the best file to use for the cutout (science image with valid WCS)
    best_file = None
    best_ext = None
    
    for file_path in downloaded_files:
        try:
            with fits.open(file_path) as hdul:
                # Check if this is a standard HST drizzled file with SCI extension
                sci_ext = None
                for i, hdu in enumerate(hdul):
                    if isinstance(hdu, fits.ImageHDU) and hdu.name == 'SCI':
                        sci_ext = i
                        break
                
                # If no SCI extension found, look for any extension with valid data and WCS
                if sci_ext is None:
                    for i, hdu in enumerate(hdul):
                        if hdu.data is not None and isinstance(hdu, (fits.ImageHDU, fits.PrimaryHDU)):
                            try:
                                wcs = WCS(hdu.header)
                                if wcs.has_celestial:
                                    sci_ext = i
                                    break
                            except Exception:
                                continue
                
                if sci_ext is not None:
                    best_file = file_path
                    best_ext = sci_ext
                    print(f"[+>] Using science data from extension {sci_ext} ({hdul[sci_ext].name}) in {file_path}")
                    break
        except Exception as e:
            print(f"[Warning] Error opening {file_path}: {e}")
    
    if not best_file:
        print("[Error] No valid WCS found in any downloaded file.")
        return None
    
    # Now make the cutout with full header preservation
    try:
        with fits.open(best_file) as hdul:
            # Get primary header (we'll preserve this completely)
            primary_header = hdul[0].header.copy()
            
            # Get the science data, header, and WCS
            sci_data = hdul[best_ext].data
            sci_header = hdul[best_ext].header.copy()
            wcs = WCS(sci_header)
            
            # Calculate pixel scale
            if 'CD1_1' in sci_header:
                pixscale = abs(sci_header['CD1_1']) * 3600.0
            elif 'CDELT1' in sci_header:
                pixscale = abs(sci_header['CDELT1']) * 3600.0
            else:
                # Try to determine from the WCS
                pixscale = wcs.proj_plane_pixel_scales()[0].to(u.arcsec).value
                
            print(f"[+>] Pixel scale: {pixscale:.4f} arcsec/pixel")
            
            # Compute cutout size in pixels
            size_pixels = int(size_arcsec / pixscale)
            print(f"[+>] Cutout size: {size_pixels}x{size_pixels} pixels")
            
            # Determine the actual center position for the cutout
            if ref_coordinate is not None:
                # Convert coord (SkyCoord) to pixel coordinates in the original image
                coord_pixel_x, coord_pixel_y = wcs.world_to_pixel(coord)
                
                # Calculate the offset between coord and ref_coordinate (in pixels)
                offset_x = ref_coordinate[0] - coord_pixel_x
                offset_y = ref_coordinate[1] - coord_pixel_y
                
                # Create a new SkyCoord for the actual cutout center
                # We do this by converting ref_coordinate back to sky coordinates
                cutout_center = wcs.pixel_to_world(ref_coordinate[0], ref_coordinate[1])
                
                print(f"[+>] Original coord pixel position: ({coord_pixel_x:.2f}, {coord_pixel_y:.2f})")
                print(f"[+>] Reference coordinate: ({ref_coordinate[0]:.2f}, {ref_coordinate[1]:.2f})")
                print(f"[+>] Offset from coord: ({offset_x:.2f}, {offset_y:.2f}) pixels")
                print(f"[+>] Using reference coordinate as cutout center")
            else:
                # Use the original coord as the cutout center
                cutout_center = coord
            
            # Make the cutout using the determined center position
            cutout = Cutout2D(sci_data, position=cutout_center, size=(size_pixels, size_pixels), wcs=wcs)
            
            # Update the science header with the new WCS information
            sci_header.update(cutout.wcs.to_header())
            
            # Create a new HDUList with the original structure
            # Start with empty primary HDU
            cutout_hdul = fits.HDUList([fits.PrimaryHDU(header=primary_header)])
            
            # Add the science extension with cutout data
            if isinstance(hdul[best_ext], fits.PrimaryHDU):
                # Special case: if the science data is in the primary HDU, 
                # we need to replace the data in the primary HDU
                cutout_hdul[0].data = cutout.data
                # Update the WCS in the primary header
                cutout_hdul[0].header.update(cutout.wcs.to_header())
            else:
                # Otherwise create a new extension with the same name/ver as original
                sci_hdu = fits.ImageHDU(
                    data=cutout.data,
                    header=sci_header,
                    name=hdul[best_ext].name,
                    ver=hdul[best_ext].ver
                )
                cutout_hdul.append(sci_hdu)
            
            # Copy any other extensions (error arrays, weight maps, etc.)
            for i, hdu in enumerate(hdul):
                # Skip primary and science extensions (already handled)
                if i == 0 or i == best_ext:
                    continue
                
                # For image extensions that have matching WCS, also cut them out
                if isinstance(hdu, fits.ImageHDU) and hdu.data is not None:
                    try:
                        ext_wcs = WCS(hdu.header)
                        if ext_wcs.has_celestial:
                            ext_cutout = Cutout2D(hdu.data, position=cutout_center, 
                                                 size=(size_pixels, size_pixels), wcs=ext_wcs)
                            new_ext_header = hdu.header.copy()
                            new_ext_header.update(ext_cutout.wcs.to_header())
                            new_hdu = fits.ImageHDU(data=ext_cutout.data, header=new_ext_header, 
                                                    name=hdu.name, ver=hdu.ver)
                            cutout_hdul.append(new_hdu)
                            continue
                    except Exception as e:
                        print(f"[Warning] Could not cutout extension {i}: {e}")
                
                # For other extensions, just copy them as-is
                cutout_hdul.append(hdu.copy())
            
            # Generate output filename if not provided
            if output_filename is None:
                if source_name and band:
                    safe_name = source_name.replace(" ", "_")
                    output_filename = f"{safe_name}_{band}_{int(size_arcsec)}arcsec.fits"
                else:
                    # Extract source name from original filename if not provided
                    base_filename = os.path.basename(best_file)
                    output_filename = f"cutout_{base_filename}"
            
            # Write the cutout to disk
            cutout_hdul.writeto(output_filename, overwrite=True)
            print(f"[+>] Cutout saved as '{output_filename}' with {len(cutout_hdul)} extensions")
            return output_filename
            
    except Exception as e:
        print(f"[Error] Failed to create cutout: {e}")
        import traceback
        traceback.print_exc()
        return None


def _fuzzy_select(prompt_text, options, user_input):
    if user_input is None or user_input.upper() not in map(str.upper, options):
        print(f"[?] Available options: {sorted(set(options))}")
        match = get_close_matches(user_input or "", options, n=1, cutoff=0.3)
        if match:
            print(f"[~] Using closest match: {match[0]}")
            return match[0]
        else:
            selected = input(f"{prompt_text} (type one): ")
            return _fuzzy_select(prompt_text, options, selected)
    return user_input

def hst_cutout_mast(source_name, band=None, instrument=None, output_filename=None,
                   download_dir=None):
    """
    Download HST images from MAST and create a cutout with interactive band/instrument selection.

    Parameters
    ----------
    source_name : str
        Name of the astronomical object.
    size_arcsec : float
        Cutout size in arcseconds.
    band : str or None
        HST filter name (e.g., 'F160W'). If None or invalid, prompts user.
    instrument : str or None
        Instrument (e.g., 'WFC3', 'ACS', 'NICMOS'). If None or invalid, prompts user.
    output_filename : str or None
        Output FITS filename. If None, auto-generates.
    download_dir : str or None
        Directory to download files. If None, uses default structure.

    Returns
    -------
    tuple
        (output_filename, downloaded_files) or (None, []) on failure.
    """     
    
    # Step 1: Resolve coordinates
    try:
        result = Ned.query_object(source_name)
        ra = result["RA"][0]
        dec = result["DEC"][0]
        coord = SkyCoord(ra=ra, dec=dec, unit="deg")
    except Exception as e:
        print(f"[Error] Could not resolve source '{source_name}': {e}")
        return [], None, band, instrument

    # Step 2: Initial region query
    obs_table = Observations.query_region(coord, radius=0.02 * u.deg)
    obs_table = obs_table[(obs_table['obs_collection'] == 'HST') & (obs_table['dataproduct_type'] == 'image')]

    # Extract available filters and instruments
    available_filters = sorted(set(str(f) for f in obs_table['filters'] if f is not None))
    available_instruments = sorted(set(str(i) for i in obs_table['instrument_name'] if i is not None))


    if not available_filters or not available_instruments:
        print(f"[Error] No HST data available for {source_name}")
        return [], None, band, instrument
    
    # Step 3: Fuzzy matching or prompting
    band = _fuzzy_select("Select filter", available_filters, band)
    instrument = _fuzzy_select("Select instrument", available_instruments, instrument)

    # Step 4: Final filtering
    obs_table = obs_table[
        (obs_table['filters'] == band) &
        ([instrument.upper() in str(instr).upper() for instr in obs_table['instrument_name']])
    ]

    if len(obs_table) == 0:
        print(f"[Error] No matching HST observations found for filter '{band}' and instrument '{instrument}'")
        return [], None, band, instrument

    # Step 5: Download science files
    data_products = Observations.get_product_list(obs_table)
    filtered_products = Observations.filter_products(
        data_products,
        productSubGroupDescription=["DRZ", "DRC"],
        extension="fits",
        mrp_only=False
    )

    
    if len(filtered_products) == 0:
        print(f"[Error] No DRZ/DRC products found.")
        return [], None, band, instrument

    if download_dir is None:
        download_dir = f"hst/{source_name}/{band}"
    
    manifest = Observations.download_products(filtered_products, 
                                             download_dir=download_dir+f"/{source_name}/{band}", 
                                             mrp_only=False)
    downloaded_files = [f for f in manifest['Local Path'] if f is not None]
    downloaded_files = np.unique(downloaded_files)
    # downloaded_files = list(set(downloaded_files))  # Remove duplicates
    print(f"[+>] Downloaded {len(downloaded_files)} FITS files.")

    return downloaded_files, coord, band, instrument


import requests
def get_source_coordinates(source_name):
    """
    Get the coordinates of a source by name (NED / SIMBAD / Sesame).
    
    Parameters
    ----------
    source_name : str
        Source name.
        
    Returns
    -------
    SkyCoord object
        The coordinates of the source.
    """
    # cosmo.resolve_source_coordinates goes through NED ObjectLookup / SIMBAD
    # TAP / Sesame (~0.5 s); Ned.query_object below is the legacy CGI path and
    # is only reached if all of those fail.
    try:
        coords = resolve_source_coordinates(source_name)
        if coords is not None:
            return coords
    except Exception as e:
        print(f"Fast resolver failed for {source_name}: {e}")

    try:
        result_table = Ned.query_object(source_name)
        ra = result_table['RA'][0]
        dec = result_table['DEC'][0]
        coords = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        print(coords)
        return coords
    except Exception as e:
        print(f"Error getting coordinates for {source_name}: {e}")
        return None

def arcsec_to_pixels(size_arcsec, 
                     pixel_scale = 0.262,  # arcsec/pixel
                     redshift=None):
    """
    Convert arcseconds to pixels for Legacy Survey DR9.
    
    Parameters
    ----------
    size_arcsec : float
        Size in arcseconds.
    redshift : float, optional
        Redshift of the source. Not used for pixel conversion since
        Legacy Survey has fixed pixel scale.
        
    Returns
    -------
    int
        Size in pixels.
    """
    # Legacy Survey DR9 pixel scale is 0.262 arcsec/pixel
    
    size_pixels = int(np.ceil(size_arcsec / pixel_scale))
    return size_pixels

def legacy_survey_cutout(source_name, 
                         size_arcsec, 
                         band, 
                         pixel_scale = 0.262,
                         output_filename=None,
                         output_path=None):
    """
    Get a cutout from the Legacy Survey DR9.
    
    Parameters
    ----------
    source_name : str
        Source name.
    size_arcsec : float
        Size of the cutout in arcseconds.
    band : str
        Band to retrieve (g, r, z, etc).
    output_filename : str, optional
        Name of the output FITS file.
        
    Returns
    -------
    str
        Name of the saved file if successful, None otherwise.
    """
    # Get source coordinates
    coords = get_source_coordinates(source_name)
    if coords is None:
        print(f"Could not find coordinates for {source_name}")
        return None
    
    # Get redshift (not directly used for pixel conversion, but saved for reference)
    try:
        redshift = find_z_NED(source_name)
    except Exception as e:
        print(f"Error getting redshift for {source_name}: {e}")
        redshift = None
    
    # Calculate size in pixels
    size_pixels = arcsec_to_pixels(size_arcsec,pixel_scale=pixel_scale)
    
    # Create output filename if not provided
    if output_filename is None:
        if output_path is not None:
            if not os.path.exists(output_path):
                os.makedirs(output_path)
            # Use the source name and band to create a unique filename  
            output_filename = f"{source_name}_{band}_{size_arcsec}arcsec.fits"
            output_filename = os.path.join(output_path, output_filename)
            print(output_filename)
        else:
            output_filename = f"{source_name}_{band}_{size_arcsec}arcsec.fits"
            print(output_filename)        
    
    # Construct the URL for the cutout
    ra = coords.ra.deg
    dec = coords.dec.deg
    
    url = (f"https://www.legacysurvey.org/viewer/fits-cutout?ra={ra}&dec={dec}"
           f"&layer=ls-dr9&pixscale={pixel_scale}&bands={band}&size={size_pixels}")
    
    try:
        # Download the cutout
        response = requests.get(url)
        if response.status_code != 200:
            print(f"Error downloading cutout: HTTP status {response.status_code}")
            return None
        
        # Save the cutout
        with open(output_filename, 'wb') as f:
            f.write(response.content)
        
        print(f"Successfully saved cutout to {output_filename}")
        return output_filename
    
    except Exception as e:
        print(f"Error downloading or saving cutout: {e}")
        return None


"""
PanSTARRS cutout function with corrected WCS handling
"""

import requests
import numpy as np
import os
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astroquery.ipac.ned import Ned
import astropy.units as u


def get_source_coordinates(source_name):
    """
    Get the coordinates of a source using NED.
    
    Parameters
    ----------
    source_name : str
        Source name.
        
    Returns
    -------
    SkyCoord object
        The coordinates of the source.
    """
    try:
        result_table = Ned.query_object(source_name)
        ra = result_table['RA'][0]
        dec = result_table['DEC'][0]
        coords = SkyCoord(ra, dec, unit=(u.deg, u.deg))
        print(f"[+>] Found coordinates for {source_name}: RA={ra:.5f}, Dec={dec:.5f}")
        return coords
    except Exception as e:
        print(f"[Error] Could not get coordinates for {source_name}: {e}")
        return None

def panstarrs_cutout(source_name, 
                     size_arcsec, 
                     band='i',
                     pixel_scale=0.25,  # PanSTARRS native pixel scale in arcsec/pixel
                     output_filename=None,
                     output_path=None,
                     data_release='dr2'):
    """
    Get a cutout from PanSTARRS (PS1) using their cutout service.
    Saves as a simple FITS file with corrected WCS information.
    
    Parameters
    ----------
    source_name : str
        Source name (will be resolved using NED).
    size_arcsec : float
        Size of the cutout in arcseconds.
    band : str
        Band to retrieve ('g', 'r', 'i', 'z', 'y'). Default is 'i'.
    pixel_scale : float
        Pixel scale in arcsec/pixel. Default is 0.25 (PanSTARRS native).
    output_filename : str, optional
        Name of the output FITS file. If None, auto-generates based on source/band.
    output_path : str, optional
        Directory path for output file. Created if doesn't exist.
    data_release : str
        PanSTARRS data release ('dr1' or 'dr2'). Default is 'dr2'.
        
    Returns
    -------
    str or None
        Path to the saved cutout file if successful, None otherwise.
    """
    
    # Validate band
    valid_bands = ['g', 'r', 'i', 'z', 'y']
    if band not in valid_bands:
        print(f"[Error] Invalid band '{band}'. Must be one of: {valid_bands}")
        return None
    
    # Get source coordinates using existing function
    coords = get_source_coordinates(source_name)
    if coords is None:
        print(f"[Error] Could not find coordinates for {source_name}")
        return None
    
    # Calculate size in pixels
    size_pixels = int(np.ceil(size_arcsec / pixel_scale))
    print(f"[+>] Cutout size: {size_pixels}x{size_pixels} pixels ({size_arcsec} arcsec)")
    
    # Create output filename if not provided
    if output_filename is None:
        # Clean source name for filename
        safe_name = source_name.replace(" ", "_").replace("/", "_")
        output_filename = f"{safe_name}_ps1_{band}_{size_arcsec}arcsec.fits"
        
        if output_path is not None:
            if not os.path.exists(output_path):
                os.makedirs(output_path)
                print(f"[+>] Created output directory: {output_path}")
            output_filename = os.path.join(output_path, output_filename)
    
    # Extract RA/Dec
    ra = coords.ra.deg
    dec = coords.dec.deg
    
    # Create temporary filename for download
    temp_filename = output_filename + '.temp'
    
    # Try primary endpoint first
    print(f"[+>] Requesting PanSTARRS {data_release.upper()} {band}-band cutout...")
    
    # Primary cutout service URL
    base_url = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi"
    
    # Build query parameters - note the size parameter is diameter in pixels
    params = {
        'ra': ra,
        'dec': dec,
        'size': size_pixels,  # This is the diameter of the cutout
        'format': 'fits',
        'filters': band,
        'output_size': size_pixels  # Ensure output matches requested size
    }
    
    # For DR2, specify stack type
    if data_release == 'dr2':
        params['type'] = 'stack'
    
    try:
        # Make the request
        response = requests.get(base_url, params=params, timeout=30)
        
        # If primary endpoint fails, try alternative approach
        if response.status_code != 200:
            print(f"[Warning] Primary endpoint returned HTTP {response.status_code}")
            print("[!] Trying alternative PanSTARRS endpoint...")
            
            # Alternative approach: Get the direct image URL first
            alt_url = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
            alt_params = {
                'ra': ra,
                'dec': dec,
                'size': size_pixels,
                'format': 'fits',
                'filters': band,
                'type': 'stack'
            }
            
            # Get the filename/URL info
            filename_response = requests.get(alt_url, params=alt_params, timeout=30)
            
            if filename_response.status_code == 200:
                # Parse the response to get the actual image URL
                lines = filename_response.text.strip().split('\n')
                if len(lines) > 1:  # Skip header line
                    # The response format includes multiple fields
                    fields = lines[1].split()
                    if len(fields) >= 8:
                        # Extract the URL (usually the last field)
                        image_url = fields[7]
                        
                        # Build complete cutout URL with size parameters
                        if 'rings.v3.skycell' in image_url:
                            # Construct fitscut URL from the filename
                            base_filename = image_url.split('/')[-1]
                            cutout_url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
                                        f"red={base_filename}&ra={ra}&dec={dec}&size={size_pixels}&"
                                        f"output_size={size_pixels}&format=fits")
                        else:
                            # Direct cutout URL
                            cutout_url = (f"https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
                                        f"ra={ra}&dec={dec}&size={size_pixels}&format=fits&"
                                        f"filters={band}&output_size={size_pixels}")
                        
                        print(f"[+>] Found image, requesting cutout...")
                        response = requests.get(cutout_url, timeout=30)
                        
                        if response.status_code != 200:
                            print(f"[Error] Alternative endpoint also failed: HTTP {response.status_code}")
                            return None
                    else:
                        print(f"[Error] Could not parse filename service response")
                        return None
                else:
                    print(f"[Error] No results from filename service")
                    return None
            else:
                print(f"[Error] Filename service failed: HTTP {filename_response.status_code}")
                return None
        
        # Check if we got valid FITS data
        if len(response.content) < 1000:  # Basic size check
            print(f"[Error] Response too small ({len(response.content)} bytes), likely not valid FITS data")
            return None
        
        # Save the temporary file
        with open(temp_filename, 'wb') as f:
            f.write(response.content)
        
        # Read the downloaded file and correct the WCS
        with fits.open(temp_filename) as hdul:
            if len(hdul) == 0:
                print(f"[Error] Invalid FITS file received")
                os.remove(temp_filename)
                return None
            
            # Find the HDU with actual image data
            image_data = None
            image_header = None
            
            # Check all HDUs for image data
            for i, hdu in enumerate(hdul):
                if hdu.data is not None and len(hdu.data.shape) == 2:
                    print(f"[+>] Found image data in HDU {i} (shape: {hdu.data.shape})")
                    image_data = hdu.data
                    image_header = hdu.header.copy()
                    break
            
            # If no image data found in extensions, check primary
            if image_data is None and hdul[0].data is not None:
                image_data = hdul[0].data
                image_header = hdul[0].header.copy()
                print(f"[+>] Found image data in primary HDU (shape: {image_data.shape})")
            
            if image_data is None:
                print(f"[Error] No image data found in FITS file")
                os.remove(temp_filename)
                return None
            
            # CRITICAL FIX: Always use Cutout2D to ensure correct WCS
            # This fixes the issue where PanSTARRS returns cutouts with incorrect CRPIX values
            print(f"[+>] Correcting WCS information using Cutout2D...")
            from astropy.nddata import Cutout2D
            
            # Get WCS from the header
            try:
                wcs = WCS(image_header)
            except Exception as e:
                print(f"[Error] Could not parse WCS from header: {e}")
                os.remove(temp_filename)
                return None
            
            # Use Cutout2D centered on our requested coordinates
            # This ensures CRPIX is correctly set relative to the cutout
            try:
                cutout = Cutout2D(image_data, position=coords, 
                                size=(size_pixels, size_pixels), 
                                wcs=wcs, mode='partial', fill_value=0.0)
                
                image_data = cutout.data
                corrected_wcs = cutout.wcs
                print(f"[+>] WCS corrected, final cutout size: {image_data.shape}")
                
            except Exception as e:
                print(f"[Warning] Cutout2D failed: {e}")
                print(f"[!] Attempting manual WCS correction...")
                
                # Fallback: manually correct the WCS
                # Calculate where the requested coordinates fall in pixel space
                try:
                    x_pixel, y_pixel = wcs.world_to_pixel(coords)
                    
                    # Calculate new CRPIX values
                    # CRPIX should point to where CRVAL is located in the cutout
                    center_x = image_data.shape[1] / 2.0
                    center_y = image_data.shape[0] / 2.0
                    
                    # Offset from actual coordinate position to image center
                    dx = center_x - x_pixel
                    dy = center_y - y_pixel
                    
                    # Update CRPIX
                    if 'CRPIX1' in image_header:
                        image_header['CRPIX1'] = image_header['CRPIX1'] + dx
                    else:
                        image_header['CRPIX1'] = center_x
                        
                    if 'CRPIX2' in image_header:
                        image_header['CRPIX2'] = image_header['CRPIX2'] + dy
                    else:
                        image_header['CRPIX2'] = center_y
                    
                    # Update CRVAL to point to our requested coordinates
                    image_header['CRVAL1'] = ra
                    image_header['CRVAL2'] = dec
                    
                    corrected_wcs = WCS(image_header)
                    print(f"[+>] Manual WCS correction applied")
                    
                except Exception as e2:
                    print(f"[Error] Manual WCS correction failed: {e2}")
                    os.remove(temp_filename)
                    return None
            
            # Create a new simple FITS file with corrected WCS
            new_hdu = fits.PrimaryHDU(data=image_data)
            
            # Start with the corrected WCS header
            wcs_header = corrected_wcs.to_header()
            for key in wcs_header:
                new_hdu.header[key] = wcs_header[key]
            
            # Add metadata
            new_hdu.header['OBJECT'] = source_name
            new_hdu.header['BAND'] = (band, 'Filter band')
            new_hdu.header['FILTER'] = (band, 'Filter band')
            new_hdu.header['PIXSCALE'] = (pixel_scale, 'Pixel scale in arcsec/pixel')
            new_hdu.header['CUTSIZE'] = (size_arcsec, 'Requested cutout size in arcsec')
            new_hdu.header['DATASRC'] = f'PanSTARRS {data_release.upper()}'
            
            # Add target coordinates (the coordinates we requested)
            new_hdu.header['RA'] = (ra, 'Right ascension in degrees')
            new_hdu.header['DEC'] = (dec, 'Declination in degrees')
            new_hdu.header['EPOCH'] = (2000.0, 'Epoch of coordinates')
            new_hdu.header['EQUINOX'] = (2000.0, 'Equinox of coordinates')
            
            # Ensure we have NAXIS keywords with correct values
            new_hdu.header['NAXIS'] = 2
            new_hdu.header['NAXIS1'] = image_data.shape[1]
            new_hdu.header['NAXIS2'] = image_data.shape[0]
            
            # Add units and other useful keywords
            if 'BUNIT' in image_header:
                new_hdu.header['BUNIT'] = image_header['BUNIT']
            else:
                new_hdu.header['BUNIT'] = ('nanomaggies', 'Pixel units')
            
            new_hdu.header['ORIGIN'] = 'STScI/PanSTARRS'
            new_hdu.header['TELESCOP'] = 'PanSTARRS'
            new_hdu.header['INSTRUME'] = 'GPC1'
            
            # Try to add redshift from NED
            try:
                result_table = Ned.query_object(source_name)
                if 'Redshift' in result_table.colnames:
                    redshift = result_table['Redshift'][0]
                    if not np.isnan(redshift):
                        new_hdu.header['REDSHIFT'] = (redshift, 'Source redshift from NED')
            except:
                pass
            
            # Add history
            new_hdu.header.add_history(f'PanSTARRS {data_release.upper()} cutout')
            new_hdu.header.add_history(f'Created using ps1images.stsci.edu cutout service')
            new_hdu.header.add_history(f'WCS corrected using astropy Cutout2D')
            new_hdu.header.add_history(f'Filter: {band}, Size: {size_arcsec} arcsec')
            new_hdu.header.add_history(f'Target: RA={ra:.5f}, Dec={dec:.5f}')
            
            # Save the simple FITS file
            new_hdu.writeto(output_filename, overwrite=True)
            
            # Clean up temporary file
            os.remove(temp_filename)
            
            print(f"[+>] Successfully saved PanSTARRS cutout to {output_filename}")
            print(f"[+>] Final image size: {image_data.shape[1]}x{image_data.shape[0]} pixels")
            print(f"[+>] WCS reference pixel (CRPIX): ({new_hdu.header['CRPIX1']:.1f}, {new_hdu.header['CRPIX2']:.1f})")
            print(f"[+>] WCS reference coord (CRVAL): ({new_hdu.header['CRVAL1']:.5f}, {new_hdu.header['CRVAL2']:.5f})")
            
        return output_filename
    
    except requests.Timeout:
        print(f"[Error] Request timed out after 30 seconds")
        if os.path.exists(temp_filename):
            os.remove(temp_filename)
        return None
    except requests.RequestException as e:
        print(f"[Error] Request failed: {e}")
        if os.path.exists(temp_filename):
            os.remove(temp_filename)
        return None
    except Exception as e:
        print(f"[Error] Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        if os.path.exists(temp_filename):
            os.remove(temp_filename)
        if os.path.exists(output_filename):
            os.remove(output_filename)
        return None

#!/usr/bin/env python
"""
Enhanced overlay function that handles WCS alignment between optical and radio images.
This addresses common issues when overlaying PanSTARRS and radio data.
"""

# def overlay_radio_optical_enhanced
def overlay_radio_optical_enhanced(hst_cutout_filename, radio_filename, output_filename=None, 
                         rms_rad=None,
                         figsize=(6, 6),
                         cutout_size=(1024, 1024), 
                         optical_stretch='asinh',
                         vmin_factor=0.1, 
                         vmax_factor=0.1, 
                         vmin_opt=3.0, 
                         optical_contour_color='black', 
                         radio_color='#EE7733', 
                         optical_cmap='Greys',
                         title=None):
    """
    Create an overlay of radio contours on optical image.
    Works with HST, Legacy Survey, and PanSTARRS data.
    
    Parameters
    ----------
    hst_cutout_filename : str
        Path to the optical FITS file (HST, Legacy Survey, or PanSTARRS)
    radio_filename : str
        Path to the radio FITS file
    output_filename : str, optional
        Path where to save the figure. If None, will not save.
    rms_rad : float, optional
        Radio RMS noise level. If None, will be calculated from data.
    figsize : tuple, optional
        Size of the figure (width, height) in inches
    cutout_size : tuple, optional
        Size of the cutout in pixels (y, x)
    optical_stretch : str, optional
        Stretch function for optical image ('linear', 'log', 'sqrt', 'asinh')
    vmin_factor : float, optional
        Factor to multiply optical std for vmin
    vmax_factor : float, optional
        Factor to multiply optical max for vmax
    vmin_opt : float, optional
        Sigma level for optical contour minimum
    optical_contour_color : str, optional
        Color for optical contours
    radio_color : str, optional
        Color for radio contours
    optical_cmap : str, optional
        Colormap for optical image
    title : str, optional
        Title for the plot. If None, no title is added.
    
    Returns
    -------
    fig, ax : matplotlib Figure and Axes objects
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.visualization import simple_norm
    from astropy.stats import mad_std
    from reproject import reproject_interp
    from astropy.nddata import Cutout2D
    import warnings
    warnings.filterwarnings('ignore')

    # Load the optical FITS image
    with fits.open(hst_cutout_filename) as optical_fits:
        optical_data = optical_fits[0].data
        optical_header = optical_fits[0].header
        optical_wcs = WCS(optical_header)

    # Load the radio data
    with fits.open(radio_filename) as radio_fits:
        radio_data = radio_fits[0].data
        radio_header = radio_fits[0].header
        radio_wcs = WCS(radio_header)
        
        # Extract the 2D slice from potentially higher-dimensional radio data
        # and properly drop the extra WCS axes
        if len(radio_data.shape) == 4:
            radio_data_2d = radio_data[0, 0, :, :]
            radio_wcs_2d = radio_wcs.dropaxis(3).dropaxis(2)
        elif len(radio_data.shape) == 3:
            radio_data_2d = radio_data[0, :, :]
            radio_wcs_2d = radio_wcs.dropaxis(2)
        else:
            radio_data_2d = radio_data
            radio_wcs_2d = radio_wcs

    # Determine the center coordinates for the cutout
    # First priority: use RA and DEC from header if available
    if 'RA' in optical_header and 'DEC' in optical_header:
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        cutout_center = SkyCoord(optical_header['RA']*u.deg, optical_header['DEC']*u.deg)
        print(f"Using header coordinates: RA={optical_header['RA']:.5f}, DEC={optical_header['DEC']:.5f}")
    elif 'CRVAL1' in optical_header and 'CRVAL2' in optical_header:
        # Fall back to WCS reference coordinates
        from astropy.coordinates import SkyCoord
        import astropy.units as u
        cutout_center = SkyCoord(optical_header['CRVAL1']*u.deg, optical_header['CRVAL2']*u.deg)
        print(f"Using CRVAL coordinates: RA={optical_header['CRVAL1']:.5f}, DEC={optical_header['CRVAL2']:.5f}")
    else:
        # Last resort: use the geometric center of the optical image
        center_x = optical_data.shape[1] / 2.0
        center_y = optical_data.shape[0] / 2.0
        cutout_center = optical_wcs.pixel_to_world(center_x, center_y)
        print(f"Using geometric center of image at pixel ({center_x:.1f}, {center_y:.1f})")
        print(f"Converted to sky coordinates: RA={cutout_center.ra.deg:.5f}, DEC={cutout_center.dec.deg:.5f}")

    # Create a cutout of the optical data
    optical_cutout = Cutout2D(optical_data, cutout_center, cutout_size, wcs=optical_wcs)
    
    # Reproject radio data to match optical cutout
    radio_cutout_data, footprint = reproject_interp(
        (radio_data_2d, radio_wcs_2d), 
        optical_cutout.wcs, 
        shape_out=optical_cutout.shape
    )

    # Handle NaN values in radio data
    radio_mask = ~np.isnan(radio_cutout_data)
    radio_cutout_data[~radio_mask] = 0

    # Calculate statistics for radio contour levels
    radio_peak = np.nanmax(radio_cutout_data)
    
    if rms_rad is not None:
        radio_std = rms_rad
    else:
        radio_std = mad_std(radio_cutout_data[radio_mask], ignore_nan=True)
    
    if radio_std == 0 or np.isnan(radio_std):
        radio_std = np.nanstd(radio_cutout_data)
    
    print(f"Radio peak: {radio_peak}, Radio std: {radio_std}")
    
    # Set contour levels
    contour_levels = np.geomspace(5.0 * radio_std, 2.0 * radio_peak, 6)

    # Create a figure and axes with WCS projection for the cutout
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=optical_cutout.wcs)

    # Calculate vmin and vmax for optical image
    vmin = vmin_factor * mad_std(optical_cutout.data[optical_cutout.data > 0])
    vmax = np.nanmax(optical_cutout.data)

    # Apply the appropriate stretch to the optical data
    norm = simple_norm(optical_cutout.data, stretch=optical_stretch, 
                      asinh_a=0.015, vmin=vmin, vmax=vmax * vmax_factor)

    # Plot the optical cutout image
    ax.imshow(optical_cutout.data, cmap=optical_cmap, origin='lower', norm=norm)

    # Add optical contours
    opt_peak = np.nanmax(optical_cutout.data)
    opt_std = mad_std(optical_cutout.data[optical_cutout.data > 0])
    
    contour_levels_opt = np.geomspace(1.0 * opt_peak, vmin_opt * opt_std, 6)
    print(contour_levels)
    print(contour_levels_opt[::-1])
    
    ax.contour(optical_cutout.data, levels=contour_levels_opt[::-1], 
              colors=optical_contour_color, linewidths=1.5, alpha=0.5)

    # Overlay radio contours on the optical cutout
    ax.contour(radio_cutout_data, levels=contour_levels, 
              colors=radio_color, linewidths=1.5, alpha=0.7)

    # Set axis labels
    ax.set_xlabel('RA (J2000)', fontsize=12)
    ax.set_ylabel('Dec (J2000)', fontsize=12)
    
    # Add title if provided
    if title:
        ax.text(0.05, 0.95, title, transform=ax.transAxes, 
               color='white', fontsize=14, weight='bold', 
               ha='left', va='top', bbox=dict(facecolor='black', alpha=0.7))

    # Add grid
    ax.grid(color='black', linestyle='--', linewidth=0.5, alpha=0.3)
    
    # Tight layout
    # plt.tight_layout()
    
    # Save the figure if an output filename is provided
    if output_filename:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved as {output_filename}")
        
    return fig, ax


# def find_sci_extension(hdul):
#     """Find the science extension with valid data in a FITS file."""
#     # First, try to find SCI extension which is standard for HST
#     for i, hdu in enumerate(hdul):
#         if hasattr(hdu, 'name') and hdu.name == 'SCI' and hdu.data is not None:
#             try:
#                 wcs = WCS(hdu.header, naxis=2)
#                 if wcs.has_celestial:
#                     return i, hdu.data, wcs
#             except Exception:
#                 pass
    
#     # If no SCI extension found, try any extension with data and valid WCS
#     for i, hdu in enumerate(hdul):
#         if hdu.data is not None:
#             try:
#                 wcs = WCS(hdu.header, naxis=2)
#                 if wcs.has_celestial:
#                     return i, hdu.data, wcs
#             except Exception:
#                 pass
    
#     # Fall back to primary HDU
#     try:
#         wcs = WCS(hdul[0].header, naxis=2)
#         return 0, hdul[0].data, wcs
#     except Exception:
#         return 0, hdul[0].data, None



def overlay_radio_optical(hst_cutout_filename, radio_filename, output_filename=None, 
                         rms_rad = None,
                         rms_opt = None,
                         figsize=(6, 6),
                         cutout_size=(1024, 1024), optical_stretch='asinh',
                         vmin_factor=0.1, vmax_factor=0.1, vmin_opt=3.0, opt_c_levels=6,
                         #optical_contour_color='grey', radio_color='limegreen',optical_cmap='magma_r',
                         optical_contour_color='black', radio_color='#EE7733', optical_cmap='Greys',
                         title=None, title_position='top'):
    """
    Create an overlay of radio contours on HST optical image.
    
    Parameters
    ----------
    hst_cutout_filename : str
        Path to the HST FITS file
    radio_filename : str
        Path to the radio FITS file
    output_filename : str, optional
        Path where to save the figure. If None, will not save.
    cutout_size : tuple, optional
        Size of the cutout in pixels (y, x)
    optical_stretch : str, optional
        Stretch function for optical image ('linear', 'log', 'sqrt', 'asinh')
    radio_color : str, optional
        Color for radio contours
    optical_cmap : str, optional
        Colormap for optical image
    optical_contour_color : str, optional
        Color for optical contours
    title : str, optional
        Title for the plot. If None, no title is added.
    
    Returns
    -------
    fig, ax : matplotlib Figure and Axes objects
    """
    import numpy as np
    import matplotlib.pyplot as plt
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.visualization import ZScaleInterval, simple_norm
    from reproject import reproject_interp
    from astropy.nddata import Cutout2D
    from astropy.stats import mad_std
    # Load the optical FITS image (HST)
    with fits.open(hst_cutout_filename) as optical_fits:
        # Find the science extension in the HST data
        sci_ext, optical_data, optical_wcs = find_sci_extension(optical_fits)
        
        if optical_data is None:
            raise ValueError("Could not find valid data in the HST file")
        
        # If we found a WCS but not in the primary header, we need to merge it with primary metadata
        primary_header = optical_fits[0].header.copy()
        if sci_ext > 0:
            # We found data in a SCI extension, get all WCS keywords
            wcs_header = optical_wcs.to_header()
            # Update primary header with WCS info
            for key in wcs_header:
                primary_header[key] = wcs_header[key]
        
        optical_header = primary_header

    # Load the radio contour data
    with fits.open(radio_filename) as radio_fits:
        radio_data = radio_fits[0].data
        radio_header = radio_fits[0].header
        
        # Extract the 2D slice from potentially higher-dimensional radio data
        if len(radio_data.shape) == 4:  # [polarization, frequency, y, x]
            radio_data_2d = radio_data[0, 0, :, :]
        elif len(radio_data.shape) == 3:  # [frequency, y, x] or [polarization, y, x]
            radio_data_2d = radio_data[0, :, :]
        else:
            radio_data_2d = radio_data
        
        # Create a WCS object for the 2D radio data
        radio_wcs_2d = WCS(radio_header, naxis=2)

    # Determine the center of the cutout (use optical peak if not provided)
    optical_peak_coords = np.unravel_index(np.argmax(optical_data), optical_data.shape)
    cutout_center = optical_wcs.pixel_to_world(optical_peak_coords[1], optical_peak_coords[0])

    # Create a cutout of both optical and radio data
    optical_cutout = Cutout2D(optical_data, cutout_center, cutout_size, wcs=optical_wcs)
    
    # Reproject radio data to match optical cutout
    radio_cutout_data, footprint = reproject_interp(
        (radio_data_2d, radio_wcs_2d), 
        optical_cutout.wcs, 
        shape_out=optical_cutout.shape
    )

    # Handle NaN values in radio data
    radio_mask = ~np.isnan(radio_cutout_data)
    radio_cutout_data[~radio_mask] = 0
    # plt.figure()
    # plt.imshow(radio_cutout_data)
    # plt.show()
    # Calculate statistics for contour levels
    radio_peak = np.nanmax(radio_cutout_data)    
    if rms_rad is not None:
        radio_std = rms_rad
    else:
        radio_std = mad_std(radio_cutout_data[radio_mask],ignore_nan=True)
    if radio_std == 0 or np.isnan(radio_std):
        radio_std = np.nanstd(radio_cutout_data)
    print(f"Radio peak: {radio_peak}, Radio std: {radio_std}")
    # Set contour levels (ensure they don't include zero)
    contour_levels = np.geomspace(5.0 * radio_std, 2.0 * radio_peak, 6)
    
    # Create a figure and axes with WCS projection for the cutout
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection=optical_cutout.wcs)

    # Plot the optical cutout image
    if rms_opt is not None:
        vmin = vmin_factor * rms_opt
    else:
        vmin = vmin_factor * mad_std(optical_cutout.data[optical_cutout.data > 0])
    vmax = np.nanmax(optical_cutout.data)

    # Apply the appropriate stretch to the optical data
    norm = simple_norm(optical_cutout.data, stretch=optical_stretch, 
                      asinh_a=0.015, vmin=vmin, vmax=vmax * vmax_factor)

    ax.imshow(optical_cutout.data, cmap=optical_cmap, origin='lower', norm=norm)

    

    
    # try:
    # Add optical contours if desired
    opt_peak = np.nanmax(optical_cutout.data)
    if rms_opt is not None:
        opt_std = rms_opt
    else:
        opt_std = mad_std(optical_cutout.data[optical_cutout.data > 0])
    # opt_std = np.nanstd(optical_cutout.data[optical_cutout.data > 0])
    contour_levels_opt = np.geomspace(1.0 * opt_peak, vmin_opt * opt_std, opt_c_levels)
    print(contour_levels)
    print(contour_levels_opt[::-1])
    ax.contour(optical_cutout.data, levels=contour_levels_opt[::-1], 
            colors=optical_contour_color, linewidths=1.5, alpha=0.5,
            #   linestyles='dashdot'
            )
    # except Exception as e:
    #     print(f"Could not plot optical contours: {e}")

    # Overlay radio contours on the optical cutout
    ax.contour(radio_cutout_data, levels=contour_levels, 
              colors=radio_color, linewidths=1.5, alpha=0.7)

    # Set axis labels
    ax.set_xlabel('RA (J2000)', fontsize=12)
    ax.set_ylabel('Dec (J2000)', fontsize=12)
    
    # Add title if provided
    if title:
        if title_position == 'top':
            ax.text(0.05, 0.95, title, transform=ax.transAxes, 
               color='white', fontsize=14, weight='bold', 
               ha='left', va='top', bbox=dict(facecolor='black', alpha=0.7))
        elif title_position == 'bottom':
            ax.text(0.05, 0.08, title, transform=ax.transAxes, 
                   color='white', fontsize=14, weight='bold', 
                   ha='left', va='top', bbox=dict(facecolor='black', alpha=0.7))
        else:
            ax.text(0.05, 0.08, title, transform=ax.transAxes, 
                   color='white', fontsize=14, weight='bold', 
                   ha='left', va='top', bbox=dict(facecolor='black', alpha=0.7))

    # Add grid
    # ax.grid(color='black', linestyle='--', linewidth=0.5, alpha=0.3)
    
    # Tight layout
    # plt.tight_layout()
    
    # Save the figure if an output filename is provided
    if output_filename:
        plt.savefig(output_filename, dpi=300, bbox_inches='tight')
        print(f"Figure saved as {output_filename}")
        
    return fig, ax






"""
#Cosmo
"""
import numpy as np
from scipy.integrate import quad
from scipy.special import gamma

def stable_one_minus_exp_minus_tau_over_tau(tau):
    tau = np.array(tau)
    result = np.ones_like(tau)
    small = np.abs(tau) < 1e-3
    result[~small] = (1 - np.exp(-tau[~small])) / tau[~small]
    # Use Taylor expansion when tau is small
    result[small] = 1 - tau[small]/2 + tau[small]**2/6
    return result

def tau_nu(nu, nu_t):
    return (nu / nu_t) ** (-2.1)

# synchrotron integral function
def synchrotron_integrand(nu, nu0, nu_t, alpha_sy, S_nu0_sy):
    return S_nu0_sy * (nu / nu0) ** alpha_sy * np.exp(-tau_nu(nu, nu_t))

# free-free integral function
def freefree_integrand(nu, nu0, nu_t, S_nu0_ff):
    tau = tau_nu(nu, nu_t)
    return S_nu0_ff * ((1 - np.exp(-tau)) / tau) * (nu / nu0) ** -0.1

# def freefree_integrand(nu, nu0, nu_t, S_nu0_ff):
#     tau = tau_nu(nu, nu_t)
#     ratio = stable_one_minus_exp_minus_tau_over_tau(tau)
#     return S_nu0_ff * ratio * (nu / nu0) ** -0.1

def int_Lum_old(z,S_nu0_sy,S_nu0_ff, 
            S_nu0_sy_err,S_nu0_ff_err, 
            nu_0, 
            nu_t=0.01e9,nu_1=1.0e9,nu_2=35.0e9,alpha_sy=-0.85):
    dist_conversion_factor = 3.08567758128 * 1e24  # m #3.08567758128*(10**24) #cm
    lum_conversion_factor = 1e-23
    DL_Mpc = luminosity_distance_cosmo(z=z)
    DL_cm = DL_Mpc * dist_conversion_factor  # Convert Mpc to cm
    
    # synchrotron_luminosity, synchrotron_luminosity_err = \
    #     quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy),
    #          epsrel=1e-12, limit=10000)
    # freefree_luminosity, freefree_luminosity_err = \
    #     quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff),
    #          epsrel=1e-12, limit=10000)
    synchrotron_luminosity, _ = \
        quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy),
             epsrel=1e-6, limit=10000
            )
    freefree_luminosity, _ = \
        quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff),
             epsrel=1e-6, limit=10000
            )

    synchrotron_luminosity_err, _ = \
        quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy_err),
             epsrel=1e-6, limit=10000
            )
    freefree_luminosity_err, _ = \
        quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff_err),
             epsrel=1e-6, limit=10000
            )
    
    # print(freefree_luminosity)
    # print(synchrotron_luminosity)
    int_luminosity_err = np.sqrt(synchrotron_luminosity_err**2 + freefree_luminosity_err**2.0)
    
    LR  = (4 * np.pi * (DL_cm**2) / ((1+z)**(alpha_sy+1))) * (synchrotron_luminosity + freefree_luminosity) * lum_conversion_factor
    LR_err  = (4 * np.pi * (DL_cm**2) / ((1+z)**(alpha_sy+1))) * (int_luminosity_err) * lum_conversion_factor
    return(LR,LR_err)


def int_Lum(z,S_nu0_sy,S_nu0_ff, 
            S_nu0_sy_err=None,S_nu0_ff_err=None, 
            nu_0=6e9, 
            nu_t=0.01e9,nu_1=1.0e9,nu_2=35.0e9,alpha_sy=-0.85):
    dist_conversion_factor = 3.08567758128 * 1e24  # m #3.08567758128*(10**24) #cm
    lum_conversion_factor = 1e-23
    DL_Mpc = luminosity_distance_cosmo(z=z)
    DL_cm = DL_Mpc * dist_conversion_factor  # Convert Mpc to cm
    
    # Check if calculating spectral luminosity (nu_1 = nu_2)
    is_spectral = (nu_1 == nu_2)
    
    if is_spectral:
        # Calculate spectral luminosity at nu_1 frequency
        synchrotron_luminosity = synchrotron_integrand(nu_1, nu_0, nu_t, alpha_sy, S_nu0_sy) #/ 1e9
        freefree_luminosity = freefree_integrand(nu_1, nu_0, nu_t, S_nu0_ff) #/ 1e9
        
        synchrotron_luminosity_err = 0.0
        freefree_luminosity_err = 0.0
        
        if S_nu0_sy_err is not None:
            synchrotron_luminosity_err = synchrotron_integrand(nu_1, nu_0, nu_t, alpha_sy, S_nu0_sy_err) #/ 1e9
            freefree_luminosity_err = freefree_integrand(nu_1, nu_0, nu_t, S_nu0_ff_err) #/ 1e9
    
    else:
        # Calculate integrated luminosity (original behavior)
        synchrotron_luminosity, _ = \
            quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy),
                 epsrel=1e-6, limit=10000
                )
        freefree_luminosity, _ = \
            quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff),
                 epsrel=1e-6, limit=10000
                )

        synchrotron_luminosity_err = 0.0
        freefree_luminosity_err = 0.0
        
        if S_nu0_sy_err is not None:
            synchrotron_luminosity_err, _ = \
                quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy_err),
                     epsrel=1e-6, limit=10000
                    )
            freefree_luminosity_err, _ = \
                quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff_err),
                     epsrel=1e-6, limit=10000
                    )

    # Error calculation (common for both cases)
    int_luminosity_err = 0.0
    LR_err = 0.0
    if S_nu0_sy_err is not None:
        int_luminosity_err = np.sqrt(synchrotron_luminosity_err**2 + freefree_luminosity_err**2.0)
        LR_err = (4 * np.pi * (DL_cm**2) / ((1+z)**(alpha_sy+1))) * (int_luminosity_err) * lum_conversion_factor

    # Final luminosity calculation (common for both cases)
    LR = (4 * np.pi * (DL_cm**2) / ((1+z)**(alpha_sy+1))) * (synchrotron_luminosity + freefree_luminosity) * lum_conversion_factor
    
    return(LR, LR_err)


# Simple robust median with asymmetric confidence intervals
def d_stats_basic(data, confidence=99.7, percentile_range=(15, 85)):
    """
    Calculate robust median with asymmetric confidence intervals.
    
    Parameters:
    -----------
    data : 2D array
        Your spectral index map
    confidence : float
        Confidence level in % (default: 99.7% 3sigma)
    percentile_range : tuple
        Percentile range to exclude outliers (default: 15th-85th percentile)
    
    Returns:
    --------
    median : float
        Robust median value
    sigma_lower : float  
        Lower asymmetric error (median - lower_bound)
    sigma_upper : float
        Upper asymmetric error (upper_bound - median)
    lower_bound : float
        Lower confidence bound
    upper_bound : float
        Upper confidence bound
    """
    
    # Remove NaNs and filter outliers using percentiles
    valid_data = data[~np.isnan(data)]
    p_low, p_high = np.percentile(valid_data, percentile_range)
    filtered_data = valid_data[(valid_data >= p_low) & (valid_data <= p_high)]
    
    # Calculate median
    median = np.median(filtered_data)
    
    # Calculate asymmetric confidence intervals using percentiles
    alpha = (100 - confidence) / 2  # e.g., for 99.7%, alpha = 0.15%
    lower_percentile = alpha
    upper_percentile = 100 - alpha
    
    lower_bound = np.percentile(filtered_data, lower_percentile)
    upper_bound = np.percentile(filtered_data, upper_percentile)
    
    # Asymmetric errors
    sigma_lower = median - lower_bound  # How much below median
    sigma_upper = upper_bound - median  # How much above median
    
    return median, sigma_lower, sigma_upper, lower_bound, upper_bound




'''
Radio SED
'''


# import numpy as np
# from scipy import ndimage
# import matplotlib.pyplot as plt

def compute_flux_with_uncertainties(g, levels, beam_area_, noise_rms, 
                                  method='noise_based', n_bootstrap=100):
    """
    Compute fluxes and their uncertainties using different methods
    
    Parameters:
    -----------
    g : 2D array
        Image data
    levels : array
        Intensity levels for apertures
    beam_area_ : float
        Beam area in pixels
    noise_rms : float
        RMS noise level of the image
    method : str
        'noise_based', 'bootstrap', 'jackknife', or 'level_sensitivity'
    n_bootstrap : int
        Number of bootstrap realizations
    
    Returns:
    --------
    fluxes, flux_errors, Lgrow, Lgrow_errors, areas
    """
    
    if method == 'noise_based':
        return _noise_based_uncertainties(g, levels, beam_area_, noise_rms)
    elif method == 'bootstrap':
        return _bootstrap_uncertainties(g, levels, beam_area_, noise_rms, n_bootstrap)
    elif method == 'jackknife':
        return _jackknife_uncertainties(g, levels, beam_area_, noise_rms)
    elif method == 'level_sensitivity':
        return _level_sensitivity_uncertainties(g, levels, beam_area_, noise_rms)
    else:
        raise ValueError("Method must be 'noise_based', 'bootstrap', 'jackknife', or 'level_sensitivity'")


def _noise_based_uncertainties(g, levels, beam_area_, noise_rms):
    """
    Calculate uncertainties based on noise statistics and aperture area
    """
    fluxes = []
    flux_errors = []
    areas = []
    
    for i in range(len(levels)):
        if i == 0:
            condition = (g >= levels[i])
        else:
            condition = ((g < levels[i - 1]) & (g >= levels[i]))
            
        flux = np.nansum(g * condition) / beam_area_
        area = np.nansum(condition)
        
        # Effective number of independent beams in the aperture
        n_independent_beams = area / beam_area_
        
        # Flux uncertainty: sigma_flux = sigma_noise * sqrt(N_independent) / beam_area
        flux_error = noise_rms * np.sqrt(n_independent_beams) / beam_area_
        
        fluxes.append(flux)
        flux_errors.append(flux_error)
        areas.append(area)
    
    fluxes = np.array(fluxes)
    flux_errors = np.array(flux_errors)
    areas = np.array(areas)
    
    # Cumulative flux and its uncertainty
    Lgrow = np.nancumsum(fluxes)
    
    # For cumulative sum, uncertainties add in quadrature
    Lgrow_errors = np.sqrt(np.nancumsum(flux_errors**2))
    
    return fluxes, flux_errors, Lgrow, Lgrow_errors, areas


def _bootstrap_uncertainties(g, levels, beam_area_, noise_rms, n_bootstrap):
    """
    Calculate uncertainties using bootstrap resampling with noise
    """
    def compute_single_realization(g_noisy):
        fluxes = []
        for i in range(len(levels)):
            if i == 0:
                condition = (g_noisy >= levels[i])
            else:
                condition = ((g_noisy < levels[i - 1]) & (g_noisy >= levels[i]))
            flux = np.nansum(g_noisy * condition) / beam_area_
            fluxes.append(flux)
        return np.array(fluxes)
    
    # Generate bootstrap realizations
    bootstrap_fluxes = []
    for _ in range(n_bootstrap):
        # Add noise realization
        noise = np.random.normal(0, noise_rms, g.shape)
        g_noisy = g + noise
        boot_fluxes = compute_single_realization(g_noisy)
        bootstrap_fluxes.append(boot_fluxes)
    
    bootstrap_fluxes = np.array(bootstrap_fluxes)
    
    # Calculate mean and standard deviation
    fluxes = np.mean(bootstrap_fluxes, axis=0)
    flux_errors = np.std(bootstrap_fluxes, axis=0)
    
    # Calculate cumulative quantities
    bootstrap_Lgrow = np.nancumsum(bootstrap_fluxes, axis=1)
    Lgrow = np.mean(bootstrap_Lgrow, axis=0)
    Lgrow_errors = np.std(bootstrap_Lgrow, axis=0)
    
    # Calculate areas from original image
    areas = []
    for i in range(len(levels)):
        if i == 0:
            condition = (g >= levels[i])
        else:
            condition = ((g < levels[i - 1]) & (g >= levels[i]))
        areas.append(np.nansum(condition))
    areas = np.array(areas)
    
    return fluxes, flux_errors, Lgrow, Lgrow_errors, areas


def _jackknife_uncertainties(g, levels, beam_area_, noise_rms):
    """
    Calculate uncertainties using spatial jackknife resampling
    """
    # Divide image into blocks for jackknife
    ny, nx = g.shape
    block_size = max(int(np.sqrt(beam_area_)), 3)  # Blocks roughly beam-sized
    
    jackknife_fluxes = []
    
    for by in range(0, ny, block_size):
        for bx in range(0, nx, block_size):
            # Create copy without this block
            g_jack = g.copy()
            g_jack[by:by+block_size, bx:bx+block_size] = np.nan
            
            # Compute fluxes without this block
            fluxes_jack = []
            for i in range(len(levels)):
                if i == 0:
                    condition = (g_jack >= levels[i])
                else:
                    condition = ((g_jack < levels[i - 1]) & (g_jack >= levels[i]))
                flux = np.nansum(g_jack * condition) / beam_area_
                fluxes_jack.append(flux)
            
            jackknife_fluxes.append(fluxes_jack)
    
    jackknife_fluxes = np.array(jackknife_fluxes)
    
    # Calculate jackknife statistics
    n_jack = len(jackknife_fluxes)
    fluxes = np.nanmean(jackknife_fluxes, axis=0)
    flux_errors = np.sqrt((n_jack - 1) / n_jack * 
                         np.nansum((jackknife_fluxes - fluxes[None, :])**2, axis=0))
    
    # Cumulative quantities
    jackknife_Lgrow = np.nancumsum(jackknife_fluxes, axis=1)
    Lgrow = np.nanmean(jackknife_Lgrow, axis=0)
    Lgrow_errors = np.sqrt((n_jack - 1) / n_jack * 
                          np.nansum((jackknife_Lgrow - Lgrow[None, :])**2, axis=0))
    
    # Calculate areas
    areas = []
    for i in range(len(levels)):
        if i == 0:
            condition = (g >= levels[i])
        else:
            condition = ((g < levels[i - 1]) & (g >= levels[i]))
        areas.append(np.nansum(condition))
    areas = np.array(areas)
    
    return fluxes, flux_errors, Lgrow, Lgrow_errors, areas


def _level_sensitivity_uncertainties(g, levels, beam_area_, noise_rms):
    """
    Calculate uncertainties by varying the contour levels
    """
    # Vary levels by Â+/-noise_rms
    level_variations = []
    
    for delta in [-noise_rms, 0, noise_rms]:
        varied_levels = levels + delta
        varied_fluxes = []
        
        for i in range(len(varied_levels)):
            if i == 0:
                condition = (g >= varied_levels[i])
            else:
                condition = ((g < varied_levels[i - 1]) & (g >= varied_levels[i]))
            flux = np.nansum(g * condition) / beam_area_
            varied_fluxes.append(flux)
        
        level_variations.append(varied_fluxes)
    
    level_variations = np.array(level_variations)
    
    # Use central values and spread as uncertainty
    fluxes = level_variations[1]  # Central values (delta=0)
    flux_errors = (np.max(level_variations, axis=0) - 
                   np.min(level_variations, axis=0)) / 2
    
    # Cumulative quantities
    Lgrow = np.nancumsum(fluxes)
    Lgrow_variations = np.nancumsum(level_variations, axis=1)
    Lgrow_errors = (np.max(Lgrow_variations, axis=0) - 
                   np.min(Lgrow_variations, axis=0)) / 2
    
    # Calculate areas
    areas = []
    for i in range(len(levels)):
        if i == 0:
            condition = (g >= levels[i])
        else:
            condition = ((g < levels[i - 1]) & (g >= levels[i]))
        areas.append(np.nansum(condition))
    areas = np.array(areas)
    
    return fluxes, flux_errors, Lgrow, Lgrow_errors, areas



import numpy as np
from scipy import stats
from typing import Optional, Tuple

def kendall_tau(
    x_data: np.ndarray, 
    y_data: np.ndarray, 
    x_errors: Optional[np.ndarray] = None, 
    y_errors: Optional[np.ndarray] = None,
    n_bootstrap: int = 1000,
    seed: Optional[int] = None,
    method: Optional[str] = 'auto',
    variant: Optional[str] = 'b'
) -> Tuple[float, float, float, int]:
    """
    Calculate Kendall's tau correlation coefficient with uncertainty estimation.
    Handles NaNs and Infs, and can incorporate measurement errors if provided.
    
    Parameters:
    -----------
    x_data : array-like
        Array of x values
    y_data : array-like
        Array of y values
    x_errors : array-like, optional
        Array of measurement errors (uncertainties) for x values
    y_errors : array-like, optional
        Array of measurement errors (uncertainties) for y values
    n_bootstrap : int, default=1000
        Number of bootstrap samples for uncertainty estimation
    seed : int, optional
        Random seed for reproducibility
        
    Returns:
    --------
    tau : float
        Kendall's tau correlation coefficient
    tau_uncertainty : float
        Estimated uncertainty (standard error) of tau
    p_value : float
        Two-sided p-value for a hypothesis test with null hypothesis: tau = 0
    n_valid : int
        Number of valid (non-NaN, non-Inf) data points used in the calculation
    """
    # Convert inputs to numpy arrays
    x_data = np.asarray(x_data)
    y_data = np.asarray(y_data)
    
    if x_errors is not None:
        x_errors = np.asarray(x_errors)
    if y_errors is not None:
        y_errors = np.asarray(y_errors)
    
    # Check for consistent input dimensions
    if len(x_data) != len(y_data):
        raise ValueError("x_data and y_data must have the same length")
    
    if x_errors is not None and len(x_data) != len(x_errors):
        raise ValueError("x_data and x_errors must have the same length")
        
    if y_errors is not None and len(y_data) != len(y_errors):
        raise ValueError("y_data and y_errors must have the same length")
    
    # Create mask for valid (non-NaN, non-Inf) data points
    valid_mask = np.isfinite(x_data) & np.isfinite(y_data)
    
    if x_errors is not None:
        valid_mask &= np.isfinite(x_errors)
    if y_errors is not None:
        valid_mask &= np.isfinite(y_errors)
    
    # Apply mask to get valid data
    x_valid = x_data[valid_mask]
    y_valid = y_data[valid_mask]
    
    # Extract valid errors if provided
    x_err_valid = None if x_errors is None else x_errors[valid_mask]
    y_err_valid = None if y_errors is None else y_errors[valid_mask]
    
    # Check if we have enough valid data points
    n_valid = len(x_valid)
    if n_valid < 2:
        return np.nan, np.nan, np.nan, n_valid
    
    # Calculate Kendall's tau on valid data
    tau, p_value = stats.kendalltau(x_valid, y_valid,method=method,variant=variant)
    
    # Initialize random number generator for bootstrap
    rng = np.random.RandomState(seed)
    
    # Perform bootstrap to estimate tau uncertainty
    tau_bootstrap = np.zeros(n_bootstrap)
    
    for i in range(n_bootstrap):
        # Generate bootstrap sample indices with replacement
        indices = rng.randint(0, n_valid, size=n_valid)
        
        # Create bootstrapped data points
        x_boot = x_valid[indices].copy()
        y_boot = y_valid[indices].copy()
        
        # If errors are provided, perturb the bootstrap samples
        if x_err_valid is not None:
            x_boot += rng.normal(0, x_err_valid[indices])
        if y_err_valid is not None:
            y_boot += rng.normal(0, y_err_valid[indices])
        
        # Calculate tau for this bootstrap sample
        tau_boot, _ = stats.kendalltau(x_boot, y_boot,method=method,variant=variant)
        tau_bootstrap[i] = tau_boot
    
    # Calculate standard error from bootstrap distribution
    tau_uncertainty = np.std(tau_bootstrap)
    
    return tau, tau_uncertainty, p_value, n_valid



import numpy as np
from scipy.stats import kendalltau
from typing import Optional, Tuple


# import numpy as np
# from scipy.stats import kendalltau
# from typing import Optional, Tuple

def kendall_tau_with_uncertainty(x_data: np.ndarray,
                                 y_data: np.ndarray,
                                 x_err: Optional[np.ndarray] = None,
                                 y_err: Optional[np.ndarray] = None,
                                 n_bootstrap: int = 1000,
                                 n_monte_carlo: int = 500,
                                 seed: Optional[int] = None,
                                 method: Optional[str] = 'auto',
                                 variant: Optional[str] = 'b') -> Tuple[float, float, float, int]:
    """
    Compute Kendall's tau correlation coefficient and its uncertainty with proper error propagation.
    
    Parameters:
        x_data (np.ndarray): Array of x data values.
        y_data (np.ndarray): Array of y data values.
        x_err (Optional[np.ndarray]): Optional 1sigma errors on x data.
        y_err (Optional[np.ndarray]): Optional 1sigma errors on y data.
        n_bootstrap (int): Number of bootstrap samples to estimate uncertainty.
        n_monte_carlo (int): Number of Monte Carlo iterations for error propagation.
        seed (Optional[int]): Random seed for reproducibility.
        method (Optional[str]): Method to compute Kendall's tau ('auto', 'asymptotic', 'exact').
        variant (Optional[str]): Variant of Kendall's tau ('b' or 'c').
        
    Returns:
        tau (float): Kendall's tau coefficient with error propagation.
        tau_std (float): Standard deviation (uncertainty) of tau.
        p_value (float): p-value of the correlation test.
        n_valid (int): Number of valid data points used.
    """
    rng = np.random.default_rng(seed)

    # Clean the data: remove NaNs, Infs
    mask = np.isfinite(x_data) & np.isfinite(y_data)
    if x_err is not None:
        mask &= np.isfinite(x_err)
    if y_err is not None:
        mask &= np.isfinite(y_err)

    x = np.asarray(x_data)[mask]
    y = np.asarray(y_data)[mask]
    n_valid = len(x)

    if x_err is not None:
        x_err = np.asarray(x_err)[mask]
    else:
        x_err = np.zeros_like(x)  # No error case
        
    if y_err is not None:
        y_err = np.asarray(y_err)[mask]
    else:
        y_err = np.zeros_like(y)  # No error case

    if n_valid < 2:
        raise ValueError("Not enough valid data points after cleaning for Kendall's tau.")

    # Original correlation for p-value calculation
    orig_tau, p_value = kendalltau(x, y, variant=variant, method=method)
    
    # If no errors provided, use the traditional approach
    if np.all(x_err == 0) and np.all(y_err == 0):
        return orig_tau, 0.0, p_value, n_valid
    
    # Monte Carlo approach to propagate measurement errors into tau
    tau_mc_samples = []
    
    for _ in range(n_monte_carlo):
        # Create perturbed dataset
        x_perturbed = rng.normal(x, x_err)
        y_perturbed = rng.normal(y, y_err)
        
        # Calculate tau for this perturbed dataset
        try:
            mc_tau, _ = kendalltau(x_perturbed, y_perturbed, variant=variant, method=method)
            if np.isfinite(mc_tau):
                tau_mc_samples.append(mc_tau)
        except Exception:
            continue
    
    # Use mean of Monte Carlo samples as our best estimate of tau
    if len(tau_mc_samples) > 0:
        tau = np.mean(tau_mc_samples)
    else:
        tau = orig_tau
    
    # Bootstrap to estimate uncertainty (combining sampling variability and measurement errors)
    tau_samples = []
    for _ in range(n_bootstrap):
        # Resample with replacement
        indices = rng.choice(n_valid, n_valid, replace=True)
        
        x_sample = x[indices]
        y_sample = y[indices]
        x_err_sample = x_err[indices]
        y_err_sample = y_err[indices]
        
        # Perturb values with Gaussian noise according to measurement errors
        x_sample = rng.normal(x_sample, x_err_sample)
        y_sample = rng.normal(y_sample, y_err_sample)
        
        # Compute tau on the bootstrap sample
        try:
            tau_boot, _ = kendalltau(x_sample, y_sample, variant=variant, method=method)
            if np.isfinite(tau_boot):
                tau_samples.append(tau_boot)
        except Exception:
            continue
    
    tau_std = np.std(tau_samples) if len(tau_samples) > 1 else np.nan

    return tau, tau_std, p_value, n_valid


def expand_limits(factor=2.0):
    ax = plt.gca()
    ax.relim()
    ax.autoscale_view()

    for get, set in [(ax.get_xlim, ax.set_xlim), (ax.get_ylim, ax.set_ylim)]:
        try:
            low, high = get()
            if not np.isfinite(low) or not np.isfinite(high) or low == high:
                continue  # Skip invalid or zero-span limits
            center = (low + high) / 2
            span = (high - low) * factor / 2
            set(center - span, center + span)
        except Exception as e:
            print(f"Skipping limit adjustment due to: {e}")


def spectral_correction(data, data_frequency, target_frequency, spectral_index):
    """
    Apply spectral correction to a radio image or a flux measurement using a power-law model.

    Parameters:
    -----------
    data : numpy.ndarray
        1D or 2D array containing the radio image or flux measurement
    data_frequency : float
        Original/reference frequency of the data (in GHz)
    target_frequency : float
        Target frequency to correct the image to (in GHz)
    spectral_index : float
        Spectral index $\alpha$ for the power-law model $S_{\nu} \propto \nu^\alpha$
        
    Returns:
    --------
    numpy.ndarray
        Corrected image at the target frequency
        
    Notes:
    ------
    Uses the power-law relation: $S_{\nu} = S_{\nu_0} \times (\nu/\nu_0)^{\alpha}$
    where $S_{\nu_0}$ is flux at reference frequency $\nu_0$
    """
    # Calculate the frequency ratio
    freq_ratio = target_frequency / data_frequency
    
    # Apply the power-law correction
    # S_target = S_ref Ã- (nu_target/nu_ref)^alpha
    correction_factor = freq_ratio ** spectral_index
    
    # Apply correction
    corrected_data = data * correction_factor

    return corrected_data



def get_err_frac_v2(x, y, x_err, y_err):
    """
    Calculate the fractional error of a ratio of two values with support for 
    both symmetric and asymmetric uncertainties.
    
    Parameters:
    -----------
    x : float or array-like
        Numerator value(s)
    y : float or array-like  
        Denominator value(s)
    x_err : float, array-like, or [lower_array, upper_array]
        Error for x. Can be:
        - Single value or array for symmetric errors
        - List/tuple of two arrays: [lower_bounds, upper_bounds] for asymmetric errors
    y_err : float, array-like, or [lower_array, upper_array] 
        Error for y. Same format options as x_err
        
    Returns:
    --------
    z : float or array-like
        The ratio x/y
    sigma_z : float, array-like, or [lower_array, upper_array]
        Error in z. Returns single value/array for symmetric errors,
        list [lower_bounds, upper_bounds] for asymmetric errors
        
    Notes:
    ------
    - For asymmetric errors, pass [lower_bounds, upper_bounds] as error arguments
    - Function handles element-wise symmetric/asymmetric cases appropriately
    - NaN values in input errors are preserved in output
    """
    
    # Calculate the central value
    z = x / y
    
    # Parse error inputs
    x_lower, x_upper, x_has_asym_structure = _parse_error_input(x_err)
    y_lower, y_upper, y_has_asym_structure = _parse_error_input(y_err)
    
    # Determine if we need asymmetric error handling
    if x_has_asym_structure or y_has_asym_structure:
        # Check if any elements actually have different bounds
        x_asym_elements = np.logical_not(np.isclose(x_lower, x_upper, equal_nan=True)) if x_has_asym_structure else False
        y_asym_elements = np.logical_not(np.isclose(y_lower, y_upper, equal_nan=True)) if y_has_asym_structure else False
        
        # If any element is truly asymmetric, return asymmetric format
        if np.any(x_asym_elements) or np.any(y_asym_elements):
            sigma_z_lower, sigma_z_upper = _calculate_asymmetric_errors(
                x, y, z, x_lower, x_upper, y_lower, y_upper
            )
            return z, [sigma_z_lower, sigma_z_upper]
    
    # Use symmetric error propagation (original method)
    # Take the upper bound as the symmetric error (identical to lower for symmetric cases)
    x_err_sym = x_upper
    y_err_sym = y_upper
    
    _sigma_z = z * np.sqrt((x_err_sym / x)**2 + (y_err_sym / y)**2)
    sigma_z = np.nan_to_num(_sigma_z, nan=0, posinf=0, neginf=0)
    
    return z, sigma_z


def _parse_error_input(err):
    """
    Parse error input and return lower bounds, upper bounds, and structure flag.
    
    Returns:
    --------
    lower, upper : error bounds (positive values)
    has_asym_structure : bool indicating if input has [lower, upper] structure
    """
    if err is None:
        return 0, 0, False
    
    # Check if input is a list or tuple with exactly 2 elements
    if isinstance(err, (list, tuple)) and len(err) == 2:
        lower = np.abs(np.asarray(err[0]))
        upper = np.abs(np.asarray(err[1]))
        return lower, upper, True
    
    # Single array or scalar - symmetric error
    err_abs = np.abs(np.asarray(err))
    return err_abs, err_abs, False


def _calculate_asymmetric_errors(x, y, z, x_lower, x_upper, y_lower, y_upper):
    """
    Calculate asymmetric error propagation using linearized approximation.
    
    For each element, determines whether to use symmetric or asymmetric propagation
    based on whether the bounds actually differ.
    """
    # Ensure all inputs are arrays for consistent handling
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    x_lower = np.asarray(x_lower)
    x_upper = np.asarray(x_upper)
    y_lower = np.asarray(y_lower)
    y_upper = np.asarray(y_upper)
    
    # Initialize output arrays
    if z.ndim == 0:
        # Scalar case
        return _propagate_single_asymmetric_error(x, y, z, x_lower, x_upper, y_lower, y_upper)
    
    # Array case - handle element by element
    sigma_z_lower = np.zeros_like(z, dtype=float)
    sigma_z_upper = np.zeros_like(z, dtype=float)
    
    # Check which elements are truly asymmetric
    x_is_asym = ~np.isclose(x_lower, x_upper, equal_nan=True)
    y_is_asym = ~np.isclose(y_lower, y_upper, equal_nan=True)
    any_asym = x_is_asym | y_is_asym
    
    # Handle symmetric elements
    sym_mask = ~any_asym
    if np.any(sym_mask):
        # Use original symmetric formula
        x_err_sym = x_upper[sym_mask]
        y_err_sym = y_upper[sym_mask]
        sigma_sym = z[sym_mask] * np.sqrt((x_err_sym / x[sym_mask])**2 + (y_err_sym / y[sym_mask])**2)
        sigma_sym = np.nan_to_num(sigma_sym, nan=0, posinf=0, neginf=0)
        sigma_z_lower[sym_mask] = sigma_sym
        sigma_z_upper[sym_mask] = sigma_sym
    
    # Handle asymmetric elements
    if np.any(any_asym):
        asym_indices = np.where(any_asym)[0]
        for i in asym_indices:
            lower_err, upper_err = _propagate_single_asymmetric_error(
                x[i], y[i], z[i], x_lower[i], x_upper[i], y_lower[i], y_upper[i]
            )
            sigma_z_lower[i] = lower_err
            sigma_z_upper[i] = upper_err
    
    return sigma_z_lower, sigma_z_upper


def _propagate_single_asymmetric_error(x, y, z, x_lower, x_upper, y_lower, y_upper):
    """
    Calculate asymmetric error propagation for a single element using linearized approximation.
    """
    # Handle NaN cases
    if np.isnan(x) or np.isnan(y) or np.isnan(x_lower) or np.isnan(x_upper) or np.isnan(y_lower) or np.isnan(y_upper):
        return 0.0, 0.0
    
    # Partial derivatives
    dz_dx = 1 / y
    dz_dy = -x / (y**2)
    
    # Calculate all possible combinations
    combinations = [
        (x_upper, -y_lower),   # x high, y low -> maximizes z
        (x_upper, y_upper),    # x high, y high
        (-x_lower, -y_lower),  # x low, y low
        (-x_lower, y_upper)    # x low, y high -> minimizes z
    ]
    
    deviations = []
    for dx, dy in combinations:
        dz = dz_dx * dx + dz_dy * dy
        if not np.isnan(dz) and not np.isinf(dz):
            deviations.append(dz)
    
    if not deviations:
        return 0.0, 0.0
    
    # Upper and lower errors
    max_deviation = max(deviations)
    min_deviation = min(deviations)
    
    sigma_z_upper = max(max_deviation, 0)
    sigma_z_lower = max(-min_deviation, 0)
    
    return sigma_z_lower, sigma_z_upper




def compute_diameters_weighted(points, weights, hull=None, percentile=90):
    """
    Compute intensity-weighted major and minor diameters of a structure.
    
    Parameters
    ----------
    points : ndarray
        Array of (x, y) coordinates.
    weights : ndarray
        Intensity weights for each point.
    hull : scipy.spatial.ConvexHull, optional
        Pre-computed ConvexHull object. If None, will be computed from points.
    percentile : float, optional
        Percentile of intensity to use for defining the effective boundary (default 90).
    
    Returns
    -------
    dict
        Dictionary containing major and minor diameters and their defining points.
    """
    from itertools import combinations

    if hull is None:
        hull = ConvexHull(points)
    
    # For weighted approach, consider high-intensity regions more heavily
    # Sort points by weight and select top percentile
    sorted_indices = np.argsort(weights)[::-1]
    n_select = max(int(len(points) * (percentile / 100.0)), 3)
    high_intensity_indices = sorted_indices[:n_select]
    high_intensity_points = points[high_intensity_indices]
    
    # Compute hull of high-intensity points for more robust diameter estimation
    if len(high_intensity_points) >= 3:
        try:
            intensity_hull = ConvexHull(high_intensity_points)
            effective_points = high_intensity_points[intensity_hull.vertices]
        except:
            # Fall back to original hull if high-intensity hull fails
            effective_points = points[hull.vertices]
    else:
        effective_points = points[hull.vertices]
    
    # Calculate Major Diameter
    major_diameter = 0
    major_points = None
    for p1, p2 in combinations(effective_points, 2):
        distance = np.linalg.norm(p1 - p2)
        if distance > major_diameter:
            major_diameter = distance
            major_points = (p1, p2)
    
    # Calculate Minor Diameter using weighted approach
    minor_diameter = float('inf')
    minor_points = None
    num_points = len(effective_points)
    
    for i in range(num_points):
        p1, p2 = effective_points[i], effective_points[(i + 1) % num_points]
        edge_vector = p2 - p1
        edge_length = np.linalg.norm(edge_vector)
        
        if edge_length == 0:
            continue
        edge_normal = np.array([-edge_vector[1], edge_vector[0]]) / edge_length
        
        distances = np.abs(np.dot(effective_points - p1, edge_normal))
        max_distance = distances.max()
        if max_distance < minor_diameter:
            minor_diameter = max_distance
            projections = effective_points[np.abs(distances - max_distance) < 1e-6]
            if len(projections) >= 2:
                minor_points = (projections[0], projections[1])
    
    return {
        "major_diameter": major_diameter,
        "major_points": major_points,
        "minor_diameter": minor_diameter,
        "minor_points": minor_points
    }


def intensity_weighted_morpho(image, mask, scale=1.0, do_plot=False, weight_power=1.0):
    """
    Perform intensity-weighted morphological analysis on a galaxy structure.
    
    This function computes morphological properties using intensity-weighted moments,
    providing more accurate measurements that account for the light distribution
    rather than just the spatial extent of the emission.
    
    Parameters
    ----------
    image : 2D ndarray
        The image data containing the galaxy emission.
    mask : 2D ndarray
        Binary mask defining the structure to be analyzed.
    scale : float, optional
        Scaling factor for axis visualization (default 1.0).
    do_plot : bool, optional
        Whether to plot the results (default False).
    weight_power : float, optional
        Power to raise intensities to for weighting (default 1.0).
        Higher values give more weight to bright regions.
    
    Returns
    -------
    dict
        Dictionary containing:
        - "PA_weighted" : Position angle in degrees
        - "q_weighted" : Axis ratio (minor/major)
        - "centroid_weighted" : Intensity-weighted centroid [x, y]
        - "major_diameter_weighted" : Weighted major diameter
        - "minor_diameter_weighted" : Weighted minor diameter
        - "PA_unweighted" : Unweighted position angle for comparison
        - "q_unweighted" : Unweighted axis ratio for comparison
        - "centroid_unweighted" : Unweighted centroid for comparison
    """
    from itertools import combinations
    # Extract coordinates and intensities
    indices = np.transpose(np.nonzero(mask))
    y, x = indices[:, 0], indices[:, 1]
    points = np.column_stack((x, y))
    
    # Get intensity values at masked positions
    intensities = image[mask].flatten()
    
    # Apply weight power for emphasis on bright regions
    weights = np.power(intensities, weight_power)
    weights = weights / np.sum(weights)  # Normalize weights
    
    # Compute intensity-weighted centroid
    centroid_weighted = np.average(points, weights=weights, axis=0)
    
    # Compute unweighted centroid for comparison
    centroid_unweighted = np.mean(points, axis=0)
    
    # Center points around weighted centroid
    centered_points = points - centroid_weighted
    
    # Compute intensity-weighted covariance matrix
    # This gives us the intensity-weighted second moments
    cov_weighted = np.zeros((2, 2))
    for i in range(len(centered_points)):
        p = centered_points[i].reshape(-1, 1)
        cov_weighted += weights[i] * np.dot(p, p.T)
    
    # Eigenanalysis of weighted covariance
    eigenvalues_w, eigenvectors_w = np.linalg.eig(cov_weighted)
    order_w = np.argsort(eigenvalues_w)[::-1]
    eigenvalues_w = eigenvalues_w[order_w]
    eigenvectors_w = eigenvectors_w[:, order_w]
    
    # Also compute unweighted for comparison
    cov_unweighted = np.cov(points, rowvar=False)
    eigenvalues_u, eigenvectors_u = np.linalg.eig(cov_unweighted)
    order_u = np.argsort(eigenvalues_u)[::-1]
    eigenvalues_u = eigenvalues_u[order_u]
    eigenvectors_u = eigenvectors_u[:, order_u]
    
    # Weighted major and minor axes
    major_axis_vector_w = eigenvectors_w[:, 0]
    minor_axis_vector_w = eigenvectors_w[:, 1]
    
    # Scale for visualization
    major_axis_w = scale * np.sqrt(eigenvalues_w[0]) * major_axis_vector_w
    minor_axis_w = scale * np.sqrt(eigenvalues_w[1]) * minor_axis_vector_w
    
    # Position angle (weighted)
    position_angle_w = np.arctan2(major_axis_vector_w[1], major_axis_vector_w[0])
    position_angle_degrees_w = np.degrees(position_angle_w)
    if position_angle_degrees_w < 0:
        position_angle_degrees_w += 360
    
    # Position angle (unweighted)
    major_axis_vector_u = eigenvectors_u[:, 0]
    position_angle_u = np.arctan2(major_axis_vector_u[1], major_axis_vector_u[0])
    position_angle_degrees_u = np.degrees(position_angle_u)
    if position_angle_degrees_u < 0:
        position_angle_degrees_u += 360
    
    # Axis ratios
    axis_ratio_w = np.sqrt(eigenvalues_w[1] / eigenvalues_w[0])
    axis_ratio_u = np.sqrt(eigenvalues_u[1] / eigenvalues_u[0])
    
    # Compute ConvexHull for diameter calculations
    hull = ConvexHull(points)
    
    # Compute weighted diameters
    diameters_weighted = compute_diameters_weighted(points, intensities, hull, percentile=90)
    
    # Compute unweighted diameters for comparison
    hull_points = points[hull.vertices]
    major_diameter_u = 0
    for p1, p2 in combinations(hull_points, 2):
        distance = np.linalg.norm(p1 - p2)
        if distance > major_diameter_u:
            major_diameter_u = distance
    
    if do_plot:
        plt.figure(figsize=(10, 5))
        
        # Left panel: Weighted analysis
        plt.subplot(1, 2, 1)
        # Use simple_norm if available, otherwise use standard normalization
        try:
            # from astropy.visualization import simple_norm
            # from astropy.stats import mad_std
            norm = simple_norm(image, stretch='sqrt', asinh_a=0.02, 
                             vmin=3*mad_std(image), vmax=0.2*np.nanmax(image))
            plt.imshow(image, cmap='gray', origin='lower', norm=norm)
        except ImportError:
            plt.imshow(image, cmap='gray', origin='lower', 
                      vmin=np.percentile(image[mask], 1),
                      vmax=np.percentile(image[mask], 99.5))
        
        # Plot weighted axes
        plt.quiver(centroid_weighted[0], centroid_weighted[1], 
                  major_axis_w[0], major_axis_w[1],
                  angles='xy', scale_units='xy', scale=1, 
                  color='limegreen', width=0.003, label='Major (weighted)')
        plt.quiver(centroid_weighted[0], centroid_weighted[1], 
                  minor_axis_w[0], minor_axis_w[1],
                  angles='xy', scale_units='xy', scale=1, 
                  color='red', width=0.003, label='Minor (weighted)')
        
        # Plot centroids
        plt.scatter(centroid_weighted[0], centroid_weighted[1], 
                   color='limegreen', s=50, zorder=5, label='Weighted centroid')
        plt.scatter(centroid_unweighted[0], centroid_unweighted[1], 
                   color='cyan', s=30, zorder=5, marker='x', label='Unweighted centroid')
        
        plt.title('Intensity-Weighted Analysis')
        plt.xlabel("x [pixels]")
        plt.ylabel("y [pixels]")
        plt.legend(fontsize=8)
        plt.axis('equal')
        
        # Right panel: Intensity distribution
        plt.subplot(1, 2, 2)
        scatter = plt.scatter(points[:, 0], points[:, 1], 
                            c=intensities, cmap='viridis', 
                            s=1, alpha=0.5)
        plt.colorbar(scatter, label='Intensity')
        
        # Show high-intensity region used for diameter calculation
        sorted_indices = np.argsort(intensities)[::-1]
        n_select = max(int(len(points) * 0.9), 3)
        high_int_indices = sorted_indices[:n_select]
        plt.scatter(points[high_int_indices, 0], points[high_int_indices, 1],
                   s=0.5, color='red', alpha=0.3, label='90th percentile')
        
        plt.scatter(centroid_weighted[0], centroid_weighted[1], 
                   color='limegreen', s=50, zorder=5)
        plt.title('Intensity Distribution')
        plt.xlabel("x [pixels]")
        plt.ylabel("y [pixels]")
        plt.legend(fontsize=8)
        plt.axis('equal')
        
        plt.tight_layout()
        plt.show()
    
    # Return comprehensive report
    report = {
        "PA_weighted": position_angle_degrees_w,
        "q_weighted": axis_ratio_w,
        "centroid_weighted": centroid_weighted,
        "major_diameter_weighted": diameters_weighted["major_diameter"],
        "minor_diameter_weighted": diameters_weighted["minor_diameter"],
        "PA_unweighted": position_angle_degrees_u,
        "q_unweighted": axis_ratio_u,
        "centroid_unweighted": centroid_unweighted,
        "major_diameter_unweighted": major_diameter_u
    }
    
    return report

def _jackknife_uncertainties(g, levels, beam_area_):
    """
    Calculate uncertainties using the jackknife resampling method
    """
    ny, nx = g.shape
    block_size = max(ny // 10, nx // 10, 1) # Define block size for jackknife       
    jackknife_fluxes = []
    for by in range(0, ny, block_size):
        for bx in range(0, nx, block_size):
            # Create a copy of the data
            g_jack = g.copy()
            # Mask out the current block
            g_jack[by:by+block_size, bx:bx+block_size] = np.nan
            # Calculate fluxes for each level
            fluxes_jack = []
            for i in range(len(levels)):
                if i == 0:
                    condition = (g_jack >= levels[i])
                else:
                    condition = ((g_jack < levels[i - 1]) & (g_jack >= levels[i]))
                flux = np.nansum(g_jack * condition) / beam_area_
                fluxes_jack.append(flux)
            jackknife_fluxes.append(fluxes_jack)
    jackknife_fluxes = np.array(jackknife_fluxes)
    n_jack = jackknife_fluxes.shape[0]
    # Mean fluxes
    fluxes = np.nanmean(jackknife_fluxes, axis=0)
    # Jackknife errors
    flux_errors = np.sqrt((n_jack - 1) / n_jack *
                            np.nansum((jackknife_fluxes - fluxes[None, :])**2, axis=0))
    # Cumulative quantities
    Lgrow = np.nancumsum(fluxes)
    jackknife_Lgrow = np.nancumsum(jackknife_fluxes, axis=1
    )
    Lgrow_errors = np.sqrt((n_jack - 1) / n_jack *
                            np.nansum((jackknife_Lgrow - Lgrow[None, :])**2, axis=0))
    # Calculate areas
    areas = []
    for i in range(len(levels)):
        if i == 0:
            condition = (g >= levels[i])
        else:
            condition = ((g < levels[i - 1]) & (g >= levels[i]))
        areas.append(np.nansum(condition))
    areas = np.array(areas)
    return fluxes, flux_errors, Lgrow, Lgrow_errors, areas



def create_radial_mask(mask, centre=None, max_radius=None,
                       iterations=5,dilation_size=5):

    # Get the shape of the mask
    y, x = np.indices(mask.shape)
    if centre is None:
        centre = (mask.shape[1]/2,mask.shape[0]/2)
        print(f'Using centre {max_radius}')
    if max_radius is None:
        max_radius = mask.shape[0]/10
        print(f'Using max radius of {max_radius}')
    
    # Calculate the distance from the center
    distance_from_centre = np.sqrt((x - centre[1])**2 + (y - centre[0])**2)
    
    # Create a new mask that keeps only pixels within the max_radius
    radial_mask = distance_from_centre <= max_radius
    
    # Apply the radial mask to the original mask
    new_mask = mask * radial_mask
    _,mask_region = mask_dilation_from_mask(mask,new_mask,
                                                 iterations=iterations,dilation_size=dilation_size,
                                                 PLOT=False,show_figure=False)

    _,mask_for_fit = mask_dilation_from_mask(mask,new_mask,
                                                 iterations=iterations,dilation_size=dilation_size*2,
                                                 PLOT=False,show_figure=False)

    plt.figure()
    plt.imshow(mask_region,origin='lower',cmap='magma')
    plt.show()

    plt.figure()
    plt.imshow(mask_for_fit,origin='lower',cmap='magma')
    plt.show()

    
    return mask_region,mask_for_fit



def elliptical_radial_profile(
    image,
    x0,
    y0,
    pa_deg,
    q,
    beam_area,
    rms=None,
    delta_r=1.0,
    sigma_clip_val=3.0,
    wedge_half_angle_deg=22.5,
    mask=None,
    beam_correction=True,
):
    """
    Compute elliptical radial intensity profiles and luminosity growth curves
    with full directional decomposition along the major axis, minor axis, and
    diagonal direction.

    This function is a geometry-aware replacement for the level-based profiling
    in ``compute_image_properties``.  It returns the same named quantities
    (``fluxes``, ``fluxes_err``, ``areas``, ``agrow``, ``Lgrow``,
    ``Lgrow_err``, ``Lgrow_norm``, ``Lgrow_err_norm``, ``radii``) so it can
    be used as a drop-in, while also returning directional profiles and the raw
    semi-major-axis grid.

    The coordinate transform is::

        xl =  (x - x0)*cos(PA) + (y - y0)*sin(PA)   # major-axis frame
        yl = -(x - x0)*sin(PA) + (y - y0)*cos(PA)   # minor-axis frame
        rl =  sqrt(xl^2 + (yl/q)^2)                    # elliptical radius (= semi-major axis)
        φ  =  arctan2(yl/q, xl)                      # azimuthal angle in de-projected frame

    The four radius conventions returned are:

    * ``radii_major``  = a               (projection along major axis)
    * ``radii_minor``  = a x q           (projection along minor axis)
    * ``radii_diag``   = a x qx sqrt 2 /  sqrt (q^2+1)  (45^o between axes)
    * ``radii_circ``   = a x  sqrt q          (area-equivalent circular radius)

    The main ``radii`` key equals ``radii_circ`` so it is directly comparable
    with the level-based (isophotal) radii produced elsewhere.

    Uncertainty model
    -----------------
    Within each annulus the effective number of independent resolution elements
    is ``N_eff = N_pix / beam_area``.  The uncertainty on the mean intensity
    (IR) and on the per-annulus flux contribution are:

        IR_err  = sigma_clip /  sqrt N_eff
        flux_err = sigma_clip x  sqrt N_eff        (= IR_err x N_pix / beam_area)

    where ``sigma_clip`` is the sigma-clipped standard deviation of pixel values
    in the annulus.  This correctly accounts for beam-to-beam correlations; the
    naive independent-pixel formula (used in the old level-based code) over-
    estimates errors by a factor of  sqrt beam_area.

    The cumulative (growth-curve) error is the quadrature sum of per-annulus
    flux errors, which is exact for independent annuli.

    Parameters
    ----------
    image : 2-D ndarray
        Image data in Jy/beam (radio) or any flux-per-pixel unit (optical).
    x0, y0 : float
        Centre pixel coordinates (x = column, y = row).
    pa_deg : float
        Position angle of the major axis in degrees, measured from the
        positive x-axis (column axis) towards the positive y-axis (row axis),
        consistent with the convention in morfometryka.
    q : float
        Axis ratio b/a (0 < q <~ 1).  q = 1 gives circular apertures.
    beam_area : float
        Beam area in pixels (output of ``beam_area2``).  Set to 1 for optical
        data where one pixel is one resolution element.
    rms : float, optional
        Global noise estimate.  If None, estimated via ``mad_std`` on finite
        pixel values.  Used as a fallback std when an annulus is too sparsely
        populated for sigma-clipping.
    delta_r : float, optional
        Semi-major-axis step in pixels (default 1.0).
    sigma_clip_val : float, optional
        Rejection threshold for sigma-clipping within each annulus (default 3).
    wedge_half_angle_deg : float, optional
        Half-angular width of each directional wedge in degrees (default 22.5^o,
        giving a 45^o-wide wedge symmetric about each axis including its
        antipode).
    mask : 2-D bool ndarray, optional
        External validity mask (True = include pixel).  Combined with a
        ``np.isfinite`` check.  If None, all finite pixels are used.
    beam_correction : bool, optional
        If True (default) use the beam-correlated flux-error formula
        ``sigma x  sqrt N_eff``.  Set to False to reproduce the naive independent-pixel
        formula ``sigma x  sqrt N_pix / beam_area`` used in the old code.

    Returns
    -------
    dict
        Keys listed below.  Arrays that depend on the semi-major-axis grid
        (``sma``) all have the same length *before* duplicate removal.
        ``fluxes``, ``fluxes_err``, ``areas``, ``agrow``, ``IR*``, ``radii_*``,
        and ``sma`` all share the same length N (= len(sma)).
        ``Lgrow*``, ``Lgrow_norm*``, and ``radii`` (main) may be shorter if
        duplicate flux values were merged.

        Drop-in equivalents
        ~~~~~~~~~~~~~~~~~~~
        fluxes          per-annulus flux density [Jy]
        fluxes_err      per-annulus flux uncertainty
        areas           number of pixels per annulus
        agrow           copy of areas (alias, for backward compatibility)
        Lgrow           cumulative flux growth curve
        Lgrow_err       cumulative flux uncertainty (quadrature)
        Lgrow_norm      Lgrow / total_flux
        Lgrow_err_norm  Lgrow_err / total_flux
        radii           area-equivalent circular radius (= radii_circ,
                        after duplicate removal; matches level-based radii)

        Semi-major-axis grid
        ~~~~~~~~~~~~~~~~~~~~
        sma             semi-major axis values used [pixels]

        Four radius representations  (length N, *not* deduplicated)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        radii_major     a  [pixels]
        radii_minor     a x q
        radii_diag      a x qx sqrt 2 /  sqrt (q^2+1)
        radii_circ      a x  sqrt q  (same as radii before dedup)

        Full-annulus intensity profiles
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IR              sigma-clipped mean intensity per annulus
        IR_err          uncertainty on IR

        Directional intensity profiles  (major / minor / diagonal)
        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        IR_major, IR_major_err
        IR_minor, IR_minor_err
        IR_diag,  IR_diag_err

        Directional growth curves
        ~~~~~~~~~~~~~~~~~~~~~~~~~~
        Lgrow_major,  Lgrow_major_err,  Lgrow_major_norm
        Lgrow_minor,  Lgrow_minor_err,  Lgrow_minor_norm
        Lgrow_diag,   Lgrow_diag_err,   Lgrow_diag_norm

    Notes
    -----
    The directional growth curves integrate *only* the flux inside each
    directional wedge, so their normalised values represent how the emission
    is distributed along that axis, not a global fraction of total flux.
    To obtain a physically meaningful fraction of total flux use the main
    ``Lgrow_norm``.
    """
    from astropy.stats import sigma_clipped_stats, mad_std as _mad_std

    g = np.asarray(image, dtype=float)
    ny, nx = g.shape

    # ------------------------------------------------------------------ noise
    if rms is None:
        finite_vals = g[np.isfinite(g)]
        rms = _mad_std(finite_vals) if finite_vals.size > 0 else 1.0

    # ----------------------------------------------------------- validity mask
    if mask is None:
        valid = np.isfinite(g)
    else:
        valid = np.asarray(mask, dtype=bool) & np.isfinite(g)

    g_work = np.where(valid, g, np.nan)

    # ------------------------------------------------- elliptical coordinates
    xx, yy = np.meshgrid(np.arange(nx, dtype=float),
                         np.arange(ny, dtype=float))

    pa_rad = np.deg2rad(pa_deg)
    cos_pa, sin_pa = np.cos(pa_rad), np.sin(pa_rad)

    xl =  (xx - x0) * cos_pa + (yy - y0) * sin_pa   # along major axis
    yl = -(xx - x0) * sin_pa + (yy - y0) * cos_pa   # along minor axis

    # elliptical radius = semi-major axis equivalent
    rl = np.sqrt(xl ** 2 + (yl / q) ** 2)

    # azimuthal angle in the de-projected frame (φ=0 → major, φ=pi/2 → minor)
    phi = np.arctan2(yl / q, xl)                      # [-pi, pi]

    # ----------------------------------------------- semi-major-axis grid
    r_max = np.hypot(ny, nx)
    sma   = np.arange(delta_r, r_max, delta_r)        # semi-major axis values
    dee   = 0.5 * np.gradient(sma)                    # half-annulus width
    N     = len(sma)

    # ------------------------------------------------- pre-compute wedge maps
    # Each wedge covers both sides of the axis (the annulus half on each side),
    # giving full 360^o sampling split into three overlapping sectors.
    dphi = np.deg2rad(wedge_half_angle_deg)

    def _angular_distance(phi_arr, phi_center):
        """Unsigned angular distance, wrapped to [0, pi]."""
        d = phi_arr - phi_center
        return np.abs(np.arctan2(np.sin(d), np.cos(d)))

    # True where φ is within dphi of φ_center *or* its antipode (φ_center +/- pi)
    def _wedge(phi_center):
        d      = _angular_distance(phi, phi_center)
        d_anti = _angular_distance(phi, phi_center + np.pi)
        return (d < dphi) | (d_anti < dphi)

    wedge_major = _wedge(0.0)          # along major axis (φ=0, 180^o)
    wedge_minor = _wedge(np.pi / 2.0)  # along minor axis (φ=90, 270^o)
    wedge_diag  = _wedge(np.pi / 4.0)  # diagonal         (φ=45, 225^o)

    # -------------------------------------------- per-annulus statistics helper
    def _annulus_stats(pix_vals):
        """
        Sigma-clipped mean, std and derived flux quantities for a pixel
        sample from one annulus (or wedge).

        Returns (mean, mean_err, flux, flux_err, n_pix).
        All values are NaN when fewer than 1 valid pixel is present.
        """
        finite = pix_vals[np.isfinite(pix_vals)]
        n_pix  = finite.size

        if n_pix == 0:
            return np.nan, np.nan, np.nan, np.nan, 0

        if n_pix < 3:
            mean  = float(np.nanmean(finite))
            std   = rms
        else:
            try:
                mean, _, std = sigma_clipped_stats(finite,
                                                   sigma=sigma_clip_val,
                                                   maxiters=5)
                if not np.isfinite(std) or std <= 0.0:
                    std = rms
            except Exception:
                mean = float(np.nanmean(finite))
                std  = rms

        n_eff    = max(n_pix / beam_area, 1.0)   # independent resolution elements
        mean_err = std / np.sqrt(n_eff)           # uncertainty on the mean

        flux     = float(np.nansum(finite)) / beam_area  # Jy (or flux-unit)

        # Beam-corrected flux error (recommended):
        #   flux_err = sigma x  sqrt N_eff  (independent beams add in quadrature)
        # Naive formula (old code, over-estimates by  sqrt beam_area):
        #   flux_err = sigma x  sqrt N_pix / beam_area
        if beam_correction:
            flux_err = std * np.sqrt(n_eff)
        else:
            flux_err = std * np.sqrt(max(n_pix, 1)) / beam_area

        return mean, mean_err, flux, flux_err, n_pix

    # ---------------------------------------------------- allocate output arrays
    IR            = np.full(N, np.nan)
    IR_err        = np.full(N, np.nan)
    IR_major      = np.full(N, np.nan);  IR_major_err = np.full(N, np.nan)
    IR_minor      = np.full(N, np.nan);  IR_minor_err = np.full(N, np.nan)
    IR_diag       = np.full(N, np.nan);  IR_diag_err  = np.full(N, np.nan)

    fluxes        = np.full(N, np.nan)
    fluxes_err    = np.full(N, np.nan)
    areas         = np.zeros(N, dtype=float)

    fluxes_major  = np.full(N, np.nan);  fluxes_major_err = np.full(N, np.nan)
    fluxes_minor  = np.full(N, np.nan);  fluxes_minor_err = np.full(N, np.nan)
    fluxes_diag   = np.full(N, np.nan);  fluxes_diag_err  = np.full(N, np.nan)

    # ----------------------------------------------------------- radial loop
    for i, a in enumerate(sma):
        d   = dee[i]
        ann = valid & (rl > (a - d)) & (rl <= (a + d))

        m, m_err, f, f_err, n = _annulus_stats(g_work[ann])
        IR[i]       = m
        IR_err[i]   = m_err
        fluxes[i]   = f
        fluxes_err[i] = f_err
        areas[i]    = n

        # directional wedge profiles - pixels must also belong to the annulus
        for wmask, ir_out, ir_err_out, fx_out, fx_err_out in (
            (wedge_major, IR_major, IR_major_err, fluxes_major, fluxes_major_err),
            (wedge_minor, IR_minor, IR_minor_err, fluxes_minor, fluxes_minor_err),
            (wedge_diag,  IR_diag,  IR_diag_err,  fluxes_diag,  fluxes_diag_err),
        ):
            wann = ann & wmask
            wm, wm_err, wf, wf_err, _ = _annulus_stats(g_work[wann])
            ir_out[i]     = wm
            ir_err_out[i] = wm_err
            fx_out[i]     = wf
            fx_err_out[i] = wf_err

    # ------------------------------------------------- four radius conventions
    # At azimuthal angle φ on the ellipse at semi-major axis a, the Euclidean
    # distance from the centre to the ellipse boundary is:
    #   r(a, φ) = a x q /  sqrt (q^2xcos^2φ + sin^2φ)
    # which gives:
    #   φ = 0   → r = a                               (major axis)
    #   φ = pi/2 → r = a x q                           (minor axis)
    #   φ = pi/4 → r = a x qx sqrt 2 /  sqrt (q^2+1)             (diagonal)
    # Area-equivalent circle (ellipse area = pixa^2xq → r_circ = ax sqrt q):
    radii_major = sma.copy()
    radii_minor = sma * q
    radii_diag  = sma * q * np.sqrt(2.0) / np.sqrt(q ** 2 + 1.0)
    radii_circ  = sma * np.sqrt(q)

    # -------------------------------------------- cumulative growth curves
    def _growth_curve(fx, fx_err):
        """Build cumulative flux, error (quadrature), and normalised versions."""
        Lg       = np.nancumsum(fx)
        Lg_err   = np.sqrt(np.nancumsum(np.where(np.isfinite(fx_err),
                                                  fx_err ** 2, 0.0)))
        total    = np.nansum(fx)
        if total > 0:
            Lg_norm     = Lg / total
            Lg_norm_err = Lg_err / total
        else:
            Lg_norm     = np.full_like(Lg, np.nan)
            Lg_norm_err = np.full_like(Lg_err, np.nan)
        return Lg, Lg_err, Lg_norm, Lg_norm_err

    Lgrow,       Lgrow_err,       Lgrow_norm,       Lgrow_err_norm       = _growth_curve(fluxes,       fluxes_err)
    Lgrow_major, Lgrow_major_err, Lgrow_major_norm, _                    = _growth_curve(fluxes_major, fluxes_major_err)
    Lgrow_minor, Lgrow_minor_err, Lgrow_minor_norm, _                    = _growth_curve(fluxes_minor, fluxes_minor_err)
    Lgrow_diag,  Lgrow_diag_err,  Lgrow_diag_norm,  _                    = _growth_curve(fluxes_diag,  fluxes_diag_err)

    # ---- duplicate removal on the main growth curve (matches existing code) ---
    agrow = areas.copy()
    Lgrow, radii_main = check_flux_duplicates(Lgrow, radii_circ.copy())
    n_dedup       = len(Lgrow)
    Lgrow_err     = Lgrow_err[:n_dedup]
    Lgrow_norm    = Lgrow_norm[:n_dedup]
    Lgrow_err_norm = Lgrow_err_norm[:n_dedup]

    # ------------------------------------------------------------------ output
    return dict(
        # --- drop-in equivalents (same names as level-based code) ---
        fluxes          = fluxes,
        fluxes_err      = fluxes_err,
        areas           = areas,
        agrow           = agrow,
        Lgrow           = Lgrow,
        Lgrow_err       = Lgrow_err,
        Lgrow_norm      = Lgrow_norm,
        Lgrow_err_norm  = Lgrow_err_norm,
        radii           = radii_main,    # area-equivalent; length may differ from N
        # --- semi-major-axis grid (length N) ---
        sma             = sma,
        # --- four radius representations (length N, not deduplicated) ---
        radii_major     = radii_major,
        radii_minor     = radii_minor,
        radii_diag      = radii_diag,
        radii_circ      = radii_circ,    # same as radii before dedup
        # --- full-annulus intensity profile ---
        IR              = IR,
        IR_err          = IR_err,
        # --- directional intensity profiles ---
        IR_major        = IR_major,
        IR_major_err    = IR_major_err,
        IR_minor        = IR_minor,
        IR_minor_err    = IR_minor_err,
        IR_diag         = IR_diag,
        IR_diag_err     = IR_diag_err,
        # --- directional growth curves (length N, not deduplicated) ---
        Lgrow_major      = Lgrow_major,
        Lgrow_major_err  = Lgrow_major_err,
        Lgrow_major_norm = Lgrow_major_norm,
        Lgrow_minor      = Lgrow_minor,
        Lgrow_minor_err  = Lgrow_minor_err,
        Lgrow_minor_norm = Lgrow_minor_norm,
        Lgrow_diag       = Lgrow_diag,
        Lgrow_diag_err   = Lgrow_diag_err,
        Lgrow_diag_norm  = Lgrow_diag_norm,
    )



def test_ned_connectivity(timeout=10, verbose=True):
    """
    Quick connectivity test for the NED service.

    Queries a well-known source (M31) with a short timeout and reports
    latency and status. Useful for diagnosing network issues before
    running batch queries.

    Parameters
    ----------
    timeout : int, optional
        Seconds to wait before declaring failure. Default: 10
    verbose : bool, optional
        Print a human-readable status line. Default: True

    Returns
    -------
    reachable : bool
        True if NED responded successfully.
    latency : float or None
        Round-trip time in seconds, or None on failure.
    """
    import time
    from astroquery.ipac.ned import Ned
    from requests.exceptions import ReadTimeout, ConnectionError as ReqConnectionError

    probe_source = 'M31'
    Ned.TIMEOUT = timeout

    t0 = time.perf_counter()
    try:
        result = Ned.query_object(probe_source)
        latency = time.perf_counter() - t0
        z = result['Redshift'].data.data[0]
        if verbose:
            print(f"[+] NED reachable  |  latency={latency:.2f}s  |  "
                  f"probe={probe_source}  z={z:.4f}")
        return True, latency

    except (ReadTimeout, TimeoutError):
        latency = time.perf_counter() - t0
        if verbose:
            print(f"[!] NED TIMEOUT after {latency:.1f}s  "
                  f"(host=ned.ipac.caltech.edu, timeout={timeout}s)")
        return False, None

    except ReqConnectionError as e:
        if verbose:
            print(f"[!] NED CONNECTION ERROR: {e}")
        return False, None

    except Exception as e:
        if verbose:
            print(f"[!] NED unexpected error: {e}")
        return False, None


