"""
#Background and Source Extraction
"""

# def adaptive_box_size(data_shape, min_boxes=20):
#     """Determine appropriate box size based on image dimensions."""
#     # Aim for at least min_boxes in each dimension
#     h, w = data_shape
#     bh = max(32, h // min_boxes)
#     bw = max(32, w // min_boxes)
#     # Ensure box size is reasonable
#     return bw, bh

# def improved_masking(data, nsigma=6, min_area=100, iterations=2):
#     """Create a more robust source mask using sigma-clipping and segmentation."""
#     from photutils.segmentation import detect_sources
#     from astropy.stats import sigma_clipped_stats
    
#     # Estimate background using sigma-clipping
#     mean, median, std = sigma_clipped_stats(data, sigma=3.0)
#     threshold = median + (nsigma * std)
    
#     # Detect sources
#     segm = detect_sources(data, threshold, npixels=min_area)
#     if segm is None:
#         return np.zeros(data.shape, dtype=bool)
        
#     # Create initial mask
#     mask = segm.data > 0
    
#     # Dilate mask if requested
#     if iterations > 0:
#         from scipy.ndimage import binary_dilation
#         mask = binary_dilation(mask, iterations=iterations)
        
#     return mask


# def multi_method_background(data, mask=None):
#     """Estimate background using multiple methods and compare results."""
#     from astropy.stats import SigmaClip, biweight_location
#     from photutils import Background2D, MedianBackground, SExtractorBackground
    
#     sigma_clip = SigmaClip(sigma=3.0)
    
#     # Method 1: SEP (your current method)
#     import sep
#     bkg_sep = sep.Background(data, mask=mask, bw=64, bh=64, fw=3, fh=3)
#     bkg_sep_image = bkg_sep.back()
    
#     # Method 2: photutils with median background
#     bkg_median = Background2D(data, box_size=64, filter_size=3, 
#                              sigma_clip=sigma_clip, mask=mask,
#                              bkg_estimator=MedianBackground())
    
#     # Method 3: photutils with SExtractor-like background
#     bkg_sextractor = Background2D(data, box_size=64, filter_size=3,
#                                  sigma_clip=sigma_clip, mask=mask,
#                                  bkg_estimator=SExtractorBackground())
    
#     # Compare results
#     methods = {
#         'SEP': bkg_sep_image,
#         'Median': bkg_median.background,
#         'SExtractor': bkg_sextractor.background
#     }
    
#     # Calculate agreement metrics
#     agreements = {}
#     for name1, bkg1 in methods.items():
#         for name2, bkg2 in methods.items():
#             if name1 >= name2:
#                 continue
#             diff = np.abs(bkg1 - bkg2)
#             agreements[f"{name1}-{name2}"] = {
#                 'mean_diff': np.mean(diff),
#                 'max_diff': np.max(diff),
#                 'relative_diff': np.mean(diff) / np.mean(bkg1)
#             }
    
#     return methods, agreements


# def visualize_background_mesh(data, box_size=64):
#     """Visualize data averaged in background estimation boxes."""
#     import matplotlib.pyplot as plt
#     from matplotlib.patches import Rectangle
    
#     h, w = data.shape
#     num_boxes_y = h // box_size
#     num_boxes_x = w // box_size
    
#     mesh_values = np.zeros((num_boxes_y, num_boxes_x))
    
#     for i in range(num_boxes_y):
#         for j in range(num_boxes_x):
#             y_slice = slice(i*box_size, min((i+1)*box_size, h))
#             x_slice = slice(j*box_size, min((j+1)*box_size, w))
#             mesh_values[i, j] = np.median(data[y_slice, x_slice])
    
#     fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
#     # Original image
#     im1 = ax1.imshow(data, origin='lower', cmap='viridis')
#     plt.colorbar(im1, ax=ax1)
#     ax1.set_title('Original Image')
    
#     # Draw box grid
#     for i in range(num_boxes_y):
#         for j in range(num_boxes_x):
#             rect = Rectangle((j*box_size, i*box_size), box_size, box_size, 
#                            fill=False, edgecolor='white', linewidth=0.5)
#             ax1.add_patch(rect)
    
#     # Mesh values
#     im2 = ax2.imshow(mesh_values, origin='lower', cmap='viridis', 
#                     extent=[0, w, 0, h])
#     plt.colorbar(im2, ax=ax2)
#     ax2.set_title('Background Mesh Values')
    
#     plt.tight_layout()
#     return fig


# def improved_background(imagename, method='sep', mask=None, apply_mask=True, 
#                        show_map=False, box_size=64, filter_size=3):
#     """
#     Enhanced background estimation with multiple methods.
    
#     Parameters
#     ----------
#     imagename : str
#         Path to the image.
#     method : str
#         Method to use ('sep', 'photutils', 'iterative', 'wavelet').
#     mask : array, optional
#         Mask to be applied to the image.
#     apply_mask : bool
#         If True, calculate the mask from the image.
#     show_map : bool
#         If True, show the background map.
#     box_size : int or tuple
#         Box size for background estimation.
#     filter_size : int or tuple
#         Filter size for smoothing the background.
        
#     Returns
#     -------
#     background : ndarray
#         Estimated background.
#     background_rms : ndarray
#         Estimated background RMS.
#     """
#     # Read data
#     import fitsio
#     data = fitsio.read(imagename)
#     if len(data.shape) == 4:
#         data = data[0][0]
    
#     # Handle NaNs
#     data_no_nan = np.copy(data)
#     nan_mask = np.isnan(data)
#     if nan_mask.any():
#         data_no_nan[nan_mask] = 0
    
#     # Create or apply mask
#     if mask is None and apply_mask:
#         mask = improved_masking(data_no_nan)
    
#     # Apply selected method
#     if method == 'sep':
#         import sep
#         if isinstance(box_size, int):
#             bw = bh = box_size
#         else:
#             bw, bh = box_size
            
#         if isinstance(filter_size, int):
#             fw = fh = filter_size
#         else:
#             fw, fh = filter_size
            
#         bkg = sep.Background(data_no_nan, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
#         background = bkg.back()
#         background_rms = bkg.rms()
        
#     elif method == 'photutils':
#         from photutils import Background2D, MedianBackground
#         from astropy.stats import SigmaClip
        
#         sigma_clip = SigmaClip(sigma=3.0)
#         bkg = Background2D(data_no_nan, box_size, filter_size=filter_size,
#                           sigma_clip=sigma_clip, bkg_estimator=MedianBackground(),
#                           mask=mask, edge_method='pad')
#         background = bkg.background
#         background_rms = bkg.background_rms
        
#     elif method == 'iterative':
#         background, background_rms = iterative_background(data_no_nan)
        
#     elif method == 'wavelet':
#         background = wavelet_background(data_no_nan)
#         # For RMS, use local standard deviation
#         from scipy.ndimage import uniform_filter
#         from astropy.stats import sigma_clipped_stats
        
#         # Compute local standard deviation
#         mean, median, std = sigma_clipped_stats(data_no_nan - background, sigma=3.0)
#         background_rms = np.ones_like(background) * std
    
#     # Restore NaN values
#     if nan_mask.any():
#         background[nan_mask] = np.nan
#         background_rms[nan_mask] = np.nan
    
#     # Visualization
#     if show_map:
#         import matplotlib.pyplot as plt
        
#         fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(15, 5))
        
#         # Original data
#         im1 = ax1.imshow(np.arcsinh(data), origin='lower',cmap='magma_r')
#         plt.colorbar(im1, ax=ax1)
#         ax1.set_title('Original Data')
        
#         # Background
#         im2 = ax2.imshow(np.arcsinh(background), origin='lower',cmap='magma_r')
#         plt.colorbar(im2, ax=ax2)
#         ax2.set_title(f'Background ({method})')
        
#         # Background-subtracted
#         im3 = ax3.imshow(np.arcsinh(data - background), origin='lower',cmap='magma_r')
#         plt.colorbar(im3, ax=ax3)
#         ax3.set_title('Background Subtracted')
        
#         plt.tight_layout()
#         plt.show()
    
#     return background, background_rms




# def photutils_background(imagename, box_size=64, filter_size=3, mask=None):
#     """
#     Estimate background using photutils Background2D.
#     """
#     from astropy.io import fits
#     from astropy.stats import SigmaClip
#     from photutils import Background2D, MedianBackground
    
#     # Read data
#     data = fits.getdata(imagename)
#     if len(data.shape) == 4:
#         data = data[0][0]
    
#     # Create mask if not provided
#     if mask is None:
#         from astropy.stats import sigma_clip
#         from scipy.ndimage import binary_dilation
        
#         clipped_data = sigma_clip(data, sigma=3)
#         mask = clipped_data.mask
#         if mask.sum() > 0:  # Only dilate if mask is not empty
#             mask = binary_dilation(mask, iterations=2)
    
#     # Sigma clipping for outlier rejection
#     sigma_clip = SigmaClip(sigma=3.0)
    
#     # Estimate background
#     bkg_estimator = MedianBackground()
#     bkg = Background2D(data, box_size, filter_size=filter_size, 
#                       sigma_clip=sigma_clip, bkg_estimator=bkg_estimator,
#                       mask=mask, edge_method='pad')
    
#     return bkg


# def iterative_background(data, n_iterations=5, sigma=3.0):
#     """
#     Estimate background through iterative sigma clipping.
#     """
#     from astropy.stats import sigma_clipped_stats
    
#     # Create working copy
#     working_data = data.copy()
    
#     for i in range(n_iterations):
#         # Calculate statistics with sigma clipping
#         mean, median, std = sigma_clipped_stats(working_data, sigma=sigma)
        
#         # Create mask for values significantly above background
#         mask = working_data > (median + sigma * std)
        
#         # Replace masked values with median for next iteration
#         working_data[mask] = median
    
#     # Final background estimate
#     background = np.ones_like(data) * median
#     background_rms = np.ones_like(data) * std
    
#     return background, background_rms

# def wavelet_background(data, level=4):
#     """
#     Estimate background using wavelet decomposition.
#     """
#     import pywt
    
#     # Perform wavelet decomposition
#     coeffs = pywt.wavedec2(data, 'haar', level=level)
    
#     # Extract approximation coefficients (lowest frequency)
#     cA = coeffs[0]
    
#     # Reconstruct using only approximation coefficients
#     new_coeffs = [cA] + [None] * level
#     background = pywt.waverec2(new_coeffs, 'haar')
    
#     # Ensure same shape as input
#     background = background[:data.shape[0], :data.shape[1]]
    
#     return background








################

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
                              mask=mask, edge_method='pad')
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

def make_hst_cutout(downloaded_files, coord, size_arcsec, output_filename=None, source_name=None, band=None):
    """
    Create a cutout from downloaded HST FITS files with complete header preservation.
    
    Parameters
    ----------
    downloaded_files : list
        List of paths to downloaded HST FITS files.
    coord : SkyCoord
        Center coordinates for the cutout.
    size_arcsec : float
        Size of the cutout in arcseconds.
    output_filename : str, optional
        Name for the output file. If None, will generate based on source_name and band.
    source_name : str, optional
        Name of the source (used for automatic filename generation).
    band : str, optional
        Filter band (used for automatic filename generation).
        
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
                    print(f"[✓] Using science data from extension {sci_ext} ({hdul[sci_ext].name}) in {file_path}")
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
                
            print(f"[✓] Pixel scale: {pixscale:.4f} arcsec/pixel")
            
            # Compute cutout size in pixels
            size_pixels = int(size_arcsec / pixscale)
            print(f"[✓] Cutout size: {size_pixels}x{size_pixels} pixels")
            
            # Make the cutout
            cutout = Cutout2D(sci_data, position=coord, size=(size_pixels, size_pixels), wcs=wcs)
            
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
                            ext_cutout = Cutout2D(hdu.data, position=coord, 
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
            print(f"[✓] Cutout saved as '{output_filename}' with {len(cutout_hdul)} extensions")
            return output_filename
            
    except Exception as e:
        print(f"[Error] Failed to create cutout: {e}")
        import traceback
        traceback.print_exc()
        return None


# def find_z_NED(source_name):
#     result_table = Ned.query_object(source_name)
#     redshift_NED = result_table['Redshift'].data.data
#     if redshift_NED.shape[0] == 0:
#         return None
#     else:
#         return redshift_NED[0]

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
    print(f"[✓] Downloaded {len(downloaded_files)} FITS files.")

    return downloaded_files, coord, band, instrument


import requests
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

# def int_Lum(z,S_nu0_sy,S_nu0_ff, 
#             S_nu0_sy_err=None,S_nu0_ff_err=None, 
#             nu_0=6e9, 
#             nu_t=0.01e9,nu_1=1.0e9,nu_2=35.0e9,alpha_sy=-0.85):
#     dist_conversion_factor = 3.08567758128 * 1e24  # m #3.08567758128*(10**24) #cm
#     lum_conversion_factor = 1e-23
#     DL_Mpc = luminosity_distance_cosmo(z=z)
#     DL_cm = DL_Mpc * dist_conversion_factor  # Convert Mpc to cm
    
#     # synchrotron_luminosity, synchrotron_luminosity_err = \
#     #     quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy),
#     #          epsrel=1e-12, limit=10000)
#     # freefree_luminosity, freefree_luminosity_err = \
#     #     quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff),
#     #          epsrel=1e-12, limit=10000)
#     synchrotron_luminosity, _ = \
#         quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy),
#              epsrel=1e-6, limit=10000
#             )
#     freefree_luminosity, _ = \
#         quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff),
#              epsrel=1e-6, limit=10000
#             )

#     synchrotron_luminosity_err = 0.0
#     freefree_luminosity_err = 0.0
#     int_luminosity_err = 0.0
#     LR_err = 0.0
#     if S_nu0_sy_err is not None:
#         synchrotron_luminosity_err, _ = \
#             quad(synchrotron_integrand, nu_1, nu_2, args=(nu_0, nu_t, alpha_sy, S_nu0_sy_err),
#                  epsrel=1e-6, limit=10000
#                 )
#         freefree_luminosity_err, _ = \
#             quad(freefree_integrand, nu_1, nu_2, args=(nu_0, nu_t, S_nu0_ff_err),
#                  epsrel=1e-6, limit=10000
#                 )
        
#         # print(freefree_luminosity)
#         # print(synchrotron_luminosity)
#         int_luminosity_err = np.sqrt(synchrotron_luminosity_err**2 + freefree_luminosity_err**2.0)
#         LR_err  = (4 * np.pi * (DL_cm**2) / ((1+z)**(alpha_sy+1))) * (int_luminosity_err) * lum_conversion_factor

    
#     LR  = (4 * np.pi * (DL_cm**2) / ((1+z)**(alpha_sy+1))) * (synchrotron_luminosity + freefree_luminosity) * lum_conversion_factor
#     return(LR,LR_err)


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
# def do_fit_spec_SY_FF_map(freqs,fluxes,fluxes_err,nu0=None,
#                           fix_alpha_nt=False,
#                           verbose=0):
#     x = freqs
#     y = fluxes
#     yerr = fluxes_err
#     if nu0 is None:
#         nu0 = np.mean(x)

#     epsilon = 1e-8
    
#     def min_func(params):
#         A_sy = params['A_sy']
#         A_ff = params['A_ff']
#         alpha_nt = params['alpha_nt']
#         model = RC_function_SY_FF(x, A_sy, A_ff, alpha_nt,nu0)
#         log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
#         return (y - model) * log_weights



#     fit_params = lmfit.Parameters()
#     fit_params.add("A_sy", value=1.0, min=1.0e-6, max=1000)
#     fit_params.add("A_ff", value=0.1, min=1.0e-6, max=100)
#     if fix_alpha_nt == True:
#         fit_params.add("alpha_nt", value=-0.85, min=-2.0, max=0.0, vary=False)
#     else:
#         fit_params.add("alpha_nt", value=-0.85, min=-3.0, max=3.0)
    

#     mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
#                            nan_policy='omit', reduce_fcn='neglogcauchy')

#     result_1 = mini.minimize(method='least_squares',
#                              max_nfev=200000,  # f_scale = 1.0,
#                              loss="cauchy",
#                              tr_solver="exact",
#                              ftol=1e-12, xtol=1e-12, gtol=1e-12,
#                              verbose=verbose
#                              )
#     second_run_params = result_1.params

#     result = mini.minimize(method='least_squares',
#                            params=second_run_params,
#                            max_nfev=200000,  # f_scale = 1.0,
#                         #    loss="huber", 
#                            loss="cauchy",
#                            tr_solver="exact",
#                            ftol=1e-12, xtol=1e-12, gtol=1e-12,
#                            verbose=verbose
#                            )
    
#     return result


# import numpy as np
# from astropy.modeling import Fittable1DModel, Parameter
# from astropy.modeling.fitting import LevMarLSQFitter, FittingWithOutlierRemoval
# from astropy.stats import sigma_clip
# from scipy.optimize import least_squares
# import warnings


# class RC_SY_FF_Model(Fittable1DModel):
#     """
#     Custom astropy model for RC_function_SY_FF fitting
#     """
#     A_sy = Parameter(default=1.0, bounds=(1.0e-6, 1000))
#     A_ff = Parameter(default=0.1, bounds=(1.0e-6, 100))
#     alpha_nt = Parameter(default=-0.85, bounds=(-3.0, 3.0))
    
#     def __init__(self, nu0=None, fix_alpha_nt=False, **kwargs):
#         super().__init__(**kwargs)
#         self.nu0 = nu0
#         if fix_alpha_nt:
#             self.alpha_nt.fixed = True
#             # Adjust bounds for fixed case
#             if self.alpha_nt.value < -2.0 or self.alpha_nt.value > 0.0:
#                 self.alpha_nt = -0.85
    
#     def evaluate(self, x, A_sy, A_ff, alpha_nt):
#         # This should call your RC_function_SY_FF
#         return RC_function_SY_FF(x, A_sy, A_ff, alpha_nt, self.nu0)


# class CauchyLossFitter:
#     """
#     Custom robust fitter using Cauchy loss with proper error calculation
#     """
#     def __init__(self, max_nfev=200000, ftol=1e-12, xtol=1e-12, gtol=1e-12):
#         self.max_nfev = max_nfev
#         self.ftol = ftol
#         self.xtol = xtol
#         self.gtol = gtol
    
#     def __call__(self, model, x, y, weights=None, **kwargs):
#         # Get parameter info
#         param_names = model.param_names
#         param_values = [getattr(model, name).value for name in param_names]
#         param_bounds = [(getattr(model, name).bounds[0] if getattr(model, name).bounds else -np.inf,
#                         getattr(model, name).bounds[1] if getattr(model, name).bounds else np.inf) 
#                        for name in param_names]
#         fixed_params = {name: getattr(model, name).fixed for name in param_names}
        
#         # Create arrays for free parameters only
#         free_param_names = [name for name in param_names if not fixed_params[name]]
#         free_param_values = [getattr(model, name).value for name in free_param_names]
#         free_param_bounds = [(getattr(model, name).bounds[0] if getattr(model, name).bounds else -np.inf,
#                              getattr(model, name).bounds[1] if getattr(model, name).bounds else np.inf) 
#                             for name in free_param_names]
        
#         # Transpose bounds for scipy format
#         bounds = list(zip(*free_param_bounds)) if free_param_bounds else ([], [])
        
#         def residual_func(free_params):
#             # Reconstruct full parameter set
#             full_params = {}
#             free_idx = 0
#             for name in param_names:
#                 if fixed_params[name]:
#                     full_params[name] = getattr(model, name).value
#                 else:
#                     full_params[name] = free_params[free_idx]
#                     free_idx += 1
            
#             # Evaluate model
#             model_vals = model.evaluate(x, **full_params)
#             residuals = y - model_vals
            
#             # Apply weights
#             if weights is not None:
#                 residuals *= weights
                
#             return residuals
        
#         # Two-stage fitting as in original
#         # First stage
#         result_1 = least_squares(
#             residual_func,
#             free_param_values,
#             bounds=bounds,
#             method='trf',
#             loss='cauchy',
#             ftol=self.ftol,
#             xtol=self.xtol,
#             gtol=self.gtol,
#             max_nfev=self.max_nfev
#         )
        
#         # Second stage
#         result = least_squares(
#             residual_func,
#             result_1.x,
#             bounds=bounds,
#             method='trf',
#             loss='cauchy',
#             ftol=self.ftol,
#             xtol=self.xtol,
#             gtol=self.gtol,
#             max_nfev=self.max_nfev
#         )
        
#         # Update model with fitted parameters
#         free_idx = 0
#         for name in param_names:
#             if not fixed_params[name]:
#                 setattr(model, name, result.x[free_idx])
#                 free_idx += 1
        
#         return result, model, free_param_names


# class LMFitCompatibleResult:
#     """
#     Result wrapper that mimics lmfit interface with proper error calculation
#     """
#     def __init__(self, scipy_result, fitted_model, free_param_names, all_param_names, 
#                  fixed_params, x, y, weights=None):
        
#         self.success = scipy_result.success
#         self.message = scipy_result.message
#         self.nfev = scipy_result.nfev
#         self.cost = scipy_result.cost
#         self.optimality = scipy_result.optimality
        
#         # Calculate chi-square
#         residuals = scipy_result.fun
#         self.chisqr = np.sum(residuals**2)
        
#         # Calculate degrees of freedom
#         n_data = len(y)
#         n_free_params = len(free_param_names)
#         self.nfree = n_data - n_free_params
#         self.redchi = self.chisqr / max(self.nfree, 1)
        
#         # Create parameter object
#         self.params = LMFitParams()
        
#         # Calculate parameter uncertainties from covariance matrix
#         param_errors = {}
#         if scipy_result.jac is not None and len(free_param_names) > 0:
#             try:
#                 # Compute covariance matrix
#                 # For weighted least squares with weights w, covariance = inv(J^T W J)
#                 # where W is diagonal weight matrix
#                 jac = scipy_result.jac
                
#                 # Apply weights to Jacobian if present
#                 if weights is not None:
#                     # Expand weights to match residual dimensions if needed
#                     weight_matrix = np.diag(weights)
#                     jac_weighted = weight_matrix @ jac
#                 else:
#                     jac_weighted = jac
                
#                 # Compute covariance: inv(J^T J) * residual_variance
#                 jtj = jac_weighted.T @ jac_weighted
                
#                 # Add small regularization to diagonal for numerical stability
#                 reg_factor = 1e-14 * np.trace(jtj) / len(jtj)
#                 jtj += reg_factor * np.eye(len(jtj))
                
#                 try:
#                     cov_matrix = np.linalg.inv(jtj)
#                     # Scale by residual variance
#                     residual_var = self.chisqr / max(self.nfree, 1)
#                     cov_matrix *= residual_var
                    
#                     param_std_errors = np.sqrt(np.diag(cov_matrix))
                    
#                     for i, param_name in enumerate(free_param_names):
#                         param_errors[param_name] = param_std_errors[i]
                        
#                 except np.linalg.LinAlgError:
#                     # Fallback: try pseudo-inverse
#                     try:
#                         cov_matrix = np.linalg.pinv(jtj)
#                         residual_var = self.chisqr / max(self.nfree, 1)
#                         cov_matrix *= residual_var
#                         param_std_errors = np.sqrt(np.abs(np.diag(cov_matrix)))
                        
#                         for i, param_name in enumerate(free_param_names):
#                             param_errors[param_name] = param_std_errors[i]
#                     except:
#                         # If all else fails, set errors to None
#                         pass
                        
#             except Exception as e:
#                 # If covariance calculation fails, errors remain None
#                 warnings.warn(f"Could not calculate parameter uncertainties: {e}")
        
#         # Populate parameters
#         for name in all_param_names:
#             param_value = getattr(fitted_model, name).value
#             param_bounds = getattr(fitted_model, name).bounds
#             is_fixed = fixed_params[name]
#             stderr = param_errors.get(name, None)
            
#             self.params.add(
#                 name, 
#                 value=param_value,
#                 min=param_bounds[0] if param_bounds else None,
#                 max=param_bounds[1] if param_bounds else None,
#                 vary=not is_fixed,
#                 stderr=stderr
#             )
        
#         # Store residual
#         self.residual = scipy_result.fun
#         self.method = 'least_squares'


# class LMFitParams:
#     """Parameter container mimicking lmfit.Parameters"""
#     def __init__(self):
#         self._params = {}
    
#     def add(self, name, value=None, min=None, max=None, vary=True, stderr=None):
#         self._params[name] = LMFitParam(value, min, max, vary, stderr)
    
#     def __getitem__(self, name):
#         return self._params[name]
    
#     def __contains__(self, name):
#         return name in self._params
    
#     def items(self):
#         return self._params.items()


# class LMFitParam:
#     """Parameter object mimicking lmfit.Parameter"""
#     def __init__(self, value, min_val, max_val, vary, stderr):
#         self.value = value
#         self.min = min_val
#         self.max = max_val
#         self.vary = vary
#         self.stderr = stderr


# def do_fit_spec_SY_FF_map(freqs, fluxes, fluxes_err, nu0=None,
#                           fix_alpha_nt=False, verbose=0):
#     """
#     Robust fitting function using astropy models with proper error calculation.
#     Maintains exact interface compatibility with original lmfit version.
#     """
#     x = freqs
#     y = fluxes
#     yerr = fluxes_err
    
#     if nu0 is None:
#         nu0 = np.mean(x)

#     # Create log weights exactly as in original
#     log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
    
#     # Create and configure model
#     model = RC_SY_FF_Model(nu0=nu0, fix_alpha_nt=fix_alpha_nt)
    
#     # Set initial values and bounds
#     model.A_sy = 1.0
#     model.A_sy.bounds = (1.0e-6, 1000)
    
#     model.A_ff = 0.1  
#     model.A_ff.bounds = (1.0e-6, 100)
    
#     model.alpha_nt = -0.85
    
#     if fix_alpha_nt:
#         model.alpha_nt.fixed = True
#         model.alpha_nt.bounds = (-0.855, -0.845)
#     else:
#         model.alpha_nt.bounds = (-3.0, -0.15)
    
#     # Get parameter information
#     param_names = model.param_names
#     fixed_params = {name: getattr(model, name).fixed for name in param_names}
    
#     # Perform robust fitting
#     fitter = CauchyLossFitter(max_nfev=200000, ftol=1e-12, xtol=1e-12, gtol=1e-12)
#     scipy_result, fitted_model, free_param_names = fitter(model, x, y, weights=log_weights)
    
#     # Create compatible result
#     result = LMFitCompatibleResult(
#         scipy_result, fitted_model, free_param_names, param_names, 
#         fixed_params, x, y, weights=log_weights
#     )
    
#     if verbose > 0:
#         print(f"Fit success: {result.success}")
#         print(f"Message: {result.message}")
#         print(f"Chi-square: {result.chisqr:.6f}")
#         print(f"Reduced chi-square: {result.redchi:.6f}")
#         print(f"Number of function evaluations: {result.nfev}")
#         print("Parameters:")
#         for name, param in result.params.items():
#             stderr_str = f"{param.stderr:.6f}" if param.stderr is not None else "N/A"
#             fixed_str = " (fixed)" if not param.vary else ""
#             print(f"  {name}: {param.value:.6f} +/- {stderr_str}{fixed_str}")
    
#     return result

























# import numpy as np
# from astropy.modeling import Fittable1DModel, Parameter
# from astropy.modeling.fitting import LevMarLSQFitter, FittingWithOutlierRemoval
# from astropy.stats import sigma_clip
# import warnings


# class RC_SY_FF_Model(Fittable1DModel):
#     """
#     Astropy model for RC_function_SY_FF fitting
#     """
#     A_sy = Parameter(default=1.0, bounds=(1.0e-6, 1000))
#     A_ff = Parameter(default=0.1, bounds=(1.0e-6, 100))
#     alpha_nt = Parameter(default=-0.85, bounds=(-3.0, 3.0))
    
#     def __init__(self, nu0=None, **kwargs):
#         super().__init__(**kwargs)
#         self.nu0 = nu0
    
#     def evaluate(self, x, A_sy, A_ff, alpha_nt):
#         return RC_function_SY_FF(x, A_sy, A_ff, alpha_nt, self.nu0)


# class LMFitCompatibleResult:
#     """
#     Wrapper to provide lmfit-compatible interface for astropy fitting results
#     """
#     def __init__(self, fitted_model, fitter, x, y, weights=None):
#         # Basic fit information
#         if hasattr(fitter, 'fit_info') and fitter.fit_info is not None:
#             fit_info = fitter.fit_info
#             self.success = True  # Astropy doesn't always provide explicit success flag
#             self.message = "Converged" if hasattr(fit_info, 'ierr') and fit_info.get('ierr', 0) <= 4 else "Unknown"
#             self.nfev = fit_info.get('nfev', 0)
#         else:
#             self.success = True
#             self.message = "Converged"
#             self.nfev = 0
        
#         # Calculate residuals and chi-square
#         model_values = fitted_model(x)
#         residuals = y - model_values
        
#         if weights is not None:
#             weighted_residuals = residuals * weights
#             self.chisqr = np.sum(weighted_residuals**2)
#         else:
#             self.chisqr = np.sum(residuals**2)
        
#         # Degrees of freedom
#         n_data = len(y)
#         n_free_params = len([p for p in fitted_model.parameters if not getattr(fitted_model, fitted_model.param_names[fitted_model.parameters.tolist().index(p)]).fixed])
#         self.nfree = n_data - n_free_params
#         self.redchi = self.chisqr / max(self.nfree, 1)
        
#         # Create parameter object
#         self.params = LMFitParams()
        
#         # Extract parameter uncertainties from covariance matrix
#         param_errors = {}
#         if hasattr(fitter, 'fit_info') and fitter.fit_info is not None:
#             cov_matrix = fitter.fit_info.get('param_cov', None)
#             if cov_matrix is not None:
#                 try:
#                     # Get indices of free parameters
#                     free_param_indices = []
#                     for i, param_name in enumerate(fitted_model.param_names):
#                         if not getattr(fitted_model, param_name).fixed:
#                             free_param_indices.append(i)
                    
#                     # Extract standard errors for free parameters
#                     if len(free_param_indices) > 0 and cov_matrix.shape[0] == len(free_param_indices):
#                         param_std_errors = np.sqrt(np.diag(cov_matrix))
#                         for idx, param_idx in enumerate(free_param_indices):
#                             param_name = fitted_model.param_names[param_idx]
#                             param_errors[param_name] = param_std_errors[idx]
#                 except Exception as e:
#                     warnings.warn(f"Could not extract parameter uncertainties: {e}")
        
#         # Populate parameters
#         for param_name in fitted_model.param_names:
#             param = getattr(fitted_model, param_name)
#             stderr = param_errors.get(param_name, None)
            
#             self.params.add(
#                 param_name,
#                 value=param.value,
#                 min=param.bounds[0] if param.bounds else None,
#                 max=param.bounds[1] if param.bounds else None,
#                 vary=not param.fixed,
#                 stderr=stderr
#             )
        
#         # Store additional info
#         self.residual = residuals
#         self.method = 'astropy_levmar'
#         self.cost = self.chisqr / 2  # For compatibility


# class LMFitParams:
#     """Parameter container mimicking lmfit.Parameters"""
#     def __init__(self):
#         self._params = {}
    
#     def add(self, name, value=None, min=None, max=None, vary=True, stderr=None):
#         self._params[name] = LMFitParam(value, min, max, vary, stderr)
    
#     def __getitem__(self, name):
#         return self._params[name]
    
#     def __contains__(self, name):
#         return name in self._params
    
#     def items(self):
#         return self._params.items()
    
#     def keys(self):
#         return self._params.keys()
    
#     def values(self):
#         return self._params.values()


# class LMFitParam:
#     """Parameter object mimicking lmfit.Parameter"""
#     def __init__(self, value, min_val, max_val, vary, stderr):
#         self.value = value
#         self.min = min_val
#         self.max = max_val
#         self.vary = vary
#         self.stderr = stderr


# def do_fit_spec_SY_FF_map(freqs, fluxes, fluxes_err, nu0=None,
#                           fix_alpha_nt=False, verbose=0):
#     """
#     Robust fitting function using astropy's standard fitting framework.
#     Maintains interface compatibility with original lmfit version.
#     """
#     x = freqs
#     y = fluxes
#     yerr = fluxes_err
    
#     if nu0 is None:
#         nu0 = np.mean(x)

#     # Create weights as in original (but simplified for astropy)
#     log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
    
#     # Create model
#     model = RC_SY_FF_Model(nu0=nu0)
    
#     # Set initial parameter values
#     model.A_sy = 1.0
#     model.A_ff = 0.1
#     model.alpha_nt = -0.85
    
#     # Set bounds
#     model.A_sy.bounds = (1.0e-6, 1000)
#     model.A_ff.bounds = (1.0e-6, 100)
    
#     if fix_alpha_nt:
#         model.alpha_nt.fixed = True
#         model.alpha_nt.bounds = (-0.855, -0.845)
#     else:
#         model.alpha_nt.bounds = (-3.0, -0.15)
    
#     # Create fitter - using robust fitting with outlier removal for robustness
#     base_fitter = LevMarLSQFitter()
#     fitter = FittingWithOutlierRemoval(base_fitter, sigma_clip, niter=3, sigma=3.0)
    
#     # Perform fit
#     try:
#         fitted_model, outlier_mask = fitter(model, x, y, weights=log_weights)
        
#         # If we used outlier removal, we need to get the base fitter for error info
#         actual_fitter = fitter.fitter
        
#         if verbose > 0:
#             n_outliers = np.sum(~outlier_mask) if outlier_mask is not None else 0
#             print(f"Identified {n_outliers} outliers during fitting")
            
#     except Exception as e:
#         if verbose > 0:
#             print(f"Robust fitting failed, falling back to standard LevMar: {e}")
        
#         # Fallback to standard fitting
#         actual_fitter = LevMarLSQFitter()
#         fitted_model = actual_fitter(model, x, y, weights=log_weights)
    
#     # Create compatible result object
#     result = LMFitCompatibleResult(fitted_model, actual_fitter, x, y, weights=log_weights)
    
#     if verbose > 0:
#         print(f"Fit success: {result.success}")
#         print(f"Message: {result.message}")
#         print(f"Chi-square: {result.chisqr:.6f}")
#         print(f"Reduced chi-square: {result.redchi:.6f}")
#         print(f"Number of function evaluations: {result.nfev}")
#         print("Parameters:")
#         for name, param in result.params.items():
#             stderr_str = f"{param.stderr:.6f}" if param.stderr is not None else "N/A"
#             fixed_str = " (fixed)" if not param.vary else ""
#             print(f"  {name}: {param.value:.6f} +/- {stderr_str}{fixed_str}")
    
#     return result




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
    # Vary levels by ±noise_rms
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

# def kendall_tau_with_uncertainty(x_data: np.ndarray,
#                                   y_data: np.ndarray,
#                                   x_err: Optional[np.ndarray] = None,
#                                   y_err: Optional[np.ndarray] = None,
#                                   n_bootstrap: int = 1000,
#                                   seed: Optional[int] = None,
#                                  method: Optional[str] = 'auto',
#                                  variant: Optional[str] = 'b',
#                                   ) -> Tuple[float, float]:
#     """
#     Compute Kendall's tau correlation coefficient and its uncertainty via bootstrapping.
    
#     Parameters:
#         x_data (np.ndarray): Array of x data values.
#         y_data (np.ndarray): Array of y data values.
#         x_err (Optional[np.ndarray]): Optional 1sigma errors on x data.
#         y_err (Optional[np.ndarray]): Optional 1sigma errors on y data.
#         n_bootstrap (int): Number of bootstrap samples to estimate uncertainty.
#         seed (Optional[int]): Random seed for reproducibility.
        
#     Returns:
#         tau (float): Kendall's tau coefficient.
#         tau_std (float): Standard deviation (uncertainty) of tau.
#     """
#     rng = np.random.default_rng(seed)

#     # Clean the data: remove NaNs, Infs
#     mask = np.isfinite(x_data) & np.isfinite(y_data)
#     if x_err is not None:
#         mask &= np.isfinite(x_err)
#     if y_err is not None:
#         mask &= np.isfinite(y_err)

#     x = np.asarray(x_data)[mask]
#     y = np.asarray(y_data)[mask]

#     if x_err is not None:
#         x_err = np.asarray(x_err)[mask]
#     if y_err is not None:
#         y_err = np.asarray(y_err)[mask]

#     if len(x) < 2:
#         raise ValueError("Not enough valid data points after cleaning for Kendall's tau.")

#     # Compute the main Kendall's tau
#     tau, _ = kendalltau(x, y,variant=variant,method=method)

#     # Bootstrap sampling
#     tau_samples = []
#     for _ in range(n_bootstrap):
#         # Resample with replacement
#         indices = rng.choice(len(x), len(x), replace=True)

#         x_sample = x[indices]
#         y_sample = y[indices]

#         # If errors provided, perturb values with Gaussian noise
#         if x_err is not None:
#             x_sample = rng.normal(x_sample, x_err[indices])
#         if y_err is not None:
#             y_sample = rng.normal(y_sample, y_err[indices])

#         # Compute tau on the bootstrap sample
#         try:
#             tau_boot, _ = kendalltau(x_sample, y_sample)
#             if np.isfinite(tau_boot):
#                 tau_samples.append(tau_boot)
#         except Exception:
#             continue  # skip if kendalltau fails for some resample

#     tau_std = np.std(tau_samples) if len(tau_samples) > 1 else np.nan

#     return tau, tau_std



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
    # S_target = S_ref × (nu_target/nu_ref)^alpha
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
                             min_cut=3*mad_std(image), max_cut=0.2*np.nanmax(image))
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