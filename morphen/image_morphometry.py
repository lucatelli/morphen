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






def circle_area(r):
    """
    Return the pixel area of a circle given the radius.

    Parameters
    ----------
    r : float, ndarray
        Radius of the circle.

    Returns
    -------
    float, ndarray
        Area of the circle.
    """
    return (np.pi * r * r)


def area_to_radii(A):
    return (np.sqrt(A / np.pi))

def radii_to_area(radii,radii_err=None):
    area_from_radii = np.pi*radii*radii
    if radii_err is not None:
        area_from_radii_err = 2.0*np.pi*radii*radii_err
        # return(area_from_radii,area_from_radii_err)
    else:
        area_from_radii_err = 0.0
    return(area_from_radii,area_from_radii_err)

def pix_area_to_kpc_area(area_pix,pix_to_pc):
    kpc_area = (area_pix * (pix_to_pc**2.0)/(1000**2.0))
    return(kpc_area)
    

def cvi(imname):
    try:
        os.system('casaviewer ' + imname)
    except:
        try:
            os.system('~/casaviewer ' + imname)
        except:
            pass


def normalise_in_log(profile):
    def normalise(x):
        return (x- x.min())/(x.max() - x.min())
    y = profile.copy()
    nu_0 = normalise(np.log(y))
    return nu_0


def shuffle_2D(image):
    height, width = image.shape

    # Reshape the image to a 1D array
    image_flat = image.copy().reshape(-1)

    # Shuffle the pixels randomly
    np.random.shuffle(image_flat)

    # Reshape the shuffled 1D array back to a 2D image
    shuffled_image = image_flat.reshape(height, width)

    # Print the shuffled image
    return(shuffled_image)


def estimate_circular_aperture(image, cellsize, std=3):
    if isinstance(image, str) == True:
        g = load_fits_data(image)
    else:
        g = image
    gmask = g[g > std * mad_std(g)]
    npix = len(gmask)
    barea = beam_area2(image, cellsize)
    nbeams = npix / barea
    circ_radii = np.sqrt(npix / np.pi)
    return (circ_radii)





def sort_masks(image, mask_array, sort_by='flux'):
    unique_labels = np.unique(mask_array)
    region_labels = unique_labels[1:]  # 0 is the zero values
    mask_areas = []
    masks = []
    if sort_by == 'area':
        for i in range(len(region_labels)):
            mask_comp = (mask_array == region_labels[i])
            masks.append(mask_comp)
            area_mask = np.sum(mask_comp)
            mask_areas.append(area_mask)
    if sort_by == 'flux':
        for i in range(len(region_labels)):
            mask_comp = (mask_array == region_labels[i])
            masks.append(mask_comp)
            area_mask = np.sum(mask_comp * image)
            mask_areas.append(area_mask)

    mask_areas = np.asarray(mask_areas)
    sorted_indices_desc = np.argsort(mask_areas)[::-1]
    sorted_arr_desc = mask_areas[sorted_indices_desc]
    return (masks, sorted_indices_desc)


def rot180(imagem, x0, y0):
    R = np.matrix([[-1, 0], [0, -1]])
    img180 = nd.affine_transform(imagem, R,
                                 offset=(2. * (y0 - 0.5), 2. * (x0 - 0.5)))
    return img180


def estimate_area(data_mask, cellsize, Omaj, Omin):
    npix = np.sum(data_mask)
    barea = beam_area(Omaj, Omin, cellsize)
    nbeams = npix / barea
    circ_radii = np.sqrt(npix / np.pi)
    # print(npix)
    return (nbeams, circ_radii, npix)


def estimate_area_nbeam(image, mask, cellsize):
    npix = np.sum(mask)
    barea = beam_area2(image, cellsize)
    nbeams = npix / barea
    circ_radii = np.sqrt(npix / np.pi)
    # print(npix)
    return (circ_radii)

def estimate_area2(image,data_mask,cellsize):
    npix = np.sum(data_mask)
    barea= beam_area2(image,cellsize)
    nbeams = npix/barea
    circ_radii = np.sqrt(npix/np.pi)
    # print(npix)
    return(nbeams,circ_radii,npix)

def get_cell_size(imagename):
    """
    Get the cell size/pixel size in arcsec from an image header wcs.
    """
    hdu = pf.open(imagename)
    ww = WCS(hdu[0].header)
    pixel_scale = (ww.pixel_scale_matrix[1,1]*3600)
    cell_size =  pixel_scale.copy()
    return(cell_size)

def get_dilation_size(image):
    omaj, omin, _, _, _ = beam_shape(image)
    dilation_size = int(
        np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
    return(dilation_size)


def mask_dilation(image, cell_size=None, sigma=6,rms=None,
                  dilation_size=None,iterations=2, dilation_type='disk',
                  PLOT=False,show_figure=False,logger=None,
                  fig_size_factor = 1,
                  special_name='',verbose=0):
    """
    Mask dilation function.

    Apply a binary dilation to a mask originated from the emission of a source.
    This was originally designed for radio images, where the dilation factor is
    proportional to the beam size. The dilation factor is computed as the geometric
    mean of the major and minor axis of the beam. The dilation factor is then
    converted to pixels using the cell size of the image. Note that the expansion
    occurs in two opposite directions, so the dilation size should be half of the
    size of the beam for one iteration, so that the final dilation will be a unit of
    the beam size. However, the exact size of the dilation will also be determinet by
    the number of iterations.


    Parameters
    ----------
    image : str
        Image name.
    cell_size : float
        Cell size of the image in arcsec.
    sigma : float
        Sigma level for the mask.
    rms : float
        RMS level for the mask.
    dilation_size : int
        Size of the dilation.
    iterations : int
        Number of iterations for the dilation.
    dilation_type : str
        Type of dilation. Options are 'disk' or 'square'.
    PLOT : bool
        Plot the mask dilation.
    show_figure : bool
        Show the figure.
    logger : logger
        Logger object.
    fig_size_factor : float
        Figure size factor.
    special_name : str
        Special name for the mask dilation.
    """

    if isinstance(image, str) == True:
        data = load_fits_data(image)
    else:
        data = image
    if rms is None:
        std = mad_std(data)
    else:
        std = rms

    if (dilation_size is None) or (dilation_size == '2x'):
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            if dilation_size == '2x':
                dilation_size = int(
                    2*np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            else:
                dilation_size = int(
                    np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            if verbose >= 1:
                if logger is not None:
                    logger.debug(f" ==>  Mask dilation size is "
                                f"{dilation_size} [px]")
                else:
                    print(f" ==>  Mask dilation size is "
                        f"{dilation_size} [px]")
        except:
            if dilation_size is None:
                dilation_size = 5

    mask = (data >= sigma * std)
    mask3 = (data >= 3 * std)
    data_mask = mask * data

    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=disk(dilation_size),
                                            iterations=iterations).astype(mask.dtype)

    if dilation_type == 'square':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=square(dilation_size),
                                            iterations=iterations).astype(mask.dtype)

    if PLOT == True:
        fig = plt.figure(figsize=(int(15*fig_size_factor), int(4*fig_size_factor)))
        ax0 = fig.add_subplot(1, 4, 1)
        ax0.imshow((mask3), origin='lower',cmap='magma')
        ax0.set_title(r'mask $>$ ' + str(3) + '$\sigma_{\mathrm{mad}}$')
        ax0.axis('off')
        ax1 = fig.add_subplot(1, 4, 2)
        #         ax1.legend(loc='lower left')
        ax1.imshow((mask), origin='lower',cmap='magma')
        ax1.set_title(r'mask $>$ ' + str(sigma) + '$\sigma_{\mathrm{mad}}$')
        ax1.axis('off')
        ax2 = fig.add_subplot(1, 4, 3)
        ax2.imshow(data_mask_d, origin='lower',cmap='magma')
        ax2.set_title(r'D(mask)'f'[{dilation_size}|{iterations}]'f'{special_name}')
        ax2.axis('off')
        ax3 = fig.add_subplot(1, 4, 4)
        ax3 = eimshow(data * data_mask_d, ax=ax3, vmin_factor=0.01,CM='magma')
        ax3.set_title(r'D(mask) $\times$ data')
        #         ax3.imshow(np.log(data*data_mask_d))
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        #         fig.tight_layout()
        ax3.axis('off')
        if show_figure == True:
            plt.show()
        else:
            plt.close()
    #         plt.savefig(image.replace('.fits','_masks.jpg'),dpi=300, bbox_inches='tight')

    # if cell_size is not None:
    #     if isinstance(image, str) == True:
    #         try:
    #             print((data * data_mask_d).sum() / beam_area2(image, cell_size))
    #             print((data * data_mask).sum() / beam_area2(image, cell_size))
    #             print((data).sum() / beam_area2(image, cell_size))
    #         except:
    #             print('Provide a cell size of the image.')
    return (mask, data_mask_d)



def t_mask_dilation(image, cell_size=None, sigma=6, rms=None,
                  dilation_size=None, iterations=2, dilation_type='disk',
                  PLOT=False, show_figure=False, logger=None,
                  fig_size_factor=1,
                  special_name='', verbose=0):
    """
    Mask dilation function with connected component filtering.
    
    Similar to original function but only keeps the largest connected component,
    typically at the center of the image.
    """
    from skimage.measure import label, regionprops
    
    if isinstance(image, str) == True:
        data = load_fits_data(image)
    else:
        data = image
    if rms is None:
        std = mad_std(data)
    else:
        std = rms

    if (dilation_size is None) or (dilation_size == '2x'):
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            if dilation_size == '2x':
                dilation_size = int(
                    2*np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            else:
                dilation_size = int(
                    np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
            if verbose >= 1:
                if logger is not None:
                    logger.debug(f" ==>  Mask dilation size is {dilation_size} [px]")
                else:
                    print(f" ==>  Mask dilation size is {dilation_size} [px]")
        except:
            if dilation_size is None:
                dilation_size = 5

    # Create initial masks
    mask = (data >= sigma * std)
    mask3 = (data >= 3 * std)
    
    # Label connected components in the mask
    labeled_mask = label(mask)
    regions = regionprops(labeled_mask)
    
    # Find the largest connected component
    if regions:
        largest_area = 0
        largest_label = None
        for region in regions:
            if region.area > largest_area:
                largest_area = region.area
                largest_label = region.label
        # Create a new mask with only the largest component
        if largest_label is not None:
            mask = labeled_mask == largest_label
    
    data_mask = mask * data

    # Perform dilation
    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=disk(dilation_size),
                                            iterations=iterations).astype(mask.dtype)
    if dilation_type == 'square':
        data_mask_d = ndimage.binary_dilation(mask,
                                            structure=square(dilation_size),
                                            iterations=iterations).astype(mask.dtype)

    # Clean up the 3-sigma mask for plotting using the same approach
    labeled_mask3 = label(mask3)
    regions3 = regionprops(labeled_mask3)
    if regions3:
        largest_area = 0
        largest_label = None
        for region in regions3:
            if region.area > largest_area:
                largest_area = region.area
                largest_label = region.label
        if largest_label is not None:
            mask3 = labeled_mask3 == largest_label

    if PLOT == True:
        fig = plt.figure(figsize=(int(15*fig_size_factor), int(4*fig_size_factor)))
        ax0 = fig.add_subplot(1, 4, 1)
        ax0.imshow((mask3), origin='lower',cmap='magma')
        ax0.set_title(r'mask $>$ ' + str(3) + '$\sigma_{\mathrm{mad}}$')
        ax0.axis('off')
        ax1 = fig.add_subplot(1, 4, 2)
        ax1.imshow((mask), origin='lower',cmap='magma')
        ax1.set_title(r'mask $>$ ' + str(sigma) + '$\sigma_{\mathrm{mad}}$')
        ax1.axis('off')
        ax2 = fig.add_subplot(1, 4, 3)
        ax2.imshow(data_mask_d, origin='lower',cmap='magma')
        ax2.set_title(r'D(mask)'f'[{dilation_size}|{iterations}]'f'{special_name}')
        ax2.axis('off')
        ax3 = fig.add_subplot(1, 4, 4)
        ax3 = eimshow(data * data_mask_d, ax=ax3, vmin_factor=0.01,CM='magma')
        ax3.set_title(r'D(mask) $\times$ data')
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        ax3.axis('off')
        if show_figure == True:
            plt.show()
        else:
            plt.close()

    return (mask, data_mask_d)



def mask_dilation_from_list(images, residuals, sigma=6, rms=None, dilation_size=None, iterations=2, dilation_type='disk',
                  PLOT=False, show_figure=True, logger=None, fig_size_factor=1, special_name='', verbose=0):
    """
    Mask dilation function for multi-frequency images.
    
    This function applies mask dilation across multiple images and computes a common mask
    that considers emission across all images. It avoids setting a smaller mask as the common mask.
    
    Parameters
    ----------
    images : list of str or ndarray
        List of image names or image data arrays.
    sigma : float
        Sigma level for the mask.
    rms : float
        RMS level for the mask.
    dilation_size : int
        Size of the dilation.
    iterations : int
        Number of iterations for the dilation.
    dilation_type : str
        Type of dilation. Options are 'disk' or 'square'.
    PLOT : bool
        Plot the mask dilation.
    show_figure : bool
        Show the figure.
    logger : logger
        Logger object.
    fig_size_factor : float
        Figure size factor.
    special_name : str
        Special name for the mask dilation.
    """

    # Function to compute mask for a single image
    def compute_mask(image, sigma, rms, dilation_size, iterations, dilation_type, verbose):
        if isinstance(image, str):
            data = load_fits_data(image)
        else:
            data = image
        if rms is None:
            std = mad_std(data)
        else:
            std = rms
        
        if (dilation_size is None) or (dilation_size == '2x'):
            try:
                omaj, omin, _, _, _ = beam_shape(image)
                if dilation_size == '2x':
                    dilation_size = int(2 * np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
                else:
                    dilation_size = int(np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
                if verbose >= 1:
                    if logger is not None:
                        logger.debug(f" ==> Mask dilation size is {dilation_size} [px]")
                    else:
                        print(f" ==> Mask dilation size is {dilation_size} [px]")
            except:
                dilation_size = 5

        mask = (data >= sigma * std)

        if dilation_type == 'disk':
            data_mask_d = ndimage.binary_dilation(mask, structure=disk(dilation_size), iterations=iterations).astype(mask.dtype)
        elif dilation_type == 'square':
            data_mask_d = ndimage.binary_dilation(mask, structure=square(dilation_size), iterations=iterations).astype(mask.dtype)

        return mask, data_mask_d

    # Compute individual masks for each image
    mask_list = []
    for kk in range(len(images)):
        rms = mad_std(load_fits_data(residuals[kk]))
        mask, dilated_mask = compute_mask(images[kk], sigma, rms, dilation_size, iterations, dilation_type, verbose)
        mask_list.append(dilated_mask)

    # Compute the common mask
    combined_mask = np.sum(mask_list, axis=0) >= len(mask_list) // 2  # Majority vote across masks

    # Optional plot of the common mask
    if PLOT:
        fig = plt.figure(figsize=(int(15 * fig_size_factor), int(4 * fig_size_factor)))
        ax = fig.add_subplot(1, 1, 1)
        ax.imshow(combined_mask, origin='lower', cmap='magma')
        ax.set_title('Common Mask')
        ax.axis('off')
        if show_figure:
            plt.show()
        else:
            plt.close()

    return combined_mask




def mask_dilation_from_mask(image, mask_init, cell_size=None, sigma=3,rms=None,
                  dilation_size=None,iterations=2, dilation_type='disk',
                  PLOT=True,show_figure=True):
    """
    Apply a dilation to an existing mask.
    """
    from scipy import ndimage
    from scipy.ndimage import morphology
    from skimage.morphology import disk, square
    from skimage.morphology import dilation

    if isinstance(image, str) == True:
        data = load_fits_data(image)
    else:
        data = image
    if rms is None:
        std = mad_std(data)
    else:
        std = rms

    if dilation_size is None:
        try:
            omaj, omin, _, _, _ = beam_shape(image)
            dilation_size = int(
                np.sqrt(omaj * omin) / (2 * get_cell_size(image)))
        except:
            if dilation_size is None:
                dilation_size = 5

    data_init = data * mask_init
    mask3 = (data >= 3 * std)
    # std = mad_std(data[mask_init])
    mask = (data_init >= sigma * std)
    data_mask = mask * data

    if dilation_type == 'disk':
        data_mask_d = ndimage.binary_dilation(mask_init,
                                            structure=disk(dilation_size),
                                            iterations=iterations).\
                                                astype(mask_init.dtype)

    if dilation_type == 'square':
        data_mask_d = ndimage.binary_dilation(mask_init,
                                            structure=square(dilation_size),
                                            iterations=iterations).\
                                                astype(mask_init.dtype)

    if PLOT == True:
        fig = plt.figure(figsize=(15, 4))
        ax0 = fig.add_subplot(1, 4, 1)
        ax0.imshow((mask3), origin='lower')
        ax0.set_title(r'Mask above' + str(3) + '$\sigma$')
        ax0.axis('off')
        ax1 = fig.add_subplot(1, 4, 2)
        #         ax1.legend(loc='lower left')
        ax1.imshow((mask), origin='lower')
        ax1.set_title(r'Mask above' + str(sigma) + '$\sigma$')
        ax1.axis('off')
        ax2 = fig.add_subplot(1, 4, 3)
        ax2.imshow(data_mask_d, origin='lower')
        ax2.set_title(r'Dilated mask')
        ax2.axis('off')
        ax3 = fig.add_subplot(1, 4, 4)
        ax3 = fast_plot2(data * data_mask_d, ax=ax3, vmin_factor=0.1)
        ax3.set_title(r'Dilated mask $\times$ data')
        #         ax3.imshow(np.log(data*data_mask_d))
        plt.subplots_adjust(wspace=0.1, hspace=0.1)
        #         fig.tight_layout()
        ax3.axis('off')
        if show_figure == True:
            plt.show()
        else:
            plt.close()

    if cell_size is not None:
        if isinstance(image, str) == True:
            try:
                print((data * data_mask_d).sum() / beam_area2(image, cell_size))
                print((data * data_mask).sum() / beam_area2(image, cell_size))
                print((data).sum() / beam_area2(image, cell_size))
            except:
                print('Provide a cell size of the image.')
    return (mask, data_mask_d)


def convert_seg_map_to_masks(seg_map):
    """
    Convert a segmentation map to a list of masks, where each mask corresponds to a unique segment.

    Parameters:
        seg_map (2D array): Segmentation map where different segments have different integer values,
                            and masked regions have a value of zero.

    Returns:
        list: A list of boolean masks, where each mask corresponds to a segment in the segmentation map.
    """
    # Find the unique segment labels, excluding zero (the background)
    unique_segments = np.unique(seg_map)
    unique_segments = unique_segments[unique_segments != 0]

    # Create a mask for each segment
    masks = [(seg_map == segment) for segment in unique_segments]

    return masks

def convert_masks_to_segmentation(masks_regions):
    # Assuming all masks have the same shape
    height, width = masks_regions[0].shape
    
    # Create an empty array for the segmentation image
    segmentation_image = np.zeros((height, width), dtype=int)
    
    # Loop over each mask and assign a unique integer value
    for idx, mask in enumerate(masks_regions, start=1):
        # idx+1 ensures the values start from 1 instead of 0
        segmentation_image[mask == 1] = idx
    
    return segmentation_image


def split_overlapping_masks(mask1, mask2):
    """
    # Compute the overlap between the masks

    GIven two masks, find the overlapping region between the two and
    split half of the pixels in this region attributing them to one maks and the
    other half to the other mask.

    # Check for overlap between the masks and split the overlapping region
    mask1_new, mask2_new = split_overlapping_masks(mask1, mask2)

    # Show the new masks
    print(mask1_new)
    print(mask2_new)

    """

    overlap_mask = (mask1 > 0) & (mask2 > 0)

    # If there is no overlap, return the original masks
    if not np.any(overlap_mask):
        return mask1, mask2

    # Split the overlapping region in half
    half_mask1 = mask1.copy()
    half_mask2 = mask2.copy()
    overlap_indices = np.where(overlap_mask)
    num_overlapping_pixels = len(overlap_indices[0])
    half_num_overlapping_pixels = num_overlapping_pixels // 2
    for i in range(num_overlapping_pixels):
        row, col = overlap_indices[0][i], overlap_indices[1][i]
        if i < half_num_overlapping_pixels:
            half_mask1[row, col] = mask1[row, col]
            half_mask2[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = mask2[row, col]

    # Return the split masks
    return half_mask1, half_mask2




def split_overlapping_masks3(mask1, mask2, mask3):
    """
    # Compute the overlap between the masks

    The same as `split_overlapping_masks`, but for three masks.
    """
    overlap_mask1 = (mask1 > 0) & (mask2 > 0) & (mask3 == 0)
    overlap_mask2 = (mask1 > 0) & (mask3 > 0) & (mask2 == 0)
    overlap_mask3 = (mask2 > 0) & (mask3 > 0) & (mask1 == 0)

    # If there is no overlap, return the original masks
    if not np.any(overlap_mask1) and not np.any(overlap_mask2) and not np.any(
            overlap_mask3):
        return mask1, mask2, mask3

    # Split the overlapping region in half for each overlapping pair
    num_overlapping_pixels1 = np.sum(overlap_mask1)
    num_overlapping_pixels2 = np.sum(overlap_mask2)
    num_overlapping_pixels3 = np.sum(overlap_mask3)
    half_num_overlapping_pixels1 = num_overlapping_pixels1 // 2
    half_num_overlapping_pixels2 = num_overlapping_pixels2 // 2
    half_num_overlapping_pixels3 = num_overlapping_pixels3 // 2

    half_mask1 = mask1.copy()
    half_mask2 = mask2.copy()
    half_mask3 = mask3.copy()

    overlap_indices1 = np.where(overlap_mask1)
    for i in range(num_overlapping_pixels1):
        row, col = overlap_indices1[0][i], overlap_indices1[1][i]
        if i < half_num_overlapping_pixels1:
            half_mask1[row, col] = mask1[row, col]
            half_mask2[row, col] = 0
            half_mask3[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = mask2[row, col]
            half_mask3[row, col] = 0

    overlap_indices2 = np.where(overlap_mask2)
    for i in range(num_overlapping_pixels2):
        row, col = overlap_indices2[0][i], overlap_indices2[1][i]
        if i < half_num_overlapping_pixels2:
            half_mask1[row, col] = mask1[row, col]
            half_mask2[row, col] = 0
            half_mask3[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = 0
            half_mask3[row, col] = mask3[row, col]

    overlap_indices3 = np.where(overlap_mask3)
    for i in range(num_overlapping_pixels3):
        row, col = overlap_indices3[0][i], overlap_indices3[1][i]
        if i < half_num_overlapping_pixels3:
            half_mask1[row, col] = 0
            half_mask2[row, col] = mask2[row, col]
            half_mask3[row, col] = 0
        else:
            half_mask1[row, col] = 0
            half_mask2[row, col] = 0
            half_mask3[row, col] = mask3[row, col]

    # Return the split masks
    return half_mask1, half_mask2, half_mask3


def trail_vector(vx, vy, v0=np.asarray([1, 0])):
    """
    Compute the trailing vector angle (in degrees) and its magnitude.
    The angle (PA) is measured counterclockwise from the positive x-axis.

    Parameters:
    vx, vy : float
        Components of the input vector.
    v0 : array-like, optional
        Reference vector (default is [1, 0] for the positive x-axis).

    Returns:
    angle_PA : float
        Position angle (PA) measured counterclockwise from the positive x-axis, in degrees.
    norm_vec : float
        Magnitude of the input vector.
    """
    from scipy.linalg import norm
    v = np.asarray([vx, vy])
    norm_vec = norm(v)
    
    # Handle zero vector case
    if norm_vec == 0:
        raise ValueError("Input vector (vx, vy) has zero magnitude.")
    
    # Normalize input vector
    v_hat = v / norm_vec
    
    # Compute the angle using arctan2 (counterclockwise convention)
    angle_radians = np.arctan2(v_hat[1], v_hat[0])  # atan2(vy, vx)
    angle_PA = np.degrees(angle_radians)  # Convert to degrees
    
    # Ensure the angle is in the range [0, 360)
    if angle_PA < 0:
        angle_PA += 360
    return angle_PA, norm_vec


def calculate_radii(image, mask):
    """
    Calculate circular radii of the emission in a 2D numpy image.

    Parameters:
        image (np.ndarray): The 2D numpy image containing the radio emission of a galaxy.
        sigma_level (float): The sigma level to use for the contour.
        background_std (float): The standard deviation of the background.

    Returns:
        float: The circular radii of the emission.
    """
    # Calculate the threshold level for the given sigma level
    #     threshold = sigma_level * background_std

    # Create a binary mask of the emission above the threshold level
    #     mask = image > threshold

    # Calculate the center of mass of the emission
    from scipy.ndimage import measurements
    com = measurements.center_of_mass(image * mask)

    # Calculate the distance of each pixel from the center of mass
    y, x = np.indices(mask.shape)
    distances = np.sqrt((x - com[1]) ** 2 + (y - com[0]) ** 2)

    # Calculate the median distance of the pixels above the threshold level
    median_distance = np.nanmedian(distances[mask])
    mean_distance = np.nanmean(distances[mask])
    std_distance = np.nanstd(distances[mask])
    # mad_std_distance = mad_std(distances[mask])

    return median_distance, mean_distance, std_distance



# Testing functions, not used anywhere.
def shannon_entropy_2d(arr):
    rows, cols = arr.shape
    result = np.zeros((rows, cols))
    for i in range(2, rows-2):
        for j in range(2, cols-2):
            box = arr[i-2:i+3, j-2:j+3].ravel()
            p = box / np.nansum(box)
            entropy = -np.nansum(p * np.log2(p))
            result[i, j] = entropy
    return result

def cvi(imname):
    try:
        os.system('casaviewer ' + imname)
    except:
        try:
            os.system('~/casaviewer ' + imname)
        except:
            pass





def get_image_statistics(imagename,cell_size=None,
                         mask_component=None,mask=None,
                         residual_name=None,region='', dic_data=None,
                         sigma_mask=6,apply_mask=True,
                         fracX=0.1, fracY=0.1):
    """
    Get some basic image statistics.



    """
    if dic_data is None:
        dic_data = {}
        dic_data['#imagename'] = os.path.basename(imagename)

    if cell_size is None:
        cell_size = get_cell_size(imagename)

    image_data = load_fits_data(imagename)
    if mask is not None:
        image_data = np.nan_to_num(image_data*mask,nan=0)
        apply_mask = False
    #     dic_data['imagename'] = imagename
    if apply_mask == True:
        omask, mask = mask_dilation(imagename, sigma=sigma_mask)
        image_data = np.nan_to_num(image_data*mask,nan=0)
    if mask_component is not None:
        image_data = np.nan_to_num(image_data * mask_component,nan=0)
        if mask is None:
            mask = mask_component
        # mask = mask_component
    if (mask_component is None) and (apply_mask == False) and (mask is None):
        unity_mask = np.ones(load_fits_data(imagename).shape)
        omask, mask = unity_mask, unity_mask
        image_data = np.nan_to_num(load_fits_data(imagename) * mask,nan=0)

    stats_im = imstat(imagename=imagename, region=region)

    box_edge, imhd = create_box(imagename, fracX=fracX, fracY=fracY)
    stats_box = imstat(imagename=imagename, box=box_edge)

    # determine the flux peak and positions of image
    flux_peak_im = stats_im['max'][0]
    flux_min_im = stats_im['min'][0]
    dic_data['max_im'] = flux_peak_im
    dic_data['min_im'] = flux_min_im
    """
    x0max,y0max = peak_center(image_data)
    dic_data['x0'], dic_data['y0'] = x0max,y0max
    # determine momentum centres.
    x0m, y0m, _, _ = momenta(image_data, PArad_0=None, q_0=None)
    dic_data['x0m'], dic_data['y0m'] = x0m, y0m

    #some geometrical measures
    # calculate PA and axis-ratio
    PA, q, x0col, y0col, PAm, qm, PAmi, qmi, PAmo, qmo,\
        x0median,y0median,\
        x0median_i,y0median_i,x0median_o,y0median_o = cal_PA_q(image_data)

    dic_data['PA'], dic_data['q'] = PA, q
    dic_data['PAm'], dic_data['qm'] = PAm, qm
    dic_data['PAm'], dic_data['qm'] = PAm, qm
    dic_data['PAmi'], dic_data['qmi'] = PAmi, qmi
    dic_data['PAmo'], dic_data['qmo'] = PAmo, qmo
    dic_data['x0m_i'], dic_data['y0m_i'] = x0median_i, y0median_i
    dic_data['x0m_o'], dic_data['y0m_o'] = x0median_o, y0median_o
    """

    # determine the rms and std of residual and of image
    rms_im = stats_im['rms'][0]
    rms_box = stats_box['rms'][0]
    sigma_im = stats_im['sigma'][0]
    sigma_box = stats_box['sigma'][0]

    dic_data['rms_im'] = rms_im
    dic_data['rms_box'] = rms_box
    dic_data['sigma_im'] = sigma_im
    dic_data['sigma_box'] = sigma_box

    # determine the image and residual flux
    flux_im = stats_im['flux'][0]
    flux_box = stats_box['flux'][0]
    dic_data['flux_im'] = flux_im
    dic_data['flux_box'] = flux_box
    sumsq_im = stats_im['sumsq'][0]
    sumsq_box = stats_box['sumsq'][0]

    q_sq = sumsq_im / sumsq_box
    q_flux = flux_im / flux_box
    # flux_ratio = flux_re/flux_im
    dic_data['q_sq'] = q_sq
    dic_data['q_flux'] = q_flux

    snr = flux_im / rms_box
    snr_im = flux_im / rms_im

    dr_e = []
    frac_ = np.linspace(0.05, 0.85, 10)
    frac_image = 0.10
    '''
    Each loop below run a sliding window,
    one to the x direction and the other to y-dircetion.
    This is to get a better estimate (in multiple regions)
    of the background rms and therefore SNR
    Each window has a fraction frac_image of the image size.
    '''
    for frac in frac_:
        box, _ = create_box(imagename, fracX=frac, fracY=frac_image)
        st = imstat(imagename, box=box)
        snr_tmp = flux_peak_im / st['rms'][0]
        dr_e.append(snr_tmp)

    dr_e2 = []
    for frac in frac_:
        box, _ = create_box(imagename, fracX=frac_image, fracY=frac)
        st = imstat(imagename, box=box)
        snr_tmp = flux_peak_im / st['rms'][0]
        dr_e2.append(snr_tmp)
    #average of the SNR -- DINAMIC RANGE
    DR_SNR_E = (np.mean(dr_e) + np.mean(dr_e2)) / 2

    dic_data['snr'] = snr
    dic_data['snr_im'] = snr_im
    dic_data['DR_SNR_E'] = DR_SNR_E

    DR_pk_rmsbox = flux_peak_im / rms_box
    DR_pk_rmsim = flux_peak_im / rms_im
    dic_data['DR_pk_rmsbox'] = DR_pk_rmsbox
    dic_data['DR_pk_rmsim'] = DR_pk_rmsim

    dic_data['bmajor'] = imhd['restoringbeam']['major']['value']
    dic_data['bminor'] = imhd['restoringbeam']['minor']['value']
    dic_data['positionangle'] = imhd['restoringbeam']['positionangle']['value']



    # if residual_name is not None:
    #     data_res = load_fits_data(residual_name)
    #     flux_res_error = 3 * np.sum(data_res * mask) \
    #                      / beam_area2(imagename, cell_size)
    #     # rms_res =imstat(residual_name)['flux'][0]
    #     flux_res = np.sum(load_fits_data(residual_name)) / beam_area2(imagename, cell_size)
    #
    #     res_error_rms =np.sqrt(
    #         np.sum((abs(data_res * mask -
    #                     np.mean(data_res * mask))) ** 2 * np.sum(mask))) / \
    #                    beam_area2(imagename,cell_size)
    #
    #     try:
    #         total_flux_tmp = dic_data['total_flux_mask']
    #     except:
    #         total_flux_tmp = flux_im
    #         total_flux_tmp = total_flux(image_data,imagename,mask=mask)
    #
    #     # print('Estimate #1 of flux error (based on sum of residual map): ')
    #     # print('Flux = ', total_flux_tmp * 1000, '+/-',
    #     #       abs(flux_res_error) * 1000, 'mJy')
    #     # print('Fractional error flux = ', flux_res_error / total_flux_tmp)
    #     print('-----------------------------------------------------------------')
    #     print('Estimate of flux error (based on rms of '
    #           'residual x area): ')
    #     print('Flux = ', total_flux_tmp * 1000, '+/-',
    #           abs(res_error_rms) * 1000, 'mJy')
    #     print('Fractional error flux = ', res_error_rms / total_flux_tmp)
    #     print('-----------------------------------------------------------------')
    #
    #     dic_data['max_residual'] = np.max(data_res * mask)
    #     dic_data['min_residual'] = np.min(data_res * mask)
    #     dic_data['flux_residual'] = flux_res
    #     dic_data['flux_error_res'] = abs(flux_res_error)
    #     dic_data['flux_error_res_2'] = abs(res_error_rms)
    #     dic_data['mad_std_residual'] = mad_std(data_res)
    #     dic_data['rms_residual'] = rms_estimate(data_res)

    #     print(' Flux=%.5f Jy/Beam' % flux_im)
    #     print(' Flux peak (image)=%.5f Jy' % flux_peak_im, 'Flux peak (residual)=%.5f Jy' % flux_peak_re)
    #     print(' flux_im/sigma_im=%.5f' % snr_im, 'flux_im/sigma_re=%.5f' % snr)
    #     print(' rms_im=%.5f' % rms_im, 'rms_re=%.5f' % rms_re)
    #     print(' flux_peak_im/rms_im=%.5f' % peak_im_rms, 'flux_peak_re/rms_re=%.5f' % peak_re_rms)
    #     print(' sumsq_im/sumsq_re=%.5f' % q)
    return (dic_data)


def level_statistics(img, cell_size=None, mask_component=None,
                    sigma=6, do_PLOT=False, crop=False,data_2D = None,data_res=None,
                    box_size=256, bkg_to_sub=None, apply_mask=True,
                    mask=None,rms=None,
                    results=None, dilation_size=None, iterations=2,
                    add_save_name='', SAVE=False, show_figure=False, ext='.jpg'):
    """
    Function old name: plot_values_std

    Slice the intensity values of an image into distinct regions.
    Then, compute information for each bin level of the emission.
    The implemented splitting is:

        1. Inner region: peak intensity -> 0.1 * peak intensity
        2. Mid-region: 0.1 * peak intensity -> 10 * rms
        3. Low-region: 10 * rms -> 6 * rms
        4. Uncertain region: 6 * rms -> 3 * rms

    Parameters
    ----------
    img : str
        Path to the image.
    cell_size : float, optional
        Cell size of the image. The default is None. In that case, the function
        get_cell_size will attempt to estimate it from the header of the image.
    mask_component : array, optional
        The default is None. This is designed to be used when an image is complex,
        and one would like to study multiple components of the emission separately,
        each one at a time.
    sigma : float, optional
        The default is 6. This is the number of standard deviations to be used during
        mask dilation.
    do_PLOT : bool, optional
        The default is False. If True, the function will plot and save the image.
    crop : bool, optional
        The default is False. If True, the function will crop the image to a square
        of size box_size.
    box_size : int, optional
        The default is 256. This is the size of the square to be used if crop=True.
    data_2D : array, optional
        The default is None. If not None, the function will use this array instead
        and consider header information from img to be used with the array data_2D.

    """
    if cell_size is None:
        cell_size = get_cell_size(img)
    if data_2D is not None:
        g_ = data_2D
    else:
        g_ = load_fits_data(img)
    g = g_.copy()
    if rms is None:
        std = mad_std(g_)
    else:
        std = rms

    if bkg_to_sub is not None:
        g = g - bkg_to_sub
    if mask_component is not None:
        g = g * mask_component

    beam_area_ = beam_area2(img)

    if mask is not None:
        g = g * mask
        apply_mask = False  # do not calculate the mask again, in case is True.
        g = g * mask
    if apply_mask == True:
        _, mask_dilated = mask_dilation(img, cell_size=cell_size, sigma=sigma,
                                        dilation_size=dilation_size, rms=rms,
                                        iterations=iterations,
                                        PLOT=False)
        g = g * mask_dilated


    g = np.nan_to_num(g,nan=0)
    
    if mask_component is not None:
        levels = np.geomspace(g.max(), (1 * std), 5)
        levels_top = np.geomspace(g.max(), g.max() * 0.1, 3)
        try:
            levels_mid = np.geomspace(g.max() * 0.1, (10 * std), 5)
        except:
            levels_mid = np.asarray([0])
        try:
            levels_low = np.geomspace(10 * std, (6.0  * std), 2)
            levels_uncertain = np.geomspace(6.0 * std, (3.0 * std), 3)
        except:
            levels_low = np.asarray([0])
            levels_uncertain = np.asarray([0])

    else:
        if apply_mask is not False:
            # print('asdasd', g.max(), std)
            levels = np.geomspace(g.max(), (1 * std), 5)
            levels_top = np.geomspace(g.max(), g.max() * 0.1, 3)
            try:
                levels_mid = np.geomspace(g.max() * 0.1, (10 * std), 5)
            except:
                levels_mid = np.asarray([0])
            try:
                levels_low = np.geomspace(10 * std, (6.0  * std), 2)
                levels_uncertain = np.geomspace(6 * std, (3.0 * std), 3)
            except:
                levels_low = np.asarray([0])
                levels_uncertain = np.asarray([0])
        else:
            levels = np.geomspace(g.max(), (3 * std), 5)
            levels_top = np.geomspace(g.max(), g.max() * 0.1, 3)
            # levels_mid = np.geomspace(g.max() * 0.1, (10 * std + dl), 5)
            try:
                levels_mid = np.geomspace(g.max() * 0.1, (10 * std), 5)
            except:
                levels_mid = np.asarray([0])
            try:
                levels_low = np.geomspace(10 * std, (6.0  * std), 2)
                levels_uncertain = np.geomspace(3 * std, (1.0 * std), 3)
            except:
                levels_low = np.asarray([0])
                levels_uncertain = np.asarray([0])

    # pix_inner = g[g >= levels_top[-1]]
    # pix_mid = g[np.where((g < levels_top[-1]) & (g >= levels_mid[-1]))]
    # pix_low = g[np.where((g < levels_mid[-1]) & (g >= levels_low[-1]))]
    # pix_uncertain = g[np.where((g < levels_low[-1]) & (g >= levels_uncertain[-1]))]
    pix_inner_mask = (g >= levels_top[-1])
    pix_inner = g * pix_inner_mask
    pix_mid_mask = (((g < levels_top[-1]) & (g >= levels_mid[-1])))
    pix_mid = g * pix_mid_mask
    pix_low_mask = (((g < levels_mid[-1]) & (g >= levels_low[-1])))
    pix_low = g * pix_low_mask
    pix_uncertain_mask = (((g < levels_low[-1]) & (g >= levels_uncertain[-1])))
    pix_uncertain = g * pix_uncertain_mask
    inner_flux = pix_inner.sum() / beam_area_
    mid_flux = pix_mid.sum() / beam_area_
    low_flux = pix_low.sum() / beam_area_
    uncertain_flux = pix_uncertain.sum() / beam_area_

    total_flux = low_flux + mid_flux + inner_flux + uncertain_flux
    ext_flux = low_flux + mid_flux + uncertain_flux
    pix_area = len(g[g >= 3 * std])
    number_of_beams = pix_area / beam_area_
    n_beams_inner = np.sum(pix_inner_mask) / beam_area_
    n_beams_mid = np.sum(pix_mid_mask) / beam_area_
    n_beams_low = np.sum(pix_low_mask) / beam_area_
    n_beams_uncertain = np.sum(pix_uncertain_mask) / beam_area_
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(img)

    # print('Low Flux (extended) Jy                    > ', low_flux, ' >> ratio=',
    #       low_flux / total_flux)
    # print('Mid Flux (outer core + inner extended) Jy > ', mid_flux, ' >> ratio=',
    #       mid_flux / total_flux)
    # print('Inner Flux (core) Jy                      > ', inner_flux,
    #       ' >> ratio=', inner_flux / total_flux)
    # print('Uncertain Flux (<5std)                    > ', uncertain_flux,
    #       ' >> ratio=', uncertain_flux / total_flux)
    # print('Total Flux Jy                             > ', total_flux)
    # print('Total area (in # ob beams)                > ', number_of_beams)
    # print('Total inner area (in # ob beams)          > ', n_beams_inner)
    # print('Total mid area (in # ob beams)            > ', n_beams_mid)
    # print('Total low area (in # ob beams)            > ', n_beams_low)
    # print('Total uncertain area (in # ob beams)      > ', n_beams_uncertain)
    # print('Inner Flux (core) fraction                > ',
    #       inner_flux / total_flux)
    # print('Outer Flux (ext)  fraction                > ', ext_flux / total_flux)

    results['peak_of_flux'] = np.max(g)
    results['total_flux'] = total_flux
    results['inner_flux'] = inner_flux
    results['low_flux'] = low_flux
    results['mid_flux'] = mid_flux
    results['uncertain_flux'] = uncertain_flux

    results['inner_flux_f'] = inner_flux / total_flux
    results['low_flux_f'] = low_flux / total_flux
    results['mid_flux_f'] = mid_flux / total_flux
    results['uncertain_flux_f'] = uncertain_flux / total_flux

    results['number_of_beams'] = number_of_beams
    results['n_beams_inner'] = n_beams_inner
    results['n_beams_mid'] = n_beams_mid
    results['n_beams_low'] = n_beams_low
    results['n_beams_uncertain'] = n_beams_uncertain

    if do_PLOT == True:
        fig = plt.figure(figsize=(6, 6))
        ax = fig.add_subplot()
        vmin = 1 * std
        vmax = g_.max()
        # norm = visualization.simple_norm(g, stretch='log')#, max_percent=max_percent_lowlevel)
        norm = simple_norm(g, stretch='asinh', asinh_a=0.005, min_cut=vmin,
                           max_cut=vmax)

        if crop == True:
            try:
                xin, xen, yin, yen = do_cutout(img, box_size=box_size,
                                               center=center, return_='box')
                g = g[xin:xen, yin:yen]
            except:
                try:
                    max_x, max_y = np.where(g == g.max())
                    xin = max_x[0] - box_size
                    xen = max_x[0] + box_size
                    yin = max_y[0] - box_size
                    yen = max_y[0] + box_size
                    g = g[xin:xen, yin:yen]
                except:
                    pass
        try:
            im_plot = ax.imshow(g, cmap='magma_r', norm=norm, alpha=1.0,
                                origin='lower')
            ax.contour(g, levels=levels_top[::-1], colors='lime',
                       alpha=1.0)  # cmap='Reds', linewidths=0.75)
            ax.contour(g, levels=levels_mid[::-1], colors='yellow',
                       linewidths=0.75)
            ax.contour(g, levels=levels_low[::-1],
                       colors='#56B4E9')  # cmap='Greens', linewidths=0.75)
            ax.contour(g, levels=levels_uncertain[::-1], colors='grey',
                       linewidths=0.4)
        except:
            pass
        # im_plot.colorbar()
        #         plt.subplots_adjust(wspace=0, hspace=0)
        #         fig.tight_layout()

        if SAVE is not None:
            plt.savefig(
                img.replace('.fits', '_std_levels') + add_save_name + ext,
                dpi=300,
                bbox_inches='tight')
        if show_figure == True:
            plt.show()
        else:
            plt.close()
    return (results)


def compute_asymetries(imagename,mask,mask_component=None,
                       bkg_to_sub=None,
                       centre=None,results=None):
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(imagename)

    if isinstance(imagename, str) == True:
        image_data = load_fits_data(imagename)
    else:
        image_data = imagename

    if bkg_to_sub is not None:
        image_data = image_data - bkg_to_sub
    if mask_component is not None:
        image_data = image_data * mask_component
    if centre is None:
        pass
    if centre is not None:
        x0A, y0A = centre
    if (mask is None) and (mask_component is not None):
        mask = mask_component
    else:
        unity_mask = np.ones(image_data.shape) == 1
        omask, mask = unity_mask, unity_mask

    try:
        BGrandom, BGmedian, BGmin, BGstd, BGx0, BGy0 \
            = background_asymmetry(image_data, mask, pre_clean=False)
    except:
        print('Error computing background assymetry.')
        BGrandom, BGmedian, BGmin, BGstd, BGx0, BGy0 = 0.0, 0.0, 0.0, 0.0, x0A, y0A

    x0A0fit, y0A0fit = fmin(assimetria0, (x0A, y0A),
                            args=(image_data, mask,), disp=0)
    x0A1fit, y0A1fit = fmin(assimetria1, (x0A, y0A),
                            args=(image_data, mask,), disp=0)
    A0 = assimetria0((x0A0fit, y0A0fit), image_data, mask) - BGmedian
    A1 = assimetria1((x0A1fit, y0A1fit), image_data, mask)
    results['A_BK_median'] = BGmedian
    results['A0'] = A0
    results['A1'] = A1
    results['x0A0fit'] = x0A0fit
    results['y0A0fit'] = y0A0fit
    results['x0A1fit'] = x0A1fit
    results['y0A1fit'] = y0A1fit
    return(results)


def convex_shape(mask):
    from scipy.spatial import ConvexHull
    indices = np.transpose(np.nonzero(mask))
    hull = ConvexHull(indices)
    convex_area = hull.area
    convex_perimeter = hull.volume
    return(convex_area,convex_perimeter)

def shape_measures(imagename, residualname, z, mask_component=None, sigma_mask=6,
             last_level=3.0, vmin_factor=1.0, plot_catalog=False,data_2D=None,
             npixels=128, fwhm=81, kernel_size=21, dilation_size=None,
             main_feature_index=0, results_final={}, iterations=2,
             fracX=0.10, fracY=0.10, deblend=False, bkg_sub=False,
             bkg_to_sub=None, rms=None,
             apply_mask=True, do_PLOT=False, SAVE=True, show_figure=True,
             mask=None,do_measurements='all',
             add_save_name=''):
    """
    Main function that perform other function calls responsible for
    all relevant calculations on the images.
    """
    cell_size = get_cell_size(imagename)
    """
    One beam area is one element resolution, so we avoind finding sub-components 
    that are smaller than the beam area. 
    """
    min_seg_pixels = beam_area2(imagename)

    if mask is not None:
        mask = mask
        apply_mask = False

    if apply_mask == True:
        _, mask_dilated = mask_dilation(imagename, cell_size=cell_size,
                                        sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations, rms=rms,
                                        PLOT=True)
        mask = mask_dilated
    else:
        mask = None

    if data_2D is not None:
        data_2D = data_2D
    else:
        data_2D = load_fits_data(imagename)

    results_final = None #start dict to store measurements.

    levels, fluxes, agrow, plt, \
        omask2, mask2, results_final = make_flux_vs_std(imagename,
                                                        cell_size=cell_size,
                                                        residual=residualname,
                                                        mask_component=mask_component,
                                                        last_level=last_level,
                                                        sigma_mask=sigma_mask,
                                                        apply_mask=False,
                                                        data_2D = data_2D,
                                                        mask=mask,
                                                        rms=rms,
                                                        vmin_factor=vmin_factor,
                                                        results=results_final,
                                                        show_figure=show_figure,
                                                        bkg_to_sub=bkg_to_sub,
                                                        add_save_name=add_save_name,
                                                        SAVE=SAVE)
    error_petro = False
    if z is not None:
        z = z
    else:
        z = 0.01

    results_final['error_petro'] = error_petro
    if mask_component is not None:
        r, ir = get_profile(load_fits_data(imagename) * mask_component)
    else:
        r, ir = get_profile(imagename)
    #     rpix = r / cell_size
    rpix = r.copy()
    r_list_arcsec = rpix * cell_size
    Rp_arcsec = results_final['C90radii'] * cell_size
    R50_arcsec = results_final['C50radii'] * cell_size
    pix_to_pc = pixsize_to_pc(z=z, cell_size=cell_size)
    r_list_pc = rpix * pix_to_pc
    Rp_pc = results_final['C90radii'] * pix_to_pc
    R50_pc = results_final['C50radii'] * pix_to_pc
    Rp_arcsec = results_final['C90radii'] * cell_size
    R50_arcsec = results_final['C50radii'] * cell_size
    r_list_arcsec = rpix * cell_size
    Rp_arcsec = results_final['C90radii'] * cell_size
    R50_arcsec = results_final['C50radii'] * cell_size


    results_final = cosmo_stats(imagename=imagename, z=z, results=results_final)

    results_final['pix_to_pc'] = pix_to_pc
    results_final['cell_size'] = cell_size
    # if error_petro == True:
    # results_final['area_beam_2Rp'] = area_beam[
    #     int(results_final['C90radii'])]
    # results_final['area_beam_R50'] = area_beam[
    #     int(results_final['C50radii'])]

    df = pd.DataFrame.from_dict(results_final, orient='index').T
    df.to_csv(imagename.replace('.fits', add_save_name + '_area_stats.csv'),
              header=True,
              index=False)

    return (results_final, mask)


def compute_flux_density(imagename, residualname, mask=None, sigma=6,
                         systematic_error_fraction=0.05):
    beam_area = beam_area2(imagename)
    image_data = load_fits_data(imagename)
    rms = mad_std(load_fits_data(residualname))
    if mask is None:
        _, mask = mask_dilation(imagename, show_figure=False,
                                PLOT=False, rms=rms, sigma=sigma,
                                      iterations=2)
    total_flux_density = np.nansum(image_data * mask) / beam_area

    # compute error in flux density
    data_res = load_fits_data(residualname)
    res_error_rms = np.sqrt(np.nansum(
        (abs(data_res * mask - np.nanmean(data_res * mask))) ** 2 * np.nansum(
            mask))) / beam_area

    total_flux_density_residual = np.nansum(data_res * mask) / beam_area

    # Calculate systematic error as a fraction of the total flux density of the image
    systematic_error = systematic_error_fraction * total_flux_density

    # Calculate total flux density error in quadrature
    total_flux_density_error = np.sqrt(
        systematic_error ** 2 + (rms / beam_area) ** 2 + total_flux_density_residual ** 2)

    # res_error_rms = 3 * np.nansum(data_res * mask)/beam_area

    print('-----------------------------------------------------------------')
    print('Estimate of flux error (based on rms of '
          'residual x area): ')
    print('Flux Density = ', total_flux_density * 1000, '+/-',
          total_flux_density_error * 1000, 'mJy')
    print('Fractional error flux = ', total_flux_density_error / total_flux_density)
    print('-----------------------------------------------------------------')
    # print(f"Flux Density = {total_flux_density*1000:.2f} mJy")
    return (total_flux_density, total_flux_density_error)

def measures(imagename, residualname, z, mask_component=None, sigma_mask=6,
             last_level=2.0, vmin_factor=1.0, plot_catalog=False,
             data_2D=None,data_res=None,
             npixels=128, fwhm=81, kernel_size=21, dilation_size=None,
             main_feature_index=0, results_final={}, iterations=2,
             fracX=0.10, fracY=0.10, deblend=False, bkg_sub=False,
             bkg_to_sub=None, rms=None,do_petro=True,
             crop=False, box_size=256,
             apply_mask=True, do_PLOT=False, SAVE=False, show_figure=True,
             mask=None,do_measurements='',compute_A=False,
             add_save_name='',logger=None,verbose=0):
    """
    Main function that perform other function calls responsible for
    all relevant calculations on the images.
    """
    cell_size = get_cell_size(imagename)
    """
    One beam area is one element resolution, so we avoind finding sub-components
    that are smaller than the beam area.
    """
    min_seg_pixels = beam_area2(imagename)

    if mask is not None:
        mask = mask
        apply_mask = False
        omask = None
        if verbose >= 1:
            if logger is not None:
                logger.info(f"  >> Using provided mask.")
            else:
                print('     >> INFO: Using provided mask.')

    if apply_mask == True:
        # if logger is not None:
        #     logger.info(f"  CALC >> Performing mask dilation.")
        # else:
        #     print('     >> CALC: Performing mask dilation.')
        original_mask, mask_dilated = mask_dilation(imagename, cell_size=cell_size,
                                        sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations, rms=rms,
                                        PLOT=do_PLOT,verbose=verbose)
        mask = mask_dilated
        omask = original_mask
    # else:
    #     print('     >> WARN: Not using any mask.')
    #     mask = None
    # """
    # Basic background estimation.
    # """
    # # bkg_ = sep_background(crop_image,apply_mask=False,mask=mask,
    # #                       bw=11, bh=11, fw=12, fh=12)
    # bkg_ = sep_background(imagename,apply_mask=True,mask=None,
    #                       bw=11, bh=11, fw=12, fh=12)
    # bkg_to_sub = bkg_.back()

    if data_2D is not None:
        data_2D = data_2D
    else:
        data_2D = load_fits_data(imagename)

    # if verbose >= 1:
    #     if logger is not None:
    #         logger.info(f"  CALC >> Performing level statistics.")
    #     else:
    #         print('     >> CALC: Performing level statistics.')

    # results_final = level_statistics(img=imagename, cell_size=cell_size,
    #                                 mask_component=mask_component,
    #                                 mask=mask, apply_mask=False,
    #                                 data_2D=data_2D,
    #                                 sigma=sigma_mask, do_PLOT=do_PLOT,
    #                                 results=results_final, bkg_to_sub=bkg_to_sub,
    #                                 show_figure=False,
    #                                 rms=rms,
    #                                 add_save_name=add_save_name,
    #                                 SAVE=SAVE, ext='.jpg')
    if verbose >=1:
        if logger is not None:
            logger.info(f"  CALC >> Computing image properties.")
        else:
            print('     >> CALC: Computing image properties.')

    levels, fluxes, agrow, \
        omask2, mask2, results_final = compute_image_properties(imagename,
                                                        residual=residualname,
                                                        cell_size=cell_size,
                                                        mask_component=mask_component,
                                                        last_level=last_level,
                                                        sigma_mask=sigma_mask,
                                                        apply_mask=False,
                                                        crop=crop,box_size=box_size,
                                                        mask=mask,
                                                        rms=rms,
                                                        data_2D=data_2D,
                                                        data_res=data_res,
                                                        vmin_factor=vmin_factor,
                                                        # results=results_final,
                                                        show_figure=show_figure,
                                                        bkg_to_sub=bkg_to_sub,
                                                        verbose=verbose,
                                                        add_save_name=add_save_name,
                                                        SAVE=SAVE,logger=logger)
    # r_list, area_arr, area_beam, p, \
    #     flux_arr, results_final = do_petrofit(imagename, cell_size,
    #                                           mask_component=mask_component,
    #                                           PLOT=do_PLOT,
    #                                           sigma_mask=sigma_mask,
    #                                           dilation_size=dilation_size,
    #                                           npixels=npixels, fwhm=fwhm,
    #                                           kernel_size=kernel_size,
    #                                           results=results_final,
    #                                           apply_mask=apply_mask)
    omask = omask2.copy()
    error_petro = False
    if do_petro == True:
        try:
            if verbose >= 1:
                if logger is not None:
                    logger.info(f"  CALC >> Computing Petrosian properties.")
                else:
                    print('     >> CALC: Computing Petrosian properties.')
            r_list, area_arr, area_beam, p, flux_arr, error_arr, results_final, cat, \
                segm, segm_deblend, sorted_idx_list = \
                compute_petrosian_properties(data_2D, imagename,
                                                mask_component=mask_component,
                                                global_mask=mask,
                                                source_props=results_final,
                                                apply_mask=False,
                                                # error = ,
                                                sigma_level=sigma_mask,
                                                bkg_sub=bkg_sub, bkg_to_sub=bkg_to_sub,
                                                vmin=vmin_factor, plot=do_PLOT,
                                                deblend=deblend,
                                                fwhm=fwhm, kernel_size=kernel_size,
                                                show_figure=show_figure,
                                                verbose = verbose,
                                                add_save_name=add_save_name,
                                                npixels=npixels,logger=logger)
            error_petro = False
        except Exception as e:
            if logger is not None:
                logger.warning(f"  -->> ERROR when computing Petrosian properties. "
                                f"Will flag error_petro as True.")
            else:
                print("     -->> ERROR when computing Petrosian properties. Will "
                        "flag error_petro as True.")
            error_petro = True
    else:
        error_petro = True

    results_final['error_petro'] = error_petro
    if mask_component is not None:
        centre = nd.maximum_position(data_2D * mask_component)[::-1]
        r, ir = get_profile(load_fits_data(imagename) * mask_component)
    else:
        centre = nd.maximum_position(data_2D)[::-1]
        r, ir = get_profile(data_2D,center=centre)
    #     rpix = r / cell_size
    rpix = r.copy()

    if z is not None:
        z = z
    else:
        z = 0.01
    if error_petro == False:
        r_list_arcsec = r_list * cell_size
        Rp_arcsec = results_final['Rp'] * cell_size
        R50_arcsec = results_final['R50'] * cell_size
        pix_to_pc = pixsize_to_pc(z=z, cell_size=cell_size)
        r_list_pc = r_list * pix_to_pc
        Rp_pc = results_final['Rp'] * pix_to_pc
        R50_pc = results_final['R50'] * pix_to_pc
        Rp_arcsec = results_final['Rp'] * cell_size
        R50_arcsec = results_final['R50'] * cell_size
        r_list_arcsec = r_list * cell_size
        Rp_arcsec = results_final['Rp'] * cell_size
        R50_arcsec = results_final['R50'] * cell_size
    if error_petro == True:
        r_list_arcsec = rpix * cell_size
        Rp_arcsec = results_final['C95radii'] * cell_size
        R50_arcsec = results_final['C50radii'] * cell_size
        pix_to_pc = pixsize_to_pc(z=z, cell_size=cell_size)
        r_list_pc = rpix * pix_to_pc
        Rp_pc = results_final['C95radii'] * pix_to_pc
        R50_pc = results_final['C50radii'] * pix_to_pc
        Rp_arcsec = results_final['C95radii'] * cell_size
        R50_arcsec = results_final['C50radii'] * cell_size
        r_list_arcsec = rpix * cell_size
        Rp_arcsec = results_final['C95radii'] * cell_size
        R50_arcsec = results_final['C50radii'] * cell_size

    
    if do_measurements=='all':
        try:
            x0c, y0c = results_final['x0m'], results_final['y0m']
            if compute_A == True:
                print('--==>> Computing asymetries...')
                results_final = compute_asymetries(imagename=imagename,
                                                mask=mask,
                                                bkg_to_sub=bkg_to_sub,
                                                mask_component=mask_component,
                                                centre=(x0c, y0c),
                                                results=results_final)

            # idx_R50 = np.where(flux_arr < 0.5 * results_final['flux_rp'])[0][-1]
            # idx_Rp = np.where(r_list < 2 * results_final['Rp'])[0][-1]
            # idx_Cradii = np.where(r_list < results_final['Cradii'])[0][-1]
            # idx_C50radii = np.where(r_list < results_final['C50radii'])[0][-1]

            # idx_R50 = np.where(flux_arr < 0.5 * results_final['flux_rp'])[0][-1]
            # flux_arr[idx_R50], r_list[idx_R50], results_final['R50']
            # if mask_component is None:
            if verbose >= 1:
                print('--==>> Computing image statistics...')
            results_final = get_image_statistics(imagename=imagename,
                                                mask_component=mask_component,
                                                mask=mask,
                                                sigma_mask=sigma_mask,
                                                apply_mask=False,
                                                residual_name=residualname, fracX=fracX,
                                                fracY=fracY,
                                                dic_data=results_final,
                                                cell_size=cell_size)
        except:
            pass
    

    results_final = cosmo_stats(imagename=imagename, z=z, results=results_final)

    results_final['pix_to_pc'] = pix_to_pc
    results_final['cell_size'] = cell_size
    if error_petro == False:
        results_final['area_beam_2Rp'] = area_beam[int(results_final['Rp'])]
        results_final['area_beam_R50'] = area_beam[int(results_final['R50'])]
    if error_petro == True:
        results_final['area_beam_2Rp'] = results_final['A95']/beam_area2(imagename)
        results_final['area_beam_R50'] = results_final['A50']/beam_area2(imagename)

    df = pd.DataFrame.from_dict(results_final, orient='index').T
    df.to_csv(imagename.replace('.fits', add_save_name + '_stats.csv'),
              header=True,
              index=False)
    return (results_final, mask, omask)



def check_flux_duplicates(cumulative_flux, radii):
    """
    Check if the array of fluxes contain duplicated values.
    """
    cumulative_flux_norm = cumulative_flux / cumulative_flux[-1]
    
    unique_elements = np.unique(cumulative_flux_norm)
    has_duplicates = len(unique_elements) != len(cumulative_flux_norm)
    if has_duplicates:
        # raise ValueError("Cumulative flux contains duplicates.")
        # print("Cumulative flux contains duplicates.")
        # print(".... attempting removal of duplicates.")
        unique_flux, indices = np.unique(cumulative_flux_norm, return_inverse=True)
        average_radii = np.zeros_like(unique_flux, dtype=float)
        average_cumulative_flux = np.zeros_like(unique_flux, dtype=float)
        # Calculate average of y-values corresponding to each unique x-value
        for i in range(len(unique_flux)):
            average_radii[i] = np.mean(radii[indices == i])
            average_cumulative_flux[i] = np.mean(cumulative_flux[indices == i])
            
        cumulative_flux_norm = unique_flux.copy()
        radii = average_radii.copy()
        cumulative_flux = average_cumulative_flux.copy()
        
    return(cumulative_flux,radii)


def find_fractional_radius(cumulative_flux, radii, 
                           fraction,do_plot=False):
    """
    Calculate the radius that encloses a given fraction of the total flux.
    
    Parameters:
    - cumulative_flux: numpy array of cumulative flux densities
    - radii: numpy array of radii corresponding to the cumulative flux densities
    - fraction: the fraction of the total flux for which to calculate the radius
    
    Returns:
    - radius_at_fraction: the radius that encloses the specified fraction of the total flux
    """
    from scipy.interpolate import interp1d
    # Ensure the fraction is between 0 and 1
    if not (0 < fraction < 1):
        raise ValueError("Fraction must be between 0 and 1.")

    cumulative_flux,radii = check_flux_duplicates(cumulative_flux,radii)
    cumulative_flux_norm = cumulative_flux / cumulative_flux[-1]
    # Calculate the total flux to find the target flux for the given fraction
    total_flux = cumulative_flux_norm[-1]
    total_flux_true= cumulative_flux[-1]
    fraction_flux = total_flux * fraction
    fraction_flux_true = total_flux_true * fraction
    
    # Use spline interpolation to create a smooth function of cumulative flux vs radius
    spline = interp1d(cumulative_flux_norm, radii, kind='linear')
    
    # Use the interpolated function to find the radius for the target flux
    radius_at_fraction = spline(fraction_flux)

    area_at_fraction = radii_to_area(radius_at_fraction)[0]

    if fraction >=0.95:
        dr_err_up = 0.025
    else:
        dr_err_up = 0.05
    if fraction <= 0.1:
        dr_err_lo = 0.025
    else:
        dr_err_lo = 0.05
    
    # Calculate the lower and upper fractions for the error range
    lower_fraction = max(fraction - dr_err_lo, 0)
    upper_fraction = min(fraction + dr_err_up, 1)

    # Calculate the flux values for the lower and upper fractions
    lower_fraction_flux = lower_fraction * total_flux
    upper_fraction_flux = upper_fraction * total_flux

    # Calculate the radii at the lower and upper fractions using the spline interpolation
    lower_radius = spline(lower_fraction_flux)
    upper_radius = spline(upper_fraction_flux)

    # Calculate the error range as the difference between the lower and upper radii
    radius_error_range = (abs(upper_radius-radius_at_fraction),abs(lower_radius-radius_at_fraction))
    # print(f"Fractional Radius with error range: {radius_at_fraction:.2f} +/- {radius_error_range[1]:.2f}/{radius_error_range[0]:.2f}")
    fraction_radius_error_mean = np.mean(radius_error_range)
    # Calculate the areas at the lower and upper radii
    lower_area = radii_to_area(lower_radius)[0]
    upper_area = radii_to_area(upper_radius)[0]

    # Calculate the area error range as the difference between the lower and upper areas
    area_error_range = (abs(upper_area - area_at_fraction), abs(lower_area - area_at_fraction))
    # print(f"Fractional Area with error range: {area_at_fraction:.2f} +/- {area_error_range[1]:.2f}/{area_error_range[0]:.2f}")
    # print(area_error_range[0],area_error_range[1],np.mean(area_error_range))
    fraction_area_error_mean = np.mean(area_error_range)
    if do_plot:
        # Plot the data points
        plt.plot(radii, cumulative_flux_norm, 'o', label='Data Points')
        
        # Create a smooth range of flux values for plotting the interpolation
        flux_range = np.linspace(np.min(cumulative_flux_norm), np.max(cumulative_flux_norm), 1000)
        
        # Plot the spline interpolation
        plt.plot(spline(flux_range), flux_range, '-', label='Spline Interpolation')
        
        # Plot a cross at the fractional radius
        plt.plot(radius_at_fraction, fraction_flux, 'rx', markersize=10, label=f'Fractional Radius ({fraction*100:.1f}%)')
        plt.axvline(lower_radius)
        plt.axvline(radius_at_fraction)
        plt.axvline(upper_radius)
        
        # Add labels and legend
        plt.ylabel('Cumulative Flux')
        plt.xlabel('Radius')
        plt.legend()
        
        # Show the plot
        plt.show()
    
    return(radius_at_fraction, fraction_flux_true, area_at_fraction, 
           fraction_radius_error_mean, fraction_area_error_mean)

def R50_to_fwhm(R50,R50_err,q_ratio=1,scale = 1):
    """
    Calculate the beam size in arcseconds for a given R50 value.
    (Lucatelli et al., 2024)
    """
    theta_maj = 2*R50*scale
    theta_min = theta_maj*q_ratio
    
    theta_maj_err = 2*R50_err*scale
    theta_min_err = theta_maj_err*q_ratio
    
    theta = (theta_maj*theta_min)**0.5
    theta_err = (theta_maj_err*theta_min_err)**0.5
    return(theta_maj,theta_maj_err,
           theta_min,theta_min_err,
           theta,theta_err)


def concentration_index_errors(C80radii, C20radii, 
                               C90radii, C50radii, 
                               C80radii_err, C20radii_err, 
                               C90radii_err, C50radii_err):
    """
    Given fractional radii and their errors, calculate the concentration indices
    and their errors.
    
    Parameters:
    - C80radii: the radius enclosing 80% of the total flux
    - C20radii: the radius enclosing 20% of the total flux
    - C90radii: the radius enclosing 90% of the total flux
    - C50radii: the radius enclosing 50% of the total flux
    - C80radii_err: the error in the radius enclosing 80% of the total flux
    - C20radii_err: the error in the radius enclosing 20% of the total flux
    - C90radii_err: the error in the radius enclosing 90% of the total flux
    - C50radii_err: the error in the radius enclosing 50% of the total flux
    
    Returns:
    - C1: the concentration index C1
    - C1_err: the error in the concentration index C1
    - C2: the concentration index C2
    - C2_err: the error in the concentration index C2
    """
    C1 = np.log10(C80radii / C20radii)
    C2 = np.log10(C90radii / C50radii)

    dC1_dC80radii = 1 / (C80radii * np.log(10))
    dC1_dC20radii = -1 / (C20radii * np.log(10))
    C1_err = np.sqrt((dC1_dC80radii * C80radii_err)**2 + (dC1_dC20radii * C20radii_err)**2)

    dC2_dC90radii = 1 / (C90radii * np.log(10))
    dC2_dC50radii = -1 / (C50radii * np.log(10))
    C2_err = np.sqrt((dC2_dC90radii * C90radii_err)**2 + (dC2_dC50radii * C50radii_err)**2)
    # print(f"C1 = {C1:.2f} +/- {C1_err:.2f}")
    # print(f"C2 = {C2:.2f} +/- {C2_err:.2f}")
    return C1, C1_err, C2, C2_err


# def deprecated(new_func_name):
#     def decorator(old_func):
#         def wrapper(*args, **kwargs):
#             warnings.warn(f"'{old_func.__name__}' is deprecated and will be "
#                           f"removed in a future version. "
#                           f"Use '{new_func_name}' instead.",
#                           category=DeprecationWarning, stacklevel=2)
#             return old_func(*args, **kwargs)
#         return wrapper
#     return decorator



def compute_image_properties(img, residual, cell_size=None, mask_component=None,
                             aspect=1, last_level=1.5, mask=None, 
                             bins = 64,
                             data_2D=None,data_res=None,
                             dilation_size=None, iterations=2,
                             dilation_type='disk',do_fit_ellipse=True,
                             sigma_mask=6, rms=None, results=None, bkg_to_sub=None,
                             apply_mask=True, vmin_factor=3, vmax_factor=0.5,
                             crop=False, box_size=256,
                             SAVE=True, add_save_name='', show_figure=True,
                             save_csv = False,
                             ext='.jpg',logger=None,verbose=0):
    """
    Params
    ------
    img : str
        The path to the image to be analyzed.
    residual : str
        The path to the residual image associated with the image.
    cell_size : float, optional
        The default is None. The size of the pixel in arcseconds.
        If None, `get_cell_size` will be used to get from the header.
    mask_component : 2D np array, optional
        The default is None. If not None, this is the mask for a specific component
        when performing a multi-component source analysis.
    aspect : float, optional (experimental)
        The default is 1. The aspect ratio of the image for plotting.
    last_level : float, optional
        New threshold level (as multiple of sigma_mad) to be used inside the existing mask.
    bins : int, optional
        Number of bins to slice intensity levels of the image. 
    mask : 2D np array, optional
        The default is None. If not None, the function `mask_dilation` will determine it with
        default parameters.
    data_2D : 2D np array, optional
        The default is None. This can be used to pass a 2D array directly to the function.
        For example, when calculating properties from an array without reading from a file,
        e.g. a model image, you can use this. But, to obtain meaningful physical units,
        you must provide the corresponding image file where this array was derived from.
    dilation_size : int, optional
        The default is None. The size of the dilation to be used in the mask dilation.
        If None, the default value will be the size of the restoring beam.
    iterations : int, optional
        The default is 2. The number of iterations to be used in the mask dilation.
        If signs of over-dilation are present, you can set to 1.
    dilation_type : str, optional
        The default is 'disk'. The type of dilation to be used in the mask dilation.
    sigma_mask : float, optional
        The default is 6. The sigma level to be used in the mask dilation.
    rms : float, optional
        The default is None. The rms value to be used in the mask dilation.
        If None, the function `mad_std` will be used to calculate it from the residual image,
        if provided. If the residual image is not provided, the function will use the image itself.
        But, in that case, the result may not be accurate (overestimated) if the image size is
        comparable in size to the size of the source structure.
    results : dict, optional
        The default is None. A dictionary to store the results.
        You can pass an existing external dictionary, so the results will be appended to it.
    bkg_to_sub : 2D np array, optional (EXPERIMENTAL)
        The default is None. The background to be subtracted from the image.
    apply_mask : bool, optional
        The default is True. If True, the mask dilation will be calcualted from the image.
    vmin_factor : float, optional
        The default is 3. The factor (as a multiple of sigma_mad) to be used in the vmin
        calculation for the image plot.
    vmax_factor : float, optional
        The default is 0.5. The factor (as a multiple of peak brightness) to be used in the vmax
        calculation for the image plot.
    crop : bool, optional
        The default is False. If True, the image will be cropped to a box_size.
    box_size : int, optional
        The default is 256. The size of the box to be used in the image cropping.
    SAVE : bool, optional
        The default is True. If True, the image plot will be saved to a file.
    add_save_name : str, optional
        The default is ''. A string to be added to the image plot file name.
    show_figure : bool, optional
        The default is True. If True, the image plot will be shown.
    ext : str, optional
        The default is '.jpg'. The file extension to be used in the image plot file name.
    logger : logging.Logger, optional
        The default is None. A logger object to be used to log messages.

    """
    if results is None:
        results = {}
        results['#imagename'] = os.path.basename(img)

    from skimage.draw import disk
    if data_2D is not None:
        g_ = data_2D
    else:
        g_ = load_fits_data(img)
    # res_ = load_fits_data(residual)
    
    g = np.nan_to_num(g_.copy(),nan=0)
    # res = res_.copy()

    if bkg_to_sub is not None:
        g = g - bkg_to_sub


    scale_units = 'arcsec'
    if cell_size is None:
        try:
            cell_size = get_cell_size(img)
        except:
            scale_units = 'px'
            cell_size = 1.0

    g_hd = imhead(img)
    freq = g_hd['refval'][2] / 1e9
    # print(freq)
    omaj = g_hd['restoringbeam']['major']['value']
    omin = g_hd['restoringbeam']['minor']['value']
    # beam_area_ = beam_area(omaj, omin, cellsize=cell_size)
    beam_area_ = beam_area2(img)

    if rms is not None:
        std = rms
    else:
        if residual is not None:
            std = mad_std(np.nan_to_num(load_fits_data(residual), nan=0))
        else:
            std = mad_std(np.nan_to_num(g_.copy(),nan=0))
    if mask is not None:
        mask = mask
        omask = mask
        g = g * mask  # *(g_>3*mad_std(g_)) + 0.1*mad_std(g_)
        # res = res * mask  # *(g_>3*mad_std(g_)) + 0.1*mad_std(g_)
        # total_flux = np.sum(g_ * (g_ > 3 * std)) / beam_area_
        total_flux = np.sum(g) / beam_area_
        total_flux_nomask = np.sum(g_) / beam_area_
        apply_mask = False

    if apply_mask == True:
        if data_res is None:
            omask, mask = mask_dilation(img, sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations,
                                        dilation_type=dilation_type,
                                        show_figure=show_figure,
                                        verbose=verbose)
        else:
            omask, mask = mask_dilation(data_2D, sigma=sigma_mask,
                                        dilation_size=dilation_size,
                                        iterations=iterations,
                                        dilation_type=dilation_type,
                                        show_figure=show_figure,
                                        verbose=verbose)
    

        # eimshow(mask)
        g = g * mask  # *(g_>3*mad_std(g_)) + 0.1*mad_std(g_)
        # res = res * mask
        # total_flux = np.sum(g_ * (g_ > 3 * std)) / beam_area_
        total_flux = np.nansum(g) / beam_area_
        total_flux_nomask = np.nansum(g_) / beam_area_
        # g = g_

    if mask_component is not None:
        #
        g = g * mask_component

        # res = res * mask_component

        total_flux_nomask = np.nansum(g) / beam_area_
        # total_flux = np.sum(g * (g > 3 * std)) / beam_area_
        total_flux = np.nansum(g) / beam_area_
        if mask is not None:
            mask = mask * mask_component
            omask = mask * mask_component
        else:
            mask = mask_component
            omask = mask_component

    # if (apply_mask is None) and  (mask is None):
    #     mask = mask_component
        # g = g_
    if (mask_component is None) and (apply_mask == False) and (mask is None):
        # g = g_
        total_flux = np.nansum(g * (g > 3 * std)) / beam_area_
        total_flux_nomask = np.nansum(g) / beam_area_

        unity_mask = np.ones(g.shape) == 1
        omask, mask = unity_mask, unity_mask

    # total_flux_nomask = np.sum(g_) / beam_area_
    # if (mask is None) and (mask_component is not None):
    #     mask = mask_component


    
    results['total_flux_nomask'] = total_flux_nomask
    results['total_flux_mask'] = total_flux
    results['peak_of_flux'] = np.nanmax(g)
    results['bmajor'] = omaj
    results['bminor'] = omin
    results['freq'] = freq*1e9
    
    
    # flux_mc, flux_error_m = mc_flux_error(img, g,
    #               res,
    #               num_threads=6, n_samples=1000)
    #
    # results['flux_mc'] = flux_mc
    # results['flux_error_mc'] = flux_error_m

    if residual is not None:
        systematic_error_fraction = 0.05
        if data_res is None: 
            data_res = load_fits_data(residual)
        else:
            data_res = data_res
        
        data_res = np.nan_to_num(data_res,nan=0)
        flux_res_error = 3 * np.sum(data_res * mask) / beam_area_
        # rms_res =imstat(residual_name)['flux'][0]

        total_flux_density_residual = np.nansum(data_res * mask) / beam_area_
        res_error_rms =np.sqrt(
            np.nansum((abs(data_res * mask -
                        np.nanmean(data_res * mask))) ** 2 * np.nansum(mask))) / beam_area_

        # Calculate systematic error as a fraction of the total flux density of the image
        systematic_error = systematic_error_fraction * results['total_flux_mask']

        # Calculate total flux density error in quadrature
        total_flux_density_error = np.sqrt(
            systematic_error ** 2 + (std / beam_area_) ** 2 + total_flux_density_residual ** 2)
        # print('total_flux_density_residual = ',total_flux_density_residual)
        # print('systematic_error = ', systematic_error)
        # print('(std / beam_area_) = ', (std / beam_area_))
        if verbose >= 1:
            print('-----------------------------------------------------------------')
            print('Flux Density and error (based on rms of residual x area): ')
            print(f"Flux Density = {results['total_flux_mask'] * 1000:.2f} "
                  f"+/- {res_error_rms * 1000:.2f} mJy")
            print(f"Fractional error = {(res_error_rms / results['total_flux_mask']):.2f}")
            print('Flux Density and error (quadrature |fract_err + res_err + rms|): ')
            print(f"Flux Density = {results['total_flux_mask'] * 1000:.2f} "
                  f"+/- {total_flux_density_error * 1000:.2f} mJy")
            print(f"Fractional error = {(total_flux_density_error / results['total_flux_mask']):.2f}")
            print('-----------------------------------------------------------------')

        results['max_residual'] = np.nanmax(data_res * mask)
        results['min_residual'] = np.nanmin(data_res * mask)
        results['flux_residual'] = total_flux_density_residual
        results['flux_error_res'] = abs(flux_res_error)
        results['flux_error_res_2'] = abs(res_error_rms)
        results['flux_error_res_3'] = abs(total_flux_density_error)
        results['mad_std_residual'] = mad_std(data_res)
        results['rms_residual'] = rms_estimate(data_res)

    try:
        # this should not be linspace, should be spaced in a logarithmic sense!!
        # print('np.nanmax(g)',np.nanmax(g))
        # print('last_level * std',last_level * std)
        # plt.figure()
        # plt.imshow(g)
        # plt.show()
        levels = np.geomspace(np.nanmax(g), last_level * std,bins)
        levels2 = np.geomspace(np.nanmax(g), last_level * std,100)
        levels_ellipse = np.geomspace(np.nanmax(g), last_level * std,32)
    except:
        levels = np.geomspace(np.nanmax(g), last_level * np.nanstd(g),bins)
        levels2 = np.geomspace(np.nanmax(g), last_level * np.nanstd(g),100)
        levels_ellipse = np.geomspace(np.nanmax(g), last_level * np.nanstd(g),32)

    fluxes = []
    areas = []
    for i in range(len(levels)):
        if i == 0:
            # condition = (g >= levels[i])
            # flux = g[np.where(condition)].sum() / beam_area_
            condition = (g >= levels[i])
            flux = np.nansum(g * (condition)) / beam_area_

            area = np.nansum(condition)
            fluxes.append(flux)
            areas.append(area)
        else:
            # condition = (g < levels[i - 1]) & (g >= levels[i])
            # flux = g[np.where(condition)].sum() / beam_area_
            condition = ((g < levels[i - 1]) & (g >= levels[i]))
            flux = np.nansum(g * (condition)) / beam_area_
            area = np.nansum(condition)
            # area = np.sum((g >= levels[i]))
            fluxes.append(flux)
            areas.append(area)
    fluxes = np.asarray(fluxes)
    areas = np.asarray(areas)
    #     print()
    #     plt.scatter(levels[:],fluxes[:])#/np.sum(fluxes))

    """
    Growth curve.
    """
    
    agrow = areas.copy()
    results['total_flux_levels'] = np.nansum(fluxes)


    # agrow_beam = agrow / beam_area_
    Lgrow = np.nancumsum(fluxes)
    # Lgrow_norm = Lgrow / fluxes.sum()
    _radii = []
    for i in range(0, len(levels)):
        mask_at_level_i = g >= levels[i]
        circular_radii = np.sqrt(np.nansum(mask_at_level_i) / np.pi)
        # print(circular_radii)
        _radii.append(circular_radii)
    radii = np.asarray(_radii)
    # print(Lgrow_norm)
    Lgrow,radii = check_flux_duplicates(Lgrow,radii)
    Lgrow_norm = Lgrow / np.nansum(fluxes)

    ####################################################################
    ## THIS NEEDS A BETTER IMPLEMENTATION USING SPLINE!!!!  ############
    ####################################################################
    mask_L20 = Lgrow_norm < 0.2
    mask_L50 = Lgrow_norm < 0.5
    mask_L80 = Lgrow_norm < 0.8
    try:
        mask_L90 = (Lgrow_norm > 0.89) & (Lgrow_norm < 0.91)
        mask_L95 = (Lgrow_norm > 0.95) & (Lgrow_norm < 0.97)
    except:
        mask_L90 = (Lgrow_norm > 0.85) & (
                    Lgrow_norm < 0.95)  # in case there not enough pixels
        mask_L95 = (Lgrow_norm > 0.92) & (
                    Lgrow_norm < 0.97)  # in case there not enough pixels

    mask_L20_idx = [i for i, x in enumerate(mask_L20) if x]
    mask_L50_idx = [i for i, x in enumerate(mask_L50) if x]
    mask_L80_idx = [i for i, x in enumerate(mask_L80) if x]
    mask_L90_idx = [i for i, x in enumerate(mask_L90) if x]
    mask_L95_idx = [i for i, x in enumerate(mask_L95) if x]
    """
    The tries and exceptions bellow are temporary solutions to avoid 
    errors when there are not enough pixels or no signal at all, when computing 
    the growth curve.    
    """
    try:
        sigma_20 = levels[mask_L20_idx[-1]]
        flag20 = False
    except:
        flag20 = True
        try:
            sigma_20 = levels[mask_L50_idx[-1]]
        except:
            sigma_20 = last_level * std

    try:
        sigma_50 = levels[mask_L50_idx[-1]]
        flag50 = False
    except:
        flag50 = True
        sigma_50 = last_level * std
    try:
        sigma_80 = levels[mask_L80_idx[-1]]
        sigma_90 = levels[mask_L90_idx[-1]]
        sigma_95 = levels[mask_L95_idx[-1]]
        flag9095 = False
    except:
        sigma_80 = last_level * std
        sigma_90 = last_level * std
        sigma_95 = last_level * std
        flag9095 = True

    # up_low_std = (np.max(gal) / (3 * mad_std(gal))) * 0.7
    # inner_shell_mask = ((gal<(up_low_std+5)*mad_std(gal)) & (gal>up_low_std*mad_std(gal)))
    # inner_shell_mask = ((g <sigma_50) & (g >sigma_50*0.95))
    # outer_shell_mask = ((g < sigma_90) & (g > sigma_90 * 0.95))

    inner_mask = (g > (sigma_50)) * mask
    outer_mask90 = (g > (sigma_90)) * mask
    inner_perimeter = perimeter_crofton(inner_mask, 4)
    outer_perimeter90 = perimeter_crofton(outer_mask90, 4)
    outer_perimeter = perimeter_crofton(mask, 4)
    
    try:
        results['convex_error_flag'] = False
        convex_properties = convex_morpho(g, mask, do_plot=False)        
        results['PA_convex'] = convex_properties['PA_convex']
        results['q_convex'] = convex_properties['q_convex']
        results['centroid_convex'] = convex_properties['centroid_convex']
        results['convex_major_diameter'] = convex_properties['major_diameter']
        results['convex_minor_diameter'] = convex_properties['minor_diameter']
    except:
        results['convex_error_flag'] = True
        results['PA_convex'] = 0.0
        results['q_convex'] = 0.0
        results['centroid_convex'] = [0,0]
        results['convex_major_diameter'] = 0.0
        results['convex_minor_diameter'] = 0.0
    

    # Geometry
    x0max, y0max = peak_center(g * mask)
    results['x0'], results['y0'] = x0max, y0max
    # determine momentum centres.
    x0m, y0m, q_mom, PAdeg_mom = momenta(g * mask, PArad_0=None, q_0=None)
    results['x0m'], results['y0m'] = x0m, y0m
    results['q_mom'], results['PAdeg_mom'] = q_mom, PAdeg_mom

    if do_fit_ellipse:
        try:
            # some geometrical measures
            # calculate PA and axis-ratio
            # region_split = [i for i, x in enumerate(levels > sigma_50) if x][-1]
            region_split = [i for i, x in enumerate(levels_ellipse > sigma_50) if x][-1]
            
            PA, q, x0col, y0col, PAm, qm, \
                PAmi, qmi, PAmo, qmo, \
                x0median, y0median, \
                x0median_i, y0median_i, \
                x0median_o, y0median_o,profiles = cal_PA_q(g * mask, Isequence=levels_ellipse,
                                                region_split=region_split,
                                                # SAVENAME=img.replace('.fits','_ellipsefit') + ext
                                                SAVENAME=img.replace('.fits','_ellipsefit')
                                                )
        except:
                PA, q, x0col, y0col, PAm, qm = 0.0, 0.0, x0max,y0max, 0.0, 0.0
                PAmi, qmi, PAmo, qmo = 0.0, 0.0, 0.0, 0.0
                x0median, y0median = x0max,y0max
                x0median_i, y0median_i = x0max,y0max
                x0median_o, y0median_o = x0max,y0max
                profiles = None

        results['PA'], results['q'] = PA, q
        results['PAm'], results['qm'] = PAm, qm
        results['PAmi'], results['qmi'] = PAmi, qmi
        results['PAmo'], results['qmo'] = PAmo, qmo
        results['x0m_i'], results['y0m_i'] = x0median_i, y0median_i
        results['x0m_o'], results['y0m_o'] = x0median_o, y0median_o

        vx = results['x0'] - results['x0m']
        vy = results['y0'] - results['y0m']
        TvPA, Tvlenght = trail_vector(vx=vx, vy=vy, v0=np.asarray([1, 0]))
        results['TvPA'] = TvPA
        results['Tvlenght'] = Tvlenght

    try:
        L20_norm = Lgrow_norm[mask_L20_idx[-1]]  # ~ 0.2
        L20 = Lgrow[mask_L20_idx[-1]]
    except:
        try:
            L20_norm = Lgrow_norm[mask_L50_idx[-1]]  # ~ 0.2
            L20 = Lgrow[mask_L50_idx[-1]]
        except:
            L20 = 0.0
            L20_norm = 0.9999
    try:
        L50_norm = Lgrow_norm[mask_L50_idx[-1]]  # ~ 0.5
        L50 = Lgrow[mask_L50_idx[-1]]
    except:
        L50_norm = 0.9999
        L50 = 0.0
    try:
        L80_norm = Lgrow_norm[mask_L80_idx[-1]]  # ~ 0.8
        L80 = Lgrow[mask_L80_idx[-1]]
    except:
        L80_norm = 0.9999
        L80 = 0.0
    try:
        """
        Not enough pixels
        """
        L90_norm = Lgrow_norm[mask_L90_idx[-1]]  # ~ 0.9
        L90 = Lgrow[mask_L90_idx[-1]]
        L95_norm = Lgrow_norm[mask_L95_idx[-1]]  # ~ 0.9
        L95 = Lgrow[mask_L95_idx[-1]]
        flagL9095 = False

    except:
        flagL9095 = True
        try:
            try:
                L90_norm = Lgrow_norm[mask_L80_idx[-1]]  # ~ 0.9
                L90 = Lgrow[mask_L80_idx[-1]]
                L95_norm = Lgrow_norm[mask_L80_idx[-1]]  # ~ 0.9
                L95 = Lgrow[mask_L80_idx[-1]]
            except:
                L90_norm = Lgrow_norm[-1]  # ~ 0.9
                L90 = Lgrow[-1]
                L95_norm = Lgrow_norm[-1]  # ~ 0.9
                L95 = Lgrow[-1]

        except:
            L90_norm = 0.9999
            L90 = fluxes.sum()
            L95_norm = 0.9999
            L95 = fluxes.sum()


    try:
        TB20 = T_B(omaj, omin, freq, L20)
    except:
        TB20 = 0.0
    TB50 = T_B(omaj, omin, freq, L50)
    TB80 = T_B(omaj, omin, freq, L80)
    TB90 = T_B(omaj, omin, freq, L90)

    TB = T_B(omaj, omin, freq, total_flux)

    levels_20 = np.asarray([sigma_20])
    levels_50 = np.asarray([sigma_50])
    levels_80 = np.asarray([sigma_80])
    levels_90 = np.asarray([sigma_90])
    levels_95 = np.asarray([sigma_95])

    levels_3sigma = np.asarray([3 * std])

    g20 = ((g * mask) > sigma_20)
    g50 = ((g * mask) > sigma_50)
    g80 = ((g * mask) > sigma_80)
    g90 = ((g * mask) > sigma_90)
    g95 = ((g * mask) > sigma_95)

    try:
        C20radii, L20, npix20, C20radii_err, npix20_err = \
            find_fractional_radius(Lgrow, radii, fraction = 0.20)
        L20_norm = 0.20
    except:
        if logger is not None:
            logger.warning(f" !!==> Not enough pixels to calculate L20-R20."
                        f"       Using flux/size of first pixel.")
        C20radii = radii[0]
        L20 = Lgrow[0]
        L20_norm = Lgrow_norm[0]
        npix20 = radii_to_area(C20radii)[0]
        C20radii_err, npix20_err = C20radii/2, npix20/2
    
    try:
        if logger is not None:
            logger.warning(f" !!==> Not enough pixels to calculate L50-R50."
                           f"       Using flux/size of first pixel.")
        C50radii, L50, npix50, C50radii_err, npix50_err = \
            find_fractional_radius(Lgrow, radii, fraction = 0.50)
        L50_norm = 0.50
    except:
        C50radii = radii[0]
        L50 = Lgrow[0]
        L50_norm = Lgrow_norm[0]
        npix50 = radii_to_area(C50radii)[0]
        C50radii_err, npix50_err = C50radii/2, npix50/2
    
    C80radii, L80, npix80, C80radii_err, npix80_err = \
        find_fractional_radius(Lgrow, radii, fraction = 0.80)
    C90radii, L90, npix90, C90radii_err, npix90_err = \
        find_fractional_radius(Lgrow, radii, fraction = 0.90)
    C95radii, L95, npix95, C95radii_err, npix95_err = \
        find_fractional_radius(Lgrow, radii, fraction = 0.95)
    C99radii, L99, npix99, C99radii_err, npix99_err = \
        find_fractional_radius(Lgrow, radii, fraction = 0.99)

    # _frac_fluxes = np.asarray([0.05, 0.1 , 0.15, 0.2 , 0.25, 0.3 , 0.35, 0.4 , 0.45, 0.5 ,
    #                            0.55, 0.6 , 0.65, 0.7 , 0.75, 0.8 , 0.85, 0.9 , 0.95, 0.99])
    # _frac_fluxes = np.asarray([0.85, 0.9 , 0.95, 0.99])
    # frac_radii  = {}
    # frac_fluxes = {}
    # for frac_flux in _frac_fluxes:
    #     # try:
    #     radius_percent, _, _ = find_fractional_radius(Lgrow, radii, fraction = frac_flux)
    #     frac_radii[f'R{int(frac_flux*100)}'] = radius_percent
    #     frac_fluxes[f'S{int(frac_flux*100)}'] = frac_flux
        # except:
        #     pass
    # print(frac_radii.keys())
    # print(frac_fluxes)
    # print(np.sqrt(((frac_radii['R99'] - frac_radii['R95'])**2.0 + (frac_radii['R90'] - frac_radii['R95'])**2.0  + (frac_radii['R95'] - frac_radii['R85'])**2.0)/3.0))
    
    # plt.figure()
    # plt.plot(frac_fluxes,frac_radii,'o')
    # plt.show()
    # np.mean(frac_radii)

    L80_norm = 0.80
    L90_norm = 0.90
    L95_norm = 0.95

    
    A20 = npix20 / beam_area_
    A50 = npix50 / beam_area_
    A80 = npix80 / beam_area_
    A90 = npix90 / beam_area_
    A95 = npix95 / beam_area_
    A99 = npix99 / beam_area_
    A20_err = npix20_err / beam_area_
    A50_err = npix50_err / beam_area_
    A80_err = npix80_err / beam_area_
    A90_err = npix90_err / beam_area_
    A95_err = npix95_err / beam_area_
    A99_err = npix99_err / beam_area_

    # A20, C20radii, npix20 = estimate_area((g > sigma_20) * mask, cell_size, omaj,
    #                                       omin)
    # A50, C50radii, npix50 = estimate_area((g > sigma_50) * mask, cell_size, omaj,
    #                                       omin)
    # A80, C80radii, npix80 = estimate_area((g > sigma_80) * mask, cell_size, omaj,
    #                                       omin)
    # A90, C90radii, npix90 = estimate_area((g > sigma_90) * mask, cell_size, omaj,
    #                                       omin)
    # A95, C95radii, npix95 = estimate_area((g > sigma_95) * mask, cell_size, omaj,
    #                                       omin)

    # if flag20 == True:
    #     if flag50 == False:
    #         A20, C20radii, npix20 = A50 / 2, C50radii / 2, npix50 / 2
    #     else:
    #         try:
    #             A20, C20radii, npix20 = A80 / 3, C80radii / 3, npix80 / 3
    #         except:
    #             A20, C20radii, npix20 = A95 / 4, C95radii / 4, npix95 / 4

    try:
        results['conv_P20'], results['conv_A20'] = convex_shape(g20)
    except:
        results['conv_P20'], results['conv_A20'] = 5,5

    try:
        results['conv_P50'], results['conv_A50'] = convex_shape(g50)
        results['conv_P80'], results['conv_A80'] = convex_shape(g80)
        results['conv_P90'], results['conv_A90'] = convex_shape(g90)
        results['conv_P95'], results['conv_A95'] = convex_shape(g95)
        results['conv_PT'], results['conv_AT'] = convex_shape(mask)
    except:
        results['conv_P50'], results['conv_A50'] = 5,5
        results['conv_P80'], results['conv_A80'] = 5,5
        results['conv_P90'], results['conv_A90'] = 5,5
        results['conv_P95'], results['conv_A95'] = 5, 5
        results['conv_PT'], results['conv_AT'] = 5,5


    # This is a more robust computation of source size.
    try:
        R20med,R20mean,R20std = calculate_radii(g, g20)
    except:
        R20med, R20mean, R20std = 1,1,1
    try:
        R50med,R50mean,R50std = calculate_radii(g, g50)
    except:
        R50med, R50mean, R50std = 1,1,1
    try:
        R80med,R80mean,R80std = calculate_radii(g, g80)
    except:
        R80med, R80mean, R80std = 1,1,1
    try:
        R90med,R90mean,R90std = calculate_radii(g, g90)
    except:
        R90med, R90mean, R90std = 1,1,1
    try:
        R95med,R95mean,R95std = calculate_radii(g, g95)
    except:
        R95med, R95mean, R95std = 1,1,1
    try:
        RTmed,RTmean,RTstd = calculate_radii(g, mask)
    except:
        RTmed, RTmean, RTstd = 1,1,1

    results['R20med'], results['R20mean'], \
        results['R20std'] = R20med,R20mean,R20std
    results['R50med'], results['R50mean'], \
        results['R50std'] = R50med,R50mean,R50std
    results['R80med'], results['R80mean'], \
        results['R80std'] = R80med,R80mean,R80std
    results['R90med'], results['R90mean'], \
        results['R90std'] = R90med,R90mean,R90std
    results['R95med'], results['R95mean'], \
        results['R95std'] = R95med,R95mean,R95std
    results['RTmed'], results['RTmean'], \
        results['RTstd'] = RTmed,RTmean,RTstd

    # # print(C20radii, C50radii, C80radii, C90radii)
    # C1 = np.log10(C80radii / C20radii)
    # C2 = np.log10(C90radii / C50radii)
    
    C1, C1_err, C2, C2_err = concentration_index_errors(C80radii, C20radii, C90radii, C50radii, 
                                                        C80radii_err, C20radii_err, C90radii_err, C50radii_err)

    AC1 = np.log10(A80 / A20)
    AC2 = np.log10(A90 / A50)

    CAC1 = np.log10(area_to_radii(results['conv_A80']) / area_to_radii(results['conv_A20']))
    CAC2 = np.log10(area_to_radii(results['conv_A90']) / area_to_radii(results['conv_A50']))



    area_total, Cradii, npix_total = estimate_area(mask, cell_size, omaj, omin)
    o_area_total, o_Cradii, o_npix_total = estimate_area(omask, cell_size, omaj,
                                                         omin)

    mask_outer = (g < sigma_50) & (g > last_level * std) * mask
    A50_100, C50_100radii, npix50_100 = estimate_area(mask_outer, cell_size,
                                                      omaj, omin)
    # gaussianity = 1 / (L50 / g.max())
    gaussianity_L50 = g.max() / L50
    gaussianity = sigma_50 / g.max()

    mask_outer_full = (g < sigma_50) * mask
    #     plt.imshow(mask_outer_full)
    A50_full, C50_full_radii, npix50_full = estimate_area(mask_outer_full,
                                                          cell_size, omaj, omin)

    radii_ratio = C50radii / C50_100radii
    radii_ratio_full = C50radii / C50_full_radii
    area_ratio = A50 / A50_100
    area_ratio_full = A50 / A50_full

    results['max'] = np.nanmax(g)
    results['std_image'] = std
    try:
        results['std_residual'] = mad_std(load_fits_data(residual))
    except:
        pass

    
    results['L50'] = L50
    results['sigma_50'] = sigma_50
    results['sigma_90'] = sigma_90
    results['sigma_95'] = sigma_95
    results['o_area_total'] = o_area_total
    results['o_Cradii'] = o_Cradii
    results['o_npix_total'] = o_npix_total
    results['area_total'] = area_total
    results['Cradii'] = Cradii
    results['npix_total'] = npix_total
    results['inner_perimeter'] = inner_perimeter
    results['outer_perimeter'] = outer_perimeter
    results['outer_perimeter90'] = outer_perimeter90

    results['TB20'] = TB20
    results['TB50'] = TB50
    results['TB80'] = TB80
    results['TB'] = TB

    #     results['nbeams_total'] = nbeams_total

    results['C1'] = C1
    results['C2'] = C2
    results['C1_err'] = C1_err
    results['C2_err'] = C2_err
    results['AC1'] = AC1
    results['AC2'] = AC2
    results['CAC1'] = CAC1
    results['CAC2'] = CAC2

    results['L20'] = L20
    results['A20'] = A20
    results['C20'] = sigma_20 / std
    results['C20radii'] = C20radii
    results['npix20'] = npix20
    results['C20radii_err'] = C20radii_err
    results['npix20_err'] = npix20_err

    results['L50'] = L50
    results['A50'] = A50
    results['C50'] = sigma_50 / std
    results['C50radii'] = C50radii
    results['C50radii_err'] = C50radii_err
    results['npix50'] = npix50
    results['npix50_err'] = npix50_err
    
    results['C50theta'] = C50radii * np.sqrt(6 * np.log(2))
    results['C50theta_err'] = C50radii_err * np.sqrt(6 * np.log(2))
    results['C50theta_dec'] = np.sqrt(results['C50theta']**2.0 - (results['bmajor']/cell_size)**2.0)
    # results['C50theta_dec_err'] = np.sqrt(results['C50theta_err']**2.0 - (results['bmajor']/cell_size)**2.0)
    results['C50theta_dec_err'] = results['C50theta_dec'] * results['C50theta_err']
    
    
    

    results['L80'] = L80
    results['A80'] = A80
    results['C80'] = sigma_80 / std
    results['C80radii'] = C80radii
    results['npix80'] = npix80
    results['C80radii_err'] = C80radii_err
    results['npix80_err'] = npix80_err

    results['L90'] = L90
    results['A90'] = A90
    results['C90'] = sigma_90 / std
    results['C90radii'] = C90radii
    results['npix90'] = npix90
    results['C90radii_err'] = C90radii_err
    results['npix90_err'] = npix90_err
    results['flag20'] = flag20
    results['flag50'] = flag50
    results['flag9095'] = flag9095
    results['flagL9095'] = flagL9095


    results['L95'] = L95
    results['A95'] = A95
    results['C95'] = sigma_95 / std
    results['C95radii'] = C95radii
    results['npix95'] = npix95
    results['C95radii_err'] = C95radii_err
    results['npix95_err'] = npix95_err

    results['L99'] = L99
    results['A99'] = A99
    results['C99radii'] = C99radii
    results['npix99'] = npix99
    results['C99radii_err'] = C99radii_err
    results['npix99_err'] = npix99_err
    
    results['A20_err'] = A20_err
    results['A50_err'] = A50_err
    results['A80_err'] = A80_err
    results['A90_err'] = A90_err
    results['A95_err'] = A95_err
    results['A99_err'] = A99_err


    results['gaussianity'] = gaussianity
    results['gaussianity_L50'] = gaussianity_L50

    results['A50_100'] = A50_100
    results['C50_100radii'] = C50_100radii
    results['npix50_100'] = npix50_100

    results['A50_full'] = A50_full
    results['C50_full_radii'] = C50_full_radii
    results['npix50_full'] = npix50_full

    results['radii_ratio'] = radii_ratio
    results['radii_ratio_full'] = radii_ratio_full
    results['area_ratio'] = area_ratio
    results['area_ratio_full'] = area_ratio_full
    results['beam_area'] = beam_area_

    if logger is not None:
        print_logger_header(title="Basic Source Properties",
                            logger=logger)
        logger.debug(f" ==>  Peak of Flux="
                     f"{results['peak_of_flux']*1000:.2f} [mJy/beam]")
        logger.debug(f" ==>  Total Flux Inside Mask='"
                     f"{results['total_flux_mask']*1000:.2f} [mJy]")
        logger.debug(f" ==>  Total Flux Image="
                     f"{results['total_flux_nomask'] * 1000:.2f} [mJy]")
        logger.debug(f" ==>  Half-Light Radii="
                     f"{results['C50radii']:.2f} [px]")
        logger.debug(f" ==>  Total Source Size="
                     f"{results['C95radii']:.2f} [px]")
        if do_fit_ellipse:
            logger.debug(f" ==>  Source Global Axis Ratio="
                        f"{results['qm']:.2f}")
            logger.debug(f" ==>  Source Global PA="
                        f"{results['PAm']:.2f} [degrees]")
            logger.debug(f" ==>  Inner Axis Ratio="
                        f"{results['qmi']:.2f}")
            logger.debug(f" ==>  Outer Axis Ratio="
                        f"{results['qmo']:.2f}")
            logger.debug(f" ==>  Inner PA="
                        f"{results['PAmi']:.2f} [degrees]")
            logger.debug(f" ==>  Outer PA="
                        f"{results['PAmo']:.2f} [degrees]")



    fig = plt.figure(figsize=(10, 4))
    ax1 = fig.add_subplot(1, 2, 1)
    ax1.scatter(radii*cell_size, Lgrow / results['total_flux_mask'],
                label='Norm Masked Flux')
    ax1.axhline(0, ls='-.', color='black')
    # ax1.axvline(g.max()*0.5,label=r'$0.5\times \max$',color='purple')
    # ax1.axvline(g.max()*0.1,label=r'$0.1\times \max$',color='#E69F00')
    ax1.axvline(C50radii*cell_size,
                label=r"$R_{50}\sim $"f"{C50radii*cell_size:0.3f}+/-"f"{C50radii_err*cell_size:0.3f}''",
                ls='-.', color='lime')
    ax1.axhline(L50_norm, ls='-.', color='lime')
    # ax1.plot(C50radii, 
    #            0.5, 'rx', markersize=10, 
    #         label=r"$R_{50}\sim $"f"{C50radii*cell_size:0.3f}''",
    #         ls='-.', color='lime')
    # ax1.axhline(L50_norm, ls='-.', color='lime')
    ax1.axvline(C95radii*cell_size,
                label=r"$R_{95}\sim $"f"{C95radii*cell_size:0.3f}+/-"f"{C95radii_err*cell_size:0.3f}''",
                # ls='--',
                color='#56B4E9')
    # ax1.axhline(0.95, ls='-.', color='#56B4E9')
    # ax1.axvline(std * 6, label=r"6.0$\times \sigma_{\mathrm{mad}}$", color='black')
    # if last_level<3:
    #     ax1.axvline(std * 3, label=r"3.0$\times \sigma_{\mathrm{mad}}$", color='brown')

    ax1.set_title("Total Integrated Flux Density \n "
                  r"($\sigma_{\mathrm{mad}}$ levels) = "
                  f"{1000*np.sum(fluxes):.2f} $\pm$ {1000*total_flux_density_error:.2f} mJy")
    #     ax1.axvline(mad_std(g)*1,label=r'$1.0\times$ std',color='gray')
    # ax1.axvline(levels[-1], label=r"Mask Dilation",
    #             color='cyan')

    ax1.set_xlabel(fr'Aperture Circular Radius $R$ [{scale_units}]')
    ax1.set_ylabel(fr"$S(\leq R)$")
    # ax1.semilogx()
    ax1.grid(alpha=0.5)
    # plt.xlim(1e-6,)
    ax1.legend(loc='lower right',prop={'size': 12})


    # fig = plt.figure(figsize=(10, 4))
    # ax1 = fig.add_subplot(1, 2, 1)
    # ax1.scatter(levels[:], np.cumsum(fluxes) / results['total_flux_mask'],
    #             label='Mask Norm Flux')

    # ax1.axhline(0, ls='-.', color='black')
    # # ax1.axvline(g.max()*0.5,label=r'$0.5\times \max$',color='purple')
    # # ax1.axvline(g.max()*0.1,label=r'$0.1\times \max$',color='#E69F00')
    # ax1.axvline(sigma_50,
    #             label=r"$R_{50}\sim $"f"{C50radii*cell_size:0.3f}''",
    #             ls='-.', color='lime')
    # ax1.axhline(L50_norm, ls='-.', color='lime')
    # ax1.axvline(sigma_95,
    #             label=r"$R_{95}\sim $"f"{C95radii*cell_size:0.3f}''",
    #             # ls='--',
    #             color='#56B4E9')
    # ax1.axvline(std * 6, label=r"6.0$\times \sigma_{\mathrm{mad}}$", color='black')
    # if last_level<3:
    #     ax1.axvline(std * 3, label=r"3.0$\times \sigma_{\mathrm{mad}}$", color='brown')

    # ax1.set_title("Total Integrated Flux Density \n "
    #               r"($\sigma_{\mathrm{mad}}$ levels) = "
    #               f"{1000*np.sum(fluxes):.2f} mJy")
    # #     ax1.axvline(mad_std(g)*1,label=r'$1.0\times$ std',color='gray')
    # ax1.axvline(levels[-1], label=r"Mask Dilation",
    #             color='cyan')

    # ax1.set_xlabel('Levels [Jy/Beam]')
    # ax1.set_ylabel("Fraction of Integrated Flux per level")
    # ax1.semilogx()
    # ax1.grid(alpha=0.5)
    # # plt.xlim(1e-6,)
    # ax1.legend(loc='lower left')

    ax2 = fig.add_subplot(1, 2, 2)

    vmin = vmin_factor * std
    #     print(g)
    vmax = vmax_factor * g.max()
    norm = simple_norm(g, stretch='asinh', asinh_a=0.05, min_cut=vmin,
                       max_cut=vmax)

    # if crop == True:
    #     try:
    #         xin, xen, yin, yen = do_cutout_2D(img, 
    #                                           box_size=box_size, 
    #                                           center=None,
    #                                           centre_mode = 'image_centre',
    #                                           return_='box')
    #         g = g[xin:xen, yin:yen]
    #     except:
    #         try:
    #             max_x, max_y = np.where(g == g.max())
    #             xin = max_x[0] - box_size
    #             xen = max_x[0] + box_size
    #             yin = max_y[0] - box_size
    #             yen = max_y[0] + box_size
    #             g = g[xin:xen, yin:yen]
    #         except:
    #             pass

    im_plot = ax2.imshow(g, cmap='magma_r', origin='lower', alpha=1.0,
                         norm=norm,
                         aspect=aspect)  # ,vmax=vmax, vmin=vmin)#norm=norm

    # levels_50 = np.asarray([-3*sigma_50,sigma_50,3*sigma_50])

    try:
        ax2.contour(g, levels=levels_50, colors='lime', linewidths=2.5,
                    alpha=1.0)  # cmap='Reds', linewidths=0.75)
        # ax2.contour(g, levels=levels_90, colors='white', linewidths=2.0,
        #             # linestyles='--',
        #             alpha=1.0)  # cmap='Reds', linewidths=0.75)
        ax2.contour(g, levels=levels_95, colors='#56B4E9', linewidths=2.0,
                    # linestyles='--',
                    alpha=1.0)  # cmap='Reds', linewidths=0.75)
        #         ax2.contour(g, levels=levels_3sigma,colors='#D55E00',
        #         linewidths=1.5,alpha=1.0)#cmap='Reds', linewidths=0.75)
        ax2.contour(g, levels=[last_level * std], colors='cyan', linewidths=0.6,
                    alpha=1.0,
                    # linestyles='--'
                    )  # cmap='Reds', linewidths=0.75)
        ax2.contour(g, levels=[6.0 * std], colors='black', linewidths=1.5,
                    alpha=0.9,
                    # linestyles='--'
                    )  # cmap='Reds', linewidths=0.75)
        ax2.contour(g, levels=[3.0 * std], colors='brown', linewidths=1.2,
                    alpha=0.9,
                    # linestyles='--'
                    )  # cmap='Reds', linewidths=0.75)
        # ax2.contour(g, levels=[5.0 * std], colors='brown', linewidths=0.6,
        #             alpha=0.3,
        #             # linestyles='--'
        #             )  # cmap='Reds', linewidths=0.75)

    except:
        print('Not plotting contours!')
    ax2.axis('off')
    plt.subplots_adjust(wspace=0, hspace=0)
    fig.tight_layout()
    if SAVE is not None:
        plt.savefig(img.replace('.fits', '_Lgrow_levels')+add_save_name + ext,
                    dpi=300,bbox_inches='tight')
    if show_figure == True:
        plt.show()
    else:
        plt.close()

    if save_csv == True:
        import csv
        # df_temp = pd.DataFrame(results).T
        # df_temp.to_csv(img.replace('.fits', '_image_properties')+add_save_name + '.csv',
        #                header=True,index=False)

        with open(img.replace('.fits', '_image_properties')+add_save_name + '.csv',
                  'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=results.keys())
            writer.writeheader()
            writer.writerow(results)
    return (levels, fluxes, agrow, omask, mask, results)


def structural_morphology(imagelist, residuallist,
                          indices, masks_deblended,
                          zd, big_mask=None, data_2D=None,sigma_mask=6.0,
                          iterations = 2,
                          do_PLOT=False, show_figure=False,
                          sigma_loop_init=6.0, do_measurements='all',
                          verbose=0):
    """
    From the emission  of a given source and its deblended components,
    run in each component the morphometry analysis.

    A list of images is accepted with a common deblended mask for all.

    This was originally intended to be used with multi-resolution images.
    First, a common-representative image is processed with a source detection and
    deblending algorithm. Then, those detected/deblended regions are used to
    run a forced morphometry on all multi-resolution images.

    << Finish documentation >>

    """
    results_conc = []
    missing_data_im = []
    missing_data_re = []
    masks = []
    for i in tqdm(range(len(imagelist))):
        try:
            crop_image = imagelist[i]
            data_2D = load_fits_data(crop_image)

            if residuallist is not None:
                crop_residual = residuallist[i]
                residual_2D = load_fits_data(crop_residual)
                std = mad_std(residual_2D)
                if verbose >= 1:
                    print('Using RMS from residual')
            else:
                crop_residual = None
                residual_2D = None
                std = mad_std(data_2D)
                if verbose >= 1:
                    print('Using RMS from data ')
            #         crop_residual = residuallist[i]
            cell_size = get_cell_size(crop_image)
            std = mad_std(data_2D)
            npixels = int(2*beam_area2(crop_image))
            #             residual_2D = load_fits_data(crop_residual)



            processing_results_source = {}  # store calculations only for source
            processing_results_source['#imagename'] = os.path.basename(
                crop_image)
            # first, run the analysis for the entire source structure
            processing_results_source, mask, _ = measures(imagename=crop_image,
                                                          residualname=crop_residual,
                                                          z=zd, deblend=False,
                                                          apply_mask=True,
                                                        #   mask = big_mask,
                                                          results_final=processing_results_source,
                                                          plot_catalog=False,
                                                          rms=std,
                                                          bkg_sub=False,
                                                          bkg_to_sub=None,
                                                          mask_component=None,
                                                          npixels=npixels, fwhm=121,
                                                          kernel_size=121,
                                                          sigma_mask=sigma_mask,
                                                          last_level=1.5,
                                                          iterations=iterations,
                                                          dilation_size=None,
                                                          do_measurements=do_measurements,
                                                          do_PLOT=do_PLOT,
                                                          show_figure=show_figure,
                                                          add_save_name='',
                                                          verbose=verbose)
            flag_subcomponent = 0
            processing_results_source['freq'] = getfreqs([crop_image])[0]
            processing_results_source['comp_ID'] = str(0)
            processing_results_source['flag_subcomponent'] = flag_subcomponent
            # ref_mask = mask*big_mask
            
            results_conc.append(processing_results_source)
            # bkg_ = sep_background(crop_image, apply_mask=True, mask=None, bw=11,
            #                       bh=11, fw=12, fh=12)

            omaj, omin, _, _, _ = beam_shape(crop_image)
            dilation_size = int(0.5*np.sqrt(omaj * omin) / (2 * get_cell_size(crop_image)))
            # print('dilation_size=', dilation_size)
            # masks_expanded = {}
            
            if len(indices) > 1:
                for j in range(len(indices)):
                    # ii = str(i+1)
                    sigma_loop = sigma_loop_init  # reset the loop
                    processing_results_components = {}  # store calculation only for individual components of the soruce
                    processing_results_components['#imagename'] = os.path.basename(crop_image)
                    

                    mask_component = masks_deblended[j]
                    data_component = mask_component*data_2D.copy()
                    add_save_name = 'comp_' + str(j+1)
                    # print('Component id ', processing_results_components['comp_ID'])
                    try:
                        # mask_new = mask_component.copy()
                        _, mask_new = \
                            mask_dilation_from_mask(data_2D,
                                                    mask_component,
                                                    sigma=sigma_loop,
                                                    PLOT=False,iterations=iterations,
                                                    dilation_size=dilation_size,
                                                    show_figure=False)

                        # dilated masks must not overlap >> non conservation of flux
                        # mask_new = mask_new * big_mask
                        for l in range(len(indices)):
                            if l != j:
                                mask_new[masks_deblended[l]] = False
                            else:
                                pass
                        # masks_expanded[f'mask_ex_{j}'] = mask_new
                        # plt.figure()
                        # plt.imshow(mask_new*big_mask)
                        # plt.show()
                        processing_results_components, mask, _ = \
                            measures(crop_image, crop_residual, z=zd,
                                        deblend=False, apply_mask=False,
                                        plot_catalog=False, bkg_sub=False,
                                        mask = big_mask,
                                        mask_component=mask_new, rms=std,
                                        iterations=iterations, npixels=50, fwhm=121,
                                        kernel_size=121, sigma_mask=sigma_loop,
                                        last_level=1.5,
                                        # bkg_to_sub = bkg_.back(),
                                        dilation_size=dilation_size,
                                        add_save_name=add_save_name,
                                        do_measurements=do_measurements,
                                        do_PLOT=do_PLOT, show_figure=show_figure,
                                        verbose=verbose,
                                        results_final=processing_results_components)
                        flag_subcomponent = 0
                        # print('Component id ', processing_results_components['comp_ID'])
                        processing_results_components['freq'] = getfreqs([crop_image])[0]
                        processing_results_components['comp_ID'] = str(j+1)
                        processing_results_components['flag_subcomponent'] = flag_subcomponent
                    except:
                        try:
                            error_mask = True
                            while error_mask and sigma_loop > 1.0:
                                try:
                                    # mask_new = mask_component.copy()
                                    _, mask_new = \
                                        mask_dilation_from_mask(load_fits_data(crop_image),
                                                                mask_component,
                                                                rms=std,
                                                                sigma=sigma_loop,
                                                                iterations=3,
                                                                PLOT=False,
                                                                dilation_size=dilation_size,
                                                                show_figure=False)

                                    # dilated masks must not overlap >> non conservation of flux
                                    for l in range(len(indices)):
                                        if l != j:
                                            mask_new[masks_deblended[l]] = False
                                        else:
                                            pass
                                    
                                    # masks_expanded[f'mask_ex_{j}'] = mask_new
                                    
                                    if sigma_loop >= 3.0:
                                        last_level = 1.5
                                    if sigma_loop < 3.0:
                                        last_level = sigma_loop - 0.5

                                    (processing_results_components, mask,
                                     _) = measures(crop_image, crop_residual,
                                                   z=zd, deblend=False,
                                                   apply_mask=False,
                                                   plot_catalog=False,
                                                   bkg_sub=False,
                                                   mask = big_mask,
                                                   mask_component=mask_new,
                                                   rms=std, iterations=3,
                                                   npixels=1000, fwhm=121,
                                                   kernel_size=121,
                                                   sigma_mask=sigma_loop,
                                                   last_level=last_level,
                                                   add_save_name=add_save_name,
                                                   dilation_size=dilation_size,
                                                   do_measurements=do_measurements,
                                                   do_PLOT=do_PLOT, 
                                                   show_figure=show_figure,
                                                   verbose=verbose,
                                                   results_final=processing_results_components)
                                    # print('Component id ', processing_results_components['comp_ID'])
                                    processing_results_components['subreg_sigma'] = sigma_loop
                                    error_mask = False
                                    flag_subcomponent = 1
                                    processing_results_components['freq'] = getfreqs([crop_image])[0]
                                    processing_results_components['comp_ID'] = str(j+1)
                                    processing_results_components['flag_subcomponent'] = flag_subcomponent
                                except Exception as e:
                                    # Handle the error, and decrease p by 0.5
                                    print(
                                        f"Error occurred with sigma={sigma_loop}: {e}")
                                    sigma_loop -= 0.5
                                    print("Reducing sigma to=", sigma_loop)

                            if not error_mask:
                                print("Function call successful with sigma_mad=",
                                      sigma_loop)
                            else:
                                print(
                                    "Unable to call function with any value of sigma_mad.")
                        except:
                            print(
                                'Last attempt to perform morphometry, '
                                'with mininum threshold allowed.')
                            processing_results_components, mask, _ = measures(
                                crop_image, crop_residual, z=zd, deblend=False,
                                apply_mask=False,mask = big_mask,
                                plot_catalog=False, bkg_sub=False,
                                # bkg_to_sub = bkg_.back(),
                                mask_component=mask_new, rms=std,
                                dilation_size=dilation_size, iterations=2,
                                npixels=int(beam_area2(crop_image)),
                                fwhm=81, kernel_size=81, sigma_mask=2.0,
                                last_level=0.5,
                                add_save_name=add_save_name,
                                do_measurements=do_measurements,
                                do_PLOT=do_PLOT, show_figure=show_figure,
                                verbose=verbose,
                                results_final=processing_results_components)
                            processing_results_components['subreg_sigma'] = 1.0
                            processing_results_components['comp_ID'] = str(j+1)
                            processing_results_components['freq'] = getfreqs([crop_image])[0]
                            flag_subcomponent = 1
                            processing_results_components['flag_subcomponent'] = flag_subcomponent

                    results_conc.append(processing_results_components)
                    # masks.append(masks_expanded)
                processing_results_source['ncomps'] = len(indices)
            else:
                processing_results_source['ncomps'] = 1
        except:
            print('Some imaging data is missing!')
            missing_data_im.append(os.path.basename(crop_image))
            missing_data_re.append(os.path.basename(crop_residual))
            pass
    return (pd.DataFrame(results_conc), processing_results_source, missing_data_im)



make_flux_vs_std = deprecated("make_flux_vs_std",
                              "compute_image_properties")(compute_image_properties)




"""
 __  __                  _                          _              
|  \/  | ___  _ __ _ __ | |__   ___  _ __ ___   ___| |_ _ __ _   _ 
| |\/| |/ _ \| '__| '_ \| '_ \ / _ \| '_ ` _ \ / _ \ __| '__| | | |
| |  | | (_) | |  | |_) | | | | (_) | | | | | |  __/ |_| |  | |_| |
|_|  |_|\___/|_|  | .__/|_| |_|\___/|_| |_| |_|\___|\__|_|   \__, |
                  |_|                                        |___/
                  
#Morphometry 
"""

def background_asymmetry(img, mask, pre_clean=False):
    """
    <<<Morfometryka-core part>>>
    """
    def measure_asymmetry_patch(pos, patch, img):
        (x0, y0) = pos
        rot_cell = rot180(patch, x0, y0)
        sub = patch - rot_cell
        return np.sum(abs(sub)) / np.sum(abs(mask * img))

    Mo,No = img.shape
    gridsize = Mo // 10  # 10% of the size of the image
    n_pix = gridsize ** 2
    xcells = Mo // gridsize
    ycells = No // gridsize
    asymmetry_grid = np.zeros((xcells, ycells))
    gal_area = mask.sum()

    for xi in range(xcells):
        for yi in range(ycells):
            cell_mask = mask[xi * gridsize:(xi + 1) * gridsize,
                        yi * gridsize:(yi + 1) * gridsize]

            if cell_mask.sum() > 0:
                asymmetry_grid[xi, yi] = 0
                continue

            cell_img = img[xi * gridsize:(xi + 1) * gridsize, yi * gridsize:(yi + 1) * gridsize]
            x0, y0 = fmin(measure_asymmetry_patch, (gridsize // 2, gridsize // 2), args=(cell_img, img), disp=0)
            asymmetry_grid[xi, yi] = (gal_area / n_pix) * measure_asymmetry_patch((x0, y0), cell_img, img)
            del cell_img

    linear = asymmetry_grid[np.where(asymmetry_grid != 0)].ravel()

    if len(linear) > 0:
        BGrandom = np.random.choice(linear, 1)[0]
        BGmedian = np.median(linear)
        BGmin = linear.min()
        BGstd = np.std(linear)
        position = np.where(asymmetry_grid == linear.min())
        x0 = position[1][0] * gridsize + gridsize // 2
        y0 = position[0][0] * gridsize + gridsize // 2

    elif pre_clean == False:
        # measure background asymmetry with original pre-clean image if it fails for the clean one
        return background_asymmetry(img, mask, pre_clean=True)
    else:
        '''
           This is a fallback for when something goes wrong with background asymmetry estimates.
           It should also appear as a QF.
        '''
        BGrandom = 0
        BGmedian = 0
        BGmin = 0
        BGstd = 0
        x0 = 0
        y0 = 0

    return BGrandom, BGmedian, BGmin, BGstd, x0, y0


def assimetria0(pos, img, mask, box=False):
    """
    <<<Morfometryka-core part>>>
    """
    # print(' @ - Computing Asymetry 0')
    (x0, y0) = pos
    # psfmask = psfmask(psfsigma, *img.shape, x0, y0)
    if (box):
        boxmask = np.zeros_like(img)
        try:
            radii_px = np.ceil(self.P.Rp * self.NRp)
        except:
            radii_px = np.ceil(self.P.Rp * 1.5)
        boxmask[int(x0 - radii_px):int(x0 + radii_px), int(y0 - radii_px):int(y0 + radii_px)] = 1
        imgorig = boxmask * img
        imgsub = boxmask * (img - rot180(img, x0, y0))
        A = np.sum(abs(imgsub)) / np.sum(abs(imgorig))
    else:
        imgorig = img * mask
        imgsub = (img - rot180(img, x0, y0)) * mask
        A = np.sum(abs(imgsub)) / np.sum(abs(imgorig))

    del imgorig, imgsub
    return A


def assimetria1(pos, img, mask,use_mask=True):
    """
    <<<Morfometryka-core part>>>
    """
    # print(' @ - Computing Asymetry 0')
    x0, y0 = pos
    A1img = np.abs(img - rot180(img, x0, y0)) / (np.sum(np.abs(img)))
    if use_mask==True:
        return np.sum(mask * A1img)
    else:
        AsySigma = 3.00
        A1mask = A1img > np.median(A1img) + AsySigma * mad_std(A1img)
        return np.sum(mask * A1mask * A1img)


def geo_mom(p, q, I, centered=True, normed=True, complex=False, verbose=False):
    """
    <<<Morfometryka-core part>>>
    return the central moment M_{p,q} of image I
    http://en.wikipedia.org/wiki/Image_moment
    F.Ferrari 2012, prior to 4th JPAS
    """


    M, N = I.shape
    x, y = np.meshgrid(np.arange(N), np.arange(M))

    M_00 = np.nansum(I)

    if centered:
        # centroids
        x_c = (1 / M_00) * np.nansum(x * I)
        y_c = (1 / M_00) * np.nansum(y * I)

        x = x - x_c
        y = y - y_c

        if verbose:
            print('centroid  at', x_c, y_c)

    if normed:
        NORM = M_00 ** (1 + (p + q) / 2.)
    else:
        NORM = 1.0

    if complex:
        XX = (x + y * 1j)
        YY = (x - y * 1j)
    else:
        XX = x
        YY = y

    M_pq = (1 / NORM) * np.nansum(XX ** p * YY ** q * I)

    return M_pq


def q_PA(image):
    """
    <<<Morfometryka-core part>>>
    Adapted version of momenta from main mfmtk.
    """
    m00 = geo_mom(0, 0, image, centered=0, normed=0)
    m10 = geo_mom(1, 0, image, centered=0, normed=0)
    m01 = geo_mom(0, 1, image, centered=0, normed=0)
    m11 = geo_mom(1, 1, image, centered=0, normed=0)
    m20 = geo_mom(2, 0, image, centered=0, normed=0)
    m02 = geo_mom(0, 2, image, centered=0, normed=0)

    mu20 = geo_mom(2, 0, image, centered=1, normed=0)
    mu02 = geo_mom(0, 2, image, centered=1, normed=0)
    mu11 = geo_mom(1, 1, image, centered=1, normed=0)

    # centroids
    x0col = m10 / m00
    y0col = m01 / m00

    # manor, minor and axis ratio
    lam1 = np.sqrt(abs((1 / 2.) * (mu20 + mu02 + np.sqrt((mu20 - mu02) ** 2 + 4 * mu11 ** 2))) / m00)
    lam2 = np.sqrt(abs((1 / 2.) * (mu20 + mu02 - np.sqrt((mu20 - mu02) ** 2 + 4 * mu11 ** 2))) / m00)
    a = max(lam1, lam2)
    b = min(lam1, lam2)

    PA = (1 / 2.) * np.arctan2(2 * mu11, (mu20 - mu02))
    if PA < 0:
        PA = PA + np.pi
    PAdeg = np.rad2deg(PA)
    a = a
    b = b
    q = b / a
    return PAdeg, b / a, x0col, y0col


def peak_center(image):
    """
    <<<Morfometryka-core part>>>
    """
    y0max, x0max = nd.maximum_position((image))
    try:
        """
        For very small images of emission, this function breaks. In that
        case, just return the the max position.
        """
        # size of peak region to consider in interpolarion for x_peak
        dp = 2
        CenterOffset = 10
        No, Mo = image.shape
        peakimage = image[y0max - dp:y0max + dp, x0max - dp:x0max + dp]
        m00 = geo_mom(0, 0, peakimage, centered=0, normed=0)
        m10 = geo_mom(1, 0, peakimage, centered=0, normed=0)
        m01 = geo_mom(0, 1, peakimage, centered=0, normed=0)

        x0peak = x0max + m10 / m00 - dp
        y0peak = y0max + m01 / m00 - dp

        # check if center is galaxy center, i.e., should be near the image center
        # otherwise apply a penalty to pixel value proportional to the center distance^2
        if np.sqrt((x0peak - No / 2.) ** 2 + (
                y0peak - Mo / 2.) ** 2) > CenterOffset:
            # define a penalty as we move from the center
            xx, yy = np.meshgrid(np.arange(No) - No / 2.,
                                 np.arange(Mo) - Mo / 2.)
            rr2 = xx ** 2 + yy ** 2
            y0peak, x0peak = nd.maximum_position((image / rr2))
        return (x0peak, y0peak)
    except:
        return (x0max, y0max)


def momenta(image, PArad_0=None, q_0=None):
    '''
    <<<Morfometryka-core part>>>
    Calculates center of mass, axis lengths and position angle
    '''

    m00 = geo_mom(0, 0, image, centered=0, normed=0)
    m10 = geo_mom(1, 0, image, centered=0, normed=0)
    m01 = geo_mom(0, 1, image, centered=0, normed=0)
    m11 = geo_mom(1, 1, image, centered=0, normed=0)
    m20 = geo_mom(2, 0, image, centered=0, normed=0)
    m02 = geo_mom(0, 2, image, centered=0, normed=0)

    mu20 = geo_mom(2, 0, image, centered=1, normed=0)
    mu02 = geo_mom(0, 2, image, centered=1, normed=0)
    mu11 = geo_mom(1, 1, image, centered=1, normed=0)

    # centroids
    x0col = m10 / m00
    y0col = m01 / m00

    # manor, minor and axis ratio
    lam1 = np.sqrt(abs((1 / 2.) *
                       (mu20 + mu02 + np.sqrt((mu20 - mu02) ** 2 +
                                              4 * mu11 ** 2))) / m00)
    lam2 = np.sqrt(abs((1 / 2.) *
                       (mu20 + mu02 - np.sqrt((mu20 - mu02) ** 2 +
                                              4 * mu11 ** 2))) / m00)
    a = max(lam1, lam2)
    b = min(lam1, lam2)

    PA = (1 / 2.) * np.arctan2(2 * mu11, (mu20 - mu02))
    if PA < 0:
        PA = PA + 2*np.pi

    # self.PArad = PA
    # self.PAdeg = np.rad2deg(self.PArad)
    # self.q = self.b/self.a
    # self.PArad = PA
    # self.PAdeg = np.rad2deg(self.PArad)
    # self.q = self.b/self.a

    # mofified by lucatelli.
    """This will force mfmtk to do photometry for the given input PA
    this can be useful when we want to study how the light profile changes
    as function of PA or q. This was indented to explore the difference between the
    profiles of elliptical and spiral galaxies, which may not change soo much for the
    former while it may for the later.
    """

    if PArad_0 is None:
        PArad = PA  # + np.pi/2
    else:
        PArad = PArad_0

    PAdeg = np.rad2deg(PArad)

    if q_0 is None:
        q = b / a
    else:
        q = q_0
    return (x0col, y0col, q, PAdeg)


def cal_PA_q(gal_image_0,Isequence = None,region_split=None,SAVENAME=None):
    '''
    <<<Morfometryka-core part>>>
    Estimates inner and outer PA nad q=(b/a)
    '''
    # mean Inner q,  mean outer q,  mean Inner PA,  mean Outer PA
    # from fitEllipse import main_test2
    
    # qmi, qmo, PAmi, PAmo, qm, PAm,\
    #     x0median,y0median,x0median_i,y0median_i,\
    #     x0median_o,y0median_o = main_test2(gal_image_0,
    #                                        Isequence = Isequence,
    #                                        region_split=region_split,
    #                                        SAVENAME=SAVENAME)
    
    from fit_ellipse import fit_ellipse_to_galaxy
    qmi, qmo, PAmi, PAmo, qm, PAm,x0median,y0median,x0median_i,y0median_i,\
        x0median_o,y0median_o,profiles = \
            fit_ellipse_to_galaxy(gal_image_0,
                                  Isequence = Isequence,
                                  region_split=region_split,
                                #   plot_results=True,
                                #   plot_profiles=True,
                                  save_name=SAVENAME)

    # global PA,  global q
    PA, q, x0col, y0col = q_PA(gal_image_0)

    # print("Initial PA and q = ", PA, q)
    # print("Median PA and q = ", PAm, qm)
    # print("Inner-Mean PA and q = ", PAmi, qmi)
    # print("Outer-Mean PA and q = ", PAmo, qmo)
    return (PA, q, x0col, y0col, PAm, qm, PAmi, qmi, PAmo, qmo,
            x0median,y0median,x0median_i,y0median_i,x0median_o,y0median_o,profiles)



def savitzky_golay_2d(z, window_size, order, derivative=None):
    """
    <<<Morfometryka-core part>>>
    http://nbviewer.ipython.org/github/pv/SciPy-CookBook/blob/master/ipython/SavitzkyGolay.ipynb
    """

    # number of terms in the polynomial expression
    n_terms = (order + 1) * (order + 2) / 2.0

    if window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size ** 2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial.
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ...
    # this line gives a list of two item tuple. Each tuple contains
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [(k - n, n) for k in range(order + 1) for n in range(k + 1)]

    # coordinates of points
    ind = np.arange(-half_size, half_size + 1, dtype=np.float64)
    dx = np.repeat(ind, window_size)
    dy = np.tile(ind, [window_size, 1]).reshape(window_size ** 2, )

    # build matrix of system of equation
    A = np.empty((window_size ** 2, len(exps)))
    for i, exp in enumerate(exps):
        A[:, i] = (dx ** exp[0]) * (dy ** exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2 * half_size, z.shape[1] + 2 * half_size
    Z = np.zeros((new_shape))
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] = band - np.abs(np.flipud(z[1:half_size + 1, :]) - band)
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band + np.abs(np.flipud(z[-half_size - 1:-1, :]) - band)
    # left band
    band = np.tile(z[:, 0].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs(np.fliplr(z[:, 1:half_size + 1]) - band)
    # right band
    band = np.tile(z[:, -1].reshape(-1, 1), [1, half_size])
    Z[half_size:-half_size, -half_size:] = band + np.abs(np.fliplr(z[:, -half_size - 1:-1]) - band)
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0, 0]
    Z[:half_size, :half_size] = band - np.abs(np.flipud(np.fliplr(z[1:half_size + 1, 1:half_size + 1])) - band)
    # bottom right corner
    band = z[-1, -1]
    Z[-half_size:, -half_size:] = band + np.abs(np.flipud(np.fliplr(z[-half_size - 1:-1, -half_size - 1:-1])) - band)

    # top right corner
    band = Z[half_size, -half_size:]
    Z[:half_size, -half_size:] = band - np.abs(np.flipud(Z[half_size + 1:2 * half_size + 1, -half_size:]) - band)
    # bottom left corner
    band = Z[-half_size:, half_size].reshape(-1, 1)
    Z[-half_size:, :half_size] = band - np.abs(np.fliplr(Z[-half_size:, half_size + 1:2 * half_size + 1]) - band)

    # solve system and convolve
    if derivative is None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')


def standartize(image, q, PArad, x0, y0):
    """
    <<<Morfometryka-core part>>>
    make a standard galaxy, id est, PA=0, q=1
    arguments are 'image' to be standartized and  its 'S' stamp and P phot classes
    """

    ##### rotate array
    R = np.array([[np.cos(PArad), np.sin(PArad)], [-np.sin(PArad), np.cos(PArad)]])

    ##### shear array
    S = np.diag([q, 1.])
    # SERSIC fit values
    # S = np.diag([self.Ss.qFit2D, 1.])

    # affine transform matrix, rotate then scale
    transform = np.dot(R, S)

    # where to transform about
    centro_i = (x0, y0)
    # contro_o: where to put center after
    centro_o = np.array(image.shape) / 2

    myoffset = centro_i - np.dot(transform, centro_o)
    bval = np.mean(image[-2:])
    stangal = nd.affine_transform(image, transform, offset=myoffset, order=2, cval=bval)

    return stangal


def polarim(image, origin=None, log=False):
    """
    <<<Morfometryka-core part>>>
    Reprojects a 2D numpy array ("image") into a polar coordinate system.
    "origin" is a tuple of (x0, y0) and defaults to the center of the image.
    http://stackoverflow.com/questions/3798333/image-information-along-a-polar-coordinate-system
    refactored by FF, 2013-2014 (see transpolar.py)
    """

    if origin is None:
        origin = np.array(image.shape) / 2.

    def cart2polar(x, y):
        r = np.sqrt(x ** 2 + y ** 2)
        #         max_radii = np.sqrt(x.max()**2+y.max()**2)
        #         rscale = x.max()/max_radii
        #         tscale = y.max()/(2*np.pi)
        theta = np.arctan2(y, x)
        return r, theta

    def polar2cart(r, theta):
        x = r * np.cos(theta)
        y = r * np.sin(theta)
        return x, y

    def cart2logpolar(x, y, M=1):
        alpha = 0.01
        r = np.sqrt(x ** 2 + y ** 2)
        rho = M * np.log(r + alpha)
        theta = np.arctan2(y, x)
        return rho, theta

    def logpolar2cart(rho, theta, M=1):
        x = np.exp(rho / M) * np.cos(theta)
        y = np.exp(rho / M) * np.sin(theta)
        return x, y

    ny, nx = image.shape
    if origin is None:
        x0, y0 = (nx // 2, ny // 2)
        origin = (x0, y0)
    else:
        x0, y0 = origin

    # Determine that the min and max r and theta coords will be...
    x, y = np.meshgrid(np.arange(nx) - x0, np.arange(ny) - y0)  # ,sparse=True )

    r, theta = cart2polar(x, y)

    # Make a regular (in polar space) grid based on the min and max r & theta
    r_i = np.linspace(r.min(), r.max(), nx)
    theta_i = np.linspace(theta.min(), theta.max(), ny)
    theta_grid, r_grid = np.meshgrid(theta_i, r_i)

    # Project the r and theta grid back into pixel coordinates
    xi, yi = polar2cart(r_grid, theta_grid)
    xi += origin[0]  # We need to shift the origin back to
    yi += origin[1]  # back to the lower-left corner...
    xi, yi = xi.flatten(), yi.flatten()
    coords = np.vstack((xi, yi))  # (map_coordinates requires a 2xn array)

    zi = nd.map_coordinates(image, coords, order=1)  # ,prefilter=False)
    galpolar = zi.reshape((nx, ny))

    r_polar = r_i
    theta_polar = theta_i

    return galpolar, r_polar, theta_polar, x, y


def Gradindex(image_data, Rp=None):
    """
    <<<Morfometryka-core part>>>
    Gradient Index
    Calculates an index based on the image gradient magnitude and orientation

    SGwindow and SGorder are Savitsky-Golay filter parameters
    F. Ferrari, 2014
    """

    def sigma_func(params):
        '''
        <<<Morfometryka-core part>>>
        calculates the sigma psi with different parameters
        called by the minimization routine '''
        (x0, y0, q, PA) = params
        #### creates standardized image
        # using segmentation geometric parameters
        # galnostarsstd = standartize(P.galnostars, S.q, S.PArad, S.y0col, S.x0col)
        # using Sersic geometric parameters
        galnostarsstd = standartize(image_data, q, PA, x0, y0)

        # creates polar imagem
        galpolar, r_polar, theta_polar, _, _ = polarim(galnostarsstd)
        #     print(galpolar)
        # print '%.5f %.5f %.5f %.5f' % (x0, y0, q, PA)

        if Rp is None:
            galpolarpetro = galpolar[0: -1, :]
        else:
            galpolarpetro = galpolar[0: 2 * int(Rp), :]
        #     galpolarpetro =  galpolar[ 0   : int(config.NRp * P.Rp), : ]

        # circular_area_radius
        #     if min(galpolarpetro.shape) <= SGwindow:
        #         SGorder -= 1

        try:
            # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
            dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, SGorder, 'both')
        except:
            SGwindow = 5
        # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
        dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, 1, 'both')

        mag = np.sqrt(dx ** 2 + dy ** 2)
        # mag = dxdy
        magmask = mag > (np.median(mag))
        ort = np.arctan2(dy, dx)
        # ort = np.arctan(dxdy)
        ortn = (ort + np.pi) % (np.pi)

        psi = circmean(ortn[magmask])
        sigma_psi = circstd(ortn[magmask])

        return sigma_psi

    def sigma_func_eval(params):
        '''
        <<<Morfometryka-core part>>>
        calculates the sigma psi with different parameters
        called by the minimization routine '''
        (x0, y0, q, PA) = params
        #### creates standardized image
        # using segmentation geometric parameters
        # galnostarsstd = standartize(P.galnostars, S.q, S.PArad, S.y0col, S.x0col)
        # using Sersic geometric parameters
        galnostarsstd = standartize(image_data, q, PA, x0, y0)

        # creates polar imagem
        galpolar, r_polar, theta_polar, _, _ = polarim(galnostarsstd)
        #     print(galpolar)
        # print '%.5f %.5f %.5f %.5f' % (x0, y0, q, PA)

        #     galpolarpetro =  galpolar[ 0   : int(config.NRp * P.Rp), : ]
        if Rp is None:
            galpolarpetro = galpolar[0: -1, :]
        else:
            galpolarpetro = galpolar[0: 2 * int(Rp), :]
        # circular_area_radius
        #     if min(galpolarpetro.shape) <= SGwindow:
        #         SGorder -= 1

        try:
            # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
            dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, SGorder, 'both')
        except:
            SGwindow = 5
        # dx,dy = np.gradient(savgol_filter(galpolarpetro, SGwindow, SGorder, 1))
        dx, dy = savitzky_golay_2d(galpolarpetro, SGwindow, 1, 'both')

        mag = np.sqrt(dx ** 2 + dy ** 2)
        # mag = dxdy
        magmask = mag > (np.median(mag))
        ort = np.arctan2(dy, dx)
        # ort = np.arctan(dxdy)
        ortn = (ort + np.pi) % (np.pi)

        psi = circmean(ortn[magmask])
        sigma_psi = circstd(ortn[magmask])

        return sigma_psi, ort, ortn, dx, dy, mag

    # SAVITSKY-GOLAY parameters
    # polynom order
    SGorder = 5

    # SG window is galaxy_size/10 and must be odd
    # CAN'T BE RELATIVE TO IMAGE SIZE... MUST BE RELATIVE TO GALAXY SIZE
    # SGwindow = int(S.Mo/10.)
    # SGwindow = int(P.Rp/2.)
    SGwindow = 5
    #     SGwindow = int(circular_area_radius/2)
    if SGwindow % 2 == 0:
        SGwindow = SGwindow + 1

    #     image_data=load_fits_data(imagelist[1])

    PA, q, x0col, y0col, PAm, qm, PAmi, qmi, PAmo, qmo,profiles = cal_PA_q(image_data)

    x0sigma, y0sigma, qsigma, PAsigma = \
        fmin(sigma_func, (x0col, y0col, qm, np.deg2rad(PAm)), ftol=0.1, xtol=1.0, disp=0)
    sigma_psi, ort, ortn, dx, dy, mag = sigma_func_eval((x0sigma, y0sigma, qsigma, PAsigma))

    #     x0sigma, y0sigma, qsigma, PAsigma = \
    #                 fmin(sigma_func, (Ss.x0Fit2D, Ss.y0Fit2D, Ss.qFit2D, np.deg2rad(Ss.PAFit2D)), ftol=0.1, xtol=1.0, disp=0)

    #     #if qsigma > 1:
    #     #   qsigma = 1/qsigma
    #     #   PAsigma = PAsigma - np.pi/2.

    #     x0sigma, y0sigma, qsigma, PAsigma = \
    #                 fmin(sigma_funcseg, (Ss.x0Fit2D, Ss.y0Fit2D, Ss.qFit2D, np.deg2rad(Ss.PAFit2D)), ftol=0.1, xtol=1.0, disp=0)
    return (sigma_psi, ort, ortn, dx, dy, mag, PAm, qm)



def evaluate_compactness(deconv_props, conv_props):
    """
    Determine if a model component of a radio structure is compact or extended.

    It uses the concentration index determined from the areas A20, A50, A80 and A90.

    """
    import pandas as pd
    ncomps = deconv_props['comp_ID'].shape[0]

    Spk_ratio = np.asarray(conv_props['peak_of_flux']) / np.asarray(deconv_props['peak_of_flux'])
    I50_ratio = np.asarray(deconv_props['I50']) / np.asarray(conv_props['I50'])
    class_criteria = {}
    for i in range(ncomps):
        class_criteria[f"comp_ID_{i + 1}"] = {}
        if Spk_ratio[i] < 0.5:
            class_criteria[f"comp_ID_{i + 1}"]['Spk_class'] = 'C'
        if Spk_ratio[i] >= 0.5:
            class_criteria[f"comp_ID_{i + 1}"]['Spk_class'] = 'D'

        if (I50_ratio[i] < 0.5) or (np.isnan(I50_ratio[i])): #nan because the component is too
            # small.
            class_criteria[f"comp_ID_{i + 1}"]['I50_class'] = 'C'
        if I50_ratio[i] >= 0.5:
            class_criteria[f"comp_ID_{i + 1}"]['I50_class'] = 'D'

        AC1_conv_check = conv_props['AC1'].iloc[i]
        AC2_conv_check = conv_props['AC2'].iloc[i]
        AC1_deconv_check = deconv_props['AC1'].iloc[i]
        AC2_deconv_check = deconv_props['AC2'].iloc[i]
        # print(AC1_conv_check)

        # if (AC1_conv_check >= 1.0) or (AC1_conv_check == np.inf):
        #     class_criteria[f"comp_ID_{i+1}"]['AC1_conv_class'] = 'C'
        # if (AC2_conv_check >= 1.0) or AC2_conv_check == np.inf:
        #     class_criteria[f"comp_ID_{i+1}"]['AC2_conv_class'] = 'C'
        if (AC1_deconv_check >= 1.0) or (np.isnan(AC1_deconv_check)):
            class_criteria[f"comp_ID_{i + 1}"]['AC1_deconv_class'] = 'C'
        if AC1_deconv_check < 1.0:
            class_criteria[f"comp_ID_{i + 1}"]['AC1_deconv_class'] = 'D'
        if AC2_deconv_check >= 0.75:
            class_criteria[f"comp_ID_{i + 1}"]['AC2_deconv_class'] = 'C'
        if AC2_deconv_check < 0.75:
            class_criteria[f"comp_ID_{i + 1}"]['AC2_deconv_class'] = 'D'

    class_results_df = pd.DataFrame(class_criteria)
    for i in range(ncomps):
        dessision_compact = np.sum(class_results_df[f"comp_ID_{i + 1}"] == 'C')
        # dessision_diffuse = np.sum(class_results_df[f"comp_ID_{i+1}"]=='D')
        if dessision_compact > 2:
            class_criteria[f"comp_ID_{i + 1}"]['final_class'] = 'C'
        if dessision_compact < 2:
            class_criteria[f"comp_ID_{i + 1}"]['final_class'] = 'D'
        if dessision_compact == 2:
            try:
                class_criteria[f"comp_ID_{i + 1}"]['final_class'] = \
                    class_criteria[(f"comp_ID_{i + 1}")]['Spk_class']
            except:
                class_criteria[f"comp_ID_{i + 1}"]['final_class'] = 'D'
    return (class_criteria)





def convex_morpho_old(image, mask, scale=1.0,do_plot=False):
    """
    Overlay semi-major and minor axes on an image and report properties.
    
    Parameters:
        image (2D array): The image data to be plotted.
        mask (2D array): mask data of the structure.
        scale (float): Scaling factor for the length of the axes vectors.
        do_plot (boolean): Plot or not the results.
    
    Returns:
        dict: Report containing basic morphometry.
    """
    # Extract y, x coordinates of the structure
    indices = np.transpose(np.nonzero(mask))
    y, x = indices[:, 0], indices[:, 1]
    points = np.column_stack((x, y))  # Convert to (x, y) coordinates
    
    # Compute ConvexHull
    hull = ConvexHull(points)

    # Compute centroid of the points
    centroid = np.mean(points, axis=0)

    # Covariance matrix and eigenvalues/vectors
    cov_matrix = np.cov(points, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

    # Order eigenvectors by eigenvalues (largest is major axis)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    # Semi-major and minor axes
    major_axis_vector = eigenvectors[:, 0]
    minor_axis_vector = eigenvectors[:, 1]

    # Scale eigenvectors by eigenvalues for visualization
    major_axis = scale * np.sqrt(eigenvalues[0]) * major_axis_vector
    minor_axis = scale * np.sqrt(eigenvalues[1]) * minor_axis_vector

    # Calculate position angle (anti-clockwise from x-axis)
    position_angle = np.arctan2(major_axis_vector[1], major_axis_vector[0])
    position_angle_degrees = np.degrees(position_angle)
    if position_angle_degrees < 0:
        position_angle_degrees += 360

    # Calculate axis ratio (minor/major)
    axis_ratio = np.sqrt(eigenvalues[1] / eigenvalues[0])

    if do_plot:
        # Plot image
        plt.figure(figsize=(5, 5))
        plt.imshow(image, cmap='gray', origin='lower')
        # plt.colorbar(label='Intensity')
        
        # Plot semi-major and minor axes
        plt.quiver(
            centroid[0], centroid[1], major_axis[0], major_axis[1],
            angles='xy', scale_units='xy', scale=1, color='r', 
            # label='Major Axis'
        )
        plt.quiver(
            centroid[0], centroid[1], minor_axis[0], minor_axis[1],
            angles='xy', scale_units='xy', scale=1, color='b', 
            # label='Minor Axis'
        )
    
        # Plot centroid
        plt.scatter(centroid[0], centroid[1], color='orange', 
                    # label='Centroid', 
                    zorder=5)
    
        # Add labels
        # plt.legend()
        # plt.title("Semi-Major and Minor Axes on Image")
        plt.xlabel("$x$ image coordinates")
        plt.ylabel("$y$ image coordinates")
        plt.axis('equal')
        plt.show()

    # Return report
    report = {
        "PA_convex": position_angle_degrees,
        "q_convex": axis_ratio,
        "centroid_convex": centroid
    }
    return report


def compute_diameters_convex(hull,points):
    """
    Compute the major and minor diameters of a structure using its ConvexHull.

    Parameters
    ----------
    hull : scipy.spatial.ConvexHull
        A ConvexHull object representing the structure, containing points
        that define the convex boundary of the structure.

    Returns
    -------
    dict
        A dictionary containing:
        - "major_diameter" : float
            The length of the major diameter (longest distance between any two 
            points on the convex hull).
        - "major_points" : tuple of ndarray
            The two points defining the major diameter.
        - "minor_diameter" : float
            The length of the minor diameter (shortest perpendicular distance 
            between parallel edges of the convex hull).
        - "minor_points" : tuple of ndarray
            The two points defining the minor diameter.

    Examples
    --------
    >>> import numpy as np
    >>> from scipy.spatial import ConvexHull
    >>> # Example points
    >>> points = np.array([[0, 0], [1, 1], [2, 2], [3, 0], [0, 3]])
    >>> # Compute the convex hull
    >>> hull = ConvexHull(points)
    >>> # Compute diameters
    >>> diameters = compute_diameters_convex(hull)
    >>> print(diameters)

    Notes
    -----
    - The major diameter is calculated as the maximum Euclidean distance
      between any two points on the convex hull.
    - The minor diameter is determined by finding the shortest perpendicular
      width across the convex hull.
    """
    from itertools import combinations

    # hull = ConvexHull(points)
    hull_points = points[hull.vertices]  # Extract points on the convex hull

    # Calculate Major Diameter (longest distance between hull points)
    major_diameter = 0
    major_points = None
    for p1, p2 in combinations(hull_points, 2):
        distance = np.linalg.norm(p1 - p2)
        if distance > major_diameter:
            major_diameter = distance
            major_points = (p1, p2)

    # Calculate Minor Diameter (shortest perpendicular distance between parallel edges)
    minor_diameter = float('inf')
    minor_points = None
    num_hull_points = len(hull_points)

    for i in range(num_hull_points):
        # Get two consecutive points forming an edge
        p1, p2 = hull_points[i], hull_points[(i + 1) % num_hull_points]
        edge_vector = p2 - p1
        edge_length = np.linalg.norm(edge_vector)
        
        # Normalize the edge vector
        if edge_length == 0:
            continue
        edge_normal = np.array([-edge_vector[1], edge_vector[0]]) / edge_length
        
        # Project all hull points onto the edge normal and find the width
        distances = np.abs(np.dot(hull_points - p1, edge_normal))
        max_distance = distances.max()
        if max_distance < minor_diameter:
            minor_diameter = max_distance
            # Points corresponding to the shortest perpendicular projection
            projections = hull_points[np.abs(distances - max_distance) < 1e-6]
            if len(projections) >= 2:
                minor_points = (projections[0], projections[1])

    return {
        "major_diameter": major_diameter,
        "major_points": major_points,
        "minor_diameter": minor_diameter,
        "minor_points": minor_points
    }


def convex_morpho(image, mask, scale=1.0, do_plot=False):
    """
    Perform morphological analysis on a structure defined by a mask and overlay results on an image.

    This function computes morphological properties of a structure, including the position angle,
    axis ratio, centroid, major diameter, and minor diameter, based on the ConvexHull of the masked region.
    Optionally, it can overlay the semi-major and semi-minor axes, centroid, and convex hull on the image.

    Parameters
    ----------
    image : 2D ndarray
        The image data where the structure is located.
    mask : 2D ndarray
        A binary mask defining the structure to be analyzed. Non-zero pixels are treated as part of the structure.
    scale : float, optional
        Scaling factor for the visualization of the axes (default is 1.0).
    do_plot : bool, optional
        Whether to plot the results overlaying the axes and convex hull on the image (default is False).

    Returns
    -------
    dict
        A dictionary containing the following morphometric properties:
        - "PA_convex" : float
            Position angle of the structure in degrees, measured counter-clockwise from the positive x-axis.
        - "q_convex" : float
            Axis ratio (minor/major) of the structure.
        - "centroid_convex" : ndarray
            Centroid coordinates of the structure as a 2-element array [x, y].
        - "major_diameter" : float
            Length of the major diameter (longest distance between any two points on the convex hull).
        - "minor_diameter" : float
            Length of the minor diameter (shortest perpendicular distance between parallel edges of the convex hull).

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> from scipy.spatial import ConvexHull

    >>> # Create a sample image and mask
    >>> image = np.random.rand(100, 100)
    >>> mask = np.zeros_like(image, dtype=bool)
    >>> mask[40:60, 45:65] = True  # Example structure

    >>> # Perform analysis
    >>> report = convex_morpho(image, mask, scale=2.0, do_plot=True)
    >>> print(report)

    Notes
    -----
    - The function uses the ConvexHull of the masked region to estimate morphological properties.
    - The position angle is computed from the eigenvector of the largest eigenvalue of the covariance matrix.
    - Axis ratio is calculated as the square root of the ratio of the smallest to largest eigenvalues.
    - Major and minor diameters are determined using the ConvexHull geometry.
    """
    # Extract y, x coordinates of the structure
    indices = np.transpose(np.nonzero(mask))
    y, x = indices[:, 0], indices[:, 1]
    points = np.column_stack((x, y))  # Convert to (x, y) coordinates
    
    # Compute ConvexHull
    hull = ConvexHull(points)

    # Compute centroid of the points
    centroid = np.mean(points, axis=0)

    # Covariance matrix and eigenvalues/vectors
    cov_matrix = np.cov(points, rowvar=False)
    eigenvalues, eigenvectors = np.linalg.eig(cov_matrix)

    # Order eigenvectors by eigenvalues (largest is major axis)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]

    # Semi-major and minor axes
    major_axis_vector = eigenvectors[:, 0]
    minor_axis_vector = eigenvectors[:, 1]

    # Scale eigenvectors by eigenvalues for visualization
    major_axis = scale * np.sqrt(eigenvalues[0]) * major_axis_vector
    minor_axis = scale * np.sqrt(eigenvalues[1]) * minor_axis_vector

    # Calculate position angle (anti-clockwise from x-axis)
    position_angle = np.arctan2(major_axis_vector[1], major_axis_vector[0])
    position_angle_degrees = np.degrees(position_angle)
    if position_angle_degrees < 0:
        position_angle_degrees += 360

    # Calculate axis ratio (minor/major)
    axis_ratio = np.sqrt(eigenvalues[1] / eigenvalues[0])

    # Compute diameters
    diameters = compute_diameters_convex(hull,points)

    if do_plot:
        # Plot image
        plt.figure(figsize=(5, 5))
        plt.imshow(image*mask, cmap='gray', origin='lower')
        
        # Plot semi-major and minor axes
        plt.quiver(
            centroid[0], centroid[1], major_axis[0], major_axis[1],
            angles='xy', scale_units='xy', scale=1, color='r'
        )
        plt.quiver(
            centroid[0], centroid[1], minor_axis[0], minor_axis[1],
            angles='xy', scale_units='xy', scale=1, color='b'
        )
    
        # Plot centroid
        plt.scatter(centroid[0], centroid[1], color='orange', zorder=5)
    
        
        plt.xlabel("$x$ image coordinates")
        plt.ylabel("$y$ image coordinates")
        plt.axis('equal')
        # plt.legend()
        plt.show()

    # Return report
    report = {
        "PA_convex": position_angle_degrees,
        "q_convex": axis_ratio,
        "centroid_convex": centroid,
        "major_diameter": diameters["major_diameter"],
        "minor_diameter": diameters["minor_diameter"]
    }
    return report
