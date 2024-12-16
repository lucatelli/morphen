
"""
 ____  _        _
/ ___|| |_ __ _| |_ ___
\___ \| __/ _` | __/ __|
 ___) | || (_| | |_\__ \
|____/ \__\__,_|\__|___/

"""

def quadrature_error(imagedata,residualdata):
    sum_squares = np.nansum(residualdata**2.0)
    mean_value = np.nanmean(imagedata)
    N_size = imagedata.size
    q_error = np.sqrt((sum_squares / (N_size - 1)) - (mean_value**2))
    return(q_error)

def rms_estimate(imagedata):
    mean_value = np.mean(imagedata)
    square_difference = (imagedata - mean_value) ** 2
    mean_squared_diff = np.mean(square_difference)
    RMS = np.sqrt(mean_squared_diff)
    return(RMS)

def mc_flux_error(imagename,image,residual, num_threads = 6, n_samples = 1000):
    """
    TEST: DO NOT USE THIS
    """
    def calculate_sum(image, residual, seed):
        # Set random seed
        np.random.seed(seed)

        # Randomly sample a subset of the data with replacement
        sample = np.random.choice(image.ravel(), size=image.size, replace=True)

        # Add residual map to sampled data
        sample_with_residual = sample + residual.ravel()

        # Calculate sum of sampled data
        sample_sum = np.sum(sample_with_residual)

        return sample_sum

    import concurrent.futures

    print(' >> Running MC flux error estimate....')
    # Create ThreadPoolExecutor
    executor = concurrent.futures.ThreadPoolExecutor(max_workers=num_threads)

    # Submit tasks to ThreadPoolExecutor
    tasks = []
    for i in range(n_samples):
        task = executor.submit(calculate_sum, image, residual, i)
        tasks.append(task)

    # Retrieve results from tasks
    sums = []
    for task in tqdm(concurrent.futures.as_completed(tasks)):
        result = task.result()
        sums.append(result)

    # Calculate mean and standard deviation of sums
    sum_mean = np.mean(sums/beam_area2(imagename))
    sum_std = np.std(sums/beam_area2(imagename))

    # Print results
    print(f"Sum: {sum_mean:.4f}")
    print(f"Error: {sum_std:.4f}")

    # Plot histogram of sums
    plt.hist(sums/beam_area2(imagename), bins=20)
    plt.axvline(sum_mean, color='r', linestyle='--', label=f"Mean: {sum_mean:.2f}")
    plt.legend()
    plt.show()
    return(sum_mean, sum_std)

def emcee_flux_error(image,residual):
    """
    TEST: DO NOT USE THIS
    """
    import numpy as np
    import emcee
    # Define log likelihood function
    def log_likelihood(flux, image, residual, sigma):
        # Calculate model image with given flux
        model = image + flux * residual

        # Calculate chi-squared statistic
        chi2 = np.sum((model - image)**2 / sigma**2)

        # Calculate log likelihood
        log_likelihood = -0.5 * chi2

        return log_likelihood

    # Define log prior function
    def log_prior(flux):
        # Uniform prior on flux between 0 and 1000
        if 0 < flux < 1000:
            return 0.0

        return -np.inf

    # Define log probability function
    def log_probability(flux, image, residual, sigma):
        # Calculate log prior
        lp = log_prior(flux)
        if not np.isfinite(lp):
            return -np.inf

        # Calculate log likelihood
        ll = log_likelihood(flux, image, residual, sigma)

        # Calculate log posterior probability
        log_prob = lp + ll

        return log_prob

    # Generate example data
    # image = np.random.normal(10.0, 1.0, size=(100, 100))
    # residual = np.random.normal(0.0, 0.1, size=(100, 100))
    sigma = 0.1

    # Define number of dimensions and number of walkers
    ndim = 1
    nwalkers = 32

    # Initialize walkers with random values
    print('# Initialize walkers with random values')
    p0 = np.random.uniform(0.0, 100.0, size=(nwalkers, ndim))

    # Initialize sampler with log probability function
    print('# Initialize sampler with log probability function')
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_probability,
                                    args=(image, residual, sigma),threads=6)

    # Burn-in phase to reach equilibrium
    print('# Burn-in phase to reach equilibrium')
    n_burn = 100
    pos, prob, state = sampler.run_mcmc(p0, n_burn, store=False)

    # Reset sampler and run production phase
    print('# Reset sampler and run production phase')
    sampler.reset()
    n_prod = 1000
    sampler.run_mcmc(pos, n_prod)

    # Extract samples from sampler
    print('# Extract samples from sampler')
    samples = sampler.get_chain()

    # Calculate mean and standard deviation of flux estimates
    print('# Calculate mean and standard deviation of flux estimates')
    flux_samples = samples[:, 0]
    flux_mean = np.mean(flux_samples)
    flux_std = np.std(flux_samples)

    # Calculate total flux estimate and error
    print('# Calculate total flux estimate and error')
    total_flux = np.sum(image + flux_mean * residual)
    total_flux_err = flux_std * np.sum(residual)

    print("Total flux estimate:", total_flux)
    print("Total flux error:", total_flux_err)
#     print("Total flux estimate:", total_flux/beam_area2(imagelist[idx]))
#     print("Total flux error:", total_flux_err/beam_area2(imagelist[idx]))
    def plot_flux_estimate(image, residual, flux_samples):
        """
        plot_flux_estimate(image,residual,flux_samples)
        """
        # Calculate total flux estimate and error
        flux_mean = np.mean(flux_samples)
        flux_std = np.std(flux_samples)
        total_flux = np.sum(image + flux_mean * residual)/beam_area2(imagelist[idx])
        total_flux_err = flux_std * np.sum(residual)/beam_area2(imagelist[idx])

        # Create plot
        fig, ax = plt.subplots(figsize=(10, 8))

    #     # Plot image
    #     ax.imshow(image, cmap='viridis', origin='lower', alpha=0.8)
    #     ax.set_xlabel('X pixel')
    #     ax.set_ylabel('Y pixel')
    #     ax.set_title('Flux estimate for 2D image')

    #     # Add colorbar
    #     cbar = plt.colorbar(ax.imshow(image, cmap='viridis', origin='lower', alpha=0.8))
    #     cbar.ax.set_ylabel('Pixel value')

        # Plot flux samples as scatter plot
        ax.scatter(flux_samples/beam_area2(imagelist[idx]),
                   total_flux -
                   (flux_samples/beam_area2(imagelist[idx])).sum(axis=1),
                   cmap='viridis', alpha=0.5)

        # Add mean and error bars
        ax.axvline(flux_mean/beam_area2(imagelist[idx]), color='red',
                   label='Flux estimate')
        ax.axhline(total_flux - np.sum(flux_samples/beam_area2(imagelist[idx])) +
                   flux_mean/beam_area2(imagelist[idx]) * len(flux_samples),
                   color='red', linestyle='--', label='Total flux estimate')
        ax.fill_betweenx([total_flux - total_flux_err, total_flux +
                          total_flux_err], flux_mean/beam_area2(imagelist[idx]) -
                         flux_std/beam_area2(imagelist[idx]),
                         flux_mean/beam_area2(imagelist[idx]) +
                         flux_std/beam_area2(imagelist[idx]),
                         alpha=0.2, color='red', label='Total flux error')

        # Add legend
        ax.legend(loc='lower right')

        # Show plot
        plt.show()


def cdiff(data):
    '''
    Manual 'fix' of numpy's diff function.
    '''
    diff_corr = np.diff(data,append=data[-1])
    diff_corr[-1] = diff_corr[-2]
    return(diff_corr)


def get_radial_profiles(imagelist, labels_profiles=None, save_fig=None):
    '''
    Given a list of images (e.g. target images, psf, etc), plot the radial
    profile of each one altogether and retur the radiis and profiles for all.
    '''
    profiles = []
    radiis = []
    for IMAGE in imagelist:
        RR, IR = get_profile(IMAGE)
        profiles.append(IR)
        radiis.append(RR)
    return (profiles, radiis)

def get_peak_pos(imagename):
    st = imstat(imagename=imagename)
    maxpos = st['maxpos'][0:2]
    # print('Peak Pos=', maxpos)
    return (maxpos)


def get_profile(imagename, center=None,binsize=1,interpnan=True,stddev=False):
    if isinstance(imagename, str) == True:
        data_2D = load_fits_data(imagename)
    else:
        data_2D = imagename

    if stddev:
        nr, radius, profile, profile_std = azimuthalAverage(data_2D,return_nr = True,
                                                interpnan=interpnan,center=center,
                                                stddev = stddev,
                                                binsize=binsize)
        return (radius, profile, profile_std)
    else:
        nr, radius, profile = azimuthalAverage(data_2D,return_nr = True,
                                                interpnan=interpnan,center=center,
                                                binsize=binsize)
        return (radius, profile)
    


def azimuthalAverage(image, center=None, stddev=False, returnradii=False, return_nr=False,
                     binsize=1.0, weights=None, steps=False, interpnan=False, left=None, right=None,
                     mask=None):
    """
    Calculate the azimuthally averaged radial profile.
    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is
             None, which then uses the center of the image (including
             fractional pixels).
    stddev - if specified, return the azimuthal standard deviation instead of the average
    returnradii - if specified, return (radii_array,radial_profile)
    return_nr   - if specified, return number of pixels per radius *and* radius
    binsize - size of the averaging bin.  Can lead to strange results if
        non-binsize factors are used to specify the center and the binsize is
        too large
    weights - can do a weighted average instead of a simple average if this keyword parameter
        is set.  weights.shape must = image.shape.  weighted stddev is undefined, so don't
        set weights and stddev.
    steps - if specified, will return a double-length bin array and radial
        profile so you can plot a step-form radial profile (which more accurately
        represents what's going on)
    interpnan - Interpolate over NAN values, i.e. bins where there is no data?
        left,right - passed to interpnan; they set the extrapolated values
    mask - can supply a mask (boolean array same size as image with True for OK and False for not)
        to average over only select data.
    If a bin contains NO DATA, it will have a NAN value because of the
    divide-by-sum-of-weights component.  I think this is a useful way to denote
    lack of data, but users let me know if an alternative is prefered...

    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if center is None:
        y0max, x0max = nd.maximum_position((image))
        # center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])
        center = np.array([x0max, y0max])

    r = np.hypot(x - center[0], y - center[1])

    if weights is None:
        weights = np.ones(image.shape)
    elif stddev:
        raise ValueError("Weighted standard deviation is not defined.")

    if mask is None:
        mask = np.ones(image.shape, dtype='bool')
    # obsolete elif len(mask.shape) > 1:
    # obsolete     mask = mask.ravel()

    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)
    nbins = int(np.round(np.nanmax(r) / binsize) + 1)
    maxbin = nbins * binsize
    bins = np.linspace(0, maxbin, nbins + 1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:] + bins[:-1]) / 2.0

    # how many per bin (i.e., histogram)?
    # there are never any in bin 0, because the lowest index returned by digitize is 1
    # nr = np.bincount(whichbin)[1:]
    nr = np.histogram(r, bins, weights=mask.astype('int'))[0]

    # recall that bins are from 1 to nbins (which is expressed in array terms by arange(nbins)+1 or range(1,nbins+1) )
    # radial_prof.shape = bin_centers.shape
    if stddev:
        # Find out which radial bin each point in the map belongs to
        whichbin = np.digitize(r.flat, bins)
        # This method is still very slow; is there a trick to do this with histograms?
        radial_prof_std = np.array([np.nanstd(image.flat[mask.flat * (whichbin == b)]) for b in range(1, nbins + 1)])
    
    radial_prof = np.histogram(r, bins, weights=(image * weights * mask))[0] / \
                    np.histogram(r, bins, weights=(mask * weights))[0]

    if interpnan:
        radial_prof = np.interp(bin_centers, bin_centers[radial_prof == radial_prof],
                                radial_prof[radial_prof == radial_prof], left=left, right=right)

    if steps:
        xarr = np.array(zip(bins[:-1], bins[1:])).ravel()
        yarr = np.array(zip(radial_prof, radial_prof)).ravel()
        return xarr, yarr
    elif returnradii:
        if stddev:
            return bin_centers, radial_prof, radial_prof_std
        else:
            return bin_centers, radial_prof
    elif return_nr:
        if stddev:
            return nr, bin_centers, radial_prof, radial_prof_std
        else:
            return nr, bin_centers, radial_prof
    else:
        return radial_prof


def azimuthalAverageBins(image, azbins, symmetric=None, center=None, **kwargs):
    """ Compute the azimuthal average over a limited range of angles
    kwargs are passed to azimuthalAverage """
    y, x = np.indices(image.shape)
    if center is None:
        center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])
    r = np.hypot(x - center[0], y - center[1])
    theta = np.arctan2(x - center[0], y - center[1])
    theta[theta < 0] += 2 * np.pi
    theta_deg = theta * 180.0 / np.pi

    if isinstance(azbins, np.ndarray):
        pass
    elif isinstance(azbins, int):
        if symmetric == 2:
            azbins = np.linspace(0, 90, azbins)
            theta_deg = theta_deg % 90
        elif symmetric == 1:
            azbins = np.linspace(0, 180, azbins)
            theta_deg = theta_deg % 180
        elif azbins == 1:
            return azbins, azimuthalAverage(image, center=center, returnradii=True, **kwargs)
        else:
            azbins = np.linspace(0, 359.9999999999999, azbins)
    else:
        raise ValueError("azbins must be an ndarray or an integer")

    azavlist = []
    for blow, bhigh in zip(azbins[:-1], azbins[1:]):
        mask = (theta_deg > (blow % 360)) * (theta_deg < (bhigh % 360))
        rr, zz = azimuthalAverage(image, center=center, mask=mask, returnradii=True, **kwargs)
        azavlist.append(zz)

    return azbins, rr, azavlist


def radialAverage(image, center=None, stddev=False, returnAz=False, return_naz=False,
                  binsize=1.0, weights=None, steps=False, interpnan=False, left=None, right=None,
                  mask=None, symmetric=None):
    """
    Calculate the radially averaged azimuthal profile.
    (this code has not been optimized; it could be speed boosted by ~20x)
    image - The 2D image
    center - The [x,y] pixel coordinates used as the center. The default is
             None, which then uses the center of the image (including
             fractional pixels).
    stddev - if specified, return the radial standard deviation instead of the average
    returnAz - if specified, return (azimuthArray,azimuthal_profile)
    return_naz   - if specified, return number of pixels per azimuth *and* azimuth
    binsize - size of the averaging bin.  Can lead to strange results if
        non-binsize factors are used to specify the center and the binsize is
        too large
    weights - can do a weighted average instead of a simple average if this keyword parameter
        is set.  weights.shape must = image.shape.  weighted stddev is undefined, so don't
        set weights and stddev.
    steps - if specified, will return a double-length bin array and azimuthal
        profile so you can plot a step-form azimuthal profile (which more accurately
        represents what's going on)
    interpnan - Interpolate over NAN values, i.e. bins where there is no data?
        left,right - passed to interpnan; they set the extrapolated values
    mask - can supply a mask (boolean array same size as image with True for OK and False for not)
        to average over only select data.
    If a bin contains NO DATA, it will have a NAN value because of the
    divide-by-sum-of-weights component.  I think this is a useful way to denote
    lack of data, but users let me know if an alternative is prefered...

    """
    # Calculate the indices from the image
    y, x = np.indices(image.shape)

    if center is None:
        center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])

    r = np.hypot(x - center[0], y - center[1])
    theta = np.arctan2(x - center[0], y - center[1])
    theta[theta < 0] += 2 * np.pi
    theta_deg = theta * 180.0 / np.pi
    maxangle = 360

    if weights is None:
        weights = np.ones(image.shape)
    elif stddev:
        raise ValueError("Weighted standard deviation is not defined.")

    if mask is None:
        # mask is only used in a flat context
        mask = np.ones(image.shape, dtype='bool').ravel()
    elif len(mask.shape) > 1:
        mask = mask.ravel()

    # allow for symmetries
    if symmetric == 2:
        theta_deg = theta_deg % 90
        maxangle = 90
    elif symmetric == 1:
        theta_deg = theta_deg % 180
        maxangle = 180

    # the 'bins' as initially defined are lower/upper bounds for each bin
    # so that values will be in [lower,upper)
    nbins = int(np.round(maxangle / binsize))
    maxbin = nbins * binsize
    bins = np.linspace(0, maxbin, nbins + 1)
    # but we're probably more interested in the bin centers than their left or right sides...
    bin_centers = (bins[1:] + bins[:-1]) / 2.0

    # Find out which azimuthal bin each point in the map belongs to
    whichbin = np.digitize(theta_deg.flat, bins)

    # how many per bin (i.e., histogram)?
    # there are never any in bin 0, because the lowest index returned by digitize is 1
    nr = np.bincount(whichbin)[1:]

    # recall that bins are from 1 to nbins (which is expressed in array terms by arange(nbins)+1 or range(1,nbins+1) )
    # azimuthal_prof.shape = bin_centers.shape
    if stddev:
        azimuthal_prof = np.array([image.flat[mask * (whichbin == b)].std() for b in range(1, nbins + 1)])
    else:
        azimuthal_prof = np.array(
            [(image * weights).flat[mask * (whichbin == b)].sum() / weights.flat[mask * (whichbin == b)].sum() for b in
             range(1, nbins + 1)])

    # import pdb; pdb.set_trace()

    if interpnan:
        azimuthal_prof = np.interp(bin_centers,
                                   bin_centers[azimuthal_prof == azimuthal_prof],
                                   azimuthal_prof[azimuthal_prof == azimuthal_prof],
                                   left=left, right=right)

    if steps:
        xarr = np.array(zip(bins[:-1], bins[1:])).ravel()
        yarr = np.array(zip(azimuthal_prof, azimuthal_prof)).ravel()
        return xarr, yarr
    elif returnAz:
        return bin_centers, azimuthal_prof
    elif return_naz:
        return nr, bin_centers, azimuthal_prof
    else:
        return azimuthal_prof


def radialAverageBins(image, radbins, corners=True, center=None, **kwargs):
    """ Compute the radial average over a limited range of radii """
    y, x = np.indices(image.shape)
    if center is None:
        center = np.array([(x.max() - x.min()) / 2.0, (y.max() - y.min()) / 2.0])
    r = np.hypot(x - center[0], y - center[1])

    if isinstance(radbins, np.ndarray):
        pass
    elif isinstance(radbins, int):
        if radbins == 1:
            return radbins, radialAverage(image, center=center, returnAz=True, **kwargs)
        elif corners:
            radbins = np.linspace(0, r.max(), radbins)
        else:
            radbins = np.linspace(0, np.max(np.abs(np.array([x - center[0], y - center[1]]))), radbins)
    else:
        raise ValueError("radbins must be an ndarray or an integer")

    radavlist = []
    for blow, bhigh in zip(radbins[:-1], radbins[1:]):
        mask = (r < bhigh) * (r > blow)
        az, zz = radialAverage(image, center=center, mask=mask, returnAz=True, **kwargs)
        radavlist.append(zz)

    return radbins, az, radavlist
