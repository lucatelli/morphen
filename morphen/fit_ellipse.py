"""
Ellipse Fitting Algorithm.

Earlier version of the code was inspired by:
http://nicky.vanforeest.com/misc/fitEllipse/fitEllipse.html

The notation has changed to that of  www.wikipedia.com

Author: Fabricio Ferrari
Co-Author: Geferson Lucatelli
v1 @ 2018 -- Morfometryka Utils

Author and Original Ideas: Geferson Lucatelli
v2 @ 2024 -- Morfometryka and Morphen Utils
Now, the code was been modified to include new features and improve the code quality.
Thanks to Claude GPT.
"""
import numpy as np
from numpy.linalg import eig, eigvals, det
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from astropy.stats import mad_std
from scipy import optimize
from dataclasses import dataclass
from typing import Tuple, List, Optional
import astropy.io.fits as fits

@dataclass
class EllipseParams:
    """Class to store ellipse parameters"""
    x0: float  # center x
    y0: float  # center y
    a: float   # semi-major axis
    b: float   # semi-minor axis
    phi: float # rotation angle in radians
    
    def to_array(self) -> np.ndarray:
        return np.array([self.x0, self.y0, self.a, self.b, self.phi])
    
    @classmethod
    def from_array(cls, params: np.ndarray) -> 'EllipseParams':
        return cls(params[0], params[1], params[2], params[3], params[4])

def ellipse_distance(params: np.ndarray, points: np.ndarray) -> np.ndarray:
    """
    Calculate the distance from points to an ellipse.
    
    Parameters:
        params: array [x0, y0, a, b, phi] - ellipse parameters
        points: array of shape (N, 2) containing the points
    
    Returns:
        Array of distances from points to ellipse
    """
    x0, y0, a, b, phi = params
    
    # Translate points to origin
    xt = points[:, 0] - x0
    yt = points[:, 1] - y0
    
    # Rotate points to align with ellipse axes
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    xr = xt * cos_phi + yt * sin_phi
    yr = -xt * sin_phi + yt * cos_phi
    
    # Calculate normalized distances
    return np.abs(xr**2/a**2 + yr**2/b**2 - 1.0)

def fit_ellipse_nonlinear(x: np.ndarray, y: np.ndarray, 
                         initial_guess: Optional[EllipseParams] = None) -> Optional[EllipseParams]:
    """
    Fit an ellipse to points using non-linear least squares.
    
    Parameters:
        x, y: arrays of point coordinates
        initial_guess: initial parameters for optimization
    
    Returns:
        EllipseParams object with fitted parameters or None if fit fails
    """
    points = np.column_stack([x, y])
    
    # If no initial guess provided, estimate one
    if initial_guess is None:
        x_mean, y_mean = np.mean(x), np.mean(y)
        x_std, y_std = np.std(x), np.std(y)
        initial_guess = EllipseParams(
            x0=x_mean,
            y0=y_mean,
            a=2*x_std,
            b=2*y_std,
            phi=0.0
        )
    
    # Define bounds for parameters
    bounds = optimize.Bounds(
        lb=[x.min(), y.min(), 0, 0, -np.pi],
        ub=[x.max(), y.max(), np.inf, np.inf, np.pi]
    )
    
    try:
        # Perform optimization
        result = optimize.minimize(
            lambda p: np.sum(ellipse_distance(p, points)**2),
            initial_guess.to_array(),
            bounds=bounds,
            method='L-BFGS-B'
        )
        
        if not result.success:
            return None
            
        # Ensure a >= b
        params = result.x
        if params[3] > params[2]:
            params[2], params[3] = params[3], params[2]
            params[4] += np.pi/2
            
        return EllipseParams.from_array(params)
    
    except:
        return None

def remove_outliers(x: np.ndarray, y: np.ndarray, n_sigma: float = 3.0) -> Tuple[np.ndarray, np.ndarray]:
    """Remove spatial outliers using sigma clipping."""
    x_clean, _, _ = sigmaclip(x, low=n_sigma, high=n_sigma)
    y_clean, _, _ = sigmaclip(y, low=n_sigma, high=n_sigma)
    
    mask_x = np.isin(x, x_clean)
    mask_y = np.isin(y, y_clean)
    mask = mask_x & mask_y
    
    return x[mask], y[mask]

def get_ellipse_points(params: EllipseParams, t: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Generate points along ellipse at parameter values t."""
    cos_t = np.cos(t)
    sin_t = np.sin(t)
    cos_phi = np.cos(params.phi)
    sin_phi = np.sin(params.phi)
    
    x = params.x0 + params.a * cos_t * cos_phi - params.b * sin_t * sin_phi
    y = params.y0 + params.a * cos_t * sin_phi + params.b * sin_t * cos_phi
    
    return x, y

def calculate_radial_profiles(fits, 
                            intensities):
    """
    Calculate radial profiles along semi-major and semi-minor axes.
    
    Parameters:
        fits: List of EllipseParams objects for each isophote
        intensities: Array of intensity values corresponding to each fit
    
    Returns:
        profiles: Dictionary containing radial measurements and corresponding quantities
        r_maj: Radial distances along semi-major axis
        r_min: Radial distances along semi-minor axis
    """
    # Extract parameters
    a_values = np.array([fit.a for fit in fits])
    b_values = np.array([fit.b for fit in fits])
    qs = b_values / a_values
    pas = np.array([np.rad2deg(fit.phi) for fit in fits])
    
    # Calculate geometric mean radius (equivalent to sqrt(a*b))
    r_geo = np.sqrt(a_values * b_values)
    
    # Calculate semi-major axis radius (a)
    r_maj = a_values
    
    # Calculate semi-minor axis radius (b)
    r_min = b_values
    
    # Create profiles dictionary
    profiles = {
        'r_geo': r_geo,          # geometric mean radius
        'r_maj': r_maj,          # semi-major axis radius
        'r_min': r_min,          # semi-minor axis radius
        'q': qs,                 # axis ratio (b/a)
        'pa': pas,               # position angle
        'intensity': intensities, # intensity values
        'mu': -2.5 * np.log10(intensities)  # surface brightness in mag/arcsec^2
    }
    
    return profiles, r_maj, r_min

def plot_radial_profiles(profiles, 
                        region_split= None,
                        save_name= None):
    """
    Plot radial profiles of various quantities.
    
    Parameters:
        profiles: Dictionary containing radial measurements
        region_split: Index to split inner/outer regions
        save_name: Filename to save the plots
    """
    fig, axs = plt.subplots(2, 2, figsize=(8, 8))
    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    # # Plot surface brightness profile
    # axs[0, 0].plot(profiles['r_maj'], profiles['mu'], 'k.', ms=3)
    # axs[0, 0].set_xlabel('Semi-major axis (pixels)')
    # axs[0, 0].set_ylabel('μ (mag/arcsec²)')
    # axs[0, 0].invert_yaxis()
    # axs[0, 0].grid(True, alpha=0.3)

    # Plot surface brightness profile
    axs[0, 0].plot(profiles['r_min'], np.asarray(profiles['intensity'])*1000, 'k.', ms=3)
    axs[0, 0].set_xlabel('Semi-minor axis (pixels)')
    axs[0, 0].set_ylabel('$\log I$ (mJy/beam)')
    # axs[0, 0].invert_yaxis()
    axs[0, 0].set_yscale('log')
    axs[0, 0].grid(True, alpha=0.3)
    
    # Plot axis ratio profile
    axs[0, 1].plot(profiles['r_maj'], profiles['q'], 'b.', ms=3)
    axs[0, 1].set_xlabel('Semi-major axis (pixels)')
    axs[0, 1].set_ylabel('Axis ratio (b/a)')
    axs[0, 1].set_ylim(0, 1)
    axs[0, 1].grid(True, alpha=0.3)
    
    
    # Plot position angle profile
    axs[1, 0].plot(profiles['r_maj'], profiles['pa'], 'r.', ms=3)
    axs[1, 0].set_xlabel('Semi-major axis (pixels)')
    axs[1, 0].set_ylabel('Position Angle (degrees)')
    axs[1, 0].set_ylim(0, 180)
    axs[1, 0].grid(True, alpha=0.3)
    
    # Plot intensity profile
    axs[1, 1].plot(profiles['r_maj'], np.asarray(profiles['intensity'])*1000, 'g.', ms=3)
    axs[1, 1].set_xlabel('Semi-major axis (pixels)')
    axs[1, 1].set_ylabel('$\log I$ (mJy/beam)')
    axs[1, 1].set_yscale('log')
    axs[1, 1].grid(True, alpha=0.3)
    
    if region_split is not None:
        for ax in axs.flat:
            ax.axvline(profiles['r_maj'][region_split], color='k', 
                      linestyle='--', alpha=0.5)
    
    plt.tight_layout()
    
    if save_name:
        plt.savefig(save_name, dpi=300, bbox_inches='tight')
        # plt.clf()
        # plt.close()
        plt.show()

# Modify the fit_isophotes function to include these new features:

def fit_isophotes(image, 
                  intensity_levels, 
                  region_split= None, 
                  plot_results = True,
                  plot_profiles = True,
                  save_name= None,
                  save_profiles= None):
    """
    Fit ellipses to isophotes in an astronomical image using non-linear least squares.
    Parameters:
        image: 2D numpy array
        intensity_levels: array of intensity values to fit
        region_split: index to split inner/outer regions
        plot_results: whether to plot the results
        save_name: filename to save the plot
        plot_profiles: whether to plot radial profiles
        save_profiles: filename to save the profile plots
    
    Returns:
        Tuple containing:
        (q_inner, q_outer, pa_inner, pa_outer, q_median, pa_median,
         x0_median, y0_median, x0_inner, y0_inner, x0_outer, y0_outer,
         profiles) 
    """
    if region_split is None:
        region_split = len(intensity_levels) // 2
    
    results = []
    previous_fit = None  # Store previous fit to use as initial guess
    
    for intensity in intensity_levels:
        # Find pixels near the isophote
        delta = 0.5 * intensity
        y, x = np.where((image > intensity - delta) & (image < intensity + delta))
        
        if len(x) < 6:
            continue
        
        # Remove outliers
        x, y = remove_outliers(x, y)
        
        if len(x) < 6:
            continue
        
        # Fit ellipse using previous fit as initial guess
        fit = fit_ellipse_nonlinear(x, y, initial_guess=previous_fit)
        if fit is None:
            continue
            
        # Store fit for next iteration
        previous_fit = fit
        results.append((fit, intensity))
    
    if not results:
        raise ValueError("No valid ellipses could be fitted")
    
    # Extract parameters
    fits, intensities = zip(*results)
    
    # Calculate statistics
    qs = np.array([fit.b/fit.a for fit in fits])
    pas = np.array([np.rad2deg(fit.phi) for fit in fits])
    x0s = np.array([fit.x0 for fit in fits])
    y0s = np.array([fit.y0 for fit in fits])
    
    # Plot if requested
    if plot_results:
        plt.figure(figsize=(5, 5))
        plt.imshow(np.log(image), origin='lower', cmap='magma_r')
        
        t = np.linspace(0, 2*np.pi, 100)
        for i, fit in enumerate(fits):
            xe, ye = get_ellipse_points(fit, t)
            color = 'green' if i < region_split else 'yellow'
            plt.plot(xe, ye, '-', color=color, alpha=0.6, lw=0.4)
            
        if save_name:
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            # plt.clf()
            # plt.close()
            plt.show()
    
    # Calculate radial profiles
    profiles, r_maj, r_min = calculate_radial_profiles(fits, intensities)
    
    # Plot profiles if requested
    if plot_profiles:
        plot_radial_profiles(profiles, region_split, save_profiles)


    # Calculate statistics
    stats = {
        'q_inner': np.nanmedian(qs[:region_split]),
        'q_outer': np.nanmedian(qs[region_split:]),
        'pa_inner': np.nanmedian(pas[:region_split]),
        'pa_outer': np.nanmedian(pas[region_split:]),
        'q_median': np.nanmedian(qs),
        'pa_median': np.nanmedian(pas),
        'x0_median': np.nanmedian(x0s),
        'y0_median': np.nanmedian(y0s),
        'x0_inner': np.nanmedian(x0s[:region_split]),
        'y0_inner': np.nanmedian(y0s[:region_split]),
        'x0_outer': np.nanmedian(x0s[region_split:]),
        'y0_outer': np.nanmedian(y0s[region_split:])
    }
    return (stats['q_inner'], stats['q_outer'],
            stats['pa_inner'], stats['pa_outer'],
            stats['q_median'], stats['pa_median'],
            stats['x0_median'], stats['y0_median'],
            stats['x0_inner'], stats['y0_inner'],
            stats['x0_outer'], stats['y0_outer'],
            profiles)  # Added profiles to return


def fit_ellipse_to_galaxy(gal_image, Isequence=None, region_split=None, 
                          save_name=None,plot_results=False,plot_profiles=False):
    """
    Main function to fit ellipses to galaxy image.
    
    Parameters:
        gal_image: 2D numpy array containing the galaxy image
        Isequence: array of intensity values to fit. If None, will be generated automatically
        region_split: index to split inner/outer regions. If None, uses middle point
        save_name: filename to save the plots
        
    Returns:
        Tuple containing:
        (q_inner, q_outer, pa_inner, pa_outer, q_median, pa_median,
         x0_median, y0_median, x0_inner, y0_inner, x0_outer, y0_outer,
         profiles)
    """
    # Generate intensity sequence if not provided
    if Isequence is None:
        Isteps = 100
        Imin = 0.01 * np.nanstd(gal_image)
        Imax = 0.99 * np.nanmax(gal_image)
        Isequence = np.geomspace(Imax, Imin, Isteps)
    
    # Perform the fitting
    return fit_isophotes(
        image=gal_image,
        intensity_levels=Isequence,
        region_split=region_split,
        plot_results=plot_results,
        plot_profiles=plot_profiles,
        save_name=save_name+'_photo.jpg' if save_name else None,
        save_profiles=save_name+'_profiles.jpg' if save_name else None
    )