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

v3 @ 11/2025 -- Enhanced with LMFIT for robust nonlinear fitting
Help from Claude Code. 

"""
import numpy as np
from numpy.linalg import eig, eigvals, det
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from astropy.stats import mad_std
from scipy import optimize
import lmfit
from dataclasses import dataclass
from typing import Tuple, List, Optional, Union
import astropy.io.fits as fits

@dataclass
class EllipseParams:
    """Class to store ellipse parameters"""
    x0: float  # center x
    y0: float  # center y
    a: float   # semi-major axis
    b: float   # semi-minor axis
    phi: float # rotation angle in radians
    
    chi2: Optional[float] = None
    redchi: Optional[float] = None
    success: bool = True
    
    def to_array(self) -> np.ndarray:
        return np.array([self.x0, self.y0, self.a, self.b, self.phi])
    
    @classmethod
    def from_array(cls, params: np.ndarray) -> 'EllipseParams':
        return cls(params[0], params[1], params[2], params[3], params[4])
    
    @classmethod
    def from_lmfit(cls, result: lmfit.minimizer.MinimizerResult) -> 'EllipseParams':
        """Create EllipseParams from LMFIT result."""
        p = result.params
        return cls(
            x0=p['x0'].value,
            y0=p['y0'].value,
            a=p['a'].value,
            b=p['b'].value,
            phi=p['phi'].value,
            chi2=getattr(result, 'chisqr', None),
            redchi=getattr(result, 'redchi', None),
            success=result.success
        )

def ellipse_residual(params: lmfit.Parameters, 
                     points: np.ndarray) -> np.ndarray:
    """
    Calculate residuals for ellipse fitting (for LMFIT).
    
    The residual is defined as the deviation from the ellipse equation:
    (x'/a)^2 + (y'/b)^2 = 1
    
    where x', y' are coordinates rotated to align with ellipse axes.
    
    Parameters:
        params: LMFIT Parameters object with x0, y0, a, b, phi
        points: array of shape (N, 2) containing the points
    
    Returns:
        Array of residuals (deviations from ellipse)
    """
    x0 = params['x0'].value
    y0 = params['y0'].value
    a = params['a'].value
    b = params['b'].value
    phi = params['phi'].value
    
    # Translate points to origin
    xt = points[:, 0] - x0
    yt = points[:, 1] - y0
    
    # Rotate points to align with ellipse axes
    cos_phi = np.cos(phi)
    sin_phi = np.sin(phi)
    xr = xt * cos_phi + yt * sin_phi
    yr = -xt * sin_phi + yt * cos_phi
    
    # Residual: deviation from ellipse equation
    # For points on the ellipse: (xr/a)^2 + (yr/b)^2 = 1
    # Residual measures how far from 1 this value is
    return (xr**2 / a**2 + yr**2 / b**2) - 1.0


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


def create_ellipse_params(x: np.ndarray, y: np.ndarray,
                          initial_guess: Optional['EllipseParams'] = None,
                          fix_center: bool = False,
                          fix_angle: bool = False,
                          dx_dy: Optional[Union[float, Tuple[float, float]]] = None) -> lmfit.Parameters:
    """
    Create LMFIT Parameters object for ellipse fitting.
    
    Parameters:
        x, y: arrays of point coordinates
        initial_guess: initial parameters for optimization
        fix_center: if True, fix the center position
        fix_angle: if True, fix the position angle
        dx_dy: maximum allowed offset from initial center position.
               Can be a single float (same for x and y) or a tuple (dx, dy).
               If None, bounds are based on data extent.
               Ignored if fix_center=True.
    
    Returns:
        LMFIT Parameters object
    """
    # Estimate initial values if not provided
    if initial_guess is None:
        x_mean, y_mean = np.mean(x), np.mean(y)
        x_std, y_std = np.std(x), np.std(y)
        x0_init = x_mean
        y0_init = y_mean
        a_init = max(2 * x_std, 1.0)
        b_init = max(2 * y_std, 1.0)
        phi_init = 0.0
    else:
        x0_init = initial_guess.x0
        y0_init = initial_guess.y0
        a_init = initial_guess.a
        b_init = initial_guess.b
        phi_init = initial_guess.phi
    
    # Create parameters with bounds
    params = lmfit.Parameters()
    
    # Center position bounds
    if fix_center:
        # Completely fixed
        params.add('x0', value=x0_init, vary=False)
        params.add('y0', value=y0_init, vary=False)
    elif dx_dy is not None:
        # Constrained to within dx_dy of initial position
        if isinstance(dx_dy, (int, float)):
            dx, dy = float(dx_dy), float(dx_dy)
        else:
            dx, dy = dx_dy
        
        params.add('x0', value=x0_init,
                   min=x0_init - dx,
                   max=x0_init + dx,
                   vary=True)
        params.add('y0', value=y0_init,
                   min=y0_init - dy,
                   max=y0_init + dy,
                   vary=True)
    else:
        # Default: allow freedom around data extent
        x_range = x.max() - x.min()
        y_range = y.max() - y.min()
        
        params.add('x0', value=x0_init, 
                   min=x.min() - 0.5 * x_range, 
                   max=x.max() + 0.5 * x_range,
                   vary=True)
        params.add('y0', value=y0_init, 
                   min=y.min() - 0.5 * y_range, 
                   max=y.max() + 0.5 * y_range,
                   vary=True)
    
    # Semi-axes - must be positive
    # Upper bound based on data extent
    x_range = x.max() - x.min()
    y_range = y.max() - y.min()
    max_extent = max(x_range, y_range) * 2
    params.add('a', value=a_init, min=0.5, max=max_extent)
    params.add('b', value=b_init, min=0.5, max=max_extent)
    
    # Position angle - constrained to [-pi, pi]
    params.add('phi', value=phi_init, min=-np.pi, max=np.pi,
               vary=not fix_angle)
    
    return params


def fit_ellipse_nonlinear(x: np.ndarray, y: np.ndarray, 
                          initial_guess: Optional[EllipseParams] = None,
                          fix_center: bool = False,
                          fix_angle: bool = False,
                          dx_dy: Optional[Union[float, Tuple[float, float]]] = None,
                          robust: bool = True,
                          max_nfev: int = 10000,
                          verbose: int = 0) -> Optional[EllipseParams]:
    """
    Fit an ellipse to points using non-linear least squares with LMFIT.
    
    This implementation uses robust fitting with Cauchy loss function
    to handle outliers in the isophote points.
    
    Parameters:
        x, y: arrays of point coordinates
        initial_guess: initial parameters for optimization
        fix_center: if True, fix the center position during fitting
        fix_angle: if True, fix the position angle during fitting
        dx_dy: maximum allowed offset from initial center position.
               Can be a single float (same for x and y) or a tuple (dx, dy).
               If None, bounds are based on data extent.
               Ignored if fix_center=True.
        robust: if True, use Cauchy loss for robustness to outliers
        max_nfev: maximum number of function evaluations
        verbose: verbosity level (0=silent, 1=summary, 2=detailed)
    
    Returns:
        EllipseParams object with fitted parameters or None if fit fails
    """
    points = np.column_stack([x, y])
    
    fit_params = create_ellipse_params(x, y, initial_guess, 
                                        fix_center, fix_angle, dx_dy)
    mini = lmfit.Minimizer(
        ellipse_residual, 
        fit_params,
        fcn_args=(points,),
        reduce_fcn='neglogcauchy' if robust else None,
        max_nfev=max_nfev,
        nan_policy='omit'
    )
    
    try:
        result = mini.minimize(
            method='least_squares',
            max_nfev=max_nfev,
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
            jac='3-point',
            x_scale = 'jac', f_scale = 1.0,
            verbose=verbose,
            tr_solver='exact',
            loss='cauchy' if robust else 'linear',
            tr_options={'regularize': True}
        )
        
        if not result.success:
            # Try fallback to Nelder-Mead if least_squares fails
            if verbose > 0:
                print("least_squares failed, trying Nelder-Mead fallback...")
            result = mini.minimize(method='nelder', max_nfev=max_nfev)
        
        if not result.success:
            return None
        
        # Extract parameters
        params = EllipseParams.from_lmfit(result)
        
        # Ensure a >= b (semi-major >= semi-minor)
        if params.b > params.a:
            params.a, params.b = params.b, params.a
            params.phi += np.pi / 2
            # Normalize phi to [-pi, pi]
            while params.phi > np.pi:
                params.phi -= np.pi
            while params.phi < -np.pi:
                params.phi += np.pi
        
        return params
    
    except Exception as e:
        if verbose > 0:
            print(f"Ellipse fitting failed: {e}")
        return None


def fit_ellipse_nonlinear_scipy(x: np.ndarray, y: np.ndarray, 
                         initial_guess: Optional[EllipseParams] = None) -> Optional[EllipseParams]:
    """
    Fit an ellipse to points using non-linear least squares (scipy version).
    
    ALTERNATIVE: This was the version prior to LMFIT (in morphen).
    
    Parameters:
        x, y: arrays of point coordinates
        initial_guess: initial parameters for optimization
    
    Returns:
        EllipseParams object with fitted parameters or None if fit fails
    """
    points = np.column_stack([x, y])
    
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
    
    bounds = optimize.Bounds(
        lb=[x.min(), y.min(), 0, 0, -np.pi],
        ub=[x.max(), y.max(), np.inf, np.inf, np.pi]
    )
    
    try:
        result = optimize.minimize(
            lambda p: np.sum(ellipse_distance(p, points)**2),
            initial_guess.to_array(),
            bounds=bounds,
            method='L-BFGS-B'
        )
        
        if not result.success:
            return None
            
        params = result.x
        if params[3] > params[2]:
            params[2], params[3] = params[3], params[2]
            params[4] += np.pi/2
            
        return EllipseParams.from_array(params)
    
    except:
        return None

def remove_outliers(x: np.ndarray, y: np.ndarray, n_sigma: float = 2.0) -> Tuple[np.ndarray, np.ndarray]:
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

def calculate_radial_profiles(fit_results, 
                            intensities):
    """
    Calculate radial profiles along semi-major and semi-minor axes.
    
    Parameters:
        image: List of EllipseParams objects for each isophote
        intensities: Array of intensity values corresponding to each fit
    
    Returns:
        profiles: Dictionary containing radial measurements and corresponding quantities
        r_maj: Radial distances along semi-major axis
        r_min: Radial distances along semi-minor axis
    """
    a_values = np.array([fit.a for fit in fit_results])
    b_values = np.array([fit.b for fit in fit_results])
    qs = b_values / a_values
    pas = np.array([np.rad2deg(fit.phi) for fit in fit_results])
    r_geo = np.sqrt(a_values * b_values)
    r_maj = a_values
    r_min = b_values
    
    profiles = {
        'r_geo': r_geo,          # geometric mean radius
        'r_maj': r_maj,          # semi-major axis radius
        'r_min': r_min,          # semi-minor axis radius
        'q': qs,                 # axis ratio (b/a)
        'pa': pas,               # position angle
        'intensity': intensities, # intensity values
        'mu': -2.5 * np.log10(intensities)  # surface brightness in mag/arcsec^2?
    }
    
    return profiles, r_maj, r_min

def plot_ellipse_profiles(profiles, 
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
    axs[1, 0].set_ylim(-180, 180)
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
                  save_profiles= None,
                  robust: bool = True,
                  fix_center: bool = False,
                  dx_dy: Optional[Union[float, Tuple[float, float]]] = None,
                  initial_params: Optional[EllipseParams] = None,
                  verbose: int = 0):
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
        robust: if True, use Cauchy loss for robustness to outliers
        fix_center: if True, fix center position after initial fit
        dx_dy: maximum allowed offset from initial center position.
               Can be a single float (same for x and y) or a tuple (dx, dy).
               Applied after the first successful fit. Ignored if fix_center=True.
        initial_params: initial EllipseParams for the first fit. If None, 
                        parameters are estimated from the data.
        verbose: verbosity level for fitting
    
    Returns:
        Tuple containing:
        (q_inner, q_outer, pa_inner, pa_outer, q_median, pa_median,
         x0_median, y0_median, x0_inner, y0_inner, x0_outer, y0_outer,
         profiles) 
    """
    if region_split is None:
        region_split = len(intensity_levels) // 2
    
    results = []
    previous_fit = initial_params  # Use provided initial params or None
    
    for i, intensity in enumerate(intensity_levels):
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
        # Optionally fix center or constrain with dx_dy after first successful fit
        use_fix_center = fix_center and (previous_fit is not None)
        use_dx_dy = dx_dy if (previous_fit is not None and not fix_center) else None
        
        fit = fit_ellipse_nonlinear(x, y, 
                                     initial_guess=previous_fit,
                                     fix_center=use_fix_center,
                                     dx_dy=use_dx_dy,
                                     robust=robust,
                                     verbose=verbose)
        if fit is None:
            continue
            
        previous_fit = fit
        results.append((fit, intensity))
    
    if not results:
        raise ValueError("No valid ellipses could be fitted")
    
    fit_results, intensities = zip(*results)    
    qs = np.array([fit.b/fit.a for fit in fit_results])
    pas = np.array([np.rad2deg(fit.phi) for fit in fit_results])
    x0s = np.array([fit.x0 for fit in fit_results])
    y0s = np.array([fit.y0 for fit in fit_results])
    
    # Plot if requested
    if plot_results:
        from astropy.visualization import simple_norm
        plt.figure(figsize=(5, 5))
        vmax=np.arcsinh(np.nanpercentile(image, 98.5))
        vmin=np.arcsinh(np.nanpercentile(image, 0.5))
        norm = simple_norm(image, 
                           stretch='asinh', 
                           asinh_a=0.07, 
                           vmin=vmin,
                           vmax=vmax)
        # plt.imshow(np.arcsinh(image), 
        #            origin='lower', cmap='Greys',
        #         #    vmax = np.nanmax(0.1*np.arcsinh(image)),
        #         #    vmin = np.arcsinh(1*mad_std(image[image>0],ignore_nan=True)),
        #            vmax=np.arcsinh(np.nanpercentile(image, 98.5)),
        #            vmin=np.arcsinh(np.nanpercentile(image, 0.5))
        #            )
        plt.imshow(image, 
                   origin='lower', cmap='RdBu_r',
                   norm = norm)
        
        t = np.linspace(0, 2*np.pi, 100)
        for i, fit in enumerate(fit_results):
            xe, ye = get_ellipse_points(fit, t)
            color = 'green' if i < region_split else 'orange'
            plt.plot(xe, ye, '-', color=color, alpha=0.7, lw=0.8)
            
        if save_name:
            plt.savefig(save_name, dpi=300, bbox_inches='tight')
            # plt.clf()
            # plt.close()
            plt.show()
    
    profiles, r_maj, r_min = calculate_radial_profiles(fit_results, intensities)
    if plot_profiles:
        plot_ellipse_profiles(profiles, region_split, save_profiles)

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
            profiles, fit_results)


def fit_ellipse_to_galaxy(gal_image, Isequence=None, region_split=None, 
                          save_name=None, plot_results=True, plot_profiles=True,
                          robust=True, fix_center=False, 
                          dx_dy: Optional[Union[float, Tuple[float, float]]] = None,
                          initial_params: Optional[EllipseParams] = None,
                          verbose=0):
    """
    Main function to fit ellipses to galaxy image.
    
    Parameters:
        gal_image: 2D numpy array containing the galaxy image
        Isequence: array of intensity values to fit. If None, will be generated automatically
        region_split: index to split inner/outer regions. If None, uses middle point
        save_name: filename to save the plots
        robust: if True, use Cauchy loss for robustness to outliers
        fix_center: if True, fix center position after initial fit
        dx_dy: maximum allowed offset from initial center position (in pixels).
               Can be a single float (same for x and y) or a tuple (dx, dy).
               Applied after the first successful fit. Ignored if fix_center=True.
        initial_params: initial EllipseParams for the first fit. If None,
                        parameters are estimated from the data. Can be created as:
                        EllipseParams(x0=100, y0=100, a=20, b=15, phi=0.5)
        verbose: verbosity level for fitting
        
    Returns:
        Tuple containing:
        (q_inner, q_outer, pa_inner, pa_outer, q_median, pa_median,
         x0_median, y0_median, x0_inner, y0_inner, x0_outer, y0_outer,
         profiles)
    """
    if Isequence is None:
        Isteps = 64
        Imin = 0.01 * np.nanstd(gal_image)
        Imax = 0.99 * np.nanmax(gal_image)
        Isequence = np.geomspace(Imax, Imin, Isteps)
    
    return fit_isophotes(
        image=gal_image,
        intensity_levels=Isequence,
        region_split=region_split,
        plot_results=plot_results,
        plot_profiles=plot_profiles,
        save_name=save_name+'_photo.jpg' if save_name else None,
        save_profiles=save_name+'_profiles.jpg' if save_name else None,
        robust=robust,
        fix_center=fix_center,
        dx_dy=dx_dy,
        initial_params=initial_params,
        verbose=verbose
    )