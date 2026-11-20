"""
Per-pixel spectral fitting used by `specidx_map`, `specidx_map_SY_FF` and
`t_specidx_map_SY_FF` (radio_sed.py).

Standalone on purpose: it imports only numpy and lmfit, never `mlibs`.
Joblib/loky workers unpickle `fit_spectra_chunk` by reference to this
module, so each worker imports this small file instead of re-running the
whole `mlibs` preamble (CASA, JAX, banner, ...).

`linear_function`, `RC_function_SY_FF`, `do_fit_spec_map` and
`do_fit_spec_SY_FF_map` were moved here verbatim from radio_sed.py;
radio_sed.py re-imports them, so `mlibs.<name>` keeps working.
The `fit_*_chunk` functions are the parallel workers.
"""
import numpy as np
import lmfit


def linear_function(x, alpha, b,x0=1.0):
    return b*((x/x0)**alpha) #+ c*x*x


def RC_function_SY_FF(nu, A_sy, A_ff, alpha_nt, nu0):
    # return A_ff * (nu0 ** (-0.1)) * ((nu / nu0) ** (-0.1)) + A_sy * (nu0 ** (alpha_nt)) * ((nu / nu0) ** (alpha_nt))
    return A_ff * ((nu / nu0) ** (-0.1)) + A_sy * ((nu / nu0) ** (alpha_nt))



def do_fit_spec_map(freqs,fluxes,fluxes_err,nu0=1.0,verbose=0,
                    weighting='log',robust_scale=1.0):
    """
    Fit a power law b*(nu/nu0)**alpha to one spectrum (lmfit least_squares,
    Cauchy loss, two passes).

    weighting : 'log' (default) or 'chi2'
        'log': residual (y - model) * 1/(1 + log1p(yerr/median(yerr))); only
            how yerr compares across the bands of this spectrum matters, and
            the Cauchy threshold is 1 in the units of y. Errors are rescaled
            by the scatter about the fit (lmfit scale_covar=True).
        'chi2': residual (y - model)/yerr, i.e. yerr is used as the 1-sigma
            uncertainty of each point. The Cauchy threshold is robust_scale
            (in sigma units), so it does not depend on flux units or source
            brightness, and errors are absolute (scale_covar=False): they
            follow from yerr, not from the scatter about the fit.
    robust_scale : float
        Only for 'chi2': residuals beyond ~robust_scale sigma are
        progressively down-weighted by the Cauchy loss.
    """
    if weighting not in ('log', 'chi2'):
        raise ValueError(f"weighting must be 'log' or 'chi2', not {weighting!r}")

    x = freqs
    y = fluxes
    yerr = fluxes_err
    
    epsilon = 1e-10
    delta = 1.0
    # def min_func(params):
    #     alpha = params['alpha']
    #     b = params['b']
    #     model = linear_function(x, alpha, b,nu0)
    #     relative_error = yerr / (np.abs(y) + epsilon)
    #     weights = np.sqrt(1 / (relative_error + epsilon))
    #     # weightned_residual = (y - model) / yerr
    #     weightned_residual = (y - model) * (weights)
    #     # weightned_residual = (y - model) / (weights * yerr)
    #     # weightned_residual = (y - model) * np.sqrt(relative_error) #bias towards single point with large error
    #     # weightned_residual = (y - model) / ((weights)*yerr)
    #     # weightned_residual = (y - model) * ((weights)*yerr)
    #     # weightned_residual = (y - model)
    #     # weightned_residual = (y - model) / (weights)
    #     # weightned_residual = (y - model) * (weights + yerr)
    #     # weightned_residual = (y - model) * np.sqrt(weights)
    #     # weightned_residual = (y - model) / np.sqrt(weights)
    #     # weightned_residual = (y - model) / (np.sqrt(weights) + yerr)
    #     # weightned_residual = (y - model) * (np.sqrt(weights) + yerr)
    #     # weightned_residual = (y - model) / (np.sqrt(weights + yerr))
    #     # weightned_residual = (y - model) * (np.sqrt(weights + yerr))
    #     # weightned_residual = (y - model) * np.sqrt(weights * yerr)
    #     # weightned_residual = (y - model) / np.sqrt(weights * yerr)
    #     # weightned_residual = (y - model) / (np.sqrt(weights) * yerr)
    #     # weightned_residual = (y - model) * (np.sqrt(weights) * yerr)
    #     # weightned_residual = (y - model) / np.log(abs(y+yerr))
    #     # weightned_residual = (y - model) / np.log(abs(yerr))
    #     # weightned_residual = (y - model) /yerr
    #     # weightned_residual = (y - model)
    #     return weightned_residual.copy()
    
    # def min_func(params):
    #     alpha = params['alpha']
    #     b = params['b']
    #     model = linear_function(x, alpha, b, nu0)
        
    #     # Huber-like weighting
    #     relative_error = yerr / (np.abs(y) + epsilon)
    #     c = np.median(relative_error)  # adaptive threshold
    #     weights = np.where(relative_error <= c,
    #                     1.0,
    #                     np.sqrt(c / relative_error))
        
    #     return (y - model) * weights
    
    def min_func_spec_map(params):
        """
        More robust implementation for uncertainties of y (yerr).
        """
        alpha = params['alpha']
        b = params['b']
        model = linear_function(x, alpha, b, nu0)
        
        # # Log-based weighting
        # log_weights = y / (1.0 + np.log1p(yerr / np.median(yerr)))
        # return (y - model) * np.sqrt(log_weights)
        
        # Log-based weighting
        # log_weights = 1 / yerr
        # return (y - model) * np.sqrt(log_weights)
        # Log-based weighting
        # l_penalty = 0.1
        # pen_weights = yerr * (alpha**2.0 + b**2.0 - 1.0) / len(x)
        # return (y - model) / pen_weights
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights
        # return ((y - model)/y) * (yerr / model)

    def min_func_spec_map_chi2(params):
        """
        Normalised residuals: yerr is the 1-sigma uncertainty of each point.
        """
        model = linear_function(x, params['alpha'], params['b'], nu0)
        return (y - model) / yerr

    
    # def min_func(params):
    #     """
    #     More robust implementation for uncertainties of y (yerr).
    #     """
    #     alpha = params['alpha']
    #     b = params['b']
    #     model = linear_function(x, alpha, b, nu0)
        
    #     # Log-based weighting
    #     log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
    #     # Frequency correction factor
    #     freq_weights = (x / np.max(x))**0.5  # Adjust exponent as needed
    #     # Combine weights
    #     final_weights = log_weights * freq_weights
    #     return (y - model) * final_weights
    
    
    # def min_func(params):
    #     alpha = params['alpha']
    #     b = params['b']
    #     model = linear_function(x, alpha, b, nu0)
        
    #     # Add scaling factors to control relative importance
    #     y_scale = 0.5  # Adjust these values based on your needs
    #     err_scale = 1.0
        
    #     normalized_y = y / np.median(np.abs(y))
    #     normalized_yerr = yerr / np.median(yerr)
        
    #     combined_weights = 1.0 / (1.0 + np.log1p(
    #         err_scale * normalized_yerr * (1 + y_scale * np.abs(normalized_y))
    #     ))
        
    #     return (y - model) * combined_weights
    
    
    # Upper bound of b (flux density at nu0, same units as y). It used to be
    # a fixed 500, which silently broke the fit of pixels brighter than that
    # (b pinned at 500, alpha compensating). Keep 500 for anything up to
    # 5 units bright (unchanged results), and scale with the data above it.
    y_peak = np.nanmax(np.abs(y)) if np.any(np.isfinite(y)) else 0.0
    b_max = max(500.0, 100.0 * y_peak)
    fit_params = lmfit.Parameters()
    fit_params.add("alpha", value=-0.5, min=-3, max=5)
    fit_params.add("b", value=5.0, min=-1, max=b_max)
    
    if weighting == 'log':
        mini = lmfit.Minimizer(min_func_spec_map, fit_params, max_nfev=15000,
                            nan_policy='omit', reduce_fcn='neglogcauchy')
        extra_kws = {}
    else:
        mini = lmfit.Minimizer(min_func_spec_map_chi2, fit_params,
                               max_nfev=15000, nan_policy='omit',
                               scale_covar=False)
        extra_kws = {'f_scale': robust_scale}
    
    result_1 = mini.minimize(method='least_squares',
                           max_nfev=15000, #f_scale = 3.0,
                        #    loss="huber", 
                           loss="cauchy",
                           tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12, 
                           verbose=verbose,
                           **extra_kws
                           )
    
    # mini = lmfit.Minimizer(min_func_spec_map, fit_params, max_nfev=15000,
    #                     nan_policy='omit')

    # result_1 = mini.minimize(method='least_squares',
    #                     loss='soft_l1',  # Alternative to Cauchy
    #                     f_scale=0.1,     # Tune this parameter
    #                     tr_solver='exact',
    #                     ftol=1e-12, xtol=1e-12, gtol=1e-12)
    
    second_run_params = result_1.params
    
    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=15000, #f_scale = 3.0,
                        #    loss="huber", 
                           loss="cauchy", 
                           tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12, 
                           verbose=0,
                           **extra_kws
                        )
    
    # result = mini.minimize(method='least_squares',
    #                        params=second_run_params,
    #                        max_nfev=15000, f_scale = 0.1,
    #                     #    loss="huber", 
    #                        loss="soft_l1", 
    #                        tr_solver="exact",
    #                        ftol=1e-12, xtol=1e-12, gtol=1e-12, 
    #                        verbose=0
    #                     )
    
    return result


def fit_spectra_chunk(x, spectra, spectra_res, nu0, flux_sys_error_frac,
                      weighting='log', robust_scale=3.0):
    """
    Fit `do_fit_spec_map` to each row of a chunk of pixel spectra.

    Parameters
    ----------
    x : array (nfreq,)
        Frequencies in GHz.
    spectra, spectra_res : arrays (npix, nfreq)
        Image spectra and their noise term (Jy/beam) of the pixels in this
        chunk: residual pixel values for weighting='log', noise rms for
        'chi2'. Either way yerr = sqrt(res**2 + (flux_sys_error_frac*S)**2).
    nu0 : float
        Reference frequency (GHz).
    flux_sys_error_frac : float
        Fractional systematic flux error.
    weighting, robust_scale :
        Passed to `do_fit_spec_map`.

    Returns
    -------
    alpha, alpha_err, b : arrays (npix,)
        Best-fit values; a missing stderr is returned as NaN.
    """
    npix = spectra.shape[0]
    alpha = np.full(npix, np.nan)
    alpha_err = np.full(npix, np.nan)
    b = np.full(npix, np.nan)
    for k in range(npix):
        # same expressions as the former per-pixel closure in specidx_map
        y = spectra[k, :] * 1000
        yerr = np.sqrt((1*spectra_res[k, :])**2.0 + (flux_sys_error_frac * spectra[k, :])**2.0) * 1000
        results_fit = do_fit_spec_map(x, y, yerr, nu0, weighting=weighting,
                                      robust_scale=robust_scale)
        alpha[k] = results_fit.params['alpha'].value
        stderr = results_fit.params['alpha'].stderr
        alpha_err[k] = np.nan if stderr is None else stderr
        b[k] = results_fit.params['b'].value
    return alpha, alpha_err, b


def do_fit_spec_SY_FF_map(freqs,fluxes,fluxes_err,nu0=None,
                          fix_alpha_nt=False,
                          verbose=0,
                          weighting='log',robust_scale=3.0,
                          alpha_nt_fixed=-0.85,
                          alpha_nt_bounds=(-2.0, 0.5)):
    """
    Fit S(nu) = A_sy*(nu/nu0)**alpha_nt + A_ff*(nu/nu0)**(-0.1) to one
    spectrum (lmfit least_squares, Cauchy loss, two passes). A_sy and A_ff
    are the synchrotron and free-free flux densities at nu0, in the units of
    `fluxes`, and are constrained to be >= 0.

    fix_alpha_nt : bool
        Keep alpha_nt fixed at `alpha_nt_fixed`; the model is then linear in
        the two amplitudes.
    weighting : 'log' (default) or 'chi2'
        Same meaning as in `do_fit_spec_map`. 'log': mild log-weights from
        fluxes_err, Cauchy threshold of 1 in flux units, errors rescaled by
        the scatter about the fit. 'chi2': residual (y - model)/fluxes_err,
        Cauchy threshold of `robust_scale` sigma, absolute errors
        (scale_covar=False), and a starting point taken from the data.
    alpha_nt_fixed : float
        Value of alpha_nt when fix_alpha_nt=True.
    alpha_nt_bounds : (min, max)
        Bounds of alpha_nt when it is free. A maximum above -0.1 lets the
        "synchrotron" term become flatter than the free-free one, at which
        point the two components are no longer distinguishable.
    """
    if weighting not in ('log', 'chi2'):
        raise ValueError(f"weighting must be 'log' or 'chi2', not {weighting!r}")
    x = freqs
    y = fluxes
    yerr = fluxes_err
    if nu0 is None:
        nu0 = np.mean(x)

    epsilon = 1e-8
    # def min_func(params):
    #     A_sy = params['A_sy']
    #     A_ff = params['A_ff']
    #     alpha_nt = params['alpha_nt']
    #     model = RC_function_SY_FF(x, A_sy, A_ff, alpha_nt,nu0)
        
    #     relative_error = yerr / (np.abs(y) + epsilon)
    #     weights = 1 / (relative_error + epsilon)
    #     # weightned_residual = (y - model) * np.sqrt(weights)
    #     weightned_residual = (y - model) / (np.sqrt(weights)*yerr)
    #     # weightned_residual = (y - model) / (np.sqrt(weights) + yerr)
    #     # res = (y - RC_function_S2(x, A1, A2, alpha_nt,nu0))/(yerr+1)
    #     # res = (y - model) / (y+np.sqrt((yerr)**2.0+0.1))#okay 
    #     # res = (y - model) / (y+yerr)
    #     # res = (y - model) / (np.log(y+yerr))
    #     # res = data - RC_function_S2(nu, A1l, alpha_nt)
    #     return weightned_residual.copy()
    
    def min_func(params):
        A_sy = params['A_sy']
        A_ff = params['A_ff']
        alpha_nt = params['alpha_nt']
        model = RC_function_SY_FF(x, A_sy, A_ff, alpha_nt,nu0)
        # Huber-like weighting
        # relative_error = yerr / (np.abs(y) + epsilon)
        # c = np.median(relative_error)  # adaptive threshold
        # weights = np.where(relative_error <= c,
        #                 1.0,
        #                 np.sqrt(c / relative_error))
        
        # return (y - model) * weights
        # Log-based weighting
        # log_weights = 1 / (1.0 + np.log1p(yerr / np.median(yerr)))
        # return (y - model) * log_weights
        # weights = 1 / yerr
        # return (y - model) * np.sqrt(weights)
        # Log-based weighting
        # return ((y - model)/y) * (yerr / model)
        # weights = 1 / yerr
        # return (y - model) * weights
        # res = (y - model) / (np.log(y+yerr))
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights
    
    def min_func_chi2(params):
        """
        Normalised residuals: yerr is the 1-sigma uncertainty of each point.
        """
        model = RC_function_SY_FF(x, params['A_sy'], params['A_ff'],
                                  params['alpha_nt'], nu0)
        return (y - model) / yerr

    # Upper bounds of the amplitudes (flux densities at nu0, units of y).
    # They used to be fixed (1000 and 100), which pinned the amplitudes of
    # bright pixels and let the other parameters compensate. Keep those
    # values for faint spectra (unchanged results) and scale with the data
    # above them. The spectrum is extrapolated to nu0 with a steep index so
    # the bound also holds when nu0 lies below the observed band.
    y_peak = np.nanmax(np.abs(y)) if np.any(np.isfinite(y)) else 0.0
    reach = max(1.0, (np.nanmin(x) / nu0) ** 2.0) if np.size(x) else 1.0
    a_sy_max = max(1000.0, 100.0 * y_peak * reach)
    a_ff_max = max(100.0, 100.0 * y_peak * reach)

    if weighting == 'log':
        # NOTE: this start is (max - min)/2000 although y is already in mJy
        # (about 1000 times too low). Kept so 'log' maps are reproducible.
        a_sy_init = abs(np.nanmax(y) - np.nanmin(y))/ (2.0 * 1000)
        a_ff_init = a_sy_init * 0.5
        a_min = 1.0e-6
    else:
        # start from the flux density of the band closest to nu0
        good = np.isfinite(y)
        if np.any(good):
            s0 = abs(y[good][np.argmin(np.abs(np.asarray(x)[good] - nu0))])
        else:
            s0 = 0.0
        s0 = min(max(s0, 1.0e-4), 0.5 * a_ff_max)
        a_sy_init = 0.7 * s0
        a_ff_init = 0.3 * s0
        a_min = 0.0
    fit_params = lmfit.Parameters()
    fit_params.add("A_sy", value=a_sy_init, min=a_min, max=a_sy_max)
    fit_params.add("A_ff", value=a_ff_init, min=a_min, max=a_ff_max)
    # fit_params.add("A1", value=0.5, min=-10, max=500)
    # fit_params.add("A2", value=0.5, min=-10, max=5000)
    if fix_alpha_nt:
        fit_params.add("alpha_nt", value=alpha_nt_fixed,
                       min=alpha_nt_fixed - 0.005, max=alpha_nt_fixed + 0.005,
                       vary=False)
    else:
        a_lo, a_hi = alpha_nt_bounds
        fit_params.add("alpha_nt", value=min(max(-0.7, a_lo + 0.05), a_hi - 0.05),
                       min=a_lo, max=a_hi)
    # fit_params.add("alpha_nt", value=-0.85, min=-0.9, max=-0.8)

    if weighting == 'log':
        mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                               nan_policy='omit', reduce_fcn='neglogcauchy')
        extra_kws = {}
    else:
        mini = lmfit.Minimizer(min_func_chi2, fit_params, max_nfev=15000,
                               nan_policy='omit', scale_covar=False)
        extra_kws = {'f_scale': robust_scale}

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                            #  loss="huber", 
                             loss="cauchy",
                             tr_solver="exact",
                             ftol=1e-12, xtol=1e-12, gtol=1e-12,
                             verbose=verbose,
                             **extra_kws
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                        #    loss="huber", 
                           loss="cauchy",
                           tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose,
                           **extra_kws
                           )
    
    return result


def fit_spectra_SY_FF_chunk(x, spectra, spectra_res, nu0, flux_sys_error_frac,
                            fix_alpha_nt):
    """
    Fit `do_fit_spec_SY_FF_map` to each row of a chunk of pixel spectra.

    Same inputs as `fit_spectra_chunk`, plus `fix_alpha_nt`.

    Returns
    -------
    list of tuples, one per pixel:
        (alpha_nt, alpha_nt_err, A_sy, A_sy_err, A_ff, A_ff_err)
    Values are returned exactly as lmfit gives them (a stderr can be None).
    """
    out = []
    for k in range(spectra.shape[0]):
        # same expressions as the former per-pixel functions
        y = spectra[k, :] * 1000
        yerr = np.sqrt((1*spectra_res[k, :])**2.0 + (flux_sys_error_frac * spectra[k, :])**2.0) * 1000
        results_fit = do_fit_spec_SY_FF_map(x, y, yerr, nu0, fix_alpha_nt)
        p = results_fit.params
        out.append((p['alpha_nt'].value, p['alpha_nt'].stderr,
                    p['A_sy'].value, p['A_sy'].stderr,
                    p['A_ff'].value, p['A_ff'].stderr))
    return out


SY_FF_KEYS = ('alpha_nt', 'alpha_nt_err', 'A_sy', 'A_sy_err', 'A_ff',
              'A_ff_err', 'cov_sy_ff', 'cov_sy_alpha', 'cov_ff_alpha',
              'chi2_red', 'n_bands')


def spectrum_yerr(spectrum, spectrum_res, flux_sys_error_frac):
    """
    (y, yerr) in mJy of one spectrum in Jy: the same expressions the chunk
    workers use, yerr = sqrt(res**2 + (flux_sys_error_frac*S)**2).
    """
    y = spectrum * 1000
    yerr = np.sqrt((1*spectrum_res)**2.0 + (flux_sys_error_frac * spectrum)**2.0) * 1000
    return y, yerr


def sy_ff_covariance(x, yerr, A_sy, alpha_nt, nu0, fix_alpha_nt):
    """
    Covariance of (A_sy, A_ff, alpha_nt) at the best fit, C = (J^T W J)^-1,
    with J the derivatives of the model with respect to the varying
    parameters and W = 1/yerr**2 (bands with a non-finite or non-positive
    yerr are left out).

    It is evaluated in the parameters themselves, so it stays defined when
    an amplitude sits on its lower bound (0), where lmfit returns no
    uncertainties; the A_ff error is then the 1-sigma sensitivity to a
    free-free component, i.e. the scale of an upper limit.

    Returns a (3, 3) array ordered (A_sy, A_ff, alpha_nt); the alpha_nt row
    and column are 0 when it is fixed, and everything is NaN when the
    problem is singular (too few bands).
    """
    cov = np.full((3, 3), np.nan)
    good = np.isfinite(yerr) & (yerr > 0)
    nvar = 2 if fix_alpha_nt else 3
    if good.sum() < nvar:
        return cov
    r = np.asarray(x, dtype=float)[good] / nu0
    w = 1.0 / yerr[good]
    cols = [r ** alpha_nt, r ** (-0.1)]
    if not fix_alpha_nt:
        cols.append(A_sy * r ** alpha_nt * np.log(r))
    J = np.column_stack(cols) * w[:, None]
    try:
        c = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        return cov
    if not np.all(np.isfinite(c)) or np.any(np.diag(c) < 0):
        return cov
    cov[:] = 0.0
    cov[:nvar, :nvar] = c
    return cov


def fit_spectra_SY_FF_cov_chunk(x, spectra, spectra_res, nu0,
                                flux_sys_error_frac, fix_alpha_nt,
                                weighting='chi2', robust_scale=3.0,
                                alpha_nt_fixed=-0.85,
                                alpha_nt_bounds=(-2.0, 0.5)):
    """
    Fit `do_fit_spec_SY_FF_map` to each row of a chunk of spectra and return
    the parameters with their covariance, as plain arrays.

    Same inputs as `fit_spectra_chunk` (spectra in Jy/beam; they are fitted
    in mJy/beam), plus the options of `do_fit_spec_SY_FF_map`.

    Uncertainties
    -------------
    'chi2': `sy_ff_covariance` at the best fit (absolute, defined at the
        bounds as well).
    'log' : what lmfit reports (rescaled by the scatter about the fit); NaN
        where lmfit could not estimate them.

    Returns
    -------
    array (npix, len(SY_FF_KEYS)), columns in the order of `SY_FF_KEYS`:
        alpha_nt, alpha_nt_err, A_sy, A_sy_err, A_ff, A_ff_err [mJy/beam],
        cov_sy_ff, cov_sy_alpha, cov_ff_alpha, chi2_red, n_bands.
        chi2_red is sum(((y - model)/yerr)**2) / (n_bands - n_free) (NaN
        without degrees of freedom). Spectra with fewer finite bands than
        free parameters, or whose fit failed, are NaN.
    """
    x = np.asarray(x, dtype=float)
    npix = spectra.shape[0]
    out = np.full((npix, len(SY_FF_KEYS)), np.nan)
    nvar = 2 if fix_alpha_nt else 3
    for k in range(npix):
        y, yerr = spectrum_yerr(spectra[k, :], spectra_res[k, :],
                                flux_sys_error_frac)
        if weighting == 'chi2':
            # points without a usable error carry no information
            y = np.where(np.isfinite(yerr) & (yerr > 0), y, np.nan)
        good = np.isfinite(y) & np.isfinite(yerr)
        n_bands = int(good.sum())
        out[k, 10] = n_bands
        if n_bands < nvar:
            continue
        try:
            results_fit = do_fit_spec_SY_FF_map(
                x, y, yerr, nu0, fix_alpha_nt, weighting=weighting,
                robust_scale=robust_scale, alpha_nt_fixed=alpha_nt_fixed,
                alpha_nt_bounds=alpha_nt_bounds)
        except Exception:
            continue
        p = results_fit.params
        alpha_nt, A_sy, A_ff = p['alpha_nt'].value, p['A_sy'].value, p['A_ff'].value
        if weighting == 'chi2':
            cov = sy_ff_covariance(x, np.where(good, yerr, np.nan), A_sy,
                                   alpha_nt, nu0, fix_alpha_nt)
        else:
            cov = np.full((3, 3), np.nan)
            names = ['A_sy', 'A_ff', 'alpha_nt']
            for a, na in enumerate(names):
                stderr = p[na].stderr
                if not p[na].vary:
                    cov[a, :] = 0.0
                    cov[:, a] = 0.0
                elif stderr is not None:
                    cov[a, a] = stderr ** 2
            covar = getattr(results_fit, 'covar', None)
            if covar is not None:
                var_names = list(results_fit.var_names)
                for a, na in enumerate(names):
                    for b, nb in enumerate(names):
                        if na in var_names and nb in var_names:
                            cov[a, b] = covar[var_names.index(na),
                                              var_names.index(nb)]
        model = RC_function_SY_FF(x[good], A_sy, A_ff, alpha_nt, nu0)
        dof = n_bands - nvar
        with np.errstate(divide='ignore', invalid='ignore'):
            chi2_red = (np.sum(((y[good] - model) / yerr[good]) ** 2) / dof
                        if dof > 0 else np.nan)
            out[k, :10] = (alpha_nt, np.sqrt(cov[2, 2]), A_sy, np.sqrt(cov[0, 0]),
                           A_ff, np.sqrt(cov[1, 1]), cov[0, 1], cov[0, 2],
                           cov[1, 2], chi2_red)
    return out


# ---------------------------------------------------------------------------
# Synchrotron + free-free with free-free absorption (FFA)
# ---------------------------------------------------------------------------
from scipy import optimize as _optimize

FFA_GEOMETRIES = ('thesis', 'screen', 'mixed')
ALPHA_FF = -0.1

# tau_ff = FF_TAU_COEFF (Te/1e4 K)^-1.35 (nu/GHz)^-2.1 (EM / pc cm^-6)
# (Mezger & Henderson 1967)
FF_TAU_COEFF = 3.28e-7


def ffa_tau(nu, nu_t):
    """
    Free-free optical depth tau(nu) = (nu/nu_t)**-2.1, so tau(nu_t) = 1.
    nu_t <= 0 means no absorber (tau = 0). nu_t is NOT the frequency of the
    spectral peak: behind a screen the peak lies at about 1.5 nu_t.
    """
    nu, nu_t = np.broadcast_arrays(np.asarray(nu, dtype=float),
                                   np.asarray(nu_t, dtype=float))
    tau = np.zeros(nu.shape)
    on = nu_t > 0
    with np.errstate(over='ignore', divide='ignore'):
        tau[on] = (nu[on] / nu_t[on]) ** (-2.1)
    return tau


def ffa_factors(nu, nu_t, geometry='thesis', f_cov=1.0):
    """
    Attenuation factors (g_sy, g_ff) of the synchrotron and free-free terms.

    With tau = (nu/nu_t)**-2.1:
        screen(tau) = 1 - f_cov*(1 - exp(-tau))    foreground screen that
                      covers a fraction f_cov of the emitter; exp(-tau) for
                      a homogeneous screen, f_cov = 1;
        mixed(tau)  = (1 - exp(-tau))/tau          emitter and absorber
                      uniformly mixed in the same volume.

    geometry
        'thesis': g_sy = screen, g_ff = mixed. The ionised gas absorbs its
                  own free-free emission and sits in front of the
                  synchrotron emitter. This is the form of
                  RC_function_SY_FF_ffa_v2 (radio_sed.py).
        'screen': g_sy = g_ff = screen.
        'mixed' : g_sy = g_ff = mixed (f_cov is not used).
    """
    if geometry not in FFA_GEOMETRIES:
        raise ValueError(f"geometry must be one of {FFA_GEOMETRIES}, "
                         f"not {geometry!r}")
    tau = ffa_tau(nu, nu_t)
    one_minus_exp = -np.expm1(-tau)
    screen = 1.0 - f_cov * one_minus_exp
    with np.errstate(divide='ignore', invalid='ignore'):
        mixed = np.where(tau > 1e-8, one_minus_exp / tau, 1.0 - 0.5 * tau)
    if geometry == 'thesis':
        return screen, mixed
    if geometry == 'screen':
        return screen, screen
    return mixed, mixed


def RC_function_SY_FF_FFA(nu, A_sy, A_ff, alpha_nt, nu_t, nu0=10.0,
                          geometry='thesis', f_cov=1.0, components=False):
    """
    Synchrotron + free-free emission with free-free absorption:

        S(nu) = A_sy (nu/nu0)**alpha_nt g_sy(tau) + A_ff (nu/nu0)**-0.1 g_ff(tau)

    A_sy and A_ff are the UNABSORBED flux densities at nu0; tau and the
    factors g are those of `ffa_factors`. With components=True the two
    absorbed terms are returned separately.
    """
    g_sy, g_ff = ffa_factors(nu, nu_t, geometry, f_cov)
    r = np.asarray(nu, dtype=float) / nu0
    sy = A_sy * r ** alpha_nt * g_sy
    ff = A_ff * r ** ALPHA_FF * g_ff
    if components:
        return sy, ff
    return sy + ff


FFA_FIT_DEFAULTS = dict(
    geometry='thesis',        # see ffa_factors
    f_cov=1.0,                # covering factor of the screen: a number
                              # (fixed), or 'fit' (geometry 'thesis' only)
    fix_alpha_nt=False,
    alpha_nt_fixed=-0.85,
    alpha_nt_bounds=(-2.0, -0.2),
    alpha_nt_prior=None,      # (mean, sigma) Gaussian prior on alpha_nt
    nu_t_bounds=None,         # GHz; None: (nu_min/10, min(nu_max, 20))
    n_alpha=41,               # grid nodes in alpha_nt
    n_nu_t=81,                # grid nodes in log10(nu_t)
    band_groups=None,         # group id per image: bands sharing a flux scale
    flux_ind_error_frac=None, # independent fractional error per image, on top
                              # of the shared one; None: 0.03 with band_groups,
                              # 0 without
    n_iter=4,                 # reweighted amplitude solves per trial
)

FFA_COV_KEYS = ('cov_sy_sy', 'cov_sy_ff', 'cov_sy_alpha', 'cov_sy_lognut',
                'cov_ff_ff', 'cov_ff_alpha', 'cov_ff_lognut',
                'cov_alpha_alpha', 'cov_alpha_lognut', 'cov_lognut_lognut')
FFA_KEYS = (('A_sy', 'A_ff', 'alpha_nt', 'log_nu_t') + FFA_COV_KEYS
            + ('chi2', 'chi2_null', 'alpha_nt_null', 'n_bands',
               'log_nu_t_lo', 'log_nu_t_hi', 'log_nu_t_up95', 'at_bound',
               'delta_bic', 'f_cov', 'f_cov_err'))


class _FFAGrid:
    """
    Grid of the non-linear parameters (alpha_nt, nu_t) and the two basis
    spectra on it. Built once per chunk: it does not depend on the data.
    Node 0 of the nu_t axis is "no absorber" (tau = 0), the model against
    which the absorption is tested.
    """

    def __init__(self, x, nu0, opts):
        self.x = np.asarray(x, dtype=float)
        self.nu0 = float(nu0)
        self.opts = opts
        self.fix_alpha = bool(opts['fix_alpha_nt'])
        self.fit_fcov = isinstance(opts['f_cov'], str)
        if self.fit_fcov:
            if opts['f_cov'] != 'fit':
                raise ValueError("f_cov must be a number or 'fit', not "
                                 f"{opts['f_cov']!r}")
            if opts['geometry'] != 'thesis':
                raise ValueError("f_cov='fit' is only available for "
                                 "geometry='thesis'")
        if self.fix_alpha:
            self.alpha = np.array([float(opts['alpha_nt_fixed'])])
        else:
            a_lo, a_hi = opts['alpha_nt_bounds']
            self.alpha = np.linspace(a_lo, a_hi, int(opts['n_alpha']))
        self.alpha_bounds = (self.alpha[0], self.alpha[-1])
        nut_bounds = opts['nu_t_bounds']
        if nut_bounds is None:
            nut_bounds = (self.x.min() / 10.0, min(self.x.max(), 20.0))
        self.lognut_bounds = (np.log10(nut_bounds[0]), np.log10(nut_bounds[1]))
        self.lognut = np.linspace(self.lognut_bounds[0], self.lognut_bounds[1],
                                  int(opts['n_nu_t']))
        self.d_lognut = self.lognut[1] - self.lognut[0]
        # nodes: nu_t axis (0 = no absorber) x alpha axis, flattened
        nut_nodes = np.concatenate([[0.0], 10.0 ** self.lognut])
        NT, AL = np.meshgrid(nut_nodes, self.alpha, indexing='ij')
        self.shape = NT.shape
        self.B = self.basis(AL.ravel(), NT.ravel())
        self.AL = AL.ravel()
        prior = opts['alpha_nt_prior']
        self.prior = None if (prior is None or self.fix_alpha) else \
            (float(prior[0]), float(prior[1]))
        groups = opts['band_groups']
        if groups is None:
            groups = np.arange(self.x.size)
        _, inv = np.unique(np.asarray(groups), return_inverse=True)
        self.G = np.zeros((inv.max() + 1, self.x.size))
        self.G[inv, np.arange(self.x.size)] = 1.0
        self.group_of = inv

    def basis(self, alpha, nu_t):
        """
        Basis spectra, each (n, nfreq), for arrays of alpha and nu_t (n,).
        The model is their non-negative combination.

        Fixed f_cov: (synchrotron, free-free), both with their attenuation.
        f_cov='fit': (covered synchrotron, uncovered synchrotron,
        free-free). The screen term 1 - f_cov (1 - e^-tau) splits into
        f_cov e^-tau + (1 - f_cov), so the covered and the uncovered part
        are two amplitudes and f_cov = A_cov / (A_cov + A_unc); no extra
        non-linear parameter is needed.
        """
        alpha = np.atleast_1d(np.asarray(alpha, dtype=float))
        nu_t = np.atleast_1d(np.asarray(nu_t, dtype=float))
        r = self.x[None, :] / self.nu0
        if self.fit_fcov:
            g_sy, g_ff = ffa_factors(self.x[None, :], nu_t[:, None], 'thesis', 1.0)
            rp = r ** alpha[:, None]
            return rp * g_sy, rp + 0.0 * g_sy, r ** ALPHA_FF * g_ff
        g_sy, g_ff = ffa_factors(self.x[None, :], nu_t[:, None],
                                 self.opts['geometry'], self.opts['f_cov'])
        return r ** alpha[:, None] * g_sy, r ** ALPHA_FF * g_ff

    def prior_chi2(self, alpha):
        if self.prior is None:
            return 0.0
        return ((alpha - self.prior[0]) / self.prior[1]) ** 2


class _Metric:
    """
    Inverse-covariance inner product <u, v> = u^T C^-1 v for

        C = diag(sig2) + sum_g f_g**2 m_g m_g^T,

    i.e. independent noise plus one multiplicative flux-scale error per
    group of bands (m: the expected flux density, f_g: fractional error of
    group g; m_g is m restricted to the group). With one image per group
    this is the usual sig2 + (f*m)**2 on the diagonal. C is block diagonal
    with rank-one blocks, so C^-1 follows from Sherman-Morrison.
    """

    def __init__(self, sig2, m, fsys_group, G, good):
        self.winv = np.where(good, 1.0 / np.where(good, sig2, 1.0), 0.0)
        self.mw = np.where(good, m, 0.0) * self.winv
        self.G = G
        s = G @ (np.where(good, m, 0.0) * self.mw)
        self.c = fsys_group ** 2 / (1.0 + fsys_group ** 2 * s)
        # ln det C (matrix determinant lemma), over the bands with data
        self.logdet = float(np.sum(np.log(sig2[good]))
                            + np.sum(np.log1p(fsys_group ** 2 * s)))

    def dot(self, U, V):
        """<U, V>; U and V are (n, nfreq) or (nfreq,)."""
        d = (U * V) @ self.winv
        return d - (((U * self.mw) @ self.G.T) * ((V * self.mw) @ self.G.T)) @ self.c


def _nnls2(a11, a12, a22, b1, b2, yy):
    """
    Minimise yy - 2(A1 b1 + A2 b2) + A1^2 a11 + 2 A1 A2 a12 + A2^2 a22 over
    A1, A2 >= 0, in closed form (vectorised). Returns A1, A2, chi2.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        det = a11 * a22 - a12 ** 2
        ok = det > 1e-12 * a11 * a22
        A1 = np.where(ok, (b1 * a22 - b2 * a12) / det, -1.0)
        A2 = np.where(ok, (b2 * a11 - b1 * a12) / det, -1.0)
        free = ok & (A1 >= 0) & (A2 >= 0)
        A1s = np.where(a11 > 0, np.maximum(b1 / a11, 0.0), 0.0)
        A2s = np.where(a22 > 0, np.maximum(b2 / a22, 0.0), 0.0)
        chi_1 = yy - 2.0 * A1s * b1 + A1s ** 2 * a11
        chi_2 = yy - 2.0 * A2s * b2 + A2s ** 2 * a22
        only1 = chi_1 <= chi_2
        A1o = np.where(free, A1, np.where(only1, A1s, 0.0))
        A2o = np.where(free, A2, np.where(only1, 0.0, A2s))
        chi = np.where(free, yy - A1 * b1 - A2 * b2,
                       np.where(only1, chi_1, chi_2))
    return A1o, A2o, chi


def _nnls3(gram, b, yy):
    """
    Minimise yy - 2 A.b + A^T G A over A >= 0 for three amplitudes
    (vectorised). gram = (g11, g12, g13, g22, g23, g33), b = (b1, b2, b3).

    The minimum is the unconstrained solution of one of the seven non-empty
    subsets of amplitudes (the others at zero) that has no negative entry;
    all of them are tried and the lowest chi2 is taken. Returns the three
    amplitudes and chi2.
    """
    g11, g12, g13, g22, g23, g33 = gram
    n = np.shape(yy)[0]
    diag = np.stack([g11, g22, g33], axis=1)
    valid = np.isfinite(diag) & (diag > 0)
    scale = np.sqrt(np.where(valid, diag, 1.0))
    G = np.empty((n, 3, 3))
    G[:, 0, 0], G[:, 1, 1], G[:, 2, 2] = g11, g22, g33
    G[:, 0, 1] = G[:, 1, 0] = g12
    G[:, 0, 2] = G[:, 2, 0] = g13
    G[:, 1, 2] = G[:, 2, 1] = g23
    Gn = G / (scale[:, :, None] * scale[:, None, :])
    bn = np.stack(b, axis=1) / scale
    Gn = np.where(np.isfinite(Gn), Gn, 0.0)
    bn = np.where(np.isfinite(bn), bn, 0.0)
    best_chi = np.array(yy, dtype=float, copy=True)        # all amplitudes zero
    best_A = np.zeros((n, 3))
    for idx in ((0,), (1,), (2,), (0, 1), (0, 2), (1, 2), (0, 1, 2)):
        idx = list(idx)
        Gs = Gn[:, idx][:, :, idx]
        bs = bn[:, idx]
        ok = np.all(valid[:, idx], axis=1) & (np.linalg.det(Gs) > 1e-10)
        Gs = np.where(ok[:, None, None], Gs, np.eye(len(idx)))
        sol = np.linalg.solve(Gs, bs[:, :, None])[:, :, 0]
        chi = yy - np.sum(sol * bs, axis=1)
        better = ok & np.all(sol >= 0, axis=1) & (chi < best_chi)
        cand = np.zeros((n, 3))
        cand[:, idx] = sol
        best_A = np.where(better[:, None], cand, best_A)
        best_chi = np.where(better, chi, best_chi)
    best_A = best_A / scale
    return best_A[:, 0], best_A[:, 1], best_A[:, 2], best_chi


def _ffa_solve(grid, metric, y, Bsy, Bff, alpha):
    """Best non-negative amplitudes and chi2 for each row of the basis."""
    A1, A2, chi = _nnls2(metric.dot(Bsy, Bsy), metric.dot(Bsy, Bff),
                         metric.dot(Bff, Bff), metric.dot(Bsy, y),
                         metric.dot(Bff, y), metric.dot(y, y))
    return A1, A2, chi + grid.prior_chi2(alpha)


def fit_SY_FF_FFA_spectrum(x, y, sigma, nu0=10.0, flux_sys_error_frac=0.1,
                           grid=None, return_profile=False, **options):
    """
    Fit synchrotron + free-free emission with free-free absorption
    (`RC_function_SY_FF_FFA`) to one spectrum.

    Method (separable least squares)
    --------------------------------
    For every trial of the non-linear parameters (alpha_nt, nu_t) the two
    amplitudes are solved with A_sy, A_ff >= 0, so only one or two
    parameters are searched: first on a grid that spans the whole allowed
    range (the global step, immune to the starting point and to local
    minima), then with a simplex refinement from the best nodes. The grid
    includes "no absorber", so the unabsorbed synchrotron + free-free model
    is fitted to the same data in the same way and the two can be compared.

    Errors and weights
    ------------------
    C = diag(sigma**2 + (f_ind*m)**2) + flux-scale term (see `_Metric`),
    where m = sqrt(|data| * model) is the geometric mean of the observed
    and the predicted flux density. This choice is deliberate:
      - m = data bounds the cost of a model that falls far BELOW a point at
        1/f sigma, so a fit can abandon the low-frequency points cheaply;
      - m = model bounds the cost of a model far ABOVE a point in the same
        way, so an unabsorbed spectrum can overshoot them cheaply;
      - the geometric mean penalises both without bound and symmetrically
        in the ratio data/model: the residual becomes
        2 sinh(ln(data/model)/2)/f, which is ln(data/model)/f (a log-normal
        flux-scale error) for small deviations.
    With `band_groups`, images of one group share a single flux-scale error
    (`flux_sys_error_frac`; the sub-bands of one observation) instead of
    each image having its own, and `flux_ind_error_frac` is the error that
    remains independent from image to image (sub-band imaging and
    calibration differences); without it the scatter between sub-bands
    would be judged against the thermal noise alone.

    C depends on the model, so the quantity minimised is the likelihood
    itself, -2 ln L = chi2 + ln det C, with C built from each trial model's
    own prediction ("self-consistent"). For a trial (alpha_nt, nu_t) the
    amplitudes and their error model are brought into agreement by `n_iter`
    reweighted solves. Fitting with weights taken from some other model is
    not equivalent: such a fit can abandon the low-frequency points, or
    overshoot them, at little cost. The reported chi2, the values for the
    unabsorbed model, delta_bic and the limits on nu_t all come from this
    likelihood. (Checked against a brute-force global optimisation of the
    same likelihood on real spectra.)

    Parameters
    ----------
    x : array (nfreq,)        frequencies, GHz
    y : array (nfreq,)        flux densities (any unit); NaN = no data
    sigma : array (nfreq,)    1-sigma noise of each point, same unit
    nu0 : float               reference frequency, GHz
    flux_sys_error_frac : float or array (nfreq,)
        Fractional flux-scale error, per image (within a group the value
        of its first image is used).
    grid : _FFAGrid, optional (built from `options` if not given)
    return_profile : bool     also return the chi2 profile over nu_t
    options : see FFA_FIT_DEFAULTS

    Returns
    -------
    dict with
        A_sy, A_ff : unabsorbed amplitudes at nu0;
        alpha_nt, nu_t (GHz), log_nu_t;
        cov : (4, 4) covariance of (A_sy, A_ff, alpha_nt, log10 nu_t); rows
            of parameters that are fixed or not constrained are 0;
        chi2, n_bands, n_free;
        chi2_null, alpha_nt_null, A_sy_null, A_ff_null : best unabsorbed fit;
        delta_bic = (-2 ln L)_null - (-2 ln L)_fit - k ln(n_bands): evidence
            for the absorption, k = 1 extra parameter (2 with f_cov='fit');
            > 6 is strong, > 10 very strong;
        log_nu_t_lo, log_nu_t_hi : 68% interval of log10 nu_t from the
            likelihood profile (delta(-2 ln L) = 1, the other parameters
            re-optimised);
            lo = -inf if the data are compatible with no absorption, hi =
            +inf if not excluded up to the upper bound;
        log_nu_t_up95 : 95% one-sided upper limit (delta = 2.71);
        at_bound : 1 A_sy = 0, 2 A_ff = 0, 4 alpha_nt on a bound, 8 nu_t at
            its lower bound (no absorption needed), 16 nu_t at its upper
            bound, 128 fitted f_cov = 1 (values add up);
        f_cov, f_cov_err : covering factor (the fixed input with error 0,
            or the fitted value with f_cov='fit');
        neg2lnL, neg2lnL_null : -2 ln L of the fit and of the unabsorbed fit;
        [profile_log_nu_t, profile_neg2lnL if return_profile: -2 ln L on
            the nu_t grid, minimised over the alpha_nt grid, for plots;
            subtract neg2lnL to get the delta that defines the limits].
    """
    opts = dict(FFA_FIT_DEFAULTS)
    unknown = set(options) - set(opts)
    if unknown:
        raise TypeError(f"unknown FFA fit options: {sorted(unknown)}")
    opts.update(options)
    x = np.asarray(x, dtype=float)
    if grid is None:
        grid = _FFAGrid(x, nu0, opts)
    y = np.asarray(y, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    good = np.isfinite(y) & np.isfinite(sigma) & (sigma > 0)
    n_bands = int(good.sum())
    fit_fcov = grid.fit_fcov
    n_free = (3 if grid.fix_alpha else 4) + (1 if fit_fcov else 0)
    out = dict(A_sy=np.nan, A_ff=np.nan, alpha_nt=np.nan, nu_t=np.nan,
               f_cov=np.nan, f_cov_err=np.nan,
               log_nu_t=np.nan, cov=np.full((4, 4), np.nan), chi2=np.nan,
               n_bands=n_bands, n_free=n_free, chi2_null=np.nan,
               alpha_nt_null=np.nan, A_sy_null=np.nan, A_ff_null=np.nan,
               delta_bic=np.nan, log_nu_t_lo=np.nan, log_nu_t_hi=np.nan,
               log_nu_t_up95=np.nan, at_bound=np.nan, neg2lnL=np.nan,
               neg2lnL_null=np.nan)
    if n_bands < n_free:
        return out
    yz = np.where(good, y, 0.0)
    goodf = good.astype(float)
    G = grid.G
    fsys = np.broadcast_to(np.asarray(flux_sys_error_frac, dtype=float),
                           x.shape)
    # one fractional error per group: that of its first image
    first = np.array([np.argmax(grid.group_of == g)
                      for g in range(G.shape[0])])
    fsys_group = fsys[first]
    f2_group = fsys_group ** 2
    sig2 = np.where(good, sigma, 1.0) ** 2
    f_ind = opts['flux_ind_error_frac']
    if f_ind is None:
        f_ind = 0.03 if opts['band_groups'] is not None else 0.0
    f_ind = np.broadcast_to(np.asarray(f_ind, dtype=float), x.shape)
    # Flux density that the fractional errors refer to: the geometric mean
    # of data and model (the data floored at the noise).
    y_abs = np.where(good, np.maximum(np.abs(yz), np.sqrt(sig2)), 0.0)
    n_irls = max(1, int(opts['n_iter']))
    n_nut, n_al = grid.shape
    lo_a, hi_a = grid.alpha_bounds
    lo_t, hi_t = grid.lognut_bounds
    da_grid = (grid.alpha[1] - grid.alpha[0]) if n_al > 1 else 0.0

    def error_model(model):
        """Pieces of C^-1 and ln det C for each row of model (n, nfreq)."""
        m = np.sqrt(np.abs(model) * y_abs)
        d = sig2 + (f_ind * m) ** 2
        winv = goodf / d
        mw = m * winv
        s_grp = (m * mw) @ G.T
        c = f2_group / (1.0 + f2_group * s_grp)
        return d, winv, mw, c, s_grp

    def solve(B, alpha):
        """
        For each row of the basis spectra B (each (n, nfreq)): amplitudes
        >= 0 whose error model agrees with their own prediction, -2 ln L
        and chi2. Returns (amplitudes, -2 ln L, chi2).
        """
        Y = yz[None, :]
        nb = len(B)
        model = np.broadcast_to(y_abs, B[0].shape)   # first pass: the data
        for _ in range(n_irls):
            _, winv, mw, c, _ = error_model(model)

            def dot(U, V):
                return ((U * V * winv).sum(axis=1)
                        - (((U * mw) @ G.T) * ((V * mw) @ G.T) * c).sum(axis=1))
            if nb == 2:
                A = _nnls2(dot(B[0], B[0]), dot(B[0], B[1]), dot(B[1], B[1]),
                           dot(B[0], Y), dot(B[1], Y), dot(Y, Y))[:2]
            else:
                A = _nnls3((dot(B[0], B[0]), dot(B[0], B[1]), dot(B[0], B[2]),
                            dot(B[1], B[1]), dot(B[1], B[2]), dot(B[2], B[2])),
                           (dot(B[0], Y), dot(B[1], Y), dot(B[2], Y)),
                           dot(Y, Y))[:3]
            model = sum(A[i][:, None] * B[i] for i in range(nb))
        d, winv, mw, c, s_grp = error_model(model)
        r = Y - model
        chi2_ = ((r * r * winv).sum(axis=1)
                 - ((((r * mw) @ G.T) ** 2) * c).sum(axis=1))
        logdet = (np.log(d) * goodf).sum(axis=1) + np.log1p(f2_group * s_grp).sum(axis=1)
        return A, chi2_ + logdet + grid.prior_chi2(alpha), chi2_

    def at(alpha_, nu_t_):
        """(amplitudes, -2 ln L, chi2) of one trial."""
        A, cost_, chi2_ = solve(grid.basis([alpha_], [nu_t_]), np.array([alpha_]))
        return tuple(float(a[0]) for a in A), float(cost_[0]), float(chi2_[0])

    def to_phys(amps):
        """(A_sy, A_ff, f_cov) from the amplitudes of the basis."""
        if fit_fcov:
            A_sy_ = amps[0] + amps[1]
            return A_sy_, amps[2], (amps[0] / A_sy_ if A_sy_ > 0 else 1.0)
        return amps[0], amps[1], float(opts['f_cov'])

    def refine(k0, j0):
        """Simplex refinement of -2 ln L from grid node (k0 >= 1, j0)."""
        if grid.fix_alpha:
            lt0 = grid.lognut[k0 - 1]
            a = max(lo_t, lt0 - 2 * grid.d_lognut)
            b = min(hi_t, lt0 + 2 * grid.d_lognut)
            res = _optimize.minimize_scalar(
                lambda lt: at(grid.alpha[0], 10.0 ** lt)[1], bounds=(a, b),
                method='bounded', options=dict(xatol=1e-6))
            lt = float(res.x)
            if at(grid.alpha[0], 10.0 ** lt0)[1] < res.fun:
                lt = lt0
            return grid.alpha[0], lt

        def f2(p):
            return at(p[0], 10.0 ** p[1])[1]
        best_p = np.array([grid.alpha[j0], grid.lognut[k0 - 1]])
        best_f = f2(best_p)
        for _ in range(2):          # one restart from the first solution
            simplex = np.array([best_p,
                                [np.clip(best_p[0] + 0.5 * da_grid, lo_a, hi_a), best_p[1]],
                                [best_p[0], np.clip(best_p[1] + 0.5 * grid.d_lognut,
                                                    lo_t, hi_t)]])
            if simplex[1, 0] == best_p[0]:
                simplex[1, 0] = best_p[0] - 0.5 * da_grid
            if simplex[2, 1] == best_p[1]:
                simplex[2, 1] = best_p[1] - 0.5 * grid.d_lognut
            res = _optimize.minimize(
                f2, best_p, method='Nelder-Mead',
                bounds=[(lo_a, hi_a), (lo_t, hi_t)],
                options=dict(xatol=1e-5, fatol=1e-8, maxiter=400,
                             initial_simplex=simplex))
            if res.fun < best_f:
                best_p, best_f = np.asarray(res.x, dtype=float), float(res.fun)
            else:
                break
        return float(best_p[0]), float(best_p[1])

    # --- global step: -2 ln L on the whole grid
    cost = solve(grid.B, grid.AL)[1].reshape(grid.shape)
    cost = np.where(np.isfinite(cost), cost, np.inf)
    prof_grid = cost.min(axis=1)            # [0]: no absorber; [1:]: nu_t grid

    # --- refine from the lowest minima of the nu_t profile (at most two:
    # the likelihood can have separate minima, e.g. one that follows the
    # low-frequency points and one that does not)
    pk = prof_grid[1:]
    is_min = np.r_[True, pk[1:] < pk[:-1]] & np.r_[pk[:-1] <= pk[1:], True]
    starts = np.nonzero(is_min)[0]
    starts = starts[np.argsort(pk[starts])][:2]
    best = None
    for k in starts:
        a_, lt_ = refine(k + 1, int(np.argmin(cost[k + 1])))
        cand = at(a_, 10.0 ** lt_) + (a_, lt_)
        if best is None or cand[1] < best[1]:
            best = cand
    amps, cost_fit, chi2, alpha, lognut = best

    # --- unabsorbed model (node 0 of the nu_t axis), fitted the same way
    j_null = int(np.argmin(cost[0]))
    alpha_null = grid.alpha[j_null]
    if not grid.fix_alpha:
        res = _optimize.minimize_scalar(
            lambda a: at(a, 0.0)[1],
            bounds=(max(lo_a, alpha_null - da_grid), min(hi_a, alpha_null + da_grid)),
            method='bounded', options=dict(xatol=1e-6))
        if res.fun <= cost[0, j_null]:
            alpha_null = float(res.x)
    amps0, cost_null, chi2_null = at(alpha_null, 0.0)
    A_sy0, A_ff0, _ = to_phys(amps0)
    if cost_null < cost_fit:
        # the absorber does not improve the fit at all
        alpha, lognut = alpha_null, lo_t
        amps, cost_fit, chi2 = at(alpha, 10.0 ** lognut)
    # the absorbed model has one more parameter than the unabsorbed one
    # (nu_t), or two when the covering factor is fitted
    delta_bic = cost_null - cost_fit - (2 if fit_fcov else 1) * np.log(n_bands)
    A_sy, A_ff, f_cov = to_phys(amps)
    nb = len(amps)

    def model_at(al_, lt_):
        B_ = grid.basis([al_], [10.0 ** lt_])
        return sum(amps[i] * B_[i][0] for i in range(nb))

    # error model of the solution (for the covariance)
    B_sol = grid.basis([alpha], [10.0 ** lognut])
    model = model_at(alpha, lognut)
    m_sol = np.sqrt(np.abs(model) * y_abs)
    metric = _Metric(sig2 + (f_ind * m_sol) ** 2, m_sol, fsys_group, G, good)

    # --- flags
    amp_tol = 1e-6 * (abs(A_sy) + abs(A_ff)) + 1e-300
    at_bound = (1 * (A_sy <= amp_tol) + 2 * (A_ff <= amp_tol)
                + (0 if grid.fix_alpha else
                   4 * ((alpha <= lo_a + 1e-4) or (alpha >= hi_a - 1e-4)))
                + 8 * (lognut <= lo_t + 1e-3) + 16 * (lognut >= hi_t - 1e-3)
                + (128 * (amps[1] <= amp_tol) if fit_fcov else 0))

    # --- covariance at the best fit, first in the parameters of the basis
    # (amplitudes, alpha_nt, log10 nu_t)
    eps = 1e-5
    J = [B_sol[i][0] for i in range(nb)]
    J.append((model_at(alpha + eps, lognut) - model_at(alpha - eps, lognut)) / (2 * eps))
    J.append((model_at(alpha, lognut + eps) - model_at(alpha, lognut - eps)) / (2 * eps))
    J = np.array(J)
    i_al, i_lt = nb, nb + 1
    cov_b = np.zeros((nb + 2, nb + 2))
    use = list(range(nb)) + ([] if grid.fix_alpha else [i_al]) + [i_lt]
    trials = [use, [u for u in use if u != i_lt]]
    if fit_fcov:
        # with (almost) nothing covered or nothing uncovered the two
        # synchrotron amplitudes cannot both be constrained
        trials += [[u for u in t if u != (1 if amps[1] <= amps[0] else 0)]
                   for t in list(trials)]
    for trial in trials:
        F = np.array([[metric.dot(J[a], J[b]) for b in trial] for a in trial])
        if grid.prior is not None and i_al in trial:
            F[trial.index(i_al), trial.index(i_al)] += 1.0 / grid.prior[1] ** 2
        try:
            c = np.linalg.inv(F)
        except np.linalg.LinAlgError:
            continue
        dg = np.diag(c)
        if (np.all(np.isfinite(c)) and np.all(dg > 0)
                and np.all(dg * np.diag(F) < 1e10)):
            cov_b[np.ix_(trial, trial)] = c
            break
    else:
        cov_b[:] = np.nan
    # ... then in (A_sy, A_ff, alpha_nt, log10 nu_t [, f_cov])
    if fit_fcov:
        T = np.zeros((5, 5))
        T[0, 0] = T[0, 1] = 1.0                 # A_sy = A_cov + A_unc
        T[1, 2] = T[2, 3] = T[3, 4] = 1.0       # A_ff, alpha_nt, log10 nu_t
        if A_sy > 0:                            # f_cov = A_cov / A_sy
            T[4, 0], T[4, 1] = amps[1] / A_sy ** 2, -amps[0] / A_sy ** 2
        cov5 = T @ cov_b @ T.T
        cov = cov5[:4, :4]
        f_cov_err = float(np.sqrt(cov5[4, 4])) if cov5[4, 4] >= 0 else np.nan
        if amps[1] <= amp_tol or amps[0] <= amp_tol:
            # on its limit (1: nothing uncovered; 0: nothing covered) the
            # covering factor has no meaningful first-order error
            f_cov_err = np.nan
    else:
        cov = cov_b
        f_cov_err = 0.0

    # --- likelihood profile over nu_t (the other parameters re-optimised
    # at each nu_t) and the limits it gives. The crossings are searched
    # directly, stepping outwards from the best fit, so they do not depend
    # on the resolution of the grid.
    def prof_at(lt, a_start):
        if grid.fix_alpha:
            return at(grid.alpha[0], 10.0 ** lt)[2], grid.alpha[0]
        for width in (3.0 * da_grid, None):
            a, b = ((lo_a, hi_a) if width is None else
                    (max(lo_a, a_start - width), min(hi_a, a_start + width)))
            res = _optimize.minimize_scalar(lambda al: at(al, 10.0 ** lt)[2],
                                            bounds=(a, b), method='bounded',
                                            options=dict(xatol=1e-4))
            inside = (res.x - a > 1e-3 * (b - a)) and (b - res.x > 1e-3 * (b - a))
            if width is None or inside or (a == lo_a and b == hi_a):
                break
        return float(res.fun), float(res.x)

    sig_lt = np.sqrt(cov[3, 3]) if np.isfinite(cov[3, 3]) and cov[3, 3] > 0 else None
    step0 = max(sig_lt, 0.05 * grid.d_lognut) if sig_lt else grid.d_lognut

    def crossings(side, levels):
        """log10 nu_t at which the delta of -2 ln L first exceeds each level."""
        roots = []
        lt_in, a_in, step = lognut, alpha, step0
        for level in levels:
            root = None
            for _ in range(60):
                lt_out = float(np.clip(lt_in + side * step, lo_t, hi_t))
                if lt_out == lt_in:
                    break                       # reached the edge of the range
                f_out, a_out = prof_at(lt_out, a_in)
                if f_out - cost_fit > level:
                    try:
                        root = _optimize.brentq(
                            lambda t: prof_at(t, a_in)[0] - cost_fit - level,
                            min(lt_in, lt_out), max(lt_in, lt_out), xtol=1e-3)
                    except ValueError:
                        root = lt_out
                    break
                lt_in, a_in = lt_out, a_out
                step *= 2.0
            roots.append(root)
            if root is None:
                roots += [None] * (len(levels) - len(roots))
                break
            lt_in = root
        return roots

    (lo_root,) = crossings(-1, [1.0])
    hi_root, up_root = crossings(+1, [1.0, 2.71])
    if lo_root is None:
        # not excluded down to the lower bound: is "no absorber" allowed?
        lt_lo = -np.inf if cost_null - cost_fit <= 1.0 else lo_t
    else:
        lt_lo = lo_root
    lt_hi = np.inf if hi_root is None else hi_root
    lt_up = np.inf if up_root is None else up_root

    out.update(A_sy=A_sy, A_ff=A_ff, alpha_nt=alpha, nu_t=10.0 ** lognut,
               log_nu_t=lognut, cov=cov, chi2=chi2, chi2_null=chi2_null,
               alpha_nt_null=alpha_null, A_sy_null=A_sy0, A_ff_null=A_ff0,
               delta_bic=delta_bic, log_nu_t_lo=lt_lo, log_nu_t_hi=lt_hi,
               log_nu_t_up95=lt_up, at_bound=float(at_bound),
               neg2lnL=cost_fit, neg2lnL_null=cost_null,
               f_cov=f_cov, f_cov_err=f_cov_err)
    if return_profile:
        out['profile_log_nu_t'] = grid.lognut
        out['profile_neg2lnL'] = prof_grid[1:]
    return out


def fit_spectra_SY_FF_FFA_chunk(x, spectra, spectra_res, nu0,
                                flux_sys_error_frac, options):
    """
    Fit `fit_SY_FF_FFA_spectrum` to each row of a chunk of spectra.

    x in GHz; spectra in Jy/beam (fitted in mJy/beam); spectra_res is the
    1-sigma noise of each point in Jy/beam; options: dict of
    FFA_FIT_DEFAULTS entries.

    Returns an array (npix, len(FFA_KEYS)), columns in the order of
    `FFA_KEYS` (amplitudes in mJy/beam). Failed or under-determined spectra
    are NaN.
    """
    x = np.asarray(x, dtype=float)
    opts = dict(FFA_FIT_DEFAULTS)
    opts.update(options or {})
    grid = _FFAGrid(x, nu0, opts)
    npix = spectra.shape[0]
    out = np.full((npix, len(FFA_KEYS)), np.nan)
    iu = np.triu_indices(4)
    for k in range(npix):
        y = np.asarray(spectra[k, :], dtype=float) * 1000
        sig = np.asarray(spectra_res[k, :], dtype=float) * 1000
        try:
            r = fit_SY_FF_FFA_spectrum(x, y, sig, nu0, flux_sys_error_frac,
                                       grid=grid, **opts)
        except Exception:
            continue
        out[k, :4] = (r['A_sy'], r['A_ff'], r['alpha_nt'], r['log_nu_t'])
        out[k, 4:14] = r['cov'][iu]
        out[k, 14:] = (r['chi2'], r['chi2_null'], r['alpha_nt_null'],
                       r['n_bands'], r['log_nu_t_lo'], r['log_nu_t_hi'],
                       r['log_nu_t_up95'], r['at_bound'], r['delta_bic'],
                       r['f_cov'], r['f_cov_err'])
    return out
