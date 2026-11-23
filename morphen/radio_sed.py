def two_point_fth(nu2,alpha,nu1=33, alpha_nt=-0.85):
    if nu2>=nu1:
        raise ValueError(f"nu2 = {nu2} should be less than nu1 = {nu1}")
    if alpha<alpha_nt:
        print('W: ++==>> alpha is smaller than alpha_nt!')
        print('W: ++==>> Setting alpha_nt = alpha - 0.1')
        alpha_nt = alpha - 0.1
    num = (nu2/nu1)**alpha - (nu2/nu1)**alpha_nt
    den = (nu2/nu1)**(-0.1) - (nu2/nu1)**alpha_nt
    return(num/den)



def RC_function_S2(nu, A1l, A2, alpha_nt, nu0):
    return A1l * ((nu / nu0) ** (-0.1)) + A2 * (nu0 ** (alpha_nt)) * ((nu / nu0) ** (alpha_nt))


# RC_function_SY_FF, linear_function, do_fit_spec_map and do_fit_spec_SY_FF_map
# live in the standalone spectral_fit.py, so that joblib workers import only
# that small module (not mlibs) when fitting pixels in parallel.
import importlib
from joblib import effective_n_jobs
import spectral_fit
importlib.reload(spectral_fit)  # keep reload_libs() picking up edits
from spectral_fit import (RC_function_SY_FF, linear_function, do_fit_spec_map,
                          do_fit_spec_SY_FF_map, fit_spectra_chunk,
                          fit_spectra_SY_FF_chunk, fit_spectra_SY_FF_cov_chunk,
                          SY_FF_KEYS, RC_function_SY_FF_FFA, ffa_tau,
                          ffa_factors, fit_SY_FF_FFA_spectrum,
                          fit_spectra_SY_FF_FFA_chunk, FFA_KEYS,
                          FFA_FIT_DEFAULTS, FF_TAU_COEFF)
# spatial binning (e.g. contour binning) before spectral fitting
import spatial_binning
importlib.reload(spatial_binning)


def RC_function_SY_ffa(nu,Snu0,fth_nu0,alpha_nt,nu_tau_t,nu0=10,f_cov=1.0):
    tau_nu = (nu/nu_tau_t)**(-2.1)
    return((1 - fth_nu0) * Snu0 * (1-f_cov*(1-np.exp(-tau_nu))) * ((nu/nu0)**alpha_nt))

def RC_function_FF_ffa(nu,Snu0,fth_nu0,nu_tau_t,nu0=10):
    tau_nu = (nu/nu_tau_t)**(-2.1)
    return(fth_nu0 * Snu0 * ((1-np.exp(-tau_nu))/(tau_nu)) * ((nu/nu0)**(-0.1)))
    

def RC_function_SY_FF_ffa(nu, Snu0, fth_nu0, alpha_nt, nu_tau_t, nu0=10, f_cov=1.0):
    # return A_ff * (nu0 ** (-0.1)) * ((nu / nu0) ** (-0.1)) + A_sy * (nu0 ** (alpha_nt)) * ((nu / nu0) ** (alpha_nt))
    # tau_nu = (nu/nu_tau_t)**(-2.1)
    # S_ff_abs = fth_nu0 * Snu0 * ((1-np.exp(-tau_nu))/(np.exp(-tau_nu))) * ((nu/nu0)**(-0.1))
    S_ff_abs = RC_function_FF_ffa(nu,Snu0,fth_nu0,nu_tau_t,nu0)
    # S_sy_abs = (1 - fth_nu0) * Snu0 * (1-f_cov*(1-np.exp(-tau_nu))) * ((nu/nu0)**alpha_nt)
    S_sy_abs = RC_function_SY_ffa(nu,Snu0,fth_nu0,alpha_nt,nu_tau_t,nu0,f_cov)
    S_total_abs = S_ff_abs + S_sy_abs
    return S_total_abs


def RC_function_SY_ffa_v2(nu,A_sy,alpha_nt,nu_tau_t,nu0=10,f_cov=1.0):
    tau_nu = (nu/nu_tau_t)**(-2.1)
    return((10**A_sy) * (1-f_cov*(1-np.exp(-tau_nu))) * ((nu/nu0)**alpha_nt))

def RC_function_FF_ffa_v2(nu,A_ff,nu_tau_t,nu0=10):
    tau_nu = (nu/nu_tau_t)**(-2.1)
    return((10**A_ff) * ((1-np.exp(-tau_nu))/(tau_nu)) * ((nu/nu0)**(-0.1)))

def RC_function_SY_FF_ffa_v2(nu, A_sy, A_ff, alpha_nt, nu_tau_t, nu0=10,f_cov=1.0):
    # return A_ff * (nu0 ** (-0.1)) * ((nu / nu0) ** (-0.1)) + A_sy * (nu0 ** (alpha_nt)) * ((nu / nu0) ** (alpha_nt))
    # tau_nu = (nu/nu_tau_t)**(-2.1)
    # S_ff_abs = fth_nu0 * Snu0 * ((1-np.exp(-tau_nu))/(np.exp(-tau_nu))) * ((nu/nu0)**(-0.1))
    S_ff_abs = RC_function_FF_ffa_v2(nu,A_ff,nu_tau_t,nu0)
    # S_sy_abs = (1 - fth_nu0) * Snu0 * (1-f_cov*(1-np.exp(-tau_nu))) * ((nu/nu0)**alpha_nt)
    S_sy_abs = RC_function_SY_ffa_v2(nu,A_sy,alpha_nt,nu_tau_t,nu0,f_cov)
    S_total_abs = S_ff_abs + S_sy_abs
    return S_total_abs

def RC_function_linear(nu, A1l, alpha,nu0=10):
    return A1l * ((nu/nu0) ** (alpha))


def power_law_phys_model(nu, Reff, S0peak, alpha,nu0=1):
    """
    Physical motivated power law model for the radio emission. 
    
    It uses the effective radius of the emission as a parameter to be minimised.
    This assumes a Gaussian model, where the total flux density is given by 
    SG = (3 * pi * Speak * Reff ** 2)/2
    Lucatelli et al 2024. 
    """
    return 0.5 * 3 * np.pi * S0peak * (Reff ** 2.0) * ((nu/nu0) ** (alpha))


def calc_specidx_region_linear(freqs, fluxes, fluxes_err, ndim=3):
    from scipy.stats import t
    x = freqs / 1e9
    nu0 = np.mean(x)
    y = fluxes
    yerr = fluxes_err

    # plt.figure(figsize=(8, 6))
    # plt.errorbar(x, y,
    #              yerr=yerr,
    #              fmt='o', label='Observed data', color='k', ecolor='gray', alpha=0.5)
    # plt.xlabel('Frequency [GHz]')
    # plt.ylabel('Flux Density [mJy]')
    # plt.legend()
    # plt.semilogy()
    # plt.semilogx()
    # plt.show()

    def log_likelihood_linear(theta, nu, y):
        A1l, alpha = theta
        model = RC_function_linear(nu, A1l, alpha,nu0)
        # sigma2 = (np.std(y)) ** 2
        # sigma2 = yerr ** 2
        # likelihood = -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))
        sigma2 = yerr ** 2
        # sigma2 = (np.sqrt((yerr)**2.0+0.1))**2.0
        likelihood = -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))
        # likelihood = (y - model) / (y+np.sqrt((yerr)**2.0+0.1))#okay 
        return likelihood
        # return np.sum(y-model)**2.0

    # def log_likelihood_linear(theta, nu, y):
    #     dof = 10  # degrees of freedom for the t-distribution; adjust as needed
    #     A1l, alpha = theta
    #     model = RC_function_linear(nu, A1l, alpha)
    #     residuals = y - model
    #     log_likelihood = np.sum(t.logpdf(residuals / yerr, df=dof, scale=yerr))
    #     return log_likelihood

    def prior_transform_linear(utheta):
        ua1l, ualpha = utheta  # nu0 is not a parameter to be estimated, so it's removed from here
        a1l = -2.0 + ua1l * (5000.0)  # Transforming ua1l to be in the range -5 to 50
        alpha = -3.0 + ualpha * (3.0)  # Transforming ualpha_nt to be in the range -2 to 1.5
        return a1l, alpha

    ndim = 2
    dsampler = dynesty.DynamicNestedSampler(log_likelihood_linear, prior_transform_linear,
                                            bound='balls',
                                            ndim=ndim, nlive=30000, sample='rwalk',  # hslice, rwalk
                                            logl_args=(x, y))

    # dsampler.run_nested(dlogz_init=0.05, nlive_init=1500, nlive_batch=500,
    #                     # use_stop=False,
    #                 # maxiter_init=10000, maxiter_batch=1000, maxbatch=10
    #                    )
    dsampler.run_nested(nlive_init=5000, nlive_batch=1000)
    results = dsampler.results

    results.summary()

    try:
        lnz_truth = ndim * -np.log(2 * 10.)  # analytic evidence solution
        fig, axes = dyplot.runplot(results, lnz_truth=lnz_truth)

        fig, axes = dyplot.traceplot(results, truths=np.zeros(ndim),
                                     truth_color='black', show_titles=True,
                                     trace_cmap='magma', connect=True,
                                     connect_highlight=range(5))
    except:
        pass

    try:
        # # initialize figure
        # fig, axes = plt.subplots(3, 3, figsize=(8, 8))
        # axes = axes.reshape((3, 3))  # reshape axes
        #
        # fg, ax = dyplot.cornerplot(results, color='blue', truths=np.zeros(ndim),
        #                            truth_color='black', show_titles=True,
        #                            max_n_ticks=3, quantiles=None,
        #                            fig=(fig, axes[:, :3]))
        fig, axes = dyplot.cornerplot(results,
                                      show_titles=True,  # Show titles with parameter values
                                      title_kwargs={'y': 1.03},  # Adjust title position
                                      labels=[r'$A_1$', r'$\alpha$'],  # Label axes
                                      label_kwargs=dict(fontsize=16),
                                      # Adjust label font size
                                      # quantiles = [0.1, 0.5, 0.9],
                                      # title_quantiles = [0.1, 0.5, 0.9],
                                      max_n_ticks=4, use_math_text=True, verbose=True,
                                      title_fmt='.3f',
                                      fig=plt.subplots(2, 2, figsize=(5, 5)))
    except:
        pass

    # Extract the best-fit parametershttps://dynesty.readthedocs.io/en/latest/examples.html
    best_fit_indices = np.argmax(results.logl)
    best_fit_params = results.samples[best_fit_indices]
    A1l_best, alpha_nt_best = best_fit_params

    # Calculate the best-fit model
    x_resample = np.linspace(np.min(x)*0.7, np.max(x)*1.3, 100)
    best_fit_model = RC_function_linear(x_resample, A1l_best, alpha_nt_best,nu0)

    # Plotting the data
    plt.figure(figsize=(8, 6))
    plt.errorbar(x, y, yerr=yerr, fmt='o', label='Observed data', color='k', ecolor='gray',
                 alpha=0.5)
    plt.plot(x_resample, best_fit_model, label='Best-fit model', color='red')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux Density [Jy]')
    plt.semilogx()
    plt.semilogy()
    plt.legend()
    plt.show()

    try:
        # Extracting samples for the corner plot
        samples = results.samples

        # Creating a corner plot
        fig = corner.corner(samples,
                            labels=[r'$A_1$', r'$\alpha_{\rm nt}$'],
                            truths=[A1l_best, alpha_nt_best],
                            quantiles=[0.025, 0.5, 0.975], show_titles=True,
                            title_kwargs={"fontsize": 12})
        plt.show()
    except:
        pass
    return (results, dsampler)


def calc_specidx_region_S2(freqs, fluxes, fluxes_err, ndim=3):
    x = freqs / 1e9
    nu0 = np.mean(x)
    y = fluxes
    yerr = fluxes_err

    plt.figure(figsize=(8, 6))
    plt.errorbar(x, y,
                 yerr=yerr,
                 fmt='o', label='Observed data', color='k', ecolor='gray', alpha=0.5)
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Flux Density [mJy]')
    plt.legend()
    plt.semilogy()
    plt.semilogx()
    plt.show()

    def log_likelihood(theta, nu, y):
        A1l, A2, alpha_nt = theta
        model = RC_function_S2(nu, A1l, A2, alpha_nt, nu0)
        sigma2 = yerr ** 2
        return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))

    def prior_transform(utheta):
        ua1l, ua2, ualpha_nt = utheta  # nu0 is not a parameter to be estimated, so it's removed from here
        a1l = 0 + ua1l * (15)  # Transforming ua1l to be in the range -5 to 50
        a2 = 0 + ua2 * (15)  # Transforming ua2 to be in the range -5 to 50
        alpha_nt = -2.0 + ualpha_nt * (2.0)  # Transforming ualpha_nt to be in the range -2 to 1.5
        return a1l, a2, alpha_nt

    dsampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, bound='balls',
                                            ndim=ndim, nlive=15000, sample='rwalk',  # hslice, rwalk
                                            logl_args=(x, y))

    # dsampler.run_nested(dlogz_init=0.05, nlive_init=1500, nlive_batch=500,
    #                     # use_stop=False,
    #                 # maxiter_init=10000, maxiter_batch=1000, maxbatch=10
    #                    )
    dsampler.run_nested(nlive_init=2500, nlive_batch=1000)
    # dsampler.run_nested()
    results = dsampler.results

    results.summary()

    try:
        lnz_truth = ndim * -np.log(2 * 10.)  # analytic evidence solution
        fig, axes = dyplot.runplot(results, lnz_truth=lnz_truth)

        fig, axes = dyplot.traceplot(results, truths=np.zeros(ndim),
                                     truth_color='black', show_titles=True,
                                     trace_cmap='magma', connect=True,
                                     connect_highlight=range(5))
    except:
        pass

    try:
        # initialize figure
        fig, axes = dyplot.cornerplot(results,
                                      show_titles=True,  # Show titles with parameter values
                                      title_kwargs={'y': 1.03},  # Adjust title position
                                      labels=[r'$A_1$',r'$A_2$', r'$\alpha$'],  # Label axes
                                      label_kwargs=dict(fontsize=16),
                                      # Adjust label font sibze
                                      # quantiles = [0.16, 0.5, 0.84],title_quantiles = [0.16, 0.5, 0.84],
                                      max_n_ticks=4, use_math_text=True, verbose=True,
                                      title_fmt='.3f',
                                      fig=plt.subplots(3, 3, figsize=(7, 7)))
    except:
        pass

    # Extract the best-fit parametershttps://dynesty.readthedocs.io/en/latest/examples.html
    best_fit_indices = np.argmax(results.logl)
    best_fit_params = results.samples[best_fit_indices]
    A1l_best, A2_best, alpha_nt_best = best_fit_params

    # Calculate the best-fit model
    x_resample = np.linspace(np.min(x)*0.5, np.max(x)*2, 100)
    best_fit_model = RC_function_S2(x_resample, A1l_best, A2_best, alpha_nt_best, nu0)

    # Plotting the data
    plt.figure(figsize=(8, 6))
    plt.errorbar(x, y, yerr=yerr, fmt='o', label='Observed data', color='k', ecolor='gray',
                 alpha=0.5)
    plt.plot(x_resample, best_fit_model, label='Best-fit model', color='red')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux Density')
    plt.semilogx()
    plt.semilogy()
    plt.legend()
    plt.show()

    try:
        # Extracting samples for the corner plot
        samples = results.samples

        # Creating a corner plot
        fig = corner.corner(samples,
                            labels=[r'$A_1$', r'$A_2$', r'$\alpha_{\rm nt}$'],
                            truths=[A1l_best, A2_best, alpha_nt_best],
                            quantiles=[0.025, 0.5, 0.975], show_titles=True,
                            title_kwargs={"fontsize": 12})
        plt.show()
    except:
        pass
    return (results, dsampler)

def RC_function_Sq(nu,S0,alpha,q):
    return S0*(nu**(alpha))*np.exp(q * (np.log(nu))**2.0)

def calc_specidx_region_Sq(freqs, fluxes, fluxes_err, ndim=3):
    x = freqs / 1e9
    nu0 = np.mean(x)
    y = fluxes
    yerr = fluxes_err
    S0_init = 2*np.max(y)

    plt.figure(figsize=(8, 6))
    plt.errorbar(x, y,
                 yerr=yerr,
                 fmt='o', label='Observed data', color='k', ecolor='gray', alpha=0.5)
    plt.xlabel('Frequency [GHz]')
    plt.ylabel('Flux Density [mJy]')
    plt.legend()
    plt.semilogy()
    plt.semilogx()
    plt.show()

    def log_likelihood(theta, nu, y):
        S0, alpha, q = theta
        model = RC_function_Sq(nu, S0, alpha, q)
        sigma2 = yerr ** 2
        return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))

    def prior_transform(utheta):
        us0, ualpha, uq = utheta
        # from here
        s0 = 0 + us0 * (S0_init)
        alpha = -2.0 + ualpha * (2.0)
        q = -0.9 + uq * (0.9)
        return s0, alpha, q


    dsampler = dynesty.DynamicNestedSampler(log_likelihood, prior_transform, bound='balls',
                                            ndim=ndim, nlive=30000, sample='rwalk',  # hslice, rwalk
                                            logl_args=(x, y))

    # dsampler.run_nested(dlogz_init=0.05, nlive_init=1500, nlive_batch=500,
    #                     # use_stop=False,
    #                 # maxiter_init=10000, maxiter_batch=1000, maxbatch=10
    #                    )
    dsampler.run_nested(nlive_init=5000, nlive_batch=1000)
    # dsampler.run_nested()
    results = dsampler.results

    results.summary()

    try:
        lnz_truth = ndim * -np.log(2 * 10.)  # analytic evidence solution
        fig, axes = dyplot.runplot(results, lnz_truth=lnz_truth)

        fig, axes = dyplot.traceplot(results, truths=np.zeros(ndim),
                                     truth_color='black', show_titles=True,
                                     trace_cmap='magma', connect=True,
                                     connect_highlight=range(5))
    except:
        pass

    try:
        # initialize figure
        fig, axes = dyplot.cornerplot(results,
                                      show_titles=True,  # Show titles with parameter values
                                      title_kwargs={'y': 1.03},  # Adjust title position
                                      labels=[r'$S_{0}$',r'$\alpha$', r'$q$'],  # Label axes
                                      label_kwargs=dict(fontsize=16),
                                      # Adjust label font size
                                      # quantiles = [0.16, 0.5, 0.84],title_quantiles = [0.16, 0.5, 0.84],
                                      max_n_ticks=4, use_math_text=True, verbose=True,
                                      title_fmt='.3f',
                                      fig=plt.subplots(3, 3, figsize=(7, 7)))
    except:
        pass

    # Extract the best-fit parametershttps://dynesty.readthedocs.io/en/latest/examples.html
    best_fit_indices = np.argmax(results.logl)
    best_fit_params = results.samples[best_fit_indices]
    S0_best, alpha_best, q_best = best_fit_params

    # Calculate the best-fit model
    x_resample = np.linspace(x[0], x[-1], 100)
    best_fit_model = RC_function_Sq(x_resample, S0_best, alpha_best, q_best)

    # Plotting the data
    plt.figure(figsize=(8, 6))
    plt.errorbar(x, y, yerr=yerr, fmt='o', label='Observed data', color='k', ecolor='gray',
                 alpha=0.5)
    plt.plot(x_resample, best_fit_model, label='Best-fit model', color='red')
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Flux Density')
    plt.semilogx()
    plt.semilogy()
    plt.legend()
    plt.show()

    try:
        # Extracting samples for the corner plot
        samples = results.samples

        # Creating a corner plot
        fig = corner.corner(samples,
                            labels=[r'$S_{0}$',r'$\alpha$', r'$q$'],
                            truths=[S0_best, alpha_best, q_best],
                            quantiles=[0.025, 0.5, 0.975], show_titles=True,
                            title_kwargs={"fontsize": 12})
        plt.show()
    except:
        pass
    return (results, dsampler)


def general_mcmc(x_data, y_data, yerr_data, 
                 result_1, param_labels, 
                 model_func, nwalkers=None,
                 burn_in = 2500,thin=5,
                 nsteps=7500, sigma_errors=1.0, 
                 prior_sigma=15,
                 quantiles = [2.5, 50, 97.5]):
    """
    General MCMC simulation for a model with n parameters.
    
    Parameters:
    - x_data: Independent variable data
    - y_data: Dependent variable data
    - yerr_data: Error in dependent variable data
    - result_1: Minimization result object from lmfit
    - param_labels: List of parameter names
    - model_func: The model function to fit
    - nwalkers: Number of walkers for the MCMC sampler (default is 25 * number of parameters)
    - nsteps: Number of MCMC steps (default is 7500)
    - sigma_errors: Number of standard deviations of the best parameters distribution, where
      p +/- sigma * p_err 
    - prior_sigma: Number of standard deviations for prior (default is 15.0)
    
    Returns:
    - samples: MCMC samples for the parameters
    """
    # Extract best-fit parameters and their standard errors from the minimization result
    best_fit_params = np.asarray([result_1.params[label].value for label in param_labels])
    params_stderr = np.asarray([result_1.params[label].stderr for label in param_labels])
    
    mask = ~np.isnan(y_data)
    x_data = x_data[mask]
    y_data = y_data[mask]
    yerr_data = yerr_data[mask]

    # Number of dimensions (parameters)
    ndim = len(best_fit_params)
    if nwalkers is None:
        nwalkers = int(25 * ndim)
    # Define the log-probability function with priors
    # def log_prob(params, x, y, yerr, best_fit_params, params_stderr):
    #     # Prior: within +/- 10 sigma
        
    #     # if params[0] < 0:
    #     #     return -np.inf
    #     # if params[1] < 0:
    #     #     return -np.inf
        
    #     if not all(best_fit_params[i] - prior_sigma*params_stderr[i] < params[i] < best_fit_params[i] + prior_sigma*params_stderr[i] for i in range(ndim)):
    #         return -np.inf
    #     # Calculate the model predictions
    #     model_pred = model_func(params, x)
    #     # Calculate the log-likelihood
    #     sigma2 = yerr**2
    #     log_likelihood = -0.5 * np.sum((y - model_pred)**2 / sigma2 + np.log(sigma2))
    #     return log_likelihood

    # def log_prob(params, x, y, yerr, best_fit_params, params_stderr):
    #     # Check for NaNs in the input data or parameters
    #     if np.any(np.isnan(x)) or np.any(np.isnan(y)) or np.any(np.isnan(yerr)) or np.any(np.isnan(params)):
    #         return -np.inf  # Return negative infinity to ignore samples with NaNs
        
    #     # if params[0] < 0:
    #     #     return -np.inf
    #     # if params[1] < 0:
    #     #     return -np.inf
        
    #     # Prior: within +/- 10 sigma
    #     if not all(best_fit_params[i] - prior_sigma*params_stderr[i] < params[i] < best_fit_params[i] + prior_sigma*params_stderr[i] for i in range(ndim)):
    #         return -np.inf
        
    #     # Calculate the model predictions
    #     model_pred = model_func(params, x)
        
    #     # Calculate the log-likelihood
    #     sigma2 = yerr**2
    #     log_likelihood = -0.5 * np.sum((y - model_pred)**2 / sigma2 + np.log(sigma2))
        
    #     return log_likelihood
    
    def log_prob(params, x, y, yerr, best_fit_params, params_stderr):
        # Check for NaN in the parameters
        if np.any(np.isnan(params)):
            return -np.inf
        
        # if params[0] < -5:
        #     return -np.inf
        # if params[1] < -5:
        #     return -np.inf

        # Prior: within +/- prior_sigma * sigma
        if not all(best_fit_params[i] - prior_sigma * params_stderr[i] < params[i] < best_fit_params[i] + prior_sigma * params_stderr[i] for i in range(ndim)):
            return -np.inf

        # Calculate the model predictions
        model_pred = model_func(params, x)

        # Check for NaNs in model predictions
        if np.any(np.isnan(model_pred)):
            return -np.inf

        # Calculate sigma squared
        sigma2 = yerr ** 2

        # Check for zeros or NaNs in sigma2
        if np.any(sigma2 <= 0) or np.any(np.isnan(sigma2)):
            return -np.inf

        # Check for NaNs in y
        if np.any(np.isnan(y)):
            return -np.inf

        # Calculate the log-likelihood
        residuals = y - model_pred
        if np.any(np.isnan(residuals)):
            return -np.inf

        log_likelihood = -0.5 * np.sum((residuals) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))

        # Check for NaN in log_likelihood
        if np.isnan(log_likelihood):
            return -np.inf

        return log_likelihood

    # Initialize the walkers
    p0 = []
    for i in range(nwalkers):
        walker_params = []
        random_params = np.asarray(generate_random_params_normal(best_fit_params, 
                                                                 params_stderr,
                                                                 sigma_errors=sigma_errors))
        p0.append(list(random_params))
    p0 = np.asarray(p0) 

    # Set up the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(x_data, 
                                                                    y_data, 
                                                                    yerr_data, 
                                                                    best_fit_params, 
                                                                    params_stderr))

    # Run MCMC
    sampler.run_mcmc(p0, nsteps, progress=True)

    # Analyze the samples
    samples = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
    
    param_dict = {}
    for i, label in enumerate(param_labels):
        q = np.percentile(samples[:, i], [quantiles[0]*100,
                                          quantiles[1]*100,
                                          quantiles[2]*100])
        param_dict[label] = {
            'best': q[1],
            'lower_bound': q[0],
            'upper_bound': q[2],
            'lower': (q[1] - q[0]),
            'upper': (q[2] - q[1])
        }
    
    # # Plotting the samples
    # fig = corner.corner(samples, 
    #                     labels=param_labels,
    #                     show_titles=True, quantiles=[0.025, 0.5, 0.975],
    #                     truths=best_fit_params)
    # fig.show()
    
    return samples,param_dict


def general_mcmc_physical(x_data, y_data, yerr_data, 
                 result_1, param_labels, 
                 model_func, S0peak, nu0=1.0,
                 nwalkers=None,
                 burn_in = 2500,thin=5,
                 nsteps=7500, sigma_errors=1.0, prior_sigma=15):
    """
    General MCMC simulation for a model with n parameters.
    
    Parameters:
    - x_data: Independent variable data
    - y_data: Dependent variable data
    - yerr_data: Error in dependent variable data
    - result_1: Minimization result object from lmfit
    - param_labels: List of parameter names
    - model_func: The model function to fit
    - nwalkers: Number of walkers for the MCMC sampler (default is 200)
    - nsteps: Number of MCMC steps (default is 7500)
    - sigma_errors: Number of standard deviations of the best parameters distribution, where
      p +/- sigma * p_err 
    - prior_sigma: Number of standard deviations for prior (default is 15.0)
    
    Returns:
    - samples: MCMC samples for the parameters
    """
    # Extract best-fit parameters and their standard errors from the minimization result
    best_fit_params = np.asarray([result_1.params[label].value for label in param_labels])
    params_stderr = np.asarray([result_1.params[label].stderr for label in param_labels])
    
    mask = ~np.isnan(y_data)
    x_data = x_data[mask]
    y_data = y_data[mask]
    yerr_data = yerr_data[mask]

    # Number of dimensions (parameters)
    ndim = len(best_fit_params)
    if nwalkers is None:
        nwalkers = int(100 * ndim)
    # Define the log-probability function with priors
    def log_prob(params, x, y, yerr, best_fit_params, params_stderr):
        # Prior: within +/- 10 sigma
        if not all(best_fit_params[i] - prior_sigma*params_stderr[i] < params[i] < best_fit_params[i] + prior_sigma*params_stderr[i] for i in range(ndim)):
            return -np.inf
        # Calculate the model predictions
        model_pred = model_func(params, x, S0peak, nu0)
        # Calculate the log-likelihood
        sigma2 = yerr**2
        log_likelihood = -0.5 * np.sum((y - model_pred)**2 / sigma2 + np.log(sigma2))
        return log_likelihood

    # Initialize the walkers
    p0 = []
    for i in range(nwalkers):
        walker_params = []
        random_params = np.asarray(generate_random_params_normal(best_fit_params, 
                                                                 params_stderr,
                                                                 sigma_errors=sigma_errors))
        p0.append(list(random_params))
    p0 = np.asarray(p0) 

    # Set up the MCMC sampler
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(x_data, 
                                                                    y_data, 
                                                                    yerr_data, 
                                                                    best_fit_params, 
                                                                    params_stderr))

    # Run MCMC
    sampler.run_mcmc(p0, nsteps, progress=True)

    # Analyze the samples
    samples = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
    
    param_dict = {}
    for i, label in enumerate(param_labels):
        q = np.percentile(samples[:, i], [2.5, 50, 97.5])
        param_dict[label] = {
            'best': q[1],
            'lower': (q[1] - q[0]),
            'upper': (q[2] - q[1])
        }
    
    # # Plotting the samples
    # fig = corner.corner(samples, 
    #                     labels=param_labels,
    #                     show_titles=True, quantiles=[0.025, 0.5, 0.975],
    #                     truths=best_fit_params)
    # fig.show()
    
    return samples,param_dict


def power_law_phys_fit(freqs,
                          fluxes,
                          fluxes_err,
                          S0peak,
                          _nu0=1.0,
                          basename_save=None,log_plot=True,
                          save_name_append = None,
                          plot_errors_shade = False,
                          do_mcmc_fit = False,
                          mcmc_version = 'general',
                          title_text = None,
                          add_fit_legend = True,
                          verbose=0):
    """
    Peform a fit to the radio spectrum using a linear model.

    Parameters
    ----------
    freqs : array
        Frequency array in Hz.
    fluxes : array
        Flux density array in mJy.
    fluxes_err : array
        Flux density error array in mJy.
    basename_save : str
        Basename to save the output files.
    verbose : int
        Verbosity level.

    """
    x = freqs / 1e9
    y = fluxes
    yerr = fluxes_err
    weights = 1.0 / yerr
    if _nu0 != 1.0:
        nu0 = _nu0/1e9
    else:
        nu0 = _nu0
        
    def min_func(params):
        Reff = params['Reff']
        alpha = params['alpha']
        model = power_law_phys_model(x, Reff, S0peak, alpha, nu0)
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights


    fit_params = lmfit.Parameters()
    fit_params.add("Reff", value=10.0, min=-0.01, max=10)
    fit_params.add("alpha", value=-0.9, min=-4.0, max=4.0)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="huber", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="huber", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))
    
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 3),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # fig, ax = plt.subplots()
    x_resample = np.linspace(np.min(x)*0.7, np.max(x)*1.3, 500)
    ax1.errorbar(x, y, yerr=yerr, 
                 fmt='o',label='Data', color='k', ecolor='gray',alpha=0.5)
    
    if do_mcmc_fit == True:
        nwalkers = int(2 * 100)
        burn_in = 2500
        if mcmc_version == 'lmfit':
            """ 
            This will be removed in the future.
            """
            results_emcee = mini.emcee(burn=burn_in,
                                    steps=4000,
                                    thin=5, nwalkers=nwalkers,
                                    is_weighted=False,
                                    float_behavior='posterior',
                                    params=result.params, workers=-1)

            _samples_emcee = results_emcee.flatchain[burn_in:]
            _Reff = np.asarray(_samples_emcee['Reff'])
            _alpha = np.asarray(_samples_emcee['alpha'])
            model_samples = np.array(
                [power_law_phys_model(x_resample, _Reff[i], S0peak, _alpha[i], nu0) for i in
                range(_samples_emcee.shape[0])])
            samples_emcee = np.asarray(_samples_emcee[['Reff', 'alpha']])
            param_dict = {}
            for i, label in enumerate(['Reff', 'alpha']):
                q = np.percentile(samples_emcee[:, i], [2.5, 50, 97.5])
                param_dict[label] = {
                    'best': q[1],
                    'lower': (q[1] - q[0]),
                    'upper': (q[2] - q[1])
                }
            print(param_dict)
        if mcmc_version == 'general':
            """
            This seems to be more robust as we have more control over
            the MCMC process.
            """
            def RC_powerlaw_physical(params,x, S0peak, x0):
                return 0.5*3.0*np.pi*S0peak*(params[0]**2.0) * ((x/x0)**params[1])
                
            samples_emcee, param_dict = general_mcmc_physical(x_data = x, 
                                                    y_data = y, 
                                                    yerr_data = yerr, 
                                                    result_1 = result, 
                                                    param_labels = ["Reff", "alpha"], 
                                                    model_func = RC_powerlaw_physical,
                                                    S0peak = S0peak, 
                                                    nu0=nu0,
                                                    burn_in = burn_in,
                                                    nwalkers = nwalkers)
            _Reff = samples_emcee.T[0]
            _alpha = samples_emcee.T[1]

            model_samples = np.array(
                [power_law_phys_model(x_resample, 
                                    _Reff[i], S0peak,
                                    _alpha[i], nu0) for i in range(_Reff.shape[0])])
            print(param_dict)

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = power_law_phys_model(x_resample,
                                    result.params['Reff'].value,
                                    S0peak,
                                    result.params['alpha'].value,
                                    nu0)

    model_best = power_law_phys_model(x,
                                    result.params['Reff'].value,
                                    S0peak,
                                    result.params['alpha'].value,
                                    nu0)

    ax1.plot(x_resample,model_resample,
             color='red', ls='-.', label='Best-fit model')
    ax2.plot(x, (y-model_best)/model_best,
             color='green', ls='dotted', label='Residual')
    ax2.set_ylim(-1.0,1.0)

    if do_mcmc_fit == True:
        ax1.plot(x_resample, model_mean,
                 color='purple', linestyle=(0, (5, 10)), label='MCMC Mean')

    # plt.fill_between(x_resample, lower_bound, upper_bound, color='lightgray', alpha=0.5,
    #                  label='Parameter uncertainty')
    if plot_errors_shade:
        if do_mcmc_fit == True:
            ax1.fill_between(x_resample,
                             model_mean - 3*model_std,
                             model_mean + 3*model_std, color='lightgray',
                            alpha=0.7)
        else:
            # Define the number of Monte Carlo samples
            num_samples = 5000

            # Generate random samples from parameter distributions
            Reff_samples = np.random.normal(result.params['Reff'].value, 
                                            result.params['Reff'].stderr,
                                            num_samples)
            alpha_samples = np.random.normal(result.params['alpha'].value,
                                             result.params['alpha'].stderr,
                                             num_samples)

            # Compute model predictions for each sample
            model_predictions = np.zeros((num_samples, len(x_resample)))
            for i in range(num_samples):
                model_predictions[i] = power_law_phys_model(x_resample, Reff_samples[i],
                                                            S0peak,
                                                            alpha_samples[i],nu0)

            median_prediction = np.median(model_predictions, axis=0)
            std_prediction = np.std(model_predictions, axis=0)

            ax1.fill_between(x_resample, median_prediction - 1*std_prediction,
                             median_prediction + 1*std_prediction,
                             color='lightgray', alpha=0.5,
                             # label='Uncertainty (1-sigma)'
                             )
    # plt.ylim(1e-3,1.2*np.max(y))
    if add_fit_legend == True:
        ax1.legend(loc=(0.05, 0.05),frameon=True,prop={'size': 11})
    if np.nanmax(y) < 5:
        ax1.set_ylim(0.1*np.nanmin(y),10.0*np.nanmax(y))
    else:
        ax1.set_ylim(0.1*np.nanmin(y),3.0*np.nanmax(y))
    ax2.set_xlabel(r'$\nu$ [GHz]')
    # plt.ylabel('Integrated Flux Density [mJy]')
    ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
    text_x, text_y = 0.65, 0.37
    text = (r"$\alpha_{\rm nth}"+rf"= {(result.params['alpha'].value):.2f}\pm "
            rf"{(result.params['alpha'].stderr):.2f}$")
    
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)
    if do_mcmc_fit == True:
        if mcmc_version == 'general':
            text_x, text_y = 0.65, 0.77
            text = (r"$\alpha_{\rm nth}^{\rm MCMC}"+f"= {(param_dict['alpha']['best']):.2f}"
                    rf"_{{-{param_dict['alpha']['lower']:.2f}}}^{{+{param_dict['alpha']['upper']:.2f}}}$")
            text_bbox = plt.text(text_x, text_y, text,
                                # ha='center', va='center',
                                fontsize=12, color='black',
                                # bbox=text_bbox_props, 
                                transform=fig.transFigure)
    
    if title_text is not None:
        ax1.set_title(title_text)

    if log_plot == True:
        ax1.semilogx()
        ax1.semilogy()

    plt.subplots_adjust(hspace=0.05)
    ax2.legend(loc=(0.7, 0.55),prop={'size': 14},frameon=False)
    # legend()

    if basename_save is not None:
        if save_name_append is None:
            save_name_append = '_RC_alpha_fit_linear'
        else:
            save_name_append = save_name_append + '_RC_alpha_fit_linear'
        plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                    bbox_inches='tight')
    # plt.show()
    corner_kwargs = {
        'bins': 30,
        'hist_bin_factor': 1.0,
        'color': 'purple',  # Change the color of the scatter plots and histograms
        'hist_kwargs': {
            'color': 'green',      # Histogram color
            'edgecolor': 'black',  # Edge color of the bins
            'linewidth': 1.5       # Width of the bin edges
        },
        'scatter_kwargs': {
            'alpha': 0.6,          # Transparency of points
            's': 10,               # Size of points
            'color': 'purple',     # Color of points
            'edgecolor': 'none'    # Edge color of points
        },
        'contour_kwargs': {
            'colors': 'blue',       # Color of contour lines
            'linewidths': 1.2       # Width of contour lines
        }
    }
    if do_mcmc_fit == True:
        try:
            from scipy.stats import gaussian_kde
            # fig_c = plt.figure(figsize=(2,2))
            fig_c = plt.figure()
            _ = corner.corner(samples_emcee,
                                    labels=[r'$R_{\rm eff}$',r'$\alpha$'],
                                    truths=[result.params['Reff'].value,
                                            result.params['alpha'].value],
                                    show_titles=True,
                                    quantiles=[0.025, 0.5, 0.975],
                                    # **corner_kwargs
                                    # fig=fig_c
                                    )
            print(samples_emcee.shape)
            if basename_save is not None:
                save_name_append_corner = save_name_append + '_corner'
                plt.savefig(basename_save.replace('.fits','_')+save_name_append_corner+'.jpg', 
                            dpi=600,bbox_inches='tight')
                
            print('++==>> Parameter Results (MCMC sampling).')
            print(lmfit.fit_report(results_emcee.params))
        except:
            pass
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))
    return(mini,result)




def do_fit_spec_RC_linear(freqs,
                          fluxes,
                          fluxes_err,
                          nu0=None,
                          basename_save=None,log_plot=True,
                          save_name_append = None,
                          plot_errors_shade = False,
                          do_mcmc_fit = False,
                          sigma_shade=3,
                          sigma_errors = 3,
                          quantiles = [0.16, 0.5, 0.84],
                          mcmc_version = 'general',
                          burn_in = 1000,
                          nsteps = 5000,
                          thin = 2,
                          title_text = None,
                          ext_save = '.jpg',
                          add_fit_legend = True,
                          figsize=(6, 3),
                          verbose=0):
    """
    Peform a fit to the radio spectrum using a linear model.

    Parameters
    ----------
    freqs : array
        Frequency array in Hz.
    fluxes : array
        Flux density array in mJy.
    fluxes_err : array
        Flux density error array in mJy.
    basename_save : str
        Basename to save the output files.
    verbose : int
        Verbosity level.

    """
    
    # if sigma_errors == 1 or sigma_errors == 1.0:
    #     quantiles = [0.16, 0.5, 0.84]
    # elif sigma_errors == 2 or sigma_errors == 2.0:
    #     quantiles = [0.025, 0.5, 0.975]
    # elif sigma_errors == 3 or sigma_errors == 3.0:
    #     quantiles = [0.00015, 0.5, 0.9985]
    # else:
    #     raise ValueError("Unsupported sigma value. Use sigma = 1, 2, or 3.")
    
    x = freqs / 1e9
    y = fluxes
    if nu0 is None:
        nu0 = np.nanmean(x)
    yerr = fluxes_err
    # weights = 1.0 / yerr
    epsilon = 1e-8
    # def min_func(params):
    #     A1 = params['A1']
    #     alpha = params['alpha']
    #     model = RC_function_linear(x, A1, alpha,nu0)
    #     relative_error = yerr / (np.abs(y) + epsilon)
    #     weights = 1 / (relative_error + epsilon)
    #     # weightned_residual = (y - model) * np.sqrt(weights)
    #     # weightned_residual = (y - model) / (np.sqrt(weights + yerr))
    #     weightned_residual = (y - model) / (np.sqrt(weights)+yerr)
    #     return weightned_residual.copy()

    def min_func(params):
        A1 = params['A1']
        alpha = params['alpha']
        model = RC_function_linear(x, A1, alpha,nu0)
        # weights = 1 / yerr
        # return (y - model) * np.sqrt(weights)
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights


    fit_params = lmfit.Parameters()
    fit_params.add("A1", value=10.0, min=-0.1, max=10000)
    fit_params.add("alpha", value=-0.8, min=-4.0, max=4.0)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    # print('++==>> Parameter Results (from least-squares fit).')
    # print(lmfit.fit_report(result.params))
    # result_1 = mini.minimize(method='nelder',
    #                          options={'maxiter': 30000, 'maxfev': 30000,
    #                                   'xatol': 1e-12, 'fatol': 1e-12,
    #                                   'return_all': True,'adaptive': True,
    #                                   'disp': True}
    #                      )
    # second_run_params = result_1.params

    # result = mini.minimize(method='nelder',params=second_run_params,
    #                        options={'maxiter': 30000, 'maxfev': 30000,
    #                                 'xatol': 1e-12, 'fatol': 1e-12,
    #                                 'return_all': True,'adaptive': True,
    #                                 'disp': True}
    #                  )

    x_resample = np.linspace(np.min(x)*0.9, np.max(x)*1.1, 500)
    
    if do_mcmc_fit == True:
        nwalkers = int(2 * 25)
        print(f"  ++==>> Number of walkers: {nwalkers}")
        if mcmc_version == 'lmfit':
            """ 
            This will be removed in the future.
            """
            results_emcee = mini.emcee(burn=burn_in,
                                    steps=4000,
                                    thin=5, nwalkers=nwalkers,
                                    is_weighted=False,
                                    float_behavior='posterior',
                                    params=result.params, workers=-1)

            _samples_emcee = results_emcee.flatchain[burn_in:]
            _A1 = np.asarray(_samples_emcee['A1'])
            _alpha = np.asarray(_samples_emcee['alpha'])
            model_samples = np.array(
                [RC_function_linear(x_resample, _A1[i], _alpha[i],nu0) for i in
                range(_samples_emcee.shape[0])])
            samples_emcee = np.asarray(_samples_emcee[['A1', 'alpha']])
            param_dict = {}
            for i, label in enumerate(['A1', 'alpha']):
                q = np.percentile(samples_emcee[:, i], [2.5, 50, 97.5])
                param_dict[label] = {
                    'best': q[1],
                    'lower': (q[1] - q[0]),
                    'upper': (q[2] - q[1])
                }
            print(param_dict)
        if mcmc_version == 'general':
            """
            This seems to be more robust as we have more control over
            the MCMC process.
            """
            def RC_powerlaw(params,x):
                return params[0] * ((x/nu0)**params[1])
                
            samples_emcee, param_dict = general_mcmc(x_data = x, 
                                                    y_data = y, 
                                                    yerr_data = yerr, 
                                                    result_1 = result, 
                                                    param_labels = ["A1", "alpha"], 
                                                    burn_in = burn_in,
                                                    nwalkers = nwalkers,
                                                    nsteps = nsteps,
                                                    thin = thin,
                                                    sigma_errors = sigma_errors,
                                                    model_func = RC_powerlaw,
                                                    quantiles=quantiles)
            _A1 = samples_emcee.T[0]
            _alpha = samples_emcee.T[1]

            model_samples = np.array(
                [RC_function_linear(x_resample, 
                                    _A1[i], 
                                    _alpha[i],nu0) for i in range(_A1.shape[0])])
            print(param_dict)

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = RC_function_linear(x_resample,
                                    result.params['A1'].value,
                                    result.params['alpha'].value,
                                    nu0)

    model_best = RC_function_linear(x,
                                    result.params['A1'].value,
                                    result.params['alpha'].value,
                                    nu0)

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, 
                                #    figsize=figsize,
                                   figsize=(figsize[0], figsize[1]*1.5),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # fig, ax = plt.subplots()
    


    # nwalkers = int(len(y) * 2 * 20)

    ax1.errorbar(x, y, yerr=yerr, 
                 fmt='o',label='Data', color='k', ecolor='gray',alpha=0.5)
    ax1.plot(x_resample,model_resample,
             color='blue', ls='-', label='Fit')

    if do_mcmc_fit == True:
        ax1.plot(x_resample, model_mean,
                 color='purple', linestyle=(0, (5, 10)), label='MCMC Mean')

    # plt.fill_between(x_resample, lower_bound, upper_bound, color='lightgray', alpha=0.5,
    #                  label='Parameter uncertainty')
    if plot_errors_shade:
        if do_mcmc_fit == True:
            ax1.fill_between(x_resample,
                             model_mean - sigma_shade*model_std,
                             model_mean + sigma_shade*model_std, color='lightgray',
                             label=rf'$\pm{sigma_shade}\sigma$', alpha=0.7)
        else:
            # Define the number of Monte Carlo samples
            num_samples = 5000

            # Generate random samples from parameter distributions
            A1_samples = np.random.normal(result.params['A1'].value, result.params['A1'].stderr,
                                          num_samples)
            alpha_samples = np.random.normal(result.params['alpha'].value,
                                             result.params['alpha'].stderr,
                                             num_samples)

            # Compute model predictions for each sample
            model_predictions = np.zeros((num_samples, len(x_resample)))
            for i in range(num_samples):
                model_predictions[i] = RC_function_linear(x_resample, A1_samples[i],
                                                          alpha_samples[i],nu0)

            median_prediction = np.nanmedian(model_predictions, axis=0)
            std_prediction = np.nanstd(model_predictions, axis=0)

            ax1.fill_between(x_resample, median_prediction - sigma_shade*std_prediction,
                             median_prediction + sigma_shade*std_prediction,
                             color='grey', alpha=0.3,
                             label=rf'$\pm{sigma_shade}\sigma$'
                             )
    ax2.plot(x, (y-model_best)/model_best,
             color='green', ls='',marker='.', lw=1.5, markersize=7, label='Residual')
    ax2.set_ylim(-2.0,2.0)

    
    # plt.ylim(1e-3,1.2*np.max(y))
    if add_fit_legend == True:
        ax1.legend(loc=(0.01, 0.01),
                handlelength=1,
                handletextpad=0.5,
                columnspacing=0.5,
                borderaxespad=0.1,
                frameon=True,
                prop={'size': 13}
                )
    if np.nanmax(y) < 5:
        ax1.set_ylim(0.1*np.nanmin(y),10.0*np.nanmax(y))
    else:
        ax1.set_ylim(0.1*np.nanmin(y),3.0*np.nanmax(y))
    ax2.set_xlabel(r'Frequency $\nu_{\mathrm{obs}}$ [GHz]')
    # plt.ylabel('Integrated Flux Density [mJy]')
    ax1.set_ylabel(r'$S_{\nu}(\nu_{\mathrm{obs}})$ [mJy]')
    text_x, text_y = 0.60, 0.37
    text = (r"$\overline{\alpha}"+rf"= {(result.params['alpha'].value):.2f}\pm "
            rf"{(result.params['alpha'].stderr):.2f}$")
    
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=13, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)
    if do_mcmc_fit == True:
        if mcmc_version == 'general':
            text_x, text_y = 0.60, 0.77
            text = (r"$\overline{\alpha}^{\rm MCMC}"+f"= {(param_dict['alpha']['best']):.2f}"
                    rf"_{{-{param_dict['alpha']['lower']:.2f}}}^{{+{param_dict['alpha']['upper']:.2f}}}$")
            text_bbox = plt.text(text_x, text_y, text,
                                # ha='center', va='center',
                                fontsize=13, color='black',
                                # bbox=text_bbox_props, 
                                transform=fig.transFigure)
    
    if title_text is not None:
        ax1.set_title(title_text)

    if log_plot == True:
        ax1.semilogx()
        ax1.semilogy()

    plt.subplots_adjust(hspace=0.05)
    ax2.legend(
        # loc=(0.01, 0.55),
        handlelength=1,
        handletextpad=0.5,
        columnspacing=0.5,
        borderaxespad=0.1,
        prop={'size': 13},frameon=True)
    # legend()

    if basename_save is not None:
        if save_name_append is None:
            save_name_append = '_RC_alpha_fit_linear'
        else:
            save_name_append = save_name_append + '_RC_alpha_fit_linear'
        fig.savefig(basename_save.replace('.fits','_')+save_name_append+ext_save, dpi=600,
                    bbox_inches='tight')
    if verbose > 0:
        plt.show()
        plt.clf()
        plt.close(fig)
    else:
        plt.clf()
        plt.close(fig)
    # plt.show()
    corner_kwargs = {
        'bins': 30,
        'hist_bin_factor': 1.0,
        'color': 'purple',  # Change the color of the scatter plots and histograms
        'hist_kwargs': {
            'color': 'green',      # Histogram color
            'edgecolor': 'black',  # Edge color of the bins
            'linewidth': 1.5       # Width of the bin edges
        },
        'scatter_kwargs': {
            'alpha': 0.6,          # Transparency of points
            's': 10,               # Size of points
            'color': 'purple',     # Color of points
            'edgecolor': 'none'    # Edge color of points
        },
        'contour_kwargs': {
            'colors': 'blue',       # Color of contour lines
            'linewidths': 1.2       # Width of contour lines
        }
    }
    if do_mcmc_fit == True:
        try:
            from scipy.stats import gaussian_kde
            # fig_c = plt.figure(figsize=(2,2))
            plt.figure()
            _ = corner.corner(samples_emcee,
                                    labels=[r'$S_{\rm \nu_0}$',r'$\alpha$'],
                                    truths=[result.params['A1'].value,
                                            result.params['alpha'].value],
                                    show_titles=True,
                                    quantiles=quantiles,
                                    # **corner_kwargs
                                    # fig=fig_c
                                    )
            print(samples_emcee.shape)
            if basename_save is not None:
                save_name_append_corner = save_name_append + '_corner'
                plt.savefig(basename_save.replace('.fits','_')+save_name_append_corner+'.jpg',
                            dpi=600,bbox_inches='tight')
            
            if verbose > 0:
                plt.show()
                plt.clf()
                plt.close()
            else:
                plt.clf()
                plt.close()
                
            print('++==>> Parameter Results (MCMC sampling).')
            print(lmfit.fit_report(results_emcee.params))
        except:
            pass
    else:
        samples_emcee = None
        param_dict = None
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))
    return(mini,result,param_dict)


def do_fit_spec_RC_curv(freqs,
                          fluxes,
                          fluxes_err,
                          basename_save=None,log_plot=True,
                          save_name_append = None,
                          plot_errors_shade = False,
                          do_mcmc_fit = False,
                          title_text = None,
                          verbose=0):
    """
    Peform a fit to the radio spectrum using a curved model.

    Parameters
    ----------
    freqs : array
        Frequency array in Hz.
    fluxes : array
        Flux density array in mJy.
    fluxes_err : array
        Flux density error array in mJy.
    basename_save : str
        Basename to save the output files.
    verbose : int
        Verbosity level.

    """
    x = freqs / 1e9
    nu0 = np.mean(x)
    y = fluxes
    yerr = fluxes_err
    S0_init = 2*np.max(y)

    def min_func(params):
        S0 = params['S0']
        alpha = params['alpha']
        q = params['q']
        model = RC_function_Sq(x, S0, alpha, q)
        # res = (y - model) / (y + np.sqrt((yerr) ** 2.0 + 0.001))  # okay
        # res = (y - model) / (y*(np.log(yerr)))  # okay
        # return res.copy()
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights

    fit_params = lmfit.Parameters()
    fit_params.add("S0", value=10.0, min=-5, max=5000)
    fit_params.add("alpha", value=-0.9, min=-5.0, max=5.0)
    fit_params.add("q", value=0.0, min=-0.9, max=0.9)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 5),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # fig = plt.figure(figsize=(8, 4))
    x_resample = np.linspace(np.min(x)*0.25, np.max(x)*1.3, 500)
    ax1.errorbar(x, y, yerr=yerr, fmt='o',
                 label='Data', color='k', ecolor='gray',
                 alpha=0.5)

    nwalkers = int(len(y) * 3 * 20)
    if do_mcmc_fit == True:
        burn_in = 500
        results_emcee = mini.emcee(burn=burn_in,
                                   steps=5000,
                                   thin=5, nwalkers=nwalkers,
                                   is_weighted=False,
                                   float_behavior='posterior',
                                   params=result.params, workers=-1)

        samples_emcee = results_emcee.flatchain[burn_in:]
        _S0 = np.asarray(samples_emcee['S0'])
        _alpha = np.asarray(samples_emcee['alpha'])
        _q = np.asarray(samples_emcee['q'])

        model_samples = np.array(
            [RC_function_Sq(x_resample, _S0[i], _alpha[i], _q[i]) for i in
             range(samples_emcee.shape[0])])

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = RC_function_Sq(x_resample,
                                    result.params['S0'].value,
                                    result.params['alpha'].value,
                                    result.params['q'].value)

    model_best = RC_function_Sq(x,
                                    result.params['S0'].value,
                                    result.params['alpha'].value,
                                    result.params['q'].value)


    ax1.plot(x_resample, model_resample,
             color='red', ls='-.', label='Best-fit model')
    ax2.plot(x, (y-model_best)/model_best,
             color='green', ls='dotted', label='Residual')
    ax2.set_xlabel(r'$\nu$ [GHz]')
    # plt.fill_between(x_resample, lower_bound, upper_bound, color='lightgray', alpha=0.5,
    #                  label='Parameter uncertainty')
    if plot_errors_shade:
        if do_mcmc_fit == True:
            ax1.plot(x_resample, model_mean,
                     color='purple', linestyle=(0, (5, 10)), label='MCMC Mean')
            ax1.fill_between(x_resample,
                             model_mean - 3*model_std,
                             model_mean + 3*model_std, color='lightgray',
                            alpha=0.5)
        else:
            # Define the number of Monte Carlo samples
            num_samples = 1000

            # Generate random samples from parameter distributions
            S0_samples = np.random.normal(result.params['S0'].value, result.params['S0'].stderr,
                                          num_samples)
            alpha_samples = np.random.normal(result.params['alpha'].value,
                                             result.params['alpha'].stderr,
                                             num_samples)
            q_samples = np.random.normal(result.params['q'].value, result.params['q'].stderr,
                                         num_samples)

            # Compute model predictions for each sample
            model_predictions = np.zeros((num_samples, len(x_resample)))
            for i in range(num_samples):
                model_predictions[i] = RC_function_Sq(x_resample, S0_samples[i],
                                                      alpha_samples[i],
                                                      q_samples[i])

            median_prediction = np.median(model_predictions, axis=0)
            std_prediction = np.std(model_predictions, axis=0)

            ax1.fill_between(x_resample,
                             median_prediction - 1*std_prediction,
                             median_prediction + 1*std_prediction,
                             color='lightgray', alpha=0.5,
                             # label='Uncertainty (1-sigma)'
                             )
    # plt.ylim(1e-3,1.2*np.max(y))
    ax1.legend()
    ax1.set_ylim(0.1*np.nanmin(y),3.0*np.nanmax(y))
    # plt.ylabel('Integrated Flux Density [mJy]')
    ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
    text_x, text_y = 0.65, 0.37
    text = (rf"$\alpha = {(result.params['alpha'].value):.2f}\pm "
            rf"{(result.params['alpha'].stderr):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = ax1.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)

    text_x, text_y = 0.48, 0.37
    text = (rf"$q = {(result.params['q'].value):.2f}\pm "
            rf"{(result.params['q'].stderr):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = ax1.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)

    if title_text is not None:
        ax1.set_title(title_text)

    if log_plot == True:
        ax1.semilogx()
        ax1.semilogy()
    plt.subplots_adjust(hspace=0.05)
    ax2.legend()
    if basename_save is not None:
        if save_name_append is None:
            save_name_append = '_RC_alpha_fit_curved'
        else:
            save_name_append = save_name_append + '_RC_alpha_fit_curved'
        plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                    bbox_inches='tight')
    plt.show()

    plt.figure()

    if do_mcmc_fit == True:
        _ = corner.corner(samples_emcee[['S0','alpha','q']],
                                labels=[r'$S_{0}$',
                                        r'$\alpha$',
                                        r'$q$'],
                                truths=[result.params['S0'].value,
                                        result.params['alpha'].value,
                                        result.params['q'].value],
                                show_titles=True,
                                quantiles=[0.025, 0.5, 0.975])
        plt.show()
        print('++==>> Parameter Results (MCMC sampling).')
        print(lmfit.fit_report(results_emcee.params))
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))


    return mini,result,samples_emcee

def RC_function_turnover(nu,nut1,A,B,alpha):
    tau = (nu/nut1)**(-2.1)
    return ((1-np.exp(-tau))*(B+A*(nu/nut1)**(0.1+alpha))*((nu/nut1)**2.0))

def do_fit_spec_RC_turnover(freqs,
                          fluxes,
                          fluxes_err,
                          basename_save=None,log_plot=True,
                          save_name_append = None,
                          plot_errors_shade = False,
                          do_mcmc_fit = False,
                          title_text = None,
                          verbose=0):
    """
    Peform a fit to the radio spectrum using a curved model.

    Parameters
    ----------
    freqs : array
        Frequency array in Hz.
    fluxes : array
        Flux density array in mJy.
    fluxes_err : array
        Flux density error array in mJy.
    basename_save : str
        Basename to save the output files.
    verbose : int
        Verbosity level.

    """
    x = freqs / 1e9
    y = fluxes
    yerr = fluxes_err

    def minimisation_function(data,model,error,delta=0.1):
        sigma2 = error ** 2
        # return -0.5 * np.sum((y - model) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))
        chi2 = np.sum((model - data)**2 / sigma2)
        # Calculate log likelihood
        log_likelihood = -0.5 * chi2
        return log_likelihood


    def min_func(params):
        nut1 = params['nut1']
        A = params['A']
        B = params['B']
        alpha = params['alpha']
        model = RC_function_turnover(x, nut1, A, B, alpha)
        # res = (y - model) / (y + np.sqrt((yerr) ** 2.0 + 0.001))  # okay
        log_likelihood = (y - model) / (y*(np.log(abs(yerr))))  # okay
        # res = (y - model) / (y * (np.sqrt(abs(yerr))))  # okay
        # sigma2 = yerr ** 2
        # chi2 = ((y - model)**2 / sigma2)
        # log_likelihood = -0.5 * chi2
        
        # log_likelihood=0.5 * ((y - model) ** 2 / sigma2 + np.log(2 * np.pi * sigma2))
        
        # res = (y - model) / (y*np.log(abs(yerr) + 1.0))  # okay 8
        return log_likelihood.copy()

    fit_params = lmfit.Parameters()
    fit_params.add("nut1", value=1.0, min=0.01, max=10)
    fit_params.add("alpha", value=-0.9, min=-2.0, max=-0.2)
    fit_params.add("A", value=10.0, min=0.00001, max=5000)
    fit_params.add("B", value=10.0, min=0.00001, max=5000)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 5),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # fig = plt.figure(figsize=(8, 4))
    x_resample = np.linspace(np.min(x)*0.25, np.max(x)*1.3, 500)
    ax1.errorbar(x, y, yerr=yerr, fmt='o',
                 label='Data', color='k', ecolor='gray',
                 alpha=0.5)

    nwalkers = int(4 * 20)
    if do_mcmc_fit == True:
        burn_in = 500
        results_emcee = mini.emcee(burn=burn_in,
                                   steps=5000,
                                   thin=5, nwalkers=nwalkers,
                                   is_weighted=False,
                                   float_behavior='posterior',
                                   params=result.params, workers=-1)

        samples_emcee = results_emcee.flatchain[burn_in:]
        _nut1 = np.asarray(samples_emcee['nut1'])
        _alpha = np.asarray(samples_emcee['alpha'])
        _A = np.asarray(samples_emcee['A'])
        _B = np.asarray(samples_emcee['B'])
        # RC_function_turnover(nu, nut1, A, B, alpha)
        model_samples = np.array(
            [RC_function_turnover(x_resample, _nut1[i], _A[i], _B[i], _alpha[i]) for i in
             range(samples_emcee.shape[0])])

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = RC_function_turnover(x_resample,
                                    result.params['nut1'].value,
                                    result.params['A'].value,
                                    result.params['B'].value,
                                    result.params['alpha'].value)

    model_best = RC_function_turnover(x,
                                    result.params['nut1'].value,
                                    result.params['A'].value,
                                    result.params['B'].value,
                                    result.params['alpha'].value)


    ax1.plot(x_resample, model_resample,
             color='red', ls='-.', label='Best-fit model')
    ax2.plot(x, (y-model_best)/model_best,
             color='green', ls='dotted', label='Residual')
    ax2.set_xlabel(r'$\nu$ [GHz]')
    # plt.fill_between(x_resample, lower_bound, upper_bound, color='lightgray', alpha=0.5,
    #                  label='Parameter uncertainty')
    if plot_errors_shade:
        if do_mcmc_fit == True:
            ax1.plot(x_resample, model_mean,
                     color='purple', linestyle=(0, (5, 10)), label='MCMC Mean')
            ax1.fill_between(x_resample,
                             model_mean - 3*model_std,
                             model_mean + 3*model_std, color='lightgray',
                            alpha=0.5)
        else:
            # Define the number of Monte Carlo samples
            num_samples = 1000

            # Generate random samples from parameter distributions
            nut1_samples = np.random.normal(result.params['nut1'].value, result.params['nut1'].stderr,
                                          num_samples)
            alpha_samples = np.random.normal(result.params['alpha'].value,
                                             result.params['alpha'].stderr,
                                             num_samples)
            A_samples = np.random.normal(result.params['A'].value, result.params['A'].stderr,
                                         num_samples)

            B_samples = np.random.normal(result.params['B'].value, result.params['B'].stderr,
                                         num_samples)

            # Compute model predictions for each sample
            model_predictions = np.zeros((num_samples, len(x_resample)))
            for i in range(num_samples):
                model_predictions[i] = RC_function_turnover(x_resample,
                                                            nut1_samples[i],
                                                            A_samples[i],
                                                            B_samples[i],
                                                            alpha_samples[i])

            median_prediction = np.median(model_predictions, axis=0)
            std_prediction = np.std(model_predictions, axis=0)

            ax1.fill_between(x_resample,
                             median_prediction - 1*std_prediction,
                             median_prediction + 1*std_prediction,
                             color='lightgray', alpha=0.5,
                             # label='Uncertainty (1-sigma)'
                             )
    # plt.ylim(1e-3,1.2*np.max(y))
    ax1.legend()
    ax1.set_ylim(0.1*np.nanmin(y),3.0*np.nanmax(y))
    # plt.ylabel('Integrated Flux Density [mJy]')
    ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
    text_x, text_y = 0.65, 0.37
    text = (rf"$\alpha = {(result.params['alpha'].value):.2f}\pm "
            rf"{(result.params['alpha'].stderr):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = ax1.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)

    text_x, text_y = 0.48, 0.37
    text = (rf"$\nu_t = {(result.params['nut1'].value):.2f}\pm "
            rf"{(result.params['nut1'].stderr):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = ax1.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)

    if title_text is not None:
        ax1.set_title(title_text)

    if log_plot == True:
        ax1.semilogx()
        ax1.semilogy()
    plt.subplots_adjust(hspace=0.05)
    ax2.legend()
    if basename_save is not None:
        if save_name_append is None:
            save_name_append = '_RC_alpha_fit_turnover'
        else:
            save_name_append = save_name_append + '_RC_alpha_fit_turnover'
        plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                    bbox_inches='tight')
    plt.show()

    plt.figure()

    if do_mcmc_fit == True:
        _ = corner.corner(samples_emcee[['nut1','A','B','alpha']],
                                labels=[r'$\nu_{t,1}$',
                                        r'$A$',
                                        r'$B$',
                                        r'$\alpha$'],
                                truths=[result.params['nut1'].value,
                                        result.params['A'].value,
                                        result.params['B'].value,
                                        result.params['alpha'].value],
                                show_titles=True,
                                quantiles=[0.025, 0.5, 0.975])
        plt.show()
        print('++==>> Parameter Results (MCMC sampling).')
        print(lmfit.fit_report(results_emcee.params))
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))


    return mini,result


    
    
def huber_loss(residuals, delta=1.0):
    return huber(delta, residuals)

def huber_like_weighting(yerr, k=1.345):
    return np.where(yerr > k, k / yerr, 1)


def fit_spec_SY_FF(freqs,
                   fluxes,
                   fluxes_err,
                   nu0=None,
                   nu_th=33.0,
                   fix_alpha_nt=False,
                   fix_alpha_nt_value=-0.85,
                   basename_save=None,log_plot=True,
                   save_name_append = None,
                   plot_errors_shade = False,
                   sigma_shade=1,
                   quantiles = [0.16, 0.5, 0.84],
                   do_mcmc_fit = False,
                   mcmc_version = 'general',
                   burn_in = 2000,
                   nsteps = 5000,
                   thin = 5,
                   sigma_errors = 1.0,
                   prior_sigma = 15.0,
                   title_text = None,
                   add_fit_legend = True,
                   plot_fit_results=True,
                   verbose=0):
    """
    Peform a fit to the radio spectrum using SY + FF model.

    Parameters
    ----------
    freqs : array
        Frequency array in Hz.
    fluxes : array
        Flux density array in mJy.
    fluxes_err : array
        Flux density error array in mJy.
    basename_save : str
        Basename to save the output files.
    verbose : int
        Verbosity level.

    """
    
    # if sigma_shade == 1 or sigma_shade == 1.0:
    #     quantiles = [0.16, 0.5, 0.84]
    # elif sigma_shade == 2 or sigma_shade == 2.0:
    #     quantiles = [0.025, 0.5, 0.975]
    # elif sigma_shade == 3 or sigma_shade == 3.0:
    #     quantiles = [0.00015, 0.5, 0.9985]
    # else:
    #     raise ValueError("Unsupported sigma value. Use sigma = 1, 2, or 3.")
    
    
    x = freqs / 1e9
    y = fluxes
    yerr = fluxes_err
    epsilon = 1e-6
    # weights = 1.0 / yerr
    if nu0 is None:
        nu0 = np.mean(x)
    print(f' ++==>> Using reference frequency of {nu0} GHz.')



        
    # def min_func(params):
    #     A_sy = params['A_sy']
    #     A_ff = params['A_ff']
    #     alpha_nt = params['alpha_nt']
    #     model = RC_function_SY_FF(x, A_sy, A_ff, alpha_nt, nu0)
    #     # res = (y - model) / (np.log(yerr))
        
    #     relative_error = yerr / (np.abs(y) + epsilon)
    #     weights = 1 / (relative_error + epsilon)
    #     # weightned_residual = (y - model) * np.sqrt(weights)
    #     weightned_residual = (y - model) / (np.sqrt(weights)*yerr)
    #     # weightned_residual = (y - model) / (np.sqrt(weights) + yerr)
    #     # weightned_residual = (y - model) * np.sqrt(2 / np.log1p(yerr**2))
        
    #     return weightned_residual.copy()

    def min_func(params):
        A_sy = params['A_sy']
        A_ff = params['A_ff']
        alpha_nt = params['alpha_nt']
        model = RC_function_SY_FF(x, A_sy, A_ff, alpha_nt, nu0)
        # res = (y - model) / (np.log(yerr))
        # weights = 1 / yerr
        # return (y - model) * np.sqrt(weights)
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights


    # def min_func(params):
    #     """
    #     Adaptive weighting to balance low and high-frequency data.
    #     """
    #     A_sy = params['A_sy']
    #     A_ff = params['A_ff']
    #     alpha_nt = params['alpha_nt']
    #     model = RC_function_SY_FF(x, A_sy, A_ff, alpha_nt, nu0)
                
    #     # Log-based weighting (same as before)
    #     log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        
    #     # Frequency correction factor
    #     freq_weights = (x / np.max(x))**0.5  # Adjust exponent as needed

    #     # Combine weights
    #     final_weights = log_weights * freq_weights

    #     return (y - model) * final_weights




    fit_params = lmfit.Parameters()
    fit_params.add("A_sy", value=0.5, min=0.01, max=5000)
    fit_params.add("A_ff", value=0.5, min=0.01, max=5000)
    if fix_alpha_nt == True:
        fit_params.add("alpha_nt", value=fix_alpha_nt_value, 
                       min=-2.0, max=0.0, vary=False)
    else:
        fit_params.add("alpha_nt", value=-0.7, min=-2.5, max=2.5)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))

    x_resample = np.linspace(np.min(x)*0.7, np.max(x)*1.3, 500)
    
    
    if do_mcmc_fit == True:
        nwalkers = int(3 * 25)
        
        if mcmc_version == 'lmfit':
            """ 
            This will be removed in the future.
            """
            results_emcee = mini.emcee(burn=burn_in,
                                    steps=nsteps,
                                    thin=thin, nwalkers=nwalkers,
                                    is_weighted=False,
                                    float_behavior='posterior',
                                    params=result.params, workers=-1)

            _samples_emcee = results_emcee.flatchain[burn_in:]
            _A_sy = np.asarray(_samples_emcee['A_sy'])
            _A_ff = np.asarray(_samples_emcee['A_ff'])
            _alpha_nt = np.asarray(_samples_emcee['alpha_nt'])
            model_samples = np.array(
                [RC_function_SY_FF(x_resample, _A_sy[i], _A_ff[i], _alpha_nt[i], nu0) for i in
                range(_samples_emcee.shape[0])])
            samples_emcee = np.asarray(_samples_emcee[['A_sy', 'A_ff', 'alpha_nt']])
            param_dict = {}
            for i, label in enumerate(['A_sy', 'A_ff','alpha_nt']):
                q = np.percentile(samples_emcee[:, i], [2.5, 50, 97.5])
                param_dict[label] = {
                    'best': q[1],
                    'lower': (q[1] - q[0]),
                    'upper': (q[2] - q[1])
                }
            print(param_dict)
        if mcmc_version == 'general':
            """
            This seems to be more robust as we have more control over
            the MCMC process.
            """
            def RC_SY_FF(params,x):
                return (params[1] * ((x/nu0)**(-0.1)) + params[0] * ((x/nu0)**params[2]))
                
            samples_emcee, param_dict = general_mcmc(x_data = x, 
                                                    y_data = y, 
                                                    yerr_data = yerr, 
                                                    result_1 = result, 
                                                    param_labels = ["A_sy", "A_ff", 'alpha_nt'], 
                                                    burn_in = burn_in,
                                                    nwalkers = nwalkers,
                                                    nsteps = nsteps,
                                                    thin = thin,
                                                    prior_sigma=prior_sigma,
                                                    sigma_errors = sigma_errors,
                                                    model_func = RC_SY_FF,
                                                    quantiles=quantiles)
            _A_sy = samples_emcee.T[0]
            _A_ff = samples_emcee.T[1]
            _alpha_nt = samples_emcee.T[2]
            # print(_A_sy.shape[0])

            model_samples = np.array(
                [RC_function_SY_FF(x_resample, 
                                    _A_sy[i], _A_ff[i], 
                                    _alpha_nt[i], nu0) for i in range(_A_sy.shape[0])])
            
            Ssy_samples = np.array(
                [RC_function_SY_FF(x_resample, 
                                    _A_sy[i], _A_ff[i]*0, 
                                    _alpha_nt[i], nu0) for i in range(_A_sy.shape[0])])
            Sff_samples = np.array(
                [RC_function_SY_FF(x_resample, 
                                    _A_sy[i]*0, _A_ff[i], 
                                    _alpha_nt[i], nu0) for i in range(_A_sy.shape[0])])
            
            print(param_dict)

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)
        Ssy_model_mean = np.mean(Ssy_samples, axis=0)
        Ssy_model_std = np.std(Ssy_samples, axis=0)
        Sff_model_mean = np.mean(Sff_samples, axis=0)
        Sff_model_std = np.std(Sff_samples, axis=0)

    model_resample = RC_function_SY_FF(x_resample,
                                    result.params['A_sy'].value,
                                    result.params['A_ff'].value,
                                    result.params['alpha_nt'].value,
                                    nu0)

    model_best = RC_function_SY_FF(x,
                                    result.params['A_sy'].value,
                                    result.params['A_ff'].value,
                                    result.params['alpha_nt'].value,
                                    nu0)



    A_sy_term = RC_function_SY_FF(x_resample,
                            result.params['A_sy'].value,
                            result.params['A_ff'].value*0,
                            result.params['alpha_nt'].value,
                            nu0)
    
    
    A_ff_term = RC_function_SY_FF(x_resample,
                            result.params['A_sy'].value*0,
                            result.params['A_ff'].value,
                            result.params['alpha_nt'].value,
                            nu0)
    
    A_ff_term_pred = RC_function_SY_FF(x,
                            result.params['A_sy'].value*0,
                            result.params['A_ff'].value,
                            result.params['alpha_nt'].value,
                            nu0)
    
    
    A_sy_term_nu_th = RC_function_SY_FF(nu_th,
                            result.params['A_sy'].value,
                            result.params['A_ff'].value*0,
                            result.params['alpha_nt'].value,
                            nu0)
    
    A_ff_term_nu_th = RC_function_SY_FF(nu_th,
                            result.params['A_sy'].value*0,
                            result.params['A_ff'].value,
                            result.params['alpha_nt'].value,
                            nu0)
    
    
    thermal_fraction_nu_th = A_ff_term_nu_th / (A_sy_term_nu_th + A_ff_term_nu_th)
    term_1 = (result.params['A_sy'].value * result.params['A_ff'].stderr) ** 2.0
    term_2 = (result.params['A_ff'].value * result.params['A_sy'].stderr) ** 2.0
    term_3 = (result.params['A_ff'].value + result.params['A_sy'].value) ** 4.0
    # term_1 = (result.params['A_sy'].value * (nu_th**(-0.1)) * (nu_th**(result.params['alpha_nt'].value)) * result.params['A_ff'].stderr) ** 2.0
    # term_2 = (result.params['A_ff'].value * (nu_th**(-0.1)) * (nu_th**(result.params['alpha_nt'].value)) * result.params['A_sy'].stderr) ** 2.0
    # term_3 = (result.params['A_ff'].value * (nu_th**(-0.1)) + result.params['A_sy'].value * (nu_th**(result.params['alpha_nt'].value))) ** 4.0
    
    thermal_fraction_nu_th_err = np.sqrt((term_1 + term_2) / term_3)
    
    thermal_fraction_freq = A_ff_term_pred / model_best
    
    thermal_fraction = {"thermal_fraction_freq": thermal_fraction_freq,
                        "thermal_fraction_freq_err": thermal_fraction_nu_th_err}

    if plot_fit_results:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 4),
                                    gridspec_kw={'height_ratios': [3, 1]})
        
        ax1.errorbar(x, y, yerr=yerr, 
                    fmt='o',label='Data', color='k', ecolor='gray',alpha=0.5)

        ax1.plot(x_resample,model_resample,
                color='red', ls='-.', label='Fit')
        ax2.plot(x, (y-model_best)/model_best,
                color='green', ls='dotted', label='Residual')
        ax2.set_ylim(-1.0,1.0)

        # if do_mcmc_fit == True:
        #     ax1.plot(x_resample, Ssy_model_mean,
        #             '-', label='sy',color='violet')
            
        #     ax1.plot(x_resample, Sff_model_mean,
        #             '-', label='ff', color='orange')
        # else:
        ax1.plot(x_resample, A_sy_term,
                '-', label='sy')
        
        ax1.plot(x_resample, A_ff_term,
                '-', label='ff')

        if do_mcmc_fit == True:
            ax1.plot(x_resample, model_mean,
                    color='purple', linestyle=(0, (5, 10)), label='MCMC')

        # plt.fill_between(x_resample, lower_bound, upper_bound, color='lightgray', alpha=0.5,
        #                  label='Parameter uncertainty')
        if plot_errors_shade:
            if do_mcmc_fit == True:
                ax1.fill_between(x_resample,
                                model_mean - sigma_shade*model_std,
                                model_mean + sigma_shade*model_std, 
                                color='lightgray',
                                alpha=0.7)
                # ax1.fill_between(x_resample,
                #                 Ssy_model_mean - sigma_shade*Ssy_model_std,
                #                 Ssy_model_mean + sigma_shade*Ssy_model_std,
                #                 color='violet',
                #                 alpha=0.3)
                # ax1.fill_between(x_resample,
                #                 Sff_model_mean - sigma_shade*Sff_model_std,
                #                 Sff_model_mean + sigma_shade*Sff_model_std,
                #                 color='orange',
                #                 alpha=0.3)
                
            else:
                # Define the number of Monte Carlo samples
                num_samples = 5000

                # Generate random samples from parameter distributions
                A_sy_samples = np.random.normal(result.params['A_sy'].value, 
                                                result.params['A_sy'].stderr,
                                                num_samples)
                A_ff_samples = np.random.normal(result.params['A_ff'].value, 
                                                result.params['A_ff'].stderr,
                                                num_samples)
                alpha_nt_samples = np.random.normal(result.params['alpha_nt'].value,
                                                result.params['alpha_nt'].stderr,
                                                num_samples)

                # Compute model predictions for each sample
                model_predictions = np.zeros((num_samples, len(x_resample)))
                sy_model_predictions = np.zeros((num_samples, len(x_resample)))
                ff_model_predictions = np.zeros((num_samples, len(x_resample)))
                
                for i in range(num_samples):
                    model_predictions[i] = RC_function_SY_FF(x_resample, 
                                                            A_sy_samples[i],
                                                            A_ff_samples[i],
                                                            alpha_nt_samples[i],
                                                            nu0)
                    sy_model_predictions[i] = RC_function_SY_FF(x_resample, 
                                                            A_sy_samples[i],
                                                            A_ff_samples[i]*0,
                                                            alpha_nt_samples[i],
                                                            nu0)
                    ff_model_predictions[i] = RC_function_SY_FF(x_resample, 
                                                            A_sy_samples[i]*0,
                                                            A_ff_samples[i],
                                                            alpha_nt_samples[i],
                                                            nu0)

                median_prediction = np.median(model_predictions, axis=0)
                std_prediction = np.std(model_predictions, axis=0)
                sy_median_prediction = np.median(sy_model_predictions, axis=0)
                sy_std_prediction = np.std(sy_model_predictions, axis=0)
                ff_median_prediction = np.median(ff_model_predictions, axis=0)
                ff_std_prediction = np.std(ff_model_predictions, axis=0)
                
                ax1.fill_between(x_resample, median_prediction - sigma_shade*std_prediction,
                                median_prediction + sigma_shade*std_prediction,
                                color='lightgray', alpha=0.3,
                                # label='Uncertainty (1-sigma)'
                                )
                ax1.fill_between(x_resample, sy_median_prediction - sigma_shade*sy_std_prediction,
                                sy_median_prediction + sigma_shade*sy_std_prediction,
                                color='violet', alpha=0.3,
                                # label='Uncertainty (1-sigma)'
                                )
                ax1.fill_between(x_resample, ff_median_prediction - sigma_shade*ff_std_prediction,
                                ff_median_prediction + sigma_shade*ff_std_prediction,
                                color='orange', alpha=0.3,
                                # label='Uncertainty (1-sigma)'
                                )
                
        # plt.ylim(1e-3,1.2*np.max(y))
        if add_fit_legend == True:
            ax1.legend(loc=(0.05, 0.05),frameon=True,
                    #    prop={'size': 11},
                       ncol=2)
        if np.nanmax(y) < 5:
            ax1.set_ylim(0.1*np.nanmin(y),10.0*np.nanmax(y))
        else:
            ax1.set_ylim(0.1*np.nanmin(y),3.0*np.nanmax(y))
        ax2.set_xlabel(r'$\nu$ [GHz]')
        # plt.ylabel('Integrated Flux Density [mJy]')
        ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
        
        text_x, text_y = 0.62, 0.72
        if fix_alpha_nt == True:
            text = (r"$\alpha_{\rm sy}"+f"= {(result.params['alpha_nt'].value):.2f}$")
        else:
            text = (r"$\alpha_{\rm sy}"+rf"= {(result.params['alpha_nt'].value):.2f}\pm "
                    rf"{(result.params['alpha_nt'].stderr):.2f}$")
        text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
        text_bbox = plt.text(text_x, text_y, text,
                            # ha='center', va='center',
                            # fontsize=12, 
                            color='black',
                            bbox=text_bbox_props, transform=fig.transFigure)
        
        text_x, text_y = 0.40, 0.84
        text =r"$f_{\rm th}$"rf"$({nu_th}) = {(thermal_fraction_nu_th):.2f}\pm{(thermal_fraction_nu_th_err):.2f}$"
        text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
        text_bbox = plt.text(text_x, text_y, text,
                            # ha='center', va='center',
                            # fontsize=12, 
                            color='black',
                            # bbox=text_bbox_props, 
                            transform=fig.transFigure)
        
        
        if do_mcmc_fit == True:
            if mcmc_version == 'general':
                text_x, text_y = 0.62, 0.80
                text = (r"$\alpha_{\rm sy}^{\rm MCMC}"+f"= {(param_dict['alpha_nt']['best']):.2f}"
                        rf"_{{-{param_dict['alpha_nt']['lower']:.2f}}}^{{+{param_dict['alpha_nt']['upper']:.2f}}}$")
                text_bbox = plt.text(text_x, text_y, text,
                                    # ha='center', va='center',
                                    # fontsize=12, 
                                    color='black',
                                    # bbox=text_bbox_props, 
                                    transform=fig.transFigure)
        
        if title_text is not None:
            ax1.set_title(title_text)

        if log_plot == True:
            ax1.semilogx()
            ax1.semilogy()

        plt.subplots_adjust(hspace=0.05)
        ax2.legend(loc=(0.7, 0.55),
                #    prop={'size': 14},
                   frameon=False,ncol=2)
        # legend()

        if basename_save is not None:
            if save_name_append is None:
                save_name_append = '_RC_sy_ff_fit'
            else:
                save_name_append = save_name_append + '_RC_sy_ff_fit'
            plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                        bbox_inches='tight')
        # plt.show()
        corner_kwargs = {
            'bins': 30,
            'hist_bin_factor': 1.0,
            'color': 'purple',  # Change the color of the scatter plots and histograms
            'hist_kwargs': {
                'color': 'green',      # Histogram color
                'edgecolor': 'black',  # Edge color of the bins
                'linewidth': 1.5       # Width of the bin edges
            },
            'scatter_kwargs': {
                'alpha': 0.6,          # Transparency of points
                's': 10,               # Size of points
                'color': 'purple',     # Color of points
                'edgecolor': 'none'    # Edge color of points
            },
            'contour_kwargs': {
                'colors': 'blue',       # Color of contour lines
                'linewidths': 1.2       # Width of contour lines
            }
        }
        if do_mcmc_fit == True:
            try:
                from scipy.stats import gaussian_kde
                # fig_c = plt.figure(figsize=(2,2))
                fig_c = plt.figure()
                if fix_alpha_nt == True:
                    _ = corner.corner(samples_emcee[:,:2],
                                            labels=[r'$A_{\rm sy}$',
                                                    r'$A_{\rm ff}$'
                                                    ],
                                            truths=[result.params['A_sy'].value,
                                                    result.params['A_ff'].value
                                                    ],
                                            show_titles=True,
                                            quantiles=quantiles,
                                            # **corner_kwargs
                                            # fig=fig_c
                                            )
                else:
                    _ = corner.corner(samples_emcee,
                                            labels=[r'$A_{\rm sy}$',
                                                    r'$A_{\rm ff}$',
                                                    r'$\alpha_{\rm nt}$'],
                                            truths=[result.params['A_sy'].value,
                                                    result.params['A_ff'].value,
                                                    result.params['alpha_nt'].value],
                                            show_titles=True,
                                            quantiles=quantiles,
                                            # **corner_kwargs
                                            # fig=fig_c
                                            )

                print(samples_emcee.shape)
                if basename_save is not None:
                    save_name_append_corner = save_name_append + '_corner'
                    plt.savefig(basename_save.replace('.fits','_')+save_name_append_corner+'.jpg', 
                                dpi=600,bbox_inches='tight')
                    
                # print('++==>> Parameter Results (MCMC sampling).')
                # print(lmfit.fit_report(results_emcee.params))
            except:
                pass
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))
    if not do_mcmc_fit:
        # return(mini,result,thermal_fraction,samples_emcee, param_dict)
    # else:
        samples_emcee = None
        param_dict = None
    # return(mini,result,thermal_fraction)
    return(mini,result,thermal_fraction,samples_emcee, param_dict)



def fit_spec_SY_FF_FFA(freqs,
                   fluxes,
                   fluxes_err,
                   nu0=None,
                   basename_save=None,log_plot=True,
                   save_name_append = None,
                   plot_errors_shade = False,
                   sigma_shade=3,
                   quantiles = [0.16, 0.5, 0.84],
                   do_mcmc_fit = False,
                   mcmc_version = 'general',
                   burn_in = 1000,
                   nsteps = 5000,
                   thin = 5,
                   sigma_errors = 1.0,
                   title_text = None,
                   add_fit_legend = True,
                   plot_fit_results=True,
                   verbose=0):
    """
    Peform a fit to the radio spectrum using SY + FF model in combination with 
    free-free absorption.

    Parameters
    ----------
    freqs : array
        Frequency array in Hz.
    fluxes : array
        Flux density array in mJy.
    fluxes_err : array
        Flux density error array in mJy.
    basename_save : str
        Basename to save the output files.
    verbose : int
        Verbosity level.

    """
    
    # if sigma_shade == 1 or sigma_shade == 1.0:
    #     quantiles = [0.16, 0.5, 0.84]
    # elif sigma_shade == 2 or sigma_shade == 2.0:
    #     quantiles = [0.025, 0.5, 0.975]
    # elif sigma_shade == 3 or sigma_shade == 3.0:
    #     quantiles = [0.00015, 0.5, 0.9985]
    # else:
    #     raise ValueError("Unsupported sigma value. Use sigma = 1, 2, or 3.")
    
    
    x = freqs / 1e9
    y = fluxes
    yerr = fluxes_err
    epsilon = 1e-8
    if nu0 is None:
        nu0 = np.mean(x)
    print(f' ++==>> Using reference frequency of {nu0} GHz.')



        
    # def min_func(params):
    #     Snu0 = params['Snu0']
    #     fth_nu0 = params['fth_nu0']
    #     alpha_nt = params['alpha_nt']
    #     nu_tau_t = params['nu_tau_t']
    #     # f_cov = params['f_cov']
    #     model = RC_function_SY_FF_ffa(x, Snu0, fth_nu0, alpha_nt, nu_tau_t, nu0)
    #     # res = (y - model) / (np.log(yerr))
    #     relative_error = yerr / (np.abs(y) + epsilon)
    #     weights = 1 / (relative_error + epsilon)
    #     # weightned_residual = (y - model) * np.sqrt(weights)
    #     weightned_residual = (y - model) / (np.sqrt(weights + yerr))
        
    #     # res = (y - model) * np.sqrt(2 / np.log1p(yerr**2))
    #     # res = (y - model)/yerr
    #     # loss = huber_loss(res)
        
    #     # return res.copy()
    #     return weightned_residual.copy()

    def min_func(params):
        Snu0 = params['Snu0']
        fth_nu0 = params['fth_nu0']
        alpha_nt = params['alpha_nt']
        nu_tau_t = params['nu_tau_t']
        # f_cov = params['f_cov']
        model = RC_function_SY_FF_ffa(x, Snu0, fth_nu0, alpha_nt, nu_tau_t, nu0)
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights


    fit_params = lmfit.Parameters()
    
    Snu0_init = np.log10(np.nanmean(y*1.0))
    if nu0 < 2.0:
        fth0_init = 0.10
    elif 2.0 <=nu0 < 10.0:
        fth0_init = 0.20
    elif 10.0 <= nu0 < 20.0:
        fth0_init = 0.30
    elif nu0 > 20.0:
        fth0_init = 0.5
    a_nth_init = np.polyfit(np.log10(x), 
                                    np.log10(y), 1)[0]
    
    print('Init Snu0=',Snu0_init)
    print('Init fth0=',fth0_init)
    print('Init a_nth=',a_nth_init)
    
    fit_params.add("Snu0", value=Snu0_init, min=0.5, max=5000)
    fit_params.add("fth_nu0", value=fth0_init, min=0.001, max=1.0)
    fit_params.add("alpha_nt", value=a_nth_init, min=-2.5, max=2.5)
    fit_params.add("nu_tau_t", value=1.0, min=0.1, max=10.0)
    # fit_params.add("f_cov", value=0.99, min=0.0, max=1.0)
    # fit_params.add("f_cov", value=1.0, min=0.9999, max=1.0)
    # fit_params.add("f_cov", value=1.0, vary=False)
    

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))

    x_resample = np.linspace(np.min(x)*0.7, np.max(x)*1.3, 500)
    
    
    if do_mcmc_fit == True:
        nwalkers = int(5 * 25)
        
        if mcmc_version == 'lmfit':
            """ 
            This will be removed in the future.
            """
            results_emcee = mini.emcee(burn=burn_in,
                                    steps=nsteps,
                                    thin=thin, nwalkers=nwalkers,
                                    is_weighted=False,
                                    float_behavior='posterior',
                                    params=result.params, workers=-1)

            _samples_emcee = results_emcee.flatchain[burn_in:]
            _Snu0 = np.asarray(_samples_emcee['Snu0'])
            _fth_nu0 = np.asarray(_samples_emcee['fth_nu0'])
            _alpha_nt = np.asarray(_samples_emcee['alpha_nt'])
            _nu_tau_t = np.asarray(_samples_emcee['nu_tau_t'])
            # _f_cov = np.asarray(_samples_emcee['f_cov'])
            model_samples = np.array(
                [RC_function_SY_FF_ffa(x_resample, _Snu0[i], _fth_nu0[i], _alpha_nt[i], _nu_tau_t, nu0) for i in
                range(_samples_emcee.shape[0])])
            samples_emcee = np.asarray(_samples_emcee[['Snu0', 'fth_nu0', 'alpha_nt', 'nu_tau_t']])
            param_dict = {}
            for i, label in enumerate(['Snu0', 'fth_nu0','alpha_nt', 'nu_tau_t']):
                q = np.percentile(samples_emcee[:, i], [2.5, 50, 97.5])
                param_dict[label] = {
                    'best': q[1],
                    'lower': (q[1] - q[0]),
                    'upper': (q[2] - q[1])
                }
            print(param_dict)
        if mcmc_version == 'general':
            """
            This seems to be more robust as we have more control over
            the MCMC process.
            """
            def RC_SY_FF_FFA(params,x):
                tau_nu = (x/params[3])**(-2.1)
                S_ff_abs = params[1] * params[0] * ((1-np.exp(-tau_nu))/(tau_nu)) * ((x/nu0)**(-0.1))
                # S_sy_abs = (1 - params[1]) * params[0] * (1-params[4]*(1-np.exp(-tau_nu))) * ((x/nu0)**params[2])
                S_sy_abs = (1 - params[1]) * params[0] * (1-1*(1-np.exp(-tau_nu))) * ((x/nu0)**params[2])
                S_total_abs = S_ff_abs + S_sy_abs
                return S_total_abs
                
            samples_emcee, param_dict = general_mcmc(x_data = x, 
                                                    y_data = y, 
                                                    yerr_data = yerr, 
                                                    result_1 = result, 
                                                    param_labels = ["Snu0", "fth_nu0", 'alpha_nt', 'nu_tau_t'], 
                                                    burn_in = burn_in,
                                                    nwalkers = nwalkers,
                                                    nsteps = nsteps,
                                                    thin = thin,
                                                    sigma_errors = sigma_errors,
                                                    model_func = RC_SY_FF_FFA,
                                                    quantiles=quantiles)
            _Snu0 = samples_emcee.T[0]
            _fth_nu0 = samples_emcee.T[1]
            _alpha_nt = samples_emcee.T[2]
            _nu_tau_t = samples_emcee.T[3]
            # _f_cov = samples_emcee.T[4]

            model_samples = np.array(
                [RC_function_SY_FF_ffa(x_resample, 
                                    _Snu0[i], _fth_nu0[i], _alpha_nt[i], 
                                    _nu_tau_t[i], nu0) for i in range(_Snu0.shape[0])])
            print(param_dict)

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = RC_function_SY_FF_ffa(x_resample,
                                        result.params['Snu0'].value,
                                        result.params['fth_nu0'].value,
                                        result.params['alpha_nt'].value,
                                        result.params['nu_tau_t'].value,
                                        # result.params['f_cov'].value,
                                    nu0)

    model_best = RC_function_SY_FF_ffa(x,
                                        result.params['Snu0'].value,
                                        result.params['fth_nu0'].value,
                                        result.params['alpha_nt'].value,
                                        result.params['nu_tau_t'].value,
                                        # result.params['f_cov'].value,
                                    nu0)



    Sy_term = RC_function_SY_ffa(x_resample,
                            result.params['Snu0'].value,
                            result.params['fth_nu0'].value,
                            result.params['alpha_nt'].value,
                            result.params['nu_tau_t'].value,
                            # result.params['f_cov'].value,
                            nu0)
    
    
    ff_term = RC_function_FF_ffa(x_resample,
                            result.params['Snu0'].value,
                            result.params['fth_nu0'].value,
                            result.params['nu_tau_t'].value,
                            nu0)
    
    ff_term_pred = RC_function_FF_ffa(x,
                            result.params['Snu0'].value,
                            result.params['fth_nu0'].value,
                            result.params['nu_tau_t'].value,
                            nu0)
    
    thermal_fraction_freq = ff_term_pred / model_best

    if plot_fit_results:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 4),
                                    gridspec_kw={'height_ratios': [3, 1]})
        
        ax1.errorbar(x, y, yerr=yerr, 
                    fmt='o',label='Data', color='k', ecolor='gray',alpha=0.5)

        ax1.plot(x_resample,model_resample,
                color='red', ls='-.', label='Fit')
        ax2.plot(x, (y-model_best)/model_best,
                color='green', ls='dotted', label='Residual')
        ax2.set_ylim(-1.0,1.0)


        ax1.plot(x_resample, Sy_term,
                '-', label='sy')
        
        ax1.plot(x_resample, ff_term,
                '-', label='ff')

        if do_mcmc_fit == True:
            ax1.plot(x_resample, model_mean,
                    color='purple', linestyle=(0, (5, 10)), label='MCMC')

        # plt.fill_between(x_resample, lower_bound, upper_bound, color='lightgray', alpha=0.5,
        #                  label='Parameter uncertainty')
        if plot_errors_shade:
            if do_mcmc_fit == True:
                ax1.fill_between(x_resample,
                                model_mean - sigma_shade*model_std,
                                model_mean + sigma_shade*model_std, color='lightgray',
                                alpha=0.7)
            else:
                # Define the number of Monte Carlo samples
                num_samples = 5000

                # Generate random samples from parameter distributions
                Snu0_samples = np.random.normal(result.params['Snu0'].value, 
                                                result.params['Snu0'].stderr,
                                                num_samples)
                fth_nu0_samples = np.random.normal(result.params['fth_nu0'].value, 
                                                result.params['fth_nu0'].stderr,
                                                num_samples)
                alpha_nt_samples = np.random.normal(result.params['alpha_nt'].value,
                                                result.params['alpha_nt'].stderr,
                                                num_samples)
                nu_tau_t_samples = np.random.normal(result.params['nu_tau_t'].value,
                                                result.params['nu_tau_t'].stderr,
                                                num_samples)
                # f_cov_samples = np.random.normal(result.params['f_cov'].value,
                #                                 result.params['f_cov'].stderr,
                #                                 num_samples)

                # Compute model predictions for each sample
                model_predictions = np.zeros((num_samples, len(x_resample)))
                for i in range(num_samples):
                    model_predictions[i] = RC_function_SY_FF_ffa(x_resample, 
                                                            Snu0_samples[i],
                                                            fth_nu0_samples[i],
                                                            alpha_nt_samples[i],
                                                            nu_tau_t_samples[i],
                                                            # f_cov_samples[i],
                                                            nu0)

                median_prediction = np.median(model_predictions, axis=0)
                std_prediction = np.std(model_predictions, axis=0)

                ax1.fill_between(x_resample, median_prediction - sigma_shade*std_prediction,
                                median_prediction + sigma_shade*std_prediction,
                                color='lightgray', alpha=0.5,
                                # label='Uncertainty (1-sigma)'
                                )
        # plt.ylim(1e-3,1.2*np.max(y))
        if add_fit_legend == True:
            ax1.legend(loc=(0.05, 0.05),frameon=True,prop={'size': 11},ncol=2)
        if np.nanmax(y) < 5:
            ax1.set_ylim(0.1*np.nanmin(y),10.0*np.nanmax(y))
        else:
            ax1.set_ylim(0.1*np.nanmin(y),3.0*np.nanmax(y))
        ax2.set_xlabel(r'$\nu$ [GHz]')
        # plt.ylabel('Integrated Flux Density [mJy]')
        ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
        
        text_x, text_y = 0.62, 0.72
        text = (r"$\alpha_{\rm nt}"+rf"= {(result.params['alpha_nt'].value):.2f}\pm "
                rf"{(result.params['alpha_nt'].stderr):.2f}$")
        text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
        text_bbox = plt.text(text_x, text_y, text,
                            # ha='center', va='center',
                            fontsize=12, color='black',
                            bbox=text_bbox_props, transform=fig.transFigure)
        
        if do_mcmc_fit == True:
            if mcmc_version == 'general':
                text_x, text_y = 0.62, 0.82
                text = (r"$\alpha_{\rm nt}^{\rm MCMC}"+f"= {(param_dict['alpha_nt']['best']):.2f}"
                        rf"_{{-{param_dict['alpha_nt']['lower']:.2f}}}^{{+{param_dict['alpha_nt']['upper']:.2f}}}$")
                text_bbox = plt.text(text_x, text_y, text,
                                    # ha='center', va='center',
                                    fontsize=12, color='black',
                                    # bbox=text_bbox_props, 
                                    transform=fig.transFigure)
        
        if title_text is not None:
            ax1.set_title(title_text)

        if log_plot == True:
            ax1.semilogx()
            ax1.semilogy()

        plt.subplots_adjust(hspace=0.05)
        ax2.legend(loc=(0.7, 0.55),prop={'size': 14},frameon=False,ncol=2)
        # legend()

        if basename_save is not None:
            if save_name_append is None:
                save_name_append = '_RC_sy_ff_ffa_fit'
            else:
                save_name_append = save_name_append + '_RC_sy_ff_ffa_fit'
            plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                        bbox_inches='tight')
        # plt.show()
        corner_kwargs = {
            'bins': 30,
            'hist_bin_factor': 1.0,
            'color': 'purple',  # Change the color of the scatter plots and histograms
            'hist_kwargs': {
                'color': 'green',      # Histogram color
                'edgecolor': 'black',  # Edge color of the bins
                'linewidth': 1.5       # Width of the bin edges
            },
            'scatter_kwargs': {
                'alpha': 0.6,          # Transparency of points
                's': 10,               # Size of points
                'color': 'purple',     # Color of points
                'edgecolor': 'none'    # Edge color of points
            },
            'contour_kwargs': {
                'colors': 'blue',       # Color of contour lines
                'linewidths': 1.2       # Width of contour lines
            }
        }
        if do_mcmc_fit == True:
            try:
                from scipy.stats import gaussian_kde
                # fig_c = plt.figure(figsize=(2,2))
                fig_c = plt.figure()
                _ = corner.corner(samples_emcee,
                                        labels=[r'$S_{\nu_{0}}$',
                                                r'$f_{\rm th}^{\nu_0}}$',
                                                r'$\alpha_{\rm nt}$',
                                                r'$\nu_{\rm t}$',
                                                r'$\Xi$'
                                                ],
                                        truths=[result.params['Snu0'].value,
                                                result.params['fth_nu0'].value,
                                                result.params['alpha_nt'].value,
                                                result.params['nu_tau_t'].value,
                                                # result.params['f_cov'].value
                                                ],
                                        show_titles=True,
                                        quantiles=quantiles,
                                        # **corner_kwargs
                                        # fig=fig_c
                                        )
                print(samples_emcee.shape)
                if basename_save is not None:
                    save_name_append_corner = save_name_append + '_corner'
                    plt.savefig(basename_save.replace('.fits','_')+save_name_append_corner+'.jpg', 
                                dpi=600,bbox_inches='tight')
                    
                print('++==>> Parameter Results (MCMC sampling).')
                print(lmfit.fit_report(results_emcee.params))
            except:
                pass
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))
    if do_mcmc_fit:
        return(mini,result,thermal_fraction_freq,samples_emcee, param_dict)
    else:
        return(mini,result,thermal_fraction_freq)



def fit_spec_SY_FF_FFA_v2(freqs,
                   fluxes,
                   fluxes_err,
                   nu0=None,
                   basename_save=None,log_plot=True,
                   save_name_append = None,
                   plot_errors_shade = False,
                   sigma_shade=3,
                   quantiles = [0.16, 0.5, 0.84],
                   do_mcmc_fit = False,
                   mcmc_version = 'general',
                   burn_in = 5000,
                   nsteps = 15000,
                   thin = 5,
                   sigma_errors = 1.0,
                   title_text = None,
                   add_fit_legend = True,
                   plot_fit_results=True,
                   verbose=0):
    """
    Peform a fit to the radio spectrum using SY + FF model in combination with 
    free-free absorption.

    Parameters
    ----------
    freqs : array
        Frequency array in Hz.
    fluxes : array
        Flux density array in mJy.
    fluxes_err : array
        Flux density error array in mJy.
    basename_save : str
        Basename to save the output files.
    verbose : int
        Verbosity level.

    """
    
    # if sigma_shade == 1 or sigma_shade == 1.0:
    #     quantiles = [0.16, 0.5, 0.84]
    # elif sigma_shade == 2 or sigma_shade == 2.0:
    #     quantiles = [0.025, 0.5, 0.975]
    # elif sigma_shade == 3 or sigma_shade == 3.0:
    #     quantiles = [0.00015, 0.5, 0.9985]
    # else:
    #     raise ValueError("Unsupported sigma value. Use sigma = 1, 2, or 3.")
    
    
    x = freqs / 1e9
    y = fluxes
    yerr = fluxes_err
    epsilon = 1e-8
    if nu0 is None:
        nu0 = np.mean(x)
    print(f' ++==>> Using reference frequency of {nu0} GHz.')



        
    # def min_func(params):
    #     A_sy = params['A_sy']
    #     A_ff = params['A_ff']
    #     alpha_nt = params['alpha_nt']
    #     nu_tau_t = params['nu_tau_t']
    #     # f_cov = params['f_cov']
    #     # model = RC_function_SY_FF_ffa_v2(x, A_sy, A_ff, alpha_nt, nu_tau_t, nu0, f_cov)
    #     model = RC_function_SY_FF_ffa_v2(x, A_sy, A_ff, alpha_nt, nu_tau_t, nu0)
    #     # res = (y - model) / (np.log(yerr))
        
    #     relative_error = yerr / (np.abs(y) + epsilon)
    #     weights = 1 / (relative_error + epsilon)
    #     # weightned_residual = (y - model) * np.sqrt(weights)
    #     weightned_residual = (y - model) / (np.sqrt(weights + yerr))
        
    #     # res = (y - model) * np.sqrt(2 / np.log1p(yerr**2))
    #     # res = (y - model)/yerr
    #     # loss = huber_loss(res)
        
    #     # return res.copy()
    #     return weightned_residual.copy()


    def min_func(params):
        A_sy = params['A_sy']
        A_ff = params['A_ff']
        alpha_nt = params['alpha_nt']
        nu_tau_t = params['nu_tau_t']
        # f_cov = params['f_cov']
        # model = RC_function_SY_FF_ffa_v2(x, A_sy, A_ff, alpha_nt, nu_tau_t, nu0, f_cov)
        model = RC_function_SY_FF_ffa_v2(x, A_sy, A_ff, alpha_nt, nu_tau_t, nu0)
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights


    fit_params = lmfit.Parameters()
    A_sy_init = np.log10(np.nanmean(y)*0.6)
    A_ff_init = np.log10(np.nanmean(y)*0.3)
    # a_nth_init = np.polyfit(np.log10(x), 
    #                                 np.log10(y), 1)[0]
    a_nth_init = -0.8
    
    print('Init log A_sy=',A_sy_init)
    print('Init log A_ff=',A_ff_init)
    print('Init A_sy=',10**A_sy_init)
    print('Init A_ff=',10**A_ff_init)
    print('Init a_nth=',a_nth_init)
    
    fit_params.add("A_sy", value=A_sy_init, min=-1, max=A_sy_init*20)
    fit_params.add("A_ff", value=A_ff_init, min=-1, max=A_ff_init*20)
    fit_params.add("alpha_nt", value=a_nth_init, min=-2.0, max=2.0)
    fit_params.add("nu_tau_t", value=1.0, min=0.1, max=10.0)
    # fit_params.add("f_cov", value=0.99, min=0.0, max=1.0)
    # fit_params.add("f_cov", value=1.0, min=0.9999, max=1.0,vary=False)
    # fit_params.add("f_cov", value=1.0, vary=False)
    

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))

    x_resample = np.linspace(np.min(x)*0.7, np.max(x)*1.3, 500)
    
    
    if do_mcmc_fit == True:
        nwalkers = int(4 * 25)
        
        if mcmc_version == 'lmfit':
            """ 
            This will be removed in the future.
            """
            results_emcee = mini.emcee(burn=burn_in,
                                    steps=nsteps,
                                    thin=thin, nwalkers=nwalkers,
                                    is_weighted=False,
                                    float_behavior='posterior',
                                    params=result.params, workers=-1)

            _samples_emcee = results_emcee.flatchain[burn_in:]
            _A_sy = np.asarray(_samples_emcee['A_sy'])
            _A_ff = np.asarray(_samples_emcee['A_ff'])
            _alpha_nt = np.asarray(_samples_emcee['alpha_nt'])
            _nu_tau_t = np.asarray(_samples_emcee['nu_tau_t'])
            _f_cov = np.asarray(_samples_emcee['f_cov'])
            model_samples = np.array(
                [RC_function_SY_FF_ffa_v2(x_resample, _A_sy[i], _A_ff[i], _alpha_nt[i], _nu_tau_t, nu0, _f_cov) for i in
                range(_samples_emcee.shape[0])])
            # samples_emcee = np.asarray(_samples_emcee[['A_sy', 'A_ff', 'alpha_nt', 'nu_tau_t', 'f_cov']])
            samples_emcee = np.asarray(_samples_emcee[['A_sy', 'A_ff', 'alpha_nt', 'nu_tau_t', 'f_cov']])
            param_dict = {}
            for i, label in enumerate(['A_sy', 'A_ff','alpha_nt', 'nu_tau_t', 'f_cov']):
                q = np.percentile(samples_emcee[:, i], [2.5, 50, 97.5])
                param_dict[label] = {
                    'best': q[1],
                    'lower': (q[1] - q[0]),
                    'upper': (q[2] - q[1])
                }
            print(param_dict)
        if mcmc_version == 'general':
            """
            This seems to be more robust as we have more control over
            the MCMC process.
            """
            def RC_SY_FF_FFA_v2(params,x):
                tau_nu = (x/params[3])**(-2.1)
                S_ff_abs = (10**params[1]) * ((1-np.exp(-tau_nu))/(tau_nu)) * ((x/nu0)**(-0.1))
                # S_sy_abs = (10**params[0]) * (1-params[4]*(1-np.exp(-tau_nu))) * ((x/nu0)**params[2])
                S_sy_abs = (10**params[0]) * (1-1*(1-np.exp(-tau_nu))) * ((x/nu0)**params[2])
                S_total_abs = S_ff_abs + S_sy_abs
                return S_total_abs
                
            samples_emcee, param_dict = general_mcmc(x_data = x, 
                                                    y_data = y, 
                                                    yerr_data = yerr, 
                                                    result_1 = result, 
                                                    param_labels = ["A_sy", "A_ff", 'alpha_nt', 'nu_tau_t'], 
                                                    # param_labels = ["A_sy", "A_ff", 'alpha_nt', 'nu_tau_t', 'f_cov'], 
                                                    burn_in = burn_in,
                                                    nwalkers = nwalkers,
                                                    nsteps = nsteps,
                                                    thin = thin,
                                                    sigma_errors = sigma_errors,
                                                    model_func = RC_SY_FF_FFA_v2,
                                                    quantiles=quantiles)
            _A_sy = samples_emcee.T[0]
            _A_ff = samples_emcee.T[1]
            _alpha_nt = samples_emcee.T[2]
            _nu_tau_t = samples_emcee.T[3]
            # _f_cov = samples_emcee.T[4]

            model_samples = np.array(
                [RC_function_SY_FF_ffa_v2(x_resample, 
                                    _A_sy[i], _A_ff[i], _alpha_nt[i], 
                                    _nu_tau_t[i], nu0) for i in range(_A_sy.shape[0])])
            print(param_dict)

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = RC_function_SY_FF_ffa_v2(x_resample,
                                        result.params['A_sy'].value,
                                        result.params['A_ff'].value,
                                        result.params['alpha_nt'].value,
                                        result.params['nu_tau_t'].value,
                                    nu0)

    model_best = RC_function_SY_FF_ffa_v2(x,
                                        result.params['A_sy'].value,
                                        result.params['A_ff'].value,
                                        result.params['alpha_nt'].value,
                                        result.params['nu_tau_t'].value,
                                    nu0)



    Sy_term = RC_function_SY_ffa_v2(x_resample,
                            result.params['A_sy'].value,
                            result.params['alpha_nt'].value,
                            result.params['nu_tau_t'].value,
                            nu0)
    
    
    ff_term = RC_function_FF_ffa_v2(x_resample,
                            result.params['A_ff'].value,
                            result.params['nu_tau_t'].value,
                            nu0)
    
    ff_term_pred = RC_function_FF_ffa_v2(x,
                            result.params['A_ff'].value,
                            result.params['nu_tau_t'].value,
                            nu0)
    
    thermal_fraction_freq = ff_term_pred / model_best

    if plot_fit_results:
        fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(7, 4),
                                    gridspec_kw={'height_ratios': [3, 1]})
        
        ax1.errorbar(x, y, yerr=yerr, 
                    fmt='o',label='Data', color='k', ecolor='gray',alpha=0.5)

        ax1.plot(x_resample,model_resample,
                color='red', ls='-.', label='Fit')
        ax2.plot(x, (y-model_best)/model_best,
                color='green', ls='dotted', label='Residual')
        ax2.set_ylim(-1.0,1.0)


        ax1.plot(x_resample, Sy_term,
                '-', label='sy')
        
        ax1.plot(x_resample, ff_term,
                '-', label='ff')

        if do_mcmc_fit == True:
            ax1.plot(x_resample, model_mean,
                    color='purple', linestyle=(0, (5, 10)), label='MCMC')

        # plt.fill_between(x_resample, lower_bound, upper_bound, color='lightgray', alpha=0.5,
        #                  label='Parameter uncertainty')
        if plot_errors_shade:
            if do_mcmc_fit == True:
                ax1.fill_between(x_resample,
                                model_mean - sigma_shade*model_std,
                                model_mean + sigma_shade*model_std, color='lightgray',
                                alpha=0.7)
            else:
                # Define the number of Monte Carlo samples
                num_samples = 5000

                # Generate random samples from parameter distributions
                A_sy_samples = np.random.normal(result.params['A_sy'].value, 
                                                result.params['A_sy'].stderr,
                                                num_samples)
                A_ff_samples = np.random.normal(result.params['A_ff'].value, 
                                                result.params['A_ff'].stderr,
                                                num_samples)
                alpha_nt_samples = np.random.normal(result.params['alpha_nt'].value,
                                                result.params['alpha_nt'].stderr,
                                                num_samples)
                nu_tau_t_samples = np.random.normal(result.params['nu_tau_t'].value,
                                                result.params['nu_tau_t'].stderr,
                                                num_samples)
                # f_cov_samples = np.random.normal(result.params['f_cov'].value,
                #                                 result.params['f_cov'].stderr,
                #                                 num_samples)

                # Compute model predictions for each sample
                model_predictions = np.zeros((num_samples, len(x_resample)))
                for i in range(num_samples):
                    model_predictions[i] = RC_function_SY_FF_ffa_v2(x_resample, 
                                                            A_sy_samples[i],
                                                            A_ff_samples[i],
                                                            alpha_nt_samples[i],
                                                            nu_tau_t_samples[i],
                                                            # f_cov_samples[i],
                                                            nu0)

                median_prediction = np.median(model_predictions, axis=0)
                std_prediction = np.std(model_predictions, axis=0)

                ax1.fill_between(x_resample, median_prediction - sigma_shade*std_prediction,
                                median_prediction + sigma_shade*std_prediction,
                                color='lightgray', alpha=0.5,
                                # label='Uncertainty (1-sigma)'
                                )
        # plt.ylim(1e-3,1.2*np.max(y))
        if add_fit_legend == True:
            ax1.legend(loc=(0.05, 0.05),frameon=True,prop={'size': 11},ncol=2)
        if np.nanmax(y) < 5:
            ax1.set_ylim(0.1*np.nanmin(y),10.0*np.nanmax(y))
        else:
            ax1.set_ylim(0.1*np.nanmin(y),3.0*np.nanmax(y))
        ax2.set_xlabel(r'$\nu$ [GHz]')
        # plt.ylabel('Integrated Flux Density [mJy]')
        ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
        
        text_x, text_y = 0.62, 0.72
        text = (r"$\alpha_{\rm nt}"+rf"= {(result.params['alpha_nt'].value):.2f}\pm "
                rf"{(result.params['alpha_nt'].stderr):.2f}$")
        text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
        text_bbox = plt.text(text_x, text_y, text,
                            # ha='center', va='center',
                            fontsize=12, color='black',
                            bbox=text_bbox_props, transform=fig.transFigure)
        
        if do_mcmc_fit == True:
            if mcmc_version == 'general':
                text_x, text_y = 0.62, 0.82
                text = (r"$\alpha_{\rm nt}^{\rm MCMC}"+f"= {(param_dict['alpha_nt']['best']):.2f}"
                        rf"_{{-{param_dict['alpha_nt']['lower']:.2f}}}^{{+{param_dict['alpha_nt']['upper']:.2f}}}$")
                text_bbox = plt.text(text_x, text_y, text,
                                    # ha='center', va='center',
                                    fontsize=12, color='black',
                                    # bbox=text_bbox_props, 
                                    transform=fig.transFigure)
        
        if title_text is not None:
            ax1.set_title(title_text)

        if log_plot == True:
            ax1.semilogx()
            ax1.semilogy()

        plt.subplots_adjust(hspace=0.05)
        ax2.legend(loc=(0.7, 0.55),prop={'size': 14},frameon=False,ncol=2)
        # legend()

        if basename_save is not None:
            if save_name_append is None:
                save_name_append = '_RC_sy_ff_ffa_v2_fit'
            else:
                save_name_append = save_name_append + '_RC_sy_ff_ffa_v2_fit'
            plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                        bbox_inches='tight')
        # plt.show()
        corner_kwargs = {
            'bins': 30,
            'hist_bin_factor': 1.0,
            'color': 'purple',  # Change the color of the scatter plots and histograms
            'hist_kwargs': {
                'color': 'green',      # Histogram color
                'edgecolor': 'black',  # Edge color of the bins
                'linewidth': 1.5       # Width of the bin edges
            },
            'scatter_kwargs': {
                'alpha': 0.6,          # Transparency of points
                's': 10,               # Size of points
                'color': 'purple',     # Color of points
                'edgecolor': 'none'    # Edge color of points
            },
            'contour_kwargs': {
                'colors': 'blue',       # Color of contour lines
                'linewidths': 1.2       # Width of contour lines
            }
        }
        if do_mcmc_fit == True:
            try:
                from scipy.stats import gaussian_kde
                # fig_c = plt.figure(figsize=(2,2))
                fig_c = plt.figure()
                _ = corner.corner(samples_emcee,
                                        labels=[r'$A_{\rm sy}$',
                                                r'$A_{\rm ff}$',
                                                r'$\alpha_{\rm nth}$',
                                                r'$\nu_{\rm t}$',
                                                r'$\Xi$'
                                                ],
                                        truths=[result.params['A_sy'].value,
                                                result.params['A_ff'].value,
                                                result.params['alpha_nt'].value,
                                                result.params['nu_tau_t'].value,
                                                # result.params['f_cov'].value
                                                ],
                                        show_titles=True,
                                        quantiles=quantiles,
                                        # **corner_kwargs
                                        # fig=fig_c
                                        )
                print(samples_emcee.shape)
                if basename_save is not None:
                    save_name_append_corner = save_name_append + '_corner'
                    plt.savefig(basename_save.replace('.fits','_')+save_name_append_corner+'.jpg', 
                                dpi=600,bbox_inches='tight')
                    
                print('++==>> Parameter Results (MCMC sampling).')
                print(lmfit.fit_report(results_emcee.params))
            except:
                pass
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))
    if do_mcmc_fit:
        return(mini,result,thermal_fraction_freq,samples_emcee, param_dict)
    else:
        return(mini,result,thermal_fraction_freq)


def do_fit_spec_RC_S2(freqs,fluxes,fluxes_err,nu0=None,
                      basename_save=None, log_plot=True,
                      save_name_append=None,
                      plot_errors_shade=False,
                      do_mcmc_fit=False,
                      title_text=None,
                      burn_in = 1000,steps=5000,
                      verbose=0):
    x = freqs / 1e9
    y = fluxes
    yerr = fluxes_err
    if nu0 is None:
        nu0 = np.mean(x)


    def min_func(params):
        A1 = params['A1']
        A2 = params['A2']
        alpha_nt = params['alpha_nt']
        # res = (y - RC_function_S2(x, A1, A2, alpha_nt,nu0))/yerr
        model = RC_function_S2(x, A1, A2, alpha_nt,nu0)
        # res = (y - model)/(y+yerr)
        # res = (y - model) / np.log(y+yerr)
        # res = (y - model) / (np.log(abs(yerr)+1.0)) # okay 9
        # res = data - RC_function_S2(nu, A1l, alpha_nt)
        # return res.copy()
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights

    fit_params = lmfit.Parameters()
    fit_params.add("A1", value=0.5, min=0, max=5000)
    fit_params.add("A2", value=0.5, min=-5, max=5000)
    fit_params.add("alpha_nt", value=-0.9, min=-2.5, max=2.5)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )


    # result_1 = mini.minimize(method='nelder',
    #                          options={'maxiter': 30000, 'maxfev': 30000,
    #                                   'xatol': 1e-12, 'fatol': 1e-12,
    #                                   'return_all': True,'adaptive': True,
    #                                   'disp': True}
    #                      )
    # second_run_params = result_1.params

    # result = mini.minimize(method='nelder',params=second_run_params,
    #                        options={'maxiter': 30000, 'maxfev': 30000,
    #                                 'xatol': 1e-12, 'fatol': 1e-12,
    #                                 'return_all': True,'adaptive': True,
    #                                 'disp': True}
    #                  )
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 5),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # fig = plt.figure(figsize=(8, 4))
    x_resample = np.linspace(np.min(x)*0.25, np.max(x)*1.3, 500)
    ax1.errorbar(x, y, yerr=yerr, fmt='o',
                 label='Data', color='k', ecolor='gray',
                 alpha=0.5)

    # nwalkers = int(len(y) * 3 * 20)
    nwalkers = int(3 * 100)

    if do_mcmc_fit == True:
        results_emcee = mini.emcee(burn=burn_in,
                                   steps=steps,
                                   thin=5, nwalkers=nwalkers,
                                   is_weighted=False,
                                   float_behavior='posterior',
                                   params=result.params, workers=-1)

        samples_emcee = results_emcee.flatchain[burn_in:]
        _A1 = np.asarray(samples_emcee['A1'])
        _A2 = np.asarray(samples_emcee['A2'])
        _alpha_nt = np.asarray(samples_emcee['alpha_nt'])

        model_samples = np.array(
            [RC_function_S2(x_resample, _A1[i], _A2[i], _alpha_nt[i], nu0) for i in
             range(samples_emcee.shape[0])])

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = RC_function_S2(x_resample,
                                 result.params['A1'].value,
                                 result.params['A2'].value,
                                 result.params['alpha_nt'].value,
                                 nu0)

    model_best = RC_function_S2(x,
                                 result.params['A1'].value,
                                 result.params['A2'].value,
                                 result.params['alpha_nt'].value,
                                 nu0)

    ax1.plot(x_resample, model_resample,
             color='red', ls='-.', label='Best-fit model')
    
    
    A1term = RC_function_S2(x_resample,
                            result.params['A1'].value,
                            result.params['A2'].value*0,
                            result.params['alpha_nt'].value,nu0)
    
    
    A2term = RC_function_S2(x_resample,
                            result.params['A1'].value*0,
                            result.params['A2'].value,
                            result.params['alpha_nt'].value,nu0)
    
    Snu0_low = RC_function_S2(1.4,
                            result.params['A1'].value,
                            result.params['A2'].value,
                            result.params['alpha_nt'].value,nu0)
    Snu0_mid = RC_function_S2(10.0,
                            result.params['A1'].value,
                            result.params['A2'].value,
                            result.params['alpha_nt'].value,nu0)
    Snu0_high = RC_function_S2(33.0,
                            result.params['A1'].value,
                            result.params['A2'].value,
                            result.params['alpha_nt'].value,nu0)
    
    thermal_fraction_low = result.params['A1'].value/Snu0_low
    thermal_fraction_mid = result.params['A1'].value/Snu0_mid
    thermal_fraction_high = result.params['A1'].value/Snu0_high
    
    ax1.plot(x_resample, A1term,
             '-', label='A1 term')
    
    ax1.plot(x_resample, A2term,
             '-', label='A2 term')

    ax2.plot(x, (y-model_best)/model_best,
             color='green', ls='dotted', label='Residual')
    ax2.set_xlabel(r'$\nu$ [GHz]')
    if plot_errors_shade:
        if do_mcmc_fit == True:
            ax1.plot(x_resample, model_mean,
                     color='purple', linestyle=(0, (5, 10)), label='MCMC Mean')
            ax1.fill_between(x_resample,
                             model_mean - 3*model_std,
                             model_mean + 3*model_std, color='lightgray',
                            alpha=0.5)
        else:
            # Define the number of Monte Carlo samples
            num_samples = 1000

            # Generate random samples from parameter distributions
            A1_samples = np.random.normal(result.params['A1'].value, result.params['A1'].stderr,
                                          num_samples)
            A2_samples = np.random.normal(result.params['A2'].value,
                                          result.params['A2'].stderr,
                                          num_samples)
            alpha_samples = np.random.normal(result.params['alpha_nt'].value,
                                             result.params['alpha_nt'].stderr,
                                             num_samples)

            # Compute model predictions for each sample
            model_predictions = np.zeros((num_samples, len(x_resample)))
            for i in range(num_samples):
                model_predictions[i] = RC_function_S2(x_resample,
                                                      A1_samples[i],
                                                      A2_samples[i],
                                                      alpha_samples[i],
                                                      nu0)

            median_prediction = np.median(model_predictions, axis=0)
            std_prediction = np.std(model_predictions, axis=0)

            ax1.fill_between(x_resample,
                             median_prediction - 1*std_prediction,
                             median_prediction + 1*std_prediction,
                             color='lightgray', alpha=0.5,
                             # label='Uncertainty (1-sigma)'
                             )
    # plt.ylim(1e-3,1.2*np.max(y))
    
    ax1.legend()
    ax1.set_ylim(0.1*np.min(y),3.0*np.max(y))
    ax1.set_xlabel(r'$\nu$ [GHz]')
    # plt.ylabel('Integrated Flux Density [mJy]')
    ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
    
    text_x, text_y = 0.72, 0.75
    text = (rf"$\alpha = {(result.params['alpha_nt'].value):.2f}\pm "
            rf"{(result.params['alpha_nt'].stderr):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)
    
    
    text_x, text_y = 0.50, 0.78
    text = (rf"$fth(10.0) = {(thermal_fraction_mid):.2f}\pm "
            rf"{(result.params['A1'].stderr/Snu0_mid):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        # bbox=text_bbox_props, 
                        transform=fig.transFigure)
    
    # text_x, text_y = 0.5, 0.35
    # text = (rf"$fth(8.0) = {(thermal_fraction_mid):.2f}\pm "
    #         rf"{(result.params['A1'].stderr/Snu0_mid):.2f}$")
    # text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    # text_bbox = plt.text(text_x, text_y, text,
    #                     # ha='center', va='center',
    #                     fontsize=12, color='black',
    #                     # bbox=text_bbox_props, 
    #                     transform=fig.transFigure)
    
    text_x, text_y = 0.50, 0.73
    text = (rf"$fth(33.0) = {(thermal_fraction_high):.2f}\pm "
            rf"{(result.params['A1'].stderr/Snu0_high):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        # bbox=text_bbox_props, 
                        transform=fig.transFigure)
    
    if title_text is not None:
        ax1.set_title(title_text)


    if log_plot == True:
        ax1.semilogx()
        ax1.semilogy()
    plt.subplots_adjust(hspace=0.05)
    ax2.legend()

    if basename_save is not None:
        if save_name_append is None:
            save_name_append = '_RC_alpha_fit_S2'
        else:
            save_name_append = save_name_append + '_RC_alpha_fit_S2'
        plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                    bbox_inches='tight')
    plt.show()
    plt.figure()

    if do_mcmc_fit == True:
        _ = corner.corner(samples_emcee[['A1','A2','alpha_nt']],
                                labels=[r'$A_1$',
                                        r'$A_2$',
                                        r'$\alpha_{\rm nt}$'],
                                truths=[result.params['A1'].value,
                                        result.params['A2'].value,
                                        result.params['alpha_nt'].value],
                                show_titles=True,
                                quantiles=[0.025, 0.5, 0.975])
        plt.show()
        print('++==>> Parameter Results (MCMC sampling).')
        print(lmfit.fit_report(results_emcee.params))
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))


    return mini,result

def RC_function_sy_ff_dust(nu, Asy, Aff, Adu, alpha_nt, alpha_du, nu0):
    dust_comp = Adu * (500 ** (alpha_du)) * ((nu / 500) ** (alpha_du))
    sy_comp = Asy * (nu0 ** (alpha_nt)) * ((nu / nu0) ** (alpha_nt))
    return Aff * ((nu / nu0) ** (-0.1)) + dust_comp + sy_comp

def do_fit_spec_RC_sy_ff_dust(freqs,fluxes,fluxes_err,nu0=None,
                      basename_save=None, log_plot=True,
                      save_name_append=None,
                      plot_errors_shade=False,
                      do_mcmc_fit=False,
                      title_text=None,
                      verbose=0):
    x = freqs / 1e9
    y = fluxes
    yerr = fluxes_err
    if nu0 is None:
        nu0 = np.mean(x)


    def min_func(params):
        Asy = params['Asy']
        Aff = params['Aff']
        Adu = params['Adu']
        alpha_nt = params['alpha_nt']
        alpha_du = params['alpha_du']
        # res = (y - RC_function_S2(x, A1, A2, alpha_nt,nu0))/yerr
        model = RC_function_sy_ff_dust(x, Asy, Aff, Adu, alpha_nt, alpha_du, nu0)
        # res = (y - model)/(y+yerr)
        # res = (y - model) / np.log(y+yerr)
        # res = (y - model) / (np.log(abs(yerr)+1.0)) # okay 9
        # res = data - RC_function_S2(nu, A1l, alpha_nt)
        # return res.copy()
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights

    fit_params = lmfit.Parameters()
    fit_params.add("Asy", value=0.5, min=0, max=5000)
    fit_params.add("Aff", value=0.5, min=0, max=5000)
    fit_params.add("Adu", value=0.5, min=0, max=5000)
    fit_params.add("alpha_nt", value=-0.9, min=-2.5, max=0.5)
    fit_params.add("alpha_du", value=3.0, min=1.8, max=4.0)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                             loss="cauchy", tr_solver="exact",
                             ftol=1e-15, xtol=1e-15, gtol=1e-15,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )


    # result_1 = mini.minimize(method='nelder',
    #                          options={'maxiter': 30000, 'maxfev': 30000,
    #                                   'xatol': 1e-12, 'fatol': 1e-12,
    #                                   'return_all': True,'adaptive': True,
    #                                   'disp': True}
    #                      )
    # second_run_params = result_1.params

    # result = mini.minimize(method='nelder',params=second_run_params,
    #                        options={'maxiter': 30000, 'maxfev': 30000,
    #                                 'xatol': 1e-12, 'fatol': 1e-12,
    #                                 'return_all': True,'adaptive': True,
    #                                 'disp': True}
    #                  )
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(8, 5),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # fig = plt.figure(figsize=(8, 4))
    x_resample = np.linspace(np.min(x)*0.25, np.max(x)*1.3, 500)
    ax1.errorbar(x, y, yerr=yerr, fmt='o',
                 label='Data', color='k', ecolor='gray',
                 alpha=0.5)

    # nwalkers = int(len(y) * 3 * 20)
    nwalkers = int(5 * 100)

    if do_mcmc_fit == True:
        burn_in = 3000
        results_emcee = mini.emcee(burn=burn_in,
                                   steps=20000,
                                   thin=5, nwalkers=nwalkers,
                                   is_weighted=False,
                                   float_behavior='posterior',
                                   params=result.params, workers=-1)

        samples_emcee = results_emcee.flatchain[burn_in:]
        _Asy = np.asarray(samples_emcee['Asy'])
        _Aff = np.asarray(samples_emcee['Aff'])
        _Adu = np.asarray(samples_emcee['Adu'])
        _alpha_nt = np.asarray(samples_emcee['alpha_nt'])
        _alpha_du = np.asarray(samples_emcee['alpha_du'])

        model_samples = np.array(
            [RC_function_sy_ff_dust(x_resample, 
                            _Asy[i], 
                            _Aff[i],
                            _Adu[i], 
                            _alpha_nt[i],
                            _alpha_du[i], 
                            nu0) for i in
             range(samples_emcee.shape[0])])

        model_mean = np.mean(model_samples, axis=0)
        model_std = np.std(model_samples, axis=0)

    model_resample = RC_function_sy_ff_dust(x_resample,
                                 result.params['Asy'].value,
                                 result.params['Aff'].value,
                                 result.params['Adu'].value,
                                 result.params['alpha_nt'].value,
                                 result.params['alpha_du'].value,
                                 nu0)

    model_best = RC_function_sy_ff_dust(x,
                                 result.params['Asy'].value,
                                 result.params['Aff'].value,
                                 result.params['Adu'].value,
                                 result.params['alpha_nt'].value,
                                 result.params['alpha_du'].value,
                                 nu0)

    ax1.plot(x_resample, model_resample,
             color='red', ls='-.', label='Best-fit model')
    
    
    Asyterm = RC_function_sy_ff_dust(x_resample,
                                 result.params['Asy'].value,
                                 result.params['Aff'].value*0,
                                 result.params['Adu'].value*0,
                                 result.params['alpha_nt'].value,
                                 result.params['alpha_du'].value,
                                 nu0)
    
    
    Affterm = RC_function_sy_ff_dust(x_resample,
                                 result.params['Asy'].value*0,
                                 result.params['Aff'].value,
                                 result.params['Adu'].value*0,
                                 result.params['alpha_nt'].value,
                                 result.params['alpha_du'].value,
                                 nu0)
    Aduterm = RC_function_sy_ff_dust(x_resample,
                                 result.params['Asy'].value*0,
                                 result.params['Aff'].value*0,
                                 result.params['Adu'].value,
                                 result.params['alpha_nt'].value,
                                 result.params['alpha_du'].value,
                                 nu0)

    Snu0_14 = RC_function_sy_ff_dust(1.4,
                                 result.params['Asy'].value,
                                 result.params['Aff'].value,
                                 result.params['Adu'].value,
                                 result.params['alpha_nt'].value,
                                 result.params['alpha_du'].value,
                            1.4)
    
    Snu0_6 = RC_function_sy_ff_dust(6.0,
                                 result.params['Asy'].value,
                                 result.params['Aff'].value,
                                 result.params['Adu'].value,
                                 result.params['alpha_nt'].value,
                                 result.params['alpha_du'].value,
                            6.0)
    
    thermal_fraction_14 = result.params['Aff'].value/Snu0_14
    
    thermal_fraction_6 = result.params['Aff'].value/Snu0_6
    
    ax1.plot(x_resample, Asyterm,
             '-', label='Sy term')
    
    ax1.plot(x_resample, Affterm,
             '-', label='FF term')
    ax1.plot(x_resample, Aduterm,
             '-', label='Dust term')

    ax2.plot(x, (y-model_best)/model_best,
             color='green', ls='dotted', label='Residual')
    ax2.set_xlabel(r'$\nu$ [GHz]')
    if plot_errors_shade:
        if do_mcmc_fit == True:
            ax1.plot(x_resample, model_mean,
                     color='purple', linestyle=(0, (5, 10)), label='MCMC Mean')
            ax1.fill_between(x_resample,
                             model_mean - 3*model_std,
                             model_mean + 3*model_std, color='lightgray',
                            alpha=0.5)
        else:
            # Define the number of Monte Carlo samples
            num_samples = 1000

            # Generate random samples from parameter distributions
            Asy_samples = np.random.normal(result.params['Asy'].value, 
                                          result.params['Asy'].stderr,
                                          num_samples)
            Aff_samples = np.random.normal(result.params['Aff'].value,
                                          result.params['Aff'].stderr,
                                          num_samples)
            Adu_samples = np.random.normal(result.params['Adu'].value,
                                          result.params['Adu'].stderr,
                                          num_samples)
            alpha_nt_samples = np.random.normal(result.params['alpha_nt'].value,
                                             result.params['alpha_nt'].stderr,
                                             num_samples)
            alpha_du_samples = np.random.normal(result.params['alpha_du'].value,
                                             result.params['alpha_du'].stderr,
                                             num_samples)

            # Compute model predictions for each sample
            model_predictions = np.zeros((num_samples, len(x_resample)))
            for i in range(num_samples):
                model_predictions[i] = RC_function_sy_ff_dust(x_resample,
                                                      Asy_samples[i],
                                                      Aff_samples[i],
                                                      Adu_samples[i],
                                                      alpha_nt_samples[i],
                                                      alpha_du_samples[i],
                                                      nu0)

            median_prediction = np.median(model_predictions, axis=0)
            std_prediction = np.std(model_predictions, axis=0)

            ax1.fill_between(x_resample,
                             median_prediction - 1*std_prediction,
                             median_prediction + 1*std_prediction,
                             color='lightgray', alpha=0.5,
                             # label='Uncertainty (1-sigma)'
                             )
    # plt.ylim(1e-3,1.2*np.max(y))
    
    ax1.legend()
    ax1.set_ylim(0.1*np.min(y),3.0*np.max(y))
    ax1.set_xlabel(r'$\nu$ [GHz]')
    # plt.ylabel('Integrated Flux Density [mJy]')
    ax1.set_ylabel(r'$S_{\nu}$ [mJy]')
    
    text_x, text_y = 0.70, 0.37
    text = (rf"$\alpha = {(result.params['alpha_nt'].value):.2f}\pm "
            rf"{(result.params['alpha_nt'].stderr):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)
    
    
    text_x, text_y = 0.5, 0.40
    text = (rf"$fth(1.4) = {(thermal_fraction_14):.2f}\pm "
            rf"{(result.params['Aff'].stderr/Snu0_14):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        # bbox=text_bbox_props, 
                        transform=fig.transFigure)
    
    text_x, text_y = 0.5, 0.35
    text = (rf"$fth(6.0) = {(thermal_fraction_6):.2f}\pm "
            rf"{(result.params['Aff'].stderr/Snu0_6):.2f}$")
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        # bbox=text_bbox_props, 
                        transform=fig.transFigure)
    
    if title_text is not None:
        ax1.set_title(title_text)


    if log_plot == True:
        ax1.semilogx()
        ax1.semilogy()
    plt.subplots_adjust(hspace=0.05)
    ax2.legend()

    if basename_save is not None:
        if save_name_append is None:
            save_name_append = '_RC_sy_ff_dust'
        else:
            save_name_append = save_name_append + '_RC_sy_ff_dust'
        plt.savefig(basename_save.replace('.fits','_')+save_name_append+'.jpg', dpi=600,
                    bbox_inches='tight')
    plt.show()
    plt.figure()

    if do_mcmc_fit == True:
        _ = corner.corner(samples_emcee[['Asy','Aff','Adu','alpha_nt','alpha_du']],
                                labels=[r'$A_{\rm sy}$',
                                        r'$A_{\rm ff}$',
                                        r'$A_{\rm du}$',
                                        r'$\alpha_{\rm nt}$',
                                        r'$\alpha_{\rm du}$'],
                                truths=[result.params['Asy'].value,
                                        result.params['Aff'].value,
                                        result.params['Adu'].value,
                                        result.params['alpha_nt'].value,
                                        result.params['alpha_du'].value],
                                show_titles=True,
                                quantiles=[0.025, 0.5, 0.975])
        plt.show()
        print('++==>> Parameter Results (MCMC sampling).')
        print(lmfit.fit_report(results_emcee.params))
    print('++==>> Parameter Results (from least-squares fit).')
    print(lmfit.fit_report(result.params))


    return mini,result



"""
#Spectral Index Maps
"""



def makecube(imagelist):
    temp = []
    # if isinstance(imagelist[0], str):
    for image in imagelist:
        temp.append(load_fits_data(image))
    # else:
    #     for image in imagelist:
    #         temp.append(image)
    return np.dstack(temp)



# def huber_loss(r):
#     """
#     Calculate the Huber loss for a residual array r.
    
#     Parameters:
#     - r (array): Residual array, where each element is (y_obs - y_pred)/yerr.
#     - delta (float): The threshold at which the loss function changes from quadratic to linear.
    
#     Returns:
#     - float: Sum of the Huber loss over all residuals.
#     """
#     delta = 1.0
#     # Calculate the condition for each residual
#     condition = np.abs(r) < delta
    
#     # Calculate the squared loss for residuals within the threshold
#     squared_loss = 0.5 * (r[condition] ** 2)
    
#     # Calculate the linear loss for residuals outside the threshold
#     linear_loss = delta * (np.abs(r[~condition]) - 0.5 * delta)
    
#     # Sum all the losses
#     total_loss = np.sum(squared_loss) + np.sum(linear_loss)
    
#     return total_loss



def huber_loss(residuals, delta=1.0):
    abs_res = np.abs(residuals)
    quadratic = 0.5 * residuals**2
    linear = delta * (abs_res - 0.5 * delta)
    return np.where(abs_res <= delta, quadratic, linear)

def calculate_weights(yerr, epsilon=0.1, w_min=0.01, w_max=100):
    """
    Calculate weights from error values.

    Args:
    yerr (numpy array): Array of measurement errors.
    epsilon (float): Small constant to prevent division by zero and soften the weight.
    w_min (float): Minimum allowable weight.
    w_max (float): Maximum allowable weight.

    Returns:
    numpy array: Weights for each measurement.
    """
    # Basic weights calculation with softening
    # weights = 1.0 / (yerr + epsilon)**2
    weights = 1/np.sqrt(yerr+epsilon)

    # Capping the weights
    # weights = np.clip(weights, w_min, w_max)

    return weights

# do_fit_spec_map now lives in spectral_fit.py (imported above).


# do_fit_spec_SY_FF_map now lives in spectral_fit.py (imported at the top).

def _fit_pixels_in_chunks(fit_chunk, idx, x, masked_cube, masked_cube_res,
                          fit_args, n_jobs=1, chunk_size=None, verbose=0,
                          label='pixel fitting', unit='pixels'):
    """
    Fit the spectra of the pixels `idx` in chunks, optionally in parallel.

    Each task receives only the spectra of its own pixels and returns plain
    values (no lmfit objects). `fit_chunk` must be a worker from
    spectral_fit.py, so joblib workers import that small module instead of
    the whole mlibs.

    Parameters
    ----------
    fit_chunk : callable
        fit_chunk(x, spectra, spectra_res, *fit_args), e.g.
        spectral_fit.fit_spectra_chunk.
    idx : array (npix, 2)
        (i, j) pixel indices to fit.
    x : array
        Frequencies in GHz.
    masked_cube, masked_cube_res : arrays (ny, nx, nfreq)
        Image and residual cubes.
    fit_args : tuple
        Extra arguments passed to fit_chunk after the spectra.
    n_jobs : int
        Number of parallel jobs (1 runs in this process, no joblib).
    chunk_size : int or None
        Pixels per task. None splits the pixels into ~16 chunks per job.
    verbose : int
        If > 0, print the number of pixels, jobs and chunks.
    label : str
        Name shown in the log message.
    unit : str
        What is being fitted, for the log message ('pixels', 'bins').

    Returns
    -------
    list
        Outputs of fit_chunk, one per chunk, in the order of `idx`.
    """
    nspec = len(idx)
    n_workers = effective_n_jobs(n_jobs)
    if chunk_size is None:
        chunk_size = int(np.clip(np.ceil(nspec / (n_workers * 16)), 1, 5000))
    chunks = [idx[k:k + chunk_size] for k in range(0, nspec, chunk_size)]

    if verbose > 0:
        print(f" > {label}: fitting {nspec} {unit} with n_jobs={n_jobs} "
              f"({n_workers} workers), in {len(chunks)} chunks of up to "
              f"{chunk_size} {unit}.")

    def _chunk_args(chunk):
        ii, jj = chunk[:, 0], chunk[:, 1]
        return (x, masked_cube[ii, jj, :], masked_cube_res[ii, jj, :],
                *fit_args)

    if n_jobs == 1:
        results = (fit_chunk(*_chunk_args(chunk)) for chunk in chunks)
    else:
        results = Parallel(n_jobs=n_jobs, return_as='generator')(
            delayed(fit_chunk)(*_chunk_args(chunk)) for chunk in chunks)

    return list(tqdm(results, total=len(chunks)))


def _noise_rms_cube(conv_cube_res, mask, noise_rms=None,
                    do_average_cube=False, bin_size=2, verbose=0):
    """
    Noise rms per pixel and band, used as sigma by weighting='chi2'.

    Parameters
    ----------
    conv_cube_res : array (ny, nx, nfreq)
        Residual cube (before masking and frequency averaging).
    mask : bool array (ny, nx)
        Source mask; by default the noise is measured outside it.
    noise_rms : None, array (nfreq,) or array (ny, nx, nfreq)
        None: per-band mad_std of the residuals outside `mask`.
        (nfreq,): one rms per band.
        (ny, nx, nfreq): rms maps (e.g. primary-beam corrected images, where
        the noise grows away from the pointing centre).
    do_average_cube, bin_size :
        Same frequency averaging as specidx_map; the rms is combined like the
        residuals there: sqrt(sum(rms**2)) / bin_size.

    Returns
    -------
    array (ny, nx, nfreq_out), read-only view when the rms is per band.
    """
    ny, nx, nf = conv_cube_res.shape
    if noise_rms is None:
        outside = ~np.asarray(mask, dtype=bool)
        rms = []
        for k in range(nf):
            band = conv_cube_res[:, :, k]
            sel = outside & np.isfinite(band)
            if sel.sum() < 100:   # mask covers (almost) all: use every pixel
                sel = np.isfinite(band)
            rms.append(mad_std(band[sel]))
        rms = np.array(rms, dtype=float)
        if verbose > 0:
            print(" > noise rms per band (outside mask): "
                  + ", ".join(f"{r:.3e}" for r in rms))
    else:
        rms = np.asarray(noise_rms, dtype=float)

    if rms.ndim == 1:
        if do_average_cube and rms.size == nf:
            rms = np.sqrt(np.sum(rms.reshape(-1, bin_size)**2, axis=1)) / bin_size
        return np.broadcast_to(rms, (ny, nx, rms.size))
    if rms.ndim == 3:
        if do_average_cube and rms.shape[2] == nf:
            rms = np.sqrt(np.sum(rms.reshape(ny, nx, -1, bin_size)**2,
                                 axis=3)) / bin_size
        return rms
    raise ValueError("noise_rms must be None, an array (nfreq,) or (ny, nx, nfreq)")


def _prepare_binned_spectra(binning, masked_cube, masked_cube_res, mask,
                            beam_area_px=1.0, binning_kwargs=None,
                            res_is_rms=False, verbose=0):
    """
    Bin the masked cube spatially and build one mean spectrum per bin.

    Parameters
    ----------
    binning : str
        Method name in spatial_binning.BINNING_METHODS (e.g. 'contour').
    masked_cube, masked_cube_res : arrays (ny, nx, nfreq)
        Image and residual cubes (NaN outside the mask).
    mask : bool array (ny, nx)
        Pixels to bin.
    beam_area_px : float
        Beam area in pixels.
    binning_kwargs : dict or None
        Options of the method. 'ref_map' (map whose contours the bins follow)
        and 'rms' (per-band noise) are optional; by default the band-mean of
        the cube and the per-band mad_std of the residuals are used.
    res_is_rms : bool
        False: masked_cube_res holds residual pixel values; the bin's
        residual is their mean (the noise of the mean spectrum).
        True: it holds the noise rms (weighting='chi2'); the bin's noise is
        the mean rms / sqrt(N_eff), N_eff = max(1, n_pix/beam_area_px).
    verbose : int
        Verbosity level.

    Returns
    -------
    labels : int array (ny, nx)
        Bin index per pixel, -1 outside the bins.
    fit_idx : int array (nbins, 2)
        (bin, 0) index of each bin in fit_cube.
    fit_cube, fit_cube_res : arrays (nbins, 1, nfreq)
        Mean image and residual spectrum of each bin.
    """
    kw = dict(binning_kwargs or {})
    ref_map = kw.pop('ref_map', None)
    rms = kw.pop('rms', None)
    if ref_map is None:
        import warnings
        with np.errstate(invalid='ignore'), warnings.catch_warnings():
            # pixels outside the mask are NaN in every band
            warnings.simplefilter('ignore', category=RuntimeWarning)
            ref_map = np.nanmean(masked_cube, axis=2)
    if rms is None:
        if res_is_rms:
            rms = np.array([np.nanmedian(band[mask])
                            for band in np.moveaxis(masked_cube_res, 2, 0)])
        else:
            rms = np.array([mad_std(band[np.isfinite(band)])
                            for band in np.moveaxis(masked_cube_res, 2, 0)])
    labels = spatial_binning.make_bins(binning, ref_map, masked_cube, mask,
                                       rms, beam_area_px=beam_area_px,
                                       verbose=verbose, **kw)
    spectra, spectra_res = spatial_binning.bin_spectra(labels, masked_cube,
                                                       masked_cube_res)
    nbins = spectra.shape[0]
    if res_is_rms:
        n_pix = np.bincount(labels[labels >= 0], minlength=nbins)
        n_eff = np.maximum(1.0, n_pix / beam_area_px)
        spectra_res = spectra_res / np.sqrt(n_eff)[:, None]
    fit_idx = np.column_stack([np.arange(nbins), np.zeros(nbins, dtype=int)])
    return labels, fit_idx, spectra[:, None, :], spectra_res[:, None, :]


def _fits_beam_cell(image):
    """
    (bmaj_arcsec, bmin_arcsec, bpa_deg), cell_arcsec and sign(CDELT1) read
    straight from a FITS header (no CASA). The cell is |CDELT2| (or |CD2_2|).
    """
    hdr = pf.getheader(image)
    beam = (hdr['BMAJ'] * 3600.0, hdr['BMIN'] * 3600.0, float(hdr.get('BPA', 0.0)))
    cdelt1 = hdr.get('CDELT1', hdr.get('CD1_1'))
    cdelt2 = hdr.get('CDELT2', hdr.get('CD2_2'))
    return beam, abs(cdelt2) * 3600.0, float(np.sign(cdelt1))


def _write_with_beam(header_source, out_name, data, beam, bunit=None):
    """
    Write `data` to `out_name` with the header of `header_source` and the
    beam `beam` = (bmaj_arcsec, bmin_arcsec, bpa_deg) in BMAJ/BMIN/BPA.
    """
    pf.writeto(out_name, np.asarray(data, dtype=np.float32), overwrite=True)
    copy_header(header_source, out_name, out_name)
    with pf.open(out_name, mode='update') as hdul:
        hdr = hdul[0].header
        hdr['BMAJ'] = beam[0] / 3600.0
        hdr['BMIN'] = beam[1] / 3600.0
        hdr['BPA'] = beam[2]
        if bunit is not None:
            hdr['BUNIT'] = bunit
        elif 'BUNIT' in hdr:
            hdr['BUNIT'] = 'JY/BEAM'


def gaussian_beam_kernel(bmaj_px, bmin_px, bpa_deg, cdelt1_sign=-1.0,
                         nsigma=20.0, max_half_size=None, normalise='peak'):
    """
    Elliptical Gaussian restoring-beam kernel on the pixel grid of a FITS
    image (north up, square pixels).

    BPA follows the FITS/WSClean convention: angle of the major axis measured
    from north towards east. East is towards -x when CDELT1 < 0 (the usual RA
    axis), so the major axis points along (dx, dy) = (cdelt1_sign*sin(BPA),
    cos(BPA)) in pixel units. This was checked against `wsclean -restore`
    (WSClean 3.6, schaapcommon RestoreImage): same peak, centroid and BPA.

    Parameters
    ----------
    bmaj_px, bmin_px : float
        FWHM of the major/minor axes in pixels.
    bpa_deg : float
        Position angle in degrees, east of north.
    cdelt1_sign : float
        Sign of CDELT1 of the image the kernel is used on (-1 for RA axes).
    nsigma : float
        Half-size of the kernel box in units of sigma_major (WSClean uses a
        40-sigma full box, i.e. 20 per side).
    max_half_size : int or None
        Clip the half-size (e.g. to the image half-size).
    normalise : 'peak' or 'sum'
        'peak': unit peak, so a model in Jy/pixel becomes Jy/beam (what
        WSClean does when restoring). 'sum': unit integral.
    """
    fwhm2sig = 1.0 / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    smaj = bmaj_px * fwhm2sig
    smin = bmin_px * fwhm2sig
    half = max(int(np.ceil(nsigma * smaj)), 1)
    if max_half_size is not None:
        half = min(half, int(max_half_size))
    y, x = np.mgrid[-half:half + 1, -half:half + 1].astype(float)
    th = np.deg2rad(bpa_deg)
    ux, uy = cdelt1_sign * np.sin(th), np.cos(th)  # major-axis unit vector
    u = x * ux + y * uy                             # along the major axis
    v = -x * uy + y * ux                            # along the minor axis
    kernel = np.exp(-0.5 * ((u / smaj) ** 2 + (v / smin) ** 2))
    if normalise == 'sum':
        kernel /= kernel.sum()
    return kernel


def restore_model_with_beam(model, bmaj_px, bmin_px, bpa_deg, residual=None,
                            cdelt1_sign=-1.0):
    """
    Restore a clean-component model (Jy/pixel) with an elliptical Gaussian
    beam: model (*) G_peak1 [+ residual], output in Jy/beam.

    This is the image-plane equivalent of WSClean's restore step
    (image = residual + model (*) G, G with unit peak). A point source of
    S Jy in the model gives a peak of S Jy/beam, and the integrated flux
    sum(output)/beam_area_px equals sum(model). The residual is added at its
    native resolution, as WSClean does.
    """
    model = np.nan_to_num(np.asarray(model, dtype=float))
    ny, nx = model.shape
    kernel = gaussian_beam_kernel(bmaj_px, bmin_px, bpa_deg,
                                  cdelt1_sign=cdelt1_sign,
                                  max_half_size=min(ny, nx) // 2)
    restored = scipy.signal.fftconvolve(model, kernel, mode='same')
    if residual is not None:
        restored = restored + np.asarray(residual, dtype=float)
    return restored


def smooth_image_to_beam(image, native_beam, target_beam, cell_arcsec,
                         cdelt1_sign=-1.0):
    """
    Smooth an image in Jy/beam from `native_beam` to `target_beam` (both
    (bmaj_arcsec, bmin_arcsec, bpa_deg)): convolve with the unit-area
    differential Gaussian and rescale by A_target/A_native so the output is
    in Jy/target_beam. Raises if the target is not broader than the native
    beam in every direction.
    """
    diff = _deconvolve_beams(target_beam, native_beam)
    if diff is None:
        raise ValueError(f"target beam {target_beam} is not broader than the "
                         f"native beam {native_beam} in every direction")
    image = np.nan_to_num(np.asarray(image, dtype=float))
    ny, nx = image.shape
    kernel = gaussian_beam_kernel(diff[0] / cell_arcsec, diff[1] / cell_arcsec,
                                  diff[2], cdelt1_sign=cdelt1_sign,
                                  max_half_size=min(ny, nx) // 2,
                                  normalise='sum')
    area_ratio = (target_beam[0] * target_beam[1]) / (native_beam[0] * native_beam[1])
    return scipy.signal.fftconvolve(image, kernel, mode='same') * area_ratio


def _prepare_spec_cube(imagelist, residuallist,
                       ref_image_conv=None,
                       freqs=None,
                       ref_image_mask=None,
                       mask=None, sigma_global_mask=6,
                       iterations=1,
                       dilation_size=2,
                       needs_convolution=False, conv_task='fft',
                       do_average_cube=False, bin_size=int(2),
                       beam_area_px=None,
                       verbose=0,
                       use_model_images=False,
                       modellist=None,
                       target_beam=None,
                       add_residual=True,
                       smooth_residual_to_target=False,
                       residual_beams=None,
                       model_out_dir=None):
    """
    Build the cubes that `specidx_map` and `t_specidx_map_SY_FF` fit: load
    the images and residuals, optionally restore the clean-component models
    (use_model_images) or convolve (needs_convolution), make or take the
    mask, mask the cubes and optionally average them in frequency.

    The arguments have the meaning documented in `specidx_map`.

    Returns
    -------
    dict with
        conv_cube, conv_cube_res : (ny, nx, nfreq) image and residual cubes
            before masking and frequency averaging;
        masked_cube, masked_cube_res : float32 cubes, NaN outside the mask,
            frequency-averaged if do_average_cube;
        mask : bool (ny, nx);
        freqs : frequencies in Hz matching masked_cube (None only for array
            inputs given without freqs);
        ref_image_conv, target_beam, beam_area_px : as resolved here (they
            can still be None).
    """
    if isinstance(imagelist[0], str):
        cube_image = makecube(imagelist)
        cube_residuals = makecube(residuallist)
    else:
        cube_image = imagelist
        cube_residuals = residuallist

    if freqs is None and isinstance(imagelist[0], str):
        freqs = getfreqs(imagelist)
    if freqs is not None:
        freqs = np.asarray(freqs, dtype=float)
        if freqs.size != cube_image.shape[2]:
            raise ValueError(f"{freqs.size} frequencies were given for "
                             f"{cube_image.shape[2]} images")

    if use_model_images:
        if needs_convolution:
            raise ValueError("use_model_images and needs_convolution are "
                             "mutually exclusive")
        if modellist is None:
            raise ValueError("use_model_images needs a modellist")
        model_files = isinstance(modellist[0], str)
        cube_model = (makecube(modellist) if model_files
                      else np.asarray(modellist, dtype=float))
        if cube_model.shape != cube_residuals.shape:
            raise ValueError("modellist must hold one model per residual, on "
                             "the same pixel grid")
        header_ref = next((lst[0] for lst in (modellist, residuallist, imagelist)
                           if isinstance(lst[0], str)), None)
        if header_ref is None:
            raise ValueError("use_model_images needs at least one list of "
                             "FITS files to read the pixel scale from")
        _, cell_model, cdelt1_sign = _fits_beam_cell(header_ref)
        if target_beam is None:
            if ref_image_conv is None:
                ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                                        return_df=False)[0][::-1][0]
            target_beam, _, _ = _fits_beam_cell(ref_image_conv)
        target_beam = tuple(float(b) for b in target_beam)
        if smooth_residual_to_target:
            if residual_beams is None or len(residual_beams) != cube_residuals.shape[2]:
                raise ValueError("smooth_residual_to_target needs one native "
                                 "beam per residual in residual_beams")
            cube_residuals = np.dstack([
                smooth_image_to_beam(cube_residuals[:, :, i], residual_beams[i],
                                     target_beam, cell_model,
                                     cdelt1_sign=cdelt1_sign)
                for i in range(cube_residuals.shape[2])])
        if verbose > 0:
            print(f" > use_model_images: restoring {cube_model.shape[2]} models "
                  f"with beam {target_beam[0]:.3f}\" x {target_beam[1]:.3f}\" "
                  f"PA {target_beam[2]:.1f} deg (cell {cell_model:.4f}\"), "
                  f"residual {'added' if add_residual else 'NOT added'}"
                  f"{' (smoothed to target)' if smooth_residual_to_target else ''}")
        cube_image = np.dstack([
            restore_model_with_beam(cube_model[:, :, i],
                                    target_beam[0] / cell_model,
                                    target_beam[1] / cell_model,
                                    target_beam[2],
                                    residual=(cube_residuals[:, :, i]
                                              if add_residual else None),
                                    cdelt1_sign=cdelt1_sign)
            for i in range(cube_model.shape[2])])
        if beam_area_px is None:
            beam_area_px = (np.pi * target_beam[0] * target_beam[1]
                            / (4.0 * np.log(2.0)) / cell_model ** 2)
        if model_out_dir is not None:
            os.makedirs(model_out_dir, exist_ok=True)
            for i in range(cube_image.shape[2]):
                src = modellist[i] if isinstance(modellist[0], str) else header_ref
                out_name = os.path.join(model_out_dir, os.path.basename(src).replace(
                    '-model', '-modelrestored' + ('' if add_residual else '-nores')))
                _write_with_beam(src, out_name, cube_image[:, :, i], target_beam)

    if mask is None:
        if not isinstance(imagelist[0], str):
            raise ValueError("a mask is required when the inputs are arrays "
                             "(it can only be built from FITS files)")
        if ref_image_mask is not None:
            ref_image_mask = ref_image_mask
        else:
            ref_image_mask = imagelist[-1]
        _,mask = mask_dilation(ref_image_mask,
                               rms=mad_std(load_fits_data(residuallist[-1])),
                                     show_figure=True,
                                     PLOT=True,
                                     iterations=iterations,
                                     dilation_size=dilation_size,
                                     sigma=sigma_global_mask)
    mask = np.asarray(mask, dtype=bool)
    if mask.shape != cube_image.shape[:2]:
        raise ValueError(f"mask shape {mask.shape} does not match the images "
                         f"{cube_image.shape[:2]}")
    inv_mask = ~mask
    inv_mask_3d = inv_mask[:, :, np.newaxis]

    num_images = cube_image.shape[2]
    if needs_convolution:
        if conv_task != 'fft':
            raise ValueError(f"conv_task={conv_task!r} is not implemented "
                             "(only 'fft')")
        if ref_image_conv is None:
            ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                                return_df = False)[0][::-1][0]
        psf_zise = int(get_beam_size_px(ref_image_conv)[0])
        psf_image_size = int(load_fits_data(ref_image_conv).shape[0])
        print(f"PSF BEAM SIZE is --> {psf_zise} px")
        print(f"PSF IMAGE SIZE is --> {psf_image_size} px")

        psf_name = tcreate_beam_psf(ref_image_conv,
                                        size=(psf_image_size, psf_image_size),
                                        aspect = 'elliptical')
        PSF_DATA = load_fits_data(psf_name)

        conv_cube = np.empty_like(cube_image)
        conv_cube_res = np.empty_like(cube_residuals)
        larger_beam_area = beam_area2(ref_image_conv)
        for i in tqdm(range(num_images)):
            conv_image_i_uncor = scipy.signal.fftconvolve(cube_image[:, :, i],
                                                        PSF_DATA, 'same')
            conv_residual_i_uncor = scipy.signal.fftconvolve(cube_residuals[:, :, i],
                                                        PSF_DATA, 'same')

            beam_area_i = beam_area2(imagelist[i])

            # factor_conv = np.max(larger_beam_image_data)/np.max(conv_image_i_uncor)
            factor_conv_i = larger_beam_area/beam_area_i

            print(f"Factor Convolution is --> {factor_conv_i}")
            # factor_conv = 1.0

            conv_cube[:, :, i] = conv_image_i_uncor * factor_conv_i
            conv_image_i = conv_cube[:, :, i]

            conv_cube_res[:, :, i] = conv_residual_i_uncor * factor_conv_i
            conv_res_i = conv_cube_res[:, :, i]

            conv_name = imagelist[i].replace('-image.','-image-conv.').replace('-image.cutout','-image-conv.cutout')
            pf.writeto(conv_name,
                            conv_image_i,
                            overwrite=True)
            copy_header(ref_image_conv,conv_name,conv_name)

            conv_res_name = residuallist[i].replace('-residual.','-residual-conv.').replace('-residual.cutout','-residual-conv.cutout')
            pf.writeto(conv_res_name,
                            conv_res_i,
                            overwrite=True)
            copy_header(ref_image_conv,conv_res_name,conv_res_name)
    else:
        conv_cube = cube_image.copy()
        conv_cube_res = cube_residuals.copy()

    masked_cube = np.where(inv_mask_3d, np.nan, conv_cube)
    masked_cube_res = np.where(inv_mask_3d, np.nan, conv_cube_res)

    if do_average_cube:
        if freqs is None:
            raise ValueError("do_average_cube needs freqs for array inputs")
        if masked_cube.shape[2] % bin_size != 0:
            raise ValueError(f"do_average_cube: the number of images "
                             f"({masked_cube.shape[2]}) is not a multiple of "
                             f"bin_size ({bin_size})")
        reshaped_cube = masked_cube.reshape(masked_cube.shape[0],
                                            masked_cube.shape[1], -1,
                                            bin_size)
        averaged_cube = reshaped_cube.mean(axis=3)

        reshaped_residual = masked_cube_res.reshape(masked_cube_res.shape[0],
                                                    masked_cube_res.shape[1], -1,
                                                    bin_size)
        averaged_cube_res = np.sqrt(np.sum(reshaped_residual**2, axis=3)) / bin_size

        print("Original cube shape:", masked_cube.shape)
        print("Averaged cube shape:", averaged_cube.shape)
        reshaped_frequencies = freqs.reshape(-1, bin_size)
        averaged_frequencies = reshaped_frequencies.mean(axis=1)
        freqs = averaged_frequencies
        masked_cube = averaged_cube.astype(np.float32)
        masked_cube_res = averaged_cube_res.astype(np.float32)
    else:
        masked_cube = masked_cube.astype(np.float32)
        masked_cube_res = masked_cube_res.astype(np.float32)

    return dict(conv_cube=conv_cube, conv_cube_res=conv_cube_res,
                masked_cube=masked_cube, masked_cube_res=masked_cube_res,
                mask=mask, freqs=freqs, ref_image_conv=ref_image_conv,
                target_beam=target_beam, beam_area_px=beam_area_px)


def _prepare_fit_spectra(imagelist, conv_cube_res, masked_cube, masked_cube_res,
                         mask, weighting='log', noise_rms=None,
                         do_average_cube=False, bin_size=int(2),
                         binning=None, binning_kwargs=None, beam_area_px=None,
                         verbose=0):
    """
    Decide what is fitted: one spectrum per pixel (binning=None), or one
    mean spectrum per spatial bin. In the binned case the bins are laid out
    as a (nbins, 1, nfreq) "cube" so the same machinery fits them.

    The noise term that goes with each spectrum is the residual pixel value
    for weighting='log' and the noise rms (`_noise_rms_cube`) for 'chi2'.

    Returns
    -------
    bin_labels : int (ny, nx) bin index per pixel (-1 outside), or None.
    fit_idx : (nspec, 2) indices of the spectra in fit_cube.
    fit_cube, fit_cube_res : spectra and their noise term.
    beam_area_px : the beam area used by the binning (as given otherwise).
    noise_term : (ny, nx, nfreq) per-pixel noise term before binning.
    """
    idx = np.column_stack(np.where(mask==True))
    if weighting == 'chi2':
        # sigma per pixel and band from the noise rms, not from the
        # residual pixel values
        noise_term = _noise_rms_cube(conv_cube_res, mask,
                                     noise_rms=noise_rms,
                                     do_average_cube=do_average_cube,
                                     bin_size=bin_size, verbose=verbose)
        if noise_term.shape != masked_cube.shape:
            raise ValueError(f"noise_rms gives a noise cube of shape "
                             f"{noise_term.shape}, but the cube to fit is "
                             f"{masked_cube.shape}")
    elif weighting == 'log':
        noise_term = masked_cube_res
    else:
        raise ValueError(f"weighting must be 'log' or 'chi2', not {weighting!r}")

    bin_labels = None
    if binning is None:
        fit_idx, fit_cube, fit_cube_res = idx, masked_cube, noise_term
    else:
        if beam_area_px is None:
            if isinstance(imagelist[0], str):
                # conservative: the largest beam among the inputs
                beam_area_px = max(beam_area2(im) for im in imagelist)
            else:
                beam_area_px = 1.0
                print(" > WARNING: beam_area_px not given for array "
                      "inputs; assuming 1 pixel per beam (bins will be "
                      "too small if pixels are correlated).")
        bin_labels, fit_idx, fit_cube, fit_cube_res = _prepare_binned_spectra(
            binning, masked_cube, noise_term, mask,
            beam_area_px=beam_area_px,
            binning_kwargs=binning_kwargs,
            res_is_rms=(weighting == 'chi2'), verbose=verbose)
    return bin_labels, fit_idx, fit_cube, fit_cube_res, beam_area_px, noise_term


def _spec_output_tag(binning=None, weighting='log', use_model_images=False,
                     add_residual=True, smooth_residual_to_target=False):
    """
    Tag added to the output FITS names, so binned, chi2 and model-restored
    maps never overwrite the default (per-pixel, 'log') ones.
    """
    tag = (('' if binning is None else f'-{binning}bin')
           + ('' if weighting == 'log' else f'-{weighting}'))
    if use_model_images:
        tag += ('-modelrestored' + ('' if add_residual else '-nores')
                + ('-ressmooth' if smooth_residual_to_target else ''))
    return tag


def _spec_output_name(ref_image_conv, product, tag='', out_dir=None):
    """
    Output FITS name of a map derived from `ref_image_conv`:
    '<...>-image<...>' becomes '<...>-<product><tag><...>', written one
    directory above the image, or in `out_dir` if given.
    """
    name = (ref_image_conv.replace('-image.', f'-{product}{tag}.').
            replace('-image.cutout', f'-{product}{tag}.cutout').
            replace('-image-pb', f'-{product}{tag}-pb'))
    if name == ref_image_conv:
        # the file name does not follow the '-image' convention: never
        # overwrite the input image
        root, ext = os.path.splitext(ref_image_conv)
        name = f'{root}-{product}{tag}{ext}'
    if out_dir is not None:
        os.makedirs(out_dir, exist_ok=True)
        return os.path.join(out_dir, os.path.basename(name))
    return os.path.dirname(os.path.dirname(name)) + '/' + os.path.basename(name)


def _report_files_written(files_written, save_fits=True, file_inputs=True):
    """Say where the FITS maps went, or why none were written."""
    if files_written:
        folder = os.path.abspath(os.path.dirname(files_written[0]))
        print(f" > {len(files_written)} FITS maps written to {folder}")
    elif not save_fits:
        print(" > no FITS maps written (save_fits=False).")
    elif not file_inputs:
        print(" > no FITS maps written: the inputs are arrays, not FITS files.")
    else:
        print(" > no FITS maps written.")


def _save_spec_maps(to_write, ref_image_conv, tag='', target_beam=None,
                    out_dir=None):
    """
    Write maps as FITS with the header of `ref_image_conv`.

    to_write : list of (product, data, bunit); the file name is
        `_spec_output_name(ref_image_conv, product, tag, out_dir)`.
    target_beam : (bmaj_arcsec, bmin_arcsec, bpa_deg) to stamp in the
        headers when the maps are not at the beam of ref_image_conv.

    Returns the list of file names written.
    """
    written = []
    for product, data, bunit in to_write:
        output_path = _spec_output_name(ref_image_conv, product, tag, out_dir)
        data = np.asarray(data)
        if np.issubdtype(data.dtype, np.floating):
            # FITS has no +-inf convention readers agree on
            data = np.where(np.isfinite(data), data, np.nan).astype(np.float32)
        pf.writeto(output_path, data, overwrite=True)
        copy_header(ref_image_conv, output_path, output_path)
        with pf.open(output_path, mode='update') as hdul:
            hdul[0].header['BUNIT'] = bunit
            if target_beam is not None:
                hdul[0].header['BMAJ'] = target_beam[0] / 3600.0
                hdul[0].header['BMIN'] = target_beam[1] / 3600.0
                hdul[0].header['BPA'] = target_beam[2]
        written.append(output_path)
    return written


def specidx_map(imagelist,residuallist,
                ref_image_conv=None,
                freqs=None,
                ref_image_mask = None,
                flux_sys_error_frac = 0.1,
                nu0=10.0,
                mask=None,sigma_global_mask=6,
                iterations=1,
                dilation_size=2,
                needs_convolution=False,conv_task='fft',
                return_only_cube=False,
                do_average_cube=False,bin_size=int(2),
                n_jobs=1,
                chunk_size=None,
                binning=None,
                binning_kwargs=None,
                beam_area_px=None,
                return_bin_map=False,
                weighting='log',
                robust_scale=3.0,
                noise_rms=None,
                verbose=0,
                use_model_images=False,
                modellist=None,
                target_beam=None,
                add_residual=True,
                smooth_residual_to_target=False,
                residual_beams=None,
                model_out_dir=None,
                out_dir=None,
                save_fits=True):

    """
    Calculates the spectral index map from a list of images.
    
    Parameters:
    ----------
    imagelist (list): List of images to calculate the spectral index map.
    residuallist (list): List of residual images to calculate the spectral index map.
    ref_image_conv (str): Reference image to use for convolution.
    freqs (list): List of frequencies corresponding to the images.
    ref_image_mask (str): Reference image to use for masking.
    flux_sys_error_frac (float): Fractional systematic error in flux.
    nu0 (float): Reference frequency.
    mask (array): Mask to use for masking.
    sigma_global_mask (float): Sigma value to use for global masking.
    iterations (int): Number of iterations to use for masking.
    dilation_size (int): Dilation size to use for masking.
    needs_convolution (bool): Whether to convolve the images.
    conv_task (str): Convolution task to use.
    return_only_cube (bool): Whether to return only the cube.
    do_average_cube (bool): Whether to average the cube.
    bin_size (int): Bin size to use for averaging.
    n_jobs (int): Number of jobs to use for parallel processing.
    chunk_size (int): Pixels fitted per parallel task. Default (None)
        splits the pixels into ~16 chunks per job.
    binning (str or None): Spatial binning before fitting. None (default):
        fit every pixel. 'contour': contour binning (Sanders 2006), bins
        follow the brightness contours and reach a target SNR in every band;
        one mean spectrum per bin is fitted and its alpha is given to all
        pixels of the bin. Methods available: spatial_binning.BINNING_METHODS.
    binning_kwargs (dict): Options of the binning method. For 'contour':
        target_snr (10), smooth_sigma_px (1.0), constraint (3.0), plus
        optional ref_map (2D map to follow; default: band-mean of the cube)
        and rms (per-band noise; default: mad_std of the residuals).
    beam_area_px (float): Beam area in pixels, used by the binning SNR.
        Default: largest beam of the input images (1.0 for array inputs).
    return_bin_map (bool): If True, also return the bin label map
        (-1 outside bins; None without binning) as a 6th output.
    weighting (str): How uncertainties enter the fit.
        'log' (default, previous behaviour): yerr = sqrt(residual_pixel**2 +
            (flux_sys_error_frac*S)**2), turned into mild log-weights; alpha
            errors are rescaled by the scatter of the points about the fit.
        'chi2': yerr = sqrt(noise_rms**2 + (flux_sys_error_frac*S)**2) is the
            1-sigma error of each point; residuals (S - model)/yerr; robust
            (Cauchy) threshold of robust_scale sigma; absolute alpha errors
            (they include flux_sys_error_frac; set it to 0 for noise only).
    robust_scale (float): 'chi2' only. Points beyond ~robust_scale sigma
        from the model are progressively down-weighted.
    noise_rms: 'chi2' only. None (default): per-band mad_std of the
        residual images outside the mask. Or an array (nfreq,) of per-band
        rms, or rms maps (ny, nx, nfreq), e.g. for primary-beam corrected
        images.
    verbose (int): Verbosity level to use.
    use_model_images (bool): Build the restored cube from the clean-component
        models instead of using imagelist: each plane is
        model (*) G_target [+ residual], with G_target a unit-peak elliptical
        Gaussian (models in Jy/pixel -> Jy/beam; a point source of S Jy gives
        a peak of S Jy/beam and the integrated flux is conserved). This is
        WSClean's own restore step, so with the same beam it reproduces
        WSClean's restored images. The residual is added at its native
        resolution (as WSClean does) and residuallist is used unchanged for
        the uncertainties, so both routes share the same error model.
        Mutually exclusive with needs_convolution. imagelist is then only
        used for headers (frequencies if freqs is None, output names).
        Output FITS names get a '-modelrestored' tag ('-nores' without the
        residual, '-ressmooth' with smooth_residual_to_target) and the
        target beam in BMAJ/BMIN/BPA; beam_area_px defaults to the target
        beam area.
    modellist (list): Model images (Jy/pixel), same order and pixel grid as
        imagelist/residuallist. Required when use_model_images=True.
    target_beam (tuple): (bmaj_arcsec, bmin_arcsec, bpa_deg) of the restoring
        beam, BPA east of north. Default: the beam in the header of
        ref_image_conv, or of the largest-beam image of imagelist.
    add_residual (bool): False leaves the residual out (convolved model
        only, a diagnostic); residuallist is still used for the errors.
    smooth_residual_to_target (bool): Also bring each residual from its
        native beam (residual_beams[i]) to target_beam, with the differential
        Gaussian and the A_target/A_native rescaling, before adding it and
        before using it for the errors. Needed when the target beam is larger
        than the native resolution of the residuals.
    residual_beams (list): (bmaj_arcsec, bmin_arcsec, bpa_deg) native beam of
        each residual, for smooth_residual_to_target.
    model_out_dir (str): If given, the restored planes are written there as
        FITS, with the target beam in BMAJ/BMIN/BPA. Nothing is written
        otherwise.
    out_dir (str): Directory for the alpha maps. Default (None): one
        directory above the input images.
    save_fits (bool): False writes no alpha maps (file inputs only write
        them when True, the default).

    Note on needs_convolution: it convolves every image with the *full*
    reference beam (not the differential kernel), so the resulting beam is
    beam_i (*) beam_ref rather than beam_ref, and the stamped header beam is
    not the real one. Kept as is for backwards compatibility.

    """
    prep = _prepare_spec_cube(
        imagelist, residuallist, ref_image_conv=ref_image_conv, freqs=freqs,
        ref_image_mask=ref_image_mask, mask=mask,
        sigma_global_mask=sigma_global_mask, iterations=iterations,
        dilation_size=dilation_size, needs_convolution=needs_convolution,
        conv_task=conv_task, do_average_cube=do_average_cube,
        bin_size=bin_size, beam_area_px=beam_area_px, verbose=verbose,
        use_model_images=use_model_images, modellist=modellist,
        target_beam=target_beam, add_residual=add_residual,
        smooth_residual_to_target=smooth_residual_to_target,
        residual_beams=residual_beams, model_out_dir=model_out_dir)
    conv_cube = prep['conv_cube']
    conv_cube_res = prep['conv_cube_res']
    masked_cube = prep['masked_cube']
    masked_cube_res = prep['masked_cube_res']
    mask = prep['mask']
    freqs = prep['freqs']
    ref_image_conv = prep['ref_image_conv']
    target_beam = prep['target_beam']
    beam_area_px = prep['beam_area_px']

    if return_only_cube:
        return(None,None,conv_cube,masked_cube_res,masked_cube)
    else:
        if freqs is None:
            raise ValueError("freqs must be given when the inputs are arrays")
        alphaimage = np.empty_like(conv_cube[:,:,0])
        alphaimage_error = np.empty_like(conv_cube[:,:,0])

        alphaimage[:] = np.nan
        alphaimage_error[:] = np.nan

        # x = np.log10(freqs/1e9)
        x = freqs/1e9

        count = 0
        pcount = 0

        # What is fitted: one spectrum per pixel (default), or one mean
        # spectrum per spatial bin (see _prepare_fit_spectra).
        (bin_labels, fit_idx, fit_cube, fit_cube_res,
         beam_area_px, _) = _prepare_fit_spectra(
            imagelist, conv_cube_res, masked_cube, masked_cube_res, mask,
            weighting=weighting, noise_rms=noise_rms,
            do_average_cube=do_average_cube, bin_size=bin_size,
            binning=binning, binning_kwargs=binning_kwargs,
            beam_area_px=beam_area_px, verbose=verbose)
        nspec = len(fit_idx)
        tenpc = int(nspec / 10.0)

        chunk_results = _fit_pixels_in_chunks(
            fit_spectra_chunk, fit_idx, x, fit_cube, fit_cube_res,
            (nu0, flux_sys_error_frac, weighting, robust_scale),
            n_jobs=n_jobs, chunk_size=chunk_size, verbose=verbose,
            label='specidx_map' if binning is None else f'specidx_map [{binning} binning]',
            unit='pixels' if binning is None else 'bins')

        alpha_values = np.empty(nspec)
        alpha_errs = np.empty(nspec)
        b_values = np.empty(nspec)
        k = 0
        for chunk_alpha, chunk_alpha_err, chunk_b in chunk_results:
            n = len(chunk_alpha)
            alpha_values[k:k + n] = chunk_alpha
            alpha_errs[k:k + n] = chunk_alpha_err
            b_values[k:k + n] = chunk_b
            k += n

        if bin_labels is not None:
            # every pixel of a bin gets the bin's fit
            in_bin = bin_labels >= 0
            alphaimage[in_bin] = alpha_values[bin_labels[in_bin]]
            alphaimage_error[in_bin] = alpha_errs[bin_labels[in_bin]]

        for p, (i, j) in enumerate(fit_idx):
            alpha_value = alpha_values[p]
            if bin_labels is None:
                alphaimage[i, j] = alpha_value
                alphaimage_error[i, j] = alpha_errs[p]

            count += 1
            if count == tenpc:
                count = 0
                pcount += 10
                if verbose>0:
                    # re-fit this pixel/bin in the parent (deterministic,
                    # same result) to get y, yerr and the report for the plot
                    b_value = b_values[p]
                    y = fit_cube[i, j, :] * 1000
                    yerr = np.sqrt((1*fit_cube_res[i, j, :])**2.0 + (flux_sys_error_frac * fit_cube[i, j, :])**2.0) * 1000
                    results_fit = do_fit_spec_map(x, y, yerr, nu0,
                                                  weighting=weighting,
                                                  robust_scale=robust_scale)
                    model_best = linear_function(x,
                                                alpha_value, 
                                                b_value,
                                                nu0)
                    plt.figure(figsize=(4, 5))

                    plt.errorbar(x, 
                                y, 
                                yerr=abs(yerr), 
                                fmt='o',
                                label='Data', color='k', ecolor='gray',
                                alpha=0.5)

                    plt.plot(x, model_best,
                            color='red', ls='-.', label='Best-fit model')

                    plt.ylim(abs(np.nanmin(y)*0.1),np.nanmax(y)*100)
                    plt.semilogy()
                    plt.semilogx()
                    plt.xlabel('Frequency [GHz]')
                    plt.ylabel('Pixel Flux Density [mJy/Beam]' if bin_labels is None
                               else f'Bin {i} mean Flux Density [mJy/Beam]')
                    plt.legend()
                    plt.show()
                    print(lmfit.fit_report(results_fit.params))

        
        if isinstance(imagelist[0], str) and save_fits:
            """
            Only save the alpha image if the input is a list of strings.
            """
            # binned and chi2 maps get their own names, never overwriting the
            # default (per-pixel, 'log') ones
            tag = _spec_output_tag(binning, weighting, use_model_images,
                                   add_residual, smooth_residual_to_target)
            if ref_image_conv is None:
                ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                                    return_df = False)[0][::-1][0]

            alphaimage_name = _spec_output_name(ref_image_conv, 'alpha', tag, out_dir)
            alphaimage_error_name = _spec_output_name(ref_image_conv, 'alphaerror', tag, out_dir)

            pf.writeto(alphaimage_name,alphaimage,overwrite=True)
            pf.writeto(alphaimage_error_name,alphaimage_error,overwrite=True)
            copy_header(ref_image_conv,alphaimage_name,alphaimage_name)
            copy_header(ref_image_conv,alphaimage_error_name,alphaimage_error_name)
            
            if bin_labels is not None:
                bins_name = _spec_output_name(ref_image_conv, 'alphabins', tag, out_dir)
                pf.writeto(bins_name,bin_labels.astype(np.int32),overwrite=True)
                copy_header(ref_image_conv,bins_name,bins_name)

            if use_model_images:
                # the maps are at the target beam, not at the header beam of
                # ref_image_conv
                for _name in [alphaimage_name, alphaimage_error_name] + (
                        [bins_name] if bin_labels is not None else []):
                    with pf.open(_name, mode='update') as hdul:
                        hdul[0].header['BMAJ'] = target_beam[0] / 3600.0
                        hdul[0].header['BMIN'] = target_beam[1] / 3600.0
                        hdul[0].header['BPA'] = target_beam[2]

        if return_bin_map:
            return(alphaimage, alphaimage_error,conv_cube,masked_cube_res,masked_cube,
                   bin_labels)
        return(alphaimage, alphaimage_error,conv_cube,masked_cube_res,masked_cube)
        
        
    
    
def specidx_map_SY_FF(imagelist,residuallist,
                      ref_image_conv=None,
                      freqs=None,
                      ref_image_mask = None,
                      flux_sys_error_frac = 0.1,
                      nu0=10.0,
                      mask=None,sigma_global_mask=6,
                      iterations=1,
                      dilation_size=2,
                      sed_model='S2',
                      fix_alpha_nt=False,
                      needs_convolution=False,conv_task='fft',
                      do_average_cube = False,bin_size=int(2),
                      n_jobs=1,
                      chunk_size=None,
                      verbose=0):
    """
    Earlier per-pixel synchrotron + free-free decomposition, kept for
    existing notebooks. `t_specidx_map_SY_FF` is the maintained version
    (spatial binning, 'chi2' weighting, errors with the parameter
    covariance, model-restored cubes).

    Known limitations of the errors returned here: they are missing (NaN)
    where an amplitude sits on its bound, and S_tot_33_err, Snu0_err and
    f_th_33_err ignore the anti-correlation between A_sy and A_ff.
    """
    if sed_model != 'S2':
        raise NotImplementedError(f"sed_model={sed_model!r}: only 'S2' "
                                  "(synchrotron + free-free) is implemented")

    if isinstance(imagelist[0], str):
        cube_image = makecube(imagelist)
        cube_residuals = makecube(residuallist)
    else:
        cube_image = imagelist
        cube_residuals = residuallist
    
    if mask is None:
        if ref_image_mask is not None:
            ref_image_mask = ref_image_mask
        else:
            ref_image_mask = imagelist[-1]
        _,mask = mask_dilation(ref_image_mask,
                               rms=mad_std(load_fits_data(residuallist[-1])),
                                     show_figure=True,
                                     PLOT=True,
                                     iterations=iterations,
                                     dilation_size=dilation_size,
                                     sigma=sigma_global_mask)
    mask_3d = mask[:, :, np.newaxis]
    inv_mask = ~mask
    inv_mask_3d = inv_mask[:, :, np.newaxis]
    
    # Reference image (largest beam by default): used for the convolution and
    # for the output names. It used to be resolved only when
    # needs_convolution=False (and a given ref_image_conv was then replaced by
    # imagelist[0]), so needs_convolution=True without ref_image_conv failed,
    # and so did array inputs.
    if ref_image_conv is None and isinstance(imagelist[0], str):
        ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                            return_df = False)[0][::-1][0]

    
    if needs_convolution and isinstance(ref_image_conv, str):
        psf_zise = int(get_beam_size_px(ref_image_conv)[0])
        psf_image_size = int(load_fits_data(ref_image_conv).shape[0])
        print(f"PSF BEAM SIZE is --> {psf_zise} px")
        print(f"PSF IMAGE SIZE is --> {psf_image_size} px")
        larger_beam_image_data = load_fits_data(ref_image_conv)
        

        psf_name = tcreate_beam_psf(ref_image_conv,
                                        size=(psf_image_size, psf_image_size),
                                        aspect = 'elliptical')
        PSF_DATA = load_fits_data(psf_name)
        
    num_images = cube_image.shape[2]
    if needs_convolution:
        conv_cube = np.empty_like(cube_image)
        conv_cube_res = np.empty_like(cube_residuals)
        if conv_task == 'fft':
            for i in tqdm(range(num_images)):
                conv_image_i_uncor = scipy.signal.fftconvolve(cube_image[:, :, i], 
                                                            PSF_DATA, 'same')
                conv_residual_i_uncor = scipy.signal.fftconvolve(cube_residuals[:, :, i], 
                                                            PSF_DATA, 'same')
                
                larger_beam_area = beam_area2(ref_image_conv)
                beam_area_i = beam_area2(imagelist[i])
                
                # factor_conv = np.max(larger_beam_image_data)/np.max(conv_image_i_uncor)
                factor_conv_i = larger_beam_area/beam_area_i
                
                print(f"Factor Convolution is --> {factor_conv_i}")
                # factor_conv = 1.0

                conv_cube[:, :, i] = conv_image_i_uncor * factor_conv_i
                conv_image_i = conv_cube[:, :, i]
                
                conv_cube_res[:, :, i] = conv_residual_i_uncor * factor_conv_i
                conv_res_i = conv_cube_res[:, :, i]
                
                conv_name = imagelist[i].replace('-image.','-image-conv.').replace('-image.cutout','-image-conv.cutout')
                pf.writeto(conv_name,
                                conv_image_i,
                                overwrite=True)
                copy_header(ref_image_conv,conv_name,conv_name)
                
                conv_res_name = residuallist[i].replace('-residual.','-residual-conv.').replace('-residual.cutout','-residual-conv.cutout')
                pf.writeto(conv_res_name,
                                conv_res_i,
                                overwrite=True)
                copy_header(ref_image_conv,conv_res_name,conv_res_name)
        # if conv_task == 'imsmooth':
        #     if isinstance(imagelist[0], str):
        #         for i in tqdm(range(num_images)):
        #             conv_name = imagelist[i].replace('-image.','-image-conv.').replace('-image.cutout','-image-conv.cutout')
        #             pf.writeto(conv_name,
        #                             cube_image[:, :, i],
        #                             overwrite=True)
        #             copy_header(imagelist[0],conv_name,conv_name)
                    
        #             convolve_2D_smooth(imagelist[0],conv_name,mode='transfer',add_prefix='')
        #             conv_cube[:,:,i] = load_fits_data(conv_name)
            
    else:
        conv_cube = cube_image.copy()
        conv_cube_res = cube_residuals.copy()
    
    # conv_cube_res = cube_residuals.copy()
    masked_cube = np.where(inv_mask_3d, np.nan, conv_cube)
    masked_cube_res = np.where(inv_mask_3d, np.nan, conv_cube_res)

    if do_average_cube:
        # bin_size = 2
        # assert masked_cube.shape[2] % bin_size == 0, "Frequency axis size must be divisible by bin size."
        reshaped_cube = masked_cube.reshape(masked_cube.shape[0], 
                                            masked_cube.shape[1], -1, 
                                            bin_size)
        averaged_cube = reshaped_cube.mean(axis=3)

        reshaped_residual = masked_cube_res.reshape(masked_cube_res.shape[0], 
                                                    masked_cube_res.shape[1], -1, 
                                                    bin_size)
        averaged_cube_res = np.sqrt(np.sum(reshaped_residual**2, axis=3)) / bin_size


        print("Original cube shape:", masked_cube.shape)
        print("Averaged cube shape:", averaged_cube.shape)
        reshaped_frequencies = freqs.reshape(-1, bin_size)
        averaged_frequencies = reshaped_frequencies.mean(axis=1)
        freqs = averaged_frequencies
        masked_cube = averaged_cube.astype(np.float32)
        masked_cube_res = averaged_cube_res.astype(np.float32)
    else:        
        masked_cube = masked_cube.astype(np.float32)
        masked_cube_res = masked_cube_res.astype(np.float32)
    
    
    # cube_image * mask_3d
    idx = np.column_stack(np.where(mask==True))
    if freqs is None:
        freqs = getfreqs(imagelist)
    alphaimage = np.empty_like(conv_cube[:,:,0])
    alphaimage_error = np.empty_like(conv_cube[:,:,0])
    A_sy_map = np.empty_like(conv_cube[:,:,0])
    A_sy_map_err = np.empty_like(conv_cube[:,:,0])
    A_ff_map = np.empty_like(conv_cube[:,:,0])
    A_ff_map_err = np.empty_like(conv_cube[:,:,0])
    Snu0 = np.empty_like(conv_cube[:,:,0])
    Snu0_err = np.empty_like(conv_cube[:,:,0])
    f_th_33 = np.empty_like(conv_cube[:,:,0])
    f_th_33_err = np.empty_like(conv_cube[:,:,0])
    S_tot_33 = np.empty_like(conv_cube[:,:,0])
    S_tot_33_err = np.empty_like(conv_cube[:,:,0])
    S_sy_33 = np.empty_like(conv_cube[:,:,0])
    S_sy_33_err = np.empty_like(conv_cube[:,:,0])
    S_ff_33 = np.empty_like(conv_cube[:,:,0])
    S_ff_33_err = np.empty_like(conv_cube[:,:,0])
    
    sy_map_33 = np.empty_like(conv_cube[:,:,0])
    sy_map_33_err = np.empty_like(conv_cube[:,:,0])
    ff_map_33 = np.empty_like(conv_cube[:,:,0])
    ff_map_33_err = np.empty_like(conv_cube[:,:,0])

    
    alphaimage[:] = np.nan
    alphaimage_error[:] = np.nan
    sy_map_33[:] = np.nan
    sy_map_33_err[:] = np.nan
    ff_map_33[:] = np.nan
    ff_map_33_err[:] = np.nan
    A_sy_map[:] = np.nan
    A_sy_map_err[:] = np.nan
    A_ff_map[:] = np.nan
    A_ff_map_err[:] = np.nan
    f_th_33[:] = np.nan
    f_th_33_err[:] = np.nan
    Snu0[:] = np.nan
    Snu0_err[:] = np.nan
    S_tot_33[:] = np.nan
    S_tot_33_err[:] = np.nan
    S_sy_33[:] = np.nan
    S_sy_33_err[:] = np.nan
    S_ff_33[:] = np.nan
    S_ff_33_err[:] = np.nan

    
#     x = np.log10(freqs)
    x = np.asarray(freqs, dtype=float)/1e9
    
    nspec = len(idx)
    tenpc = int(nspec / 10.0)
    count = 0
    pcount = 0
    
    masked_cube = masked_cube.astype(np.float32)
    masked_cube_res = masked_cube_res.astype(np.float32)
    
    
    chunk_results = _fit_pixels_in_chunks(
        fit_spectra_SY_FF_chunk, idx, x, masked_cube, masked_cube_res,
        (nu0, flux_sys_error_frac, fix_alpha_nt),
        n_jobs=n_jobs, chunk_size=chunk_size, verbose=verbose,
        label='specidx_map_SY_FF')
    # one (alpha_nt, alpha_nt_err, A_sy, A_sy_err, A_ff, A_ff_err) per pixel;
    # a missing stderr (None, e.g. a parameter on its bound) becomes NaN
    results = [(i, j, *(np.nan if v is None else v for v in pixel_fit))
               for (i, j), pixel_fit in
               zip(idx, (p for chunk in chunk_results for p in chunk))]

    if sed_model == 'S2':
        # for i, j in tqdm(idx):
        # for i, j in idx:
        for (i, j, alpha_nt_value, alpha_nt_err,
             A_sy_value,A_sy_err,
             A_ff_value,A_ff_err) in results:
            # if count == 0:
            #     print(str(pcount) + '%...')
            
            
    #         y = np.log10(masked_cube[i,j,:])
            # y = masked_cube[i,j,:]*1000
            # yerr = masked_cube_res[i,j,:]*1000
        
            # results_fit = do_fit_spec_SY_FF_map(x,y,yerr,nu0)
            
            alphaimage[i, j] = alpha_nt_value
            alphaimage_error[i, j] = alpha_nt_err
            A_sy_map[i, j] = A_sy_value
            A_sy_map_err[i, j] = A_sy_err
            A_ff_map[i, j] = A_ff_value
            A_ff_map_err[i, j] = A_ff_err
            
            sy_map_33[i,j] = A_sy_value * (33/nu0)**alpha_nt_value
            ff_map_33[i,j] = A_ff_value * (33/nu0)**(-0.1)
            sy_map_33_err[i,j] = sy_map_33[i,j] * np.sqrt((A_sy_err/A_sy_value)**2.0 + (alpha_nt_err * np.log(33/nu0))**2.0)
            # the free-free index (-0.1) is fixed: it carries no error
            ff_map_33_err[i,j] = ff_map_33[i,j] * np.sqrt((A_ff_err/A_ff_value)**2.0)
            
            Snu0[i,j] = RC_function_SY_FF(nu0, 
                                    A_sy_value, 
                                    A_ff_value, 
                                    alpha_nt_value,
                                    nu0)
            S_tot_33[i,j] = RC_function_SY_FF(33, 
                                    A_sy_value, 
                                    A_ff_value, 
                                    alpha_nt_value,
                                    nu0)
            
            S_tot_33_err[i,j] = RC_function_SY_FF(33, 
                                    A_sy_err, 
                                    A_ff_err, 
                                    alpha_nt_value,
                                    nu0)
            
            S_sy_33_px = RC_function_SY_FF(33, 
                                    A_sy_value,
                                    0.0, 
                                    alpha_nt_value,
                                    nu0)
            S_sy_33_err_px = RC_function_SY_FF(33, 
                                    A_sy_err, 
                                    0.0, 
                                    alpha_nt_value,
                                    nu0)
            
            S_ff_33_px = RC_function_SY_FF(33, 
                                    0.0, 
                                    A_ff_value, 
                                    0.0,
                                    nu0)
            S_ff_33_err_px = RC_function_SY_FF(33, 
                                    0.0, 
                                    A_ff_err, 
                                    0.0,
                                    nu0)
            
            S_sy_33[i,j] = S_sy_33_px
            S_sy_33_err[i,j] = S_sy_33_err_px
            S_ff_33[i,j] = S_ff_33_px
            S_ff_33_err[i,j] = S_ff_33_err_px
            
            f_th_33[i,j] = S_ff_33_px / S_tot_33[i,j]
            
            f_th_33_err[i,j] = f_th_33[i,j] * np.sqrt((S_ff_33_err_px/S_ff_33_px)**2.0 + (S_tot_33_err[i,j]/(S_ff_33_px+S_sy_33_px))**2.0)
            
            try:
                Snu0_err[i,j] = RC_function_SY_FF(nu0, 
                                        A_sy_err, 
                                        A_ff_err, 
                                        alpha_nt_value,
                                        nu0)
            except:
                Snu0_err[i,j] = np.nan                 
                
            count += 1
            if count == tenpc:
                count = 0
                pcount += 10
                if verbose>0:
                    # re-fit this pixel in the parent (deterministic, same
                    # result) to get y, yerr and the report for the plot
                    y = masked_cube[i, j, :] * 1000
                    yerr = np.sqrt((1*masked_cube_res[i, j, :])**2.0 + (flux_sys_error_frac * masked_cube[i, j, :])**2.0) * 1000
                    results_fit = do_fit_spec_SY_FF_map(x, y, yerr, nu0,fix_alpha_nt)
                    model_best = RC_function_SY_FF(x,
                                                A_sy_value,
                                                A_ff_value,
                                                alpha_nt_value,
                                                nu0
                                                )

                    model_best_sy = RC_function_SY_FF(x, 
                                                A_sy_value, 
                                                0.0,
                                                alpha_nt_value,
                                                nu0
                                                )
                    model_best_ff = RC_function_SY_FF(x, 
                                                0.0, 
                                                A_ff_value,
                                                alpha_nt_value,
                                                nu0
                                                )
                    
                    
                    plt.figure(figsize=(4, 5))

                    plt.errorbar(x, 
                                y, 
                                yerr=abs(yerr), 
                                fmt='o',
                                label='Data', color='k', ecolor='gray',
                                alpha=0.5)
                    plt.plot(x, model_best,
                            color='red', ls='-.', label='Best-fit model')
                    plt.plot(x, model_best_sy,
                            color='blue', ls='-.', label='SY term')
                    plt.plot(x, model_best_ff,
                            color='orange', ls='-.', label='FF term')
                    
                    
                    plt.ylim(abs(np.nanmin(y)*0.1),np.nanmax(y)*100)
                    plt.semilogy()
                    plt.semilogx()
                    plt.xlabel('Frequency [GHz]')
                    plt.ylabel('Pixel Flux Density [mJy/Beam]')
                    plt.legend()
                    plt.show()
                    print(lmfit.fit_report(results_fit.params))
        
    
    # 33 GHz decomposition maps (all entries are 2D maps)
    sed_sy_ff = {'S_sy_33':S_sy_33,
                 'S_sy_33_err':S_sy_33_err,
                 'S_ff_33':S_ff_33,
                 'S_ff_33_err':S_ff_33_err,
                 'sy_map_33':sy_map_33,
                 'sy_map_33_err':sy_map_33_err,
                 'ff_map_33':ff_map_33,
                 'ff_map_33_err':ff_map_33_err,
                 'S_tot_33':S_tot_33,
                 'S_tot_33_err':S_tot_33_err}
    
    if isinstance(imagelist[0], str):
        """
        Only save the alpha image if the input is a list of strings.
        """
        # alphaimage_name = imagelist[0].replace('-0000-image.fits','-alpha.fits')
        # alphaimage_error_name = imagelist[0].replace('-0000-image.fits','-alphaerror.fits')
        _alphaimage_name = (ref_image_conv.replace('-image.','-alpha_nt.').
                        replace('-image.cutout','-alpha_nt.cutout').
                        replace('-image-pb','-alpha_nt-pb'))
        
        _alphaimage_error_name = (ref_image_conv.replace('-image.','-alpha_nt_error.').
                                replace('-image.cutout','-alpha_nt_error.cutout').
                                replace('-image-pb','-alpha_nt_error-pb'))
        
        alphaimage_name = os.path.dirname(os.path.dirname(_alphaimage_name))+'/'+os.path.basename(_alphaimage_name)
        alphaimage_error_name = os.path.dirname(os.path.dirname(_alphaimage_error_name))+'/'+os.path.basename(_alphaimage_error_name)

        pf.writeto(alphaimage_name,alphaimage,overwrite=True)
        pf.writeto(alphaimage_error_name,alphaimage_error,overwrite=True)
        copy_header(ref_image_conv,alphaimage_name,alphaimage_name)
        copy_header(ref_image_conv,alphaimage_error_name,alphaimage_error_name)
    
    return(alphaimage, alphaimage_error,Snu0,Snu0_err,
           A_ff_map, A_ff_map_err, 
           A_sy_map, A_sy_map_err, 
           f_th_33,f_th_33_err,
           conv_cube,masked_cube_res,masked_cube,sed_sy_ff)
    







# import numpy as np
# import os
# from tqdm import tqdm
# from joblib import Parallel, delayed
# import astropy.io.fits as pf
# import scipy.signal
# import matplotlib.pyplot as plt
# import lmfit

def create_frequency_maps(A_sy, A_sy_err, A_ff, A_ff_err, alpha_nt, alpha_nt_err,
                         frequencies, nu0=10.0,
                         cov_sy_ff=None, cov_sy_alpha=None, cov_ff_alpha=None):
    """
    Calculate synchrotron and free-free maps at specified frequencies.

    With S = A_sy*(nu/nu0)**alpha_nt and F = A_ff*(nu/nu0)**(-0.1), the
    total is T = S + F and the thermal fraction is f_th = F/T. Errors are
    propagated to first order with the covariance of the fitted parameters:
    A_sy and A_ff are strongly anti-correlated, so the errors of T and f_th
    are not the quadrature sum of the errors of S and F. The free-free index
    (-0.1) is fixed and carries no error.

    Parameters:
    -----------
    A_sy, A_sy_err : array_like
        Synchrotron amplitude and error maps
    A_ff, A_ff_err : array_like
        Free-free amplitude and error maps
    alpha_nt, alpha_nt_err : array_like
        Non-thermal spectral index and error maps (error 0 when it is fixed)
    frequencies : array_like
        Target frequencies in GHz
    nu0 : float
        Reference frequency in GHz
    cov_sy_ff, cov_sy_alpha, cov_ff_alpha : array_like, optional
        Covariance maps cov(A_sy, A_ff), cov(A_sy, alpha_nt) and
        cov(A_ff, alpha_nt). None means uncorrelated.

    Returns:
    --------
    dict : Dictionary containing maps for each frequency ('<freq>GHz'):
        synchrotron, free_free, total, thermal_fraction and their '_err'.
    """
    A_sy = np.asarray(A_sy, dtype=float)
    A_ff = np.asarray(A_ff, dtype=float)
    alpha_nt = np.asarray(alpha_nt, dtype=float)
    var_sy = np.asarray(A_sy_err, dtype=float) ** 2
    var_ff = np.asarray(A_ff_err, dtype=float) ** 2
    var_al = np.asarray(alpha_nt_err, dtype=float) ** 2
    zero = np.zeros_like(A_sy)
    c_sf = zero if cov_sy_ff is None else np.asarray(cov_sy_ff, dtype=float)
    c_sa = zero if cov_sy_alpha is None else np.asarray(cov_sy_alpha, dtype=float)
    c_fa = zero if cov_ff_alpha is None else np.asarray(cov_ff_alpha, dtype=float)

    results = {}

    for freq in frequencies:
        freq_key = f'{freq:.1f}GHz'
        r = freq / nu0
        lnr = np.log(r)
        with np.errstate(divide='ignore', invalid='ignore'):
            r_nt = r ** alpha_nt
            r_ff = r ** (-0.1)

            # Synchrotron component
            sy_map = A_sy * r_nt
            var_S = (r_nt ** 2 * var_sy + (sy_map * lnr) ** 2 * var_al
                     + 2.0 * r_nt * sy_map * lnr * c_sa)

            # Free-free component (fixed spectral index = -0.1)
            ff_map = A_ff * r_ff
            var_F = r_ff ** 2 * var_ff

            # Total emission
            cov_SF = r_nt * r_ff * c_sf + sy_map * lnr * r_ff * c_fa
            total_map = sy_map + ff_map
            var_T = var_S + var_F + 2.0 * cov_SF

            # Thermal fraction f = F/(S+F): df/dS = -F/T**2, df/dF = S/T**2
            thermal_fraction = ff_map / total_map
            var_f = (ff_map ** 2 * var_S + sy_map ** 2 * var_F
                     - 2.0 * ff_map * sy_map * cov_SF) / total_map ** 4

            results[freq_key] = {
                'synchrotron': sy_map,
                'synchrotron_err': np.sqrt(np.maximum(var_S, 0.0)),
                'free_free': ff_map,
                'free_free_err': np.sqrt(np.maximum(var_F, 0.0)),
                'total': total_map,
                'total_err': np.sqrt(np.maximum(var_T, 0.0)),
                'thermal_fraction': thermal_fraction,
                'thermal_fraction_err': np.sqrt(np.maximum(var_f, 0.0))
            }

    return results

def prepare_input_data(imagelist, residuallist, freqs=None):
    """
    Prepare input cubes from image and residual lists.
    
    Parameters:
    -----------
    imagelist : list or array
        List of image filenames or pre-loaded cube
    residuallist : list or array  
        List of residual filenames or pre-loaded cube
    freqs : array_like, optional
        Frequency array
        
    Returns:
    --------
    tuple : (cube_image, cube_residuals, frequencies)
    """
    if isinstance(imagelist[0], str):
        cube_image = makecube(imagelist)
        cube_residuals = makecube(residuallist)
        if freqs is None:
            freqs = getfreqs(imagelist)
    else:
        cube_image = imagelist
        cube_residuals = residuallist
        if freqs is None:
            raise ValueError("Frequencies must be provided when using pre-loaded cubes")
    
    return cube_image, cube_residuals, freqs

def create_processing_mask(ref_image_mask, residuallist, mask=None, 
                          sigma_global_mask=6, iterations=1, dilation_size=2):
    """
    Create or validate processing mask.
    
    Parameters:
    -----------
    ref_image_mask : str or None
        Reference image for mask creation
    residuallist : list
        List of residual images
    mask : array_like, optional
        Pre-existing mask
    sigma_global_mask : float
        Sigma threshold for mask creation
    iterations : int
        Dilation iterations
    dilation_size : int
        Dilation kernel size
        
    Returns:
    --------
    tuple : (mask, inverse_mask_3d)
    """
    if mask is None:
        if ref_image_mask is None:
            ref_image_mask = residuallist[-1]
        
        _, mask = mask_dilation(
            ref_image_mask,
            rms=mad_std(load_fits_data(residuallist[-1])),
            show_figure=True,
            PLOT=True,
            iterations=iterations,
            dilation_size=dilation_size,
            sigma=sigma_global_mask
        )
    
    mask_3d = mask[:, :, np.newaxis]
    inv_mask = ~mask
    inv_mask_3d = inv_mask[:, :, np.newaxis]
    
    return mask, inv_mask_3d

def convolve_images_to_common_resolution(cube_image, cube_residuals, imagelist, 
                                       residuallist, ref_image_conv, conv_task='fft'):
    """
    Convolve all images to common resolution.
    
    Parameters:
    -----------
    cube_image, cube_residuals : array_like
        Input image and residual cubes
    imagelist, residuallist : list
        File paths for saving convolved images
    ref_image_conv : str
        Reference image for target resolution
    conv_task : str
        Convolution method ('fft')
        
    Returns:
    --------
    tuple : (convolved_cube, convolved_residuals_cube)
    """
    # Get PSF information
    psf_size = int(get_beam_size_px(ref_image_conv)[0])
    psf_image_size = int(load_fits_data(ref_image_conv).shape[0])
    larger_beam_image_data = load_fits_data(ref_image_conv)
    
    print(f"PSF BEAM SIZE is --> {psf_size} px")
    print(f"PSF IMAGE SIZE is --> {psf_image_size} px")
    
    # Create PSF
    psf_name = tcreate_beam_psf(
        ref_image_conv,
        size=(psf_image_size, psf_image_size),
        aspect='elliptical'
    )
    PSF_DATA = load_fits_data(psf_name)
    
    # Initialize output cubes
    num_images = cube_image.shape[2]
    conv_cube = np.empty_like(cube_image)
    conv_cube_res = np.empty_like(cube_residuals)
    
    if conv_task == 'fft':
        larger_beam_area = beam_area2(ref_image_conv)
        
        for i in tqdm(range(num_images)):
            # Convolve images
            conv_image_i_uncor = scipy.signal.fftconvolve(
                cube_image[:, :, i], PSF_DATA, 'same'
            )
            conv_residual_i_uncor = scipy.signal.fftconvolve(
                cube_residuals[:, :, i], PSF_DATA, 'same'
            )
            
            # Apply beam area correction
            beam_area_i = beam_area2(imagelist[i])
            factor_conv_i = larger_beam_area / beam_area_i
            print(f"Factor Convolution is --> {factor_conv_i}")
            
            conv_cube[:, :, i] = conv_image_i_uncor * factor_conv_i
            conv_cube_res[:, :, i] = conv_residual_i_uncor * factor_conv_i
            
            # Save convolved images if input was file paths
            if isinstance(imagelist[0], str):
                save_convolved_image(imagelist[i], conv_cube[:, :, i], ref_image_conv)
                save_convolved_residual(residuallist[i], conv_cube_res[:, :, i], ref_image_conv)
    
    return conv_cube, conv_cube_res

def save_convolved_image(original_path, convolved_data, ref_image):
    """Save convolved image with proper naming and header."""
    conv_name = original_path.replace('-image.', '-image-conv.').replace('-image.cutout', '-image-conv.cutout')
    pf.writeto(conv_name, convolved_data, overwrite=True)
    copy_header(ref_image, conv_name, conv_name)

def save_convolved_residual(original_path, convolved_data, ref_image):
    """Save convolved residual with proper naming and header."""
    conv_res_name = original_path.replace('-residual.', '-residual-conv.').replace('-residual.cutout', '-residual-conv.cutout')
    pf.writeto(conv_res_name, convolved_data, overwrite=True)
    copy_header(ref_image, conv_res_name, conv_res_name)

def apply_frequency_averaging(masked_cube, masked_cube_res, freqs, bin_size=2):
    """
    Apply frequency averaging to reduce noise.
    
    Parameters:
    -----------
    masked_cube, masked_cube_res : array_like
        Masked image and residual cubes
    freqs : array_like
        Frequency array
    bin_size : int
        Number of channels to average
        
    Returns:
    --------
    tuple : (averaged_cube, averaged_residuals, averaged_frequencies)
    """
    reshaped_cube = masked_cube.reshape(
        masked_cube.shape[0], masked_cube.shape[1], -1, bin_size
    )
    averaged_cube = reshaped_cube.mean(axis=3)
    
    reshaped_residual = masked_cube_res.reshape(
        masked_cube_res.shape[0], masked_cube_res.shape[1], -1, bin_size
    )
    averaged_cube_res = np.sqrt(np.sum(reshaped_residual**2, axis=3)) / bin_size
    
    reshaped_frequencies = freqs.reshape(-1, bin_size)
    averaged_frequencies = reshaped_frequencies.mean(axis=1)
    
    print("Original cube shape:", masked_cube.shape)
    print("Averaged cube shape:", averaged_cube.shape)
    
    return averaged_cube.astype(np.float32), averaged_cube_res.astype(np.float32), averaged_frequencies

def initialize_output_maps(image_shape):
    """
    Initialize all output parameter maps with NaN values.
    
    Parameters:
    -----------
    image_shape : tuple
        Shape of output maps (ny, nx)
        
    Returns:
    --------
    dict : Dictionary of initialized maps
    """
    maps = {}
    
    # Primary fit parameters
    parameter_names = [
        'alpha_nt', 'alpha_nt_err',
        'A_sy', 'A_sy_err', 
        'A_ff', 'A_ff_err',
        'S_nu0', 'S_nu0_err'
    ]
    
    for name in parameter_names:
        maps[name] = np.full(image_shape, np.nan, dtype=np.float32)
    
    return maps

def fit_single_pixel_spectrum(i, j, frequencies, masked_cube, masked_cube_res, 
                             nu0, flux_sys_error_frac, fix_alpha_nt):
    """
    Fit SED model to single pixel spectrum.
    
    Parameters:
    -----------
    i, j : int
        Pixel coordinates
    frequencies : array_like
        Frequency array in GHz
    masked_cube, masked_cube_res : array_like
        Data and residual cubes
    nu0 : float
        Reference frequency
    flux_sys_error_frac : float
        Systematic flux error fraction
    fix_alpha_nt : bool
        Whether to fix non-thermal spectral index
        
    Returns:
    --------
    tuple : Pixel coordinates and fit results
    """
    y = masked_cube[i, j, :] * 1000  # Convert to mJy
    yerr = np.sqrt(
        (1 * masked_cube_res[i, j, :])**2 + 
        (flux_sys_error_frac * masked_cube[i, j, :])**2
    ) * 1000
    
    results_fit = do_fit_spec_SY_FF_map(frequencies, y, yerr, nu0, fix_alpha_nt)
    
    return (i, j, results_fit, y, yerr)

def extract_fit_parameters(results_fit):
    """
    Extract fit parameters from lmfit results.
    
    Parameters:
    -----------
    results_fit : lmfit.ModelResult
        Fit results object
        
    Returns:
    --------
    dict : Dictionary of parameter values and errors
    """
    params = {}
    
    for param_name in ['alpha_nt', 'A_sy', 'A_ff']:
        param = results_fit.params[param_name]
        params[param_name] = param.value
        params[f'{param_name}_err'] = param.stderr if param.stderr is not None else np.nan
    
    return params

def save_parameter_maps(output_maps, ref_image_conv, imagelist, tag='',
                        frequency_maps=None, bin_labels=None, target_beam=None,
                        out_dir=None):
    """
    Save parameter maps to FITS files.

    Parameters:
    -----------
    output_maps : dict
        Dictionary of parameter maps
    ref_image_conv : str
        Reference image for header copying (and for the output names)
    imagelist : list
        List of input images; nothing is written for array inputs
    tag : str
        Tag added to every file name (binning, weighting, ...), so different
        kinds of runs do not overwrite each other.
    frequency_maps : dict, optional
        Output of `create_frequency_maps`; the thermal fraction and its
        error are saved for each frequency.
    bin_labels : int array, optional
        Bin label map (saved when given).
    target_beam : tuple, optional
        (bmaj_arcsec, bmin_arcsec, bpa_deg) to stamp in the headers when the
        maps are not at the beam of ref_image_conv (use_model_images).
    out_dir : str, optional
        Directory for the files. Default: one directory above the images.

    Returns:
    --------
    list : names of the files written
    """
    if not isinstance(imagelist[0], str) or not isinstance(ref_image_conv, str):
        return []

    # map key -> (file suffix, BUNIT)
    save_params = {
        'alpha_nt': ('alpha_nt', ''),
        'alpha_nt_err': ('alpha_nt_error', ''),
        'A_sy': ('A_sy', 'mJy/beam'),
        'A_sy_err': ('A_sy_error', 'mJy/beam'),
        'A_ff': ('A_ff', 'mJy/beam'),
        'A_ff_err': ('A_ff_error', 'mJy/beam'),
        'S_nu0': ('S_nu0', 'mJy/beam'),
        'S_nu0_err': ('S_nu0_error', 'mJy/beam'),
        'chi2_red': ('chi2red_syff', ''),
        'fit_flag': ('fitflag_syff', ''),
    }
    to_write = [(suffix, output_maps[key], bunit)
                for key, (suffix, bunit) in save_params.items()
                if key in output_maps]
    for freq_key, maps in (frequency_maps or {}).items():
        to_write += [
            (f'S_sy_{freq_key}', maps['synchrotron'], 'mJy/beam'),
            (f'S_sy_error_{freq_key}', maps['synchrotron_err'], 'mJy/beam'),
            (f'S_ff_{freq_key}', maps['free_free'], 'mJy/beam'),
            (f'S_ff_error_{freq_key}', maps['free_free_err'], 'mJy/beam'),
            (f'S_total_{freq_key}', maps['total'], 'mJy/beam'),
            (f'S_total_error_{freq_key}', maps['total_err'], 'mJy/beam'),
            (f'fth_{freq_key}', maps['thermal_fraction'], ''),
            (f'fth_error_{freq_key}', maps['thermal_fraction_err'], ''),
        ]
    if bin_labels is not None:
        to_write.append(('syffbins', bin_labels.astype(np.int32), ''))

    return _save_spec_maps(to_write, ref_image_conv, tag=tag,
                           target_beam=target_beam, out_dir=out_dir)

def t_specidx_map_SY_FF(imagelist, residuallist,
                      ref_image_conv=None,
                      freqs=None,
                      ref_image_mask=None,
                      flux_sys_error_frac=0.1,
                      nu0=10.0,
                      mask=None,
                      sigma_global_mask=6,
                      iterations=1,
                      dilation_size=2,
                      sed_model='S2',
                      fix_alpha_nt=False,
                      needs_convolution=False,
                      conv_task='fft',
                      do_average_cube=False,
                      bin_size=2,
                      n_jobs=1,
                      chunk_size=None,
                      verbose=0,
                      output_frequencies=(1.4, 6.0, 33.0),
                      binning=None,
                      binning_kwargs=None,
                      beam_area_px=None,
                      weighting='chi2',
                      robust_scale=3.0,
                      noise_rms=None,
                      alpha_nt_fixed=-0.85,
                      alpha_nt_bounds=(-2.5, -0.2),
                      use_model_images=False,
                      modellist=None,
                      target_beam=None,
                      add_residual=True,
                      smooth_residual_to_target=False,
                      residual_beams=None,
                      model_out_dir=None,
                      out_dir=None,
                      save_fits=True):
    """
    Resolved SED decomposition into synchrotron and free-free emission:

        S(nu) = A_sy*(nu/nu0)**alpha_nt + A_ff*(nu/nu0)**(-0.1)

    fitted to the spectrum of every pixel, or of every spatial bin, of a
    list of images at different frequencies. A_sy and A_ff are the
    synchrotron and free-free surface brightness at nu0 (both >= 0).

    The cube preparation (mask, convolution, model restoring, frequency
    averaging), the spatial binning and the noise model are the ones of
    `specidx_map`, so the arguments they share mean the same.

    Parameters:
    -----------
    imagelist : list or array
        Input images (FITS names) or image cube (ny, nx, nfreq), Jy/beam
    residuallist : list or array
        Residual images or residual cube
    ref_image_conv : str, optional
        Reference image for convolution and for the output names/headers.
        Default: the image with the largest beam.
    freqs : array_like, optional
        Frequency array in Hz (read from the headers if None)
    ref_image_mask : str, optional
        Reference image for mask creation
    flux_sys_error_frac : float
        Fractional flux error added in quadrature to the noise of every
        point: yerr = sqrt(noise**2 + (flux_sys_error_frac*S)**2). It is
        treated as independent per band; with 'chi2' it enters the reported
        errors directly (set it to 0 for noise-only errors).
    nu0 : float
        Reference frequency in GHz
    mask : array_like, optional
        Processing mask
    sigma_global_mask, iterations, dilation_size :
        Used to build the mask when `mask` is None
    sed_model : str
        Only 'S2' (synchrotron + free-free) is implemented
    fix_alpha_nt : bool
        Fix the non-thermal spectral index at `alpha_nt_fixed`. The model is
        then linear in A_sy and A_ff, and the decomposition is much better
        constrained: with alpha_nt free, alpha_nt and A_ff are strongly
        degenerate over a limited frequency range.
    needs_convolution, conv_task :
        As in `specidx_map`
    do_average_cube, bin_size :
        Average the cube in groups of `bin_size` consecutive frequencies
    n_jobs : int
        Number of parallel jobs
    chunk_size : int, optional
        Spectra fitted per parallel task. Default (None) splits them into
        ~16 chunks per job.
    verbose : int
        Verbosity level; > 0 also shows ten example fits
    output_frequencies : sequence
        Frequencies for derived quantity maps (GHz)
    binning : str or None
        None: fit every pixel. 'contour': contour binning (see
        `specidx_map`); one mean spectrum per bin is fitted and its
        parameters are given to all pixels of the bin.
    binning_kwargs : dict
        Options of the binning method (target_snr, smooth_sigma_px,
        constraint, ref_map, rms). The thermal fraction comes from the
        curvature of the spectrum, so it needs more signal-to-noise than a
        single spectral index, in particular when alpha_nt is free.
    beam_area_px : float
        Beam area in pixels for the binning. Default: the largest beam of
        the inputs (1 for array inputs).
    weighting : str
        'chi2' (default): yerr = sqrt(noise_rms**2 + (flux_sys_error_frac*S)**2)
            is the 1-sigma error of each point, residuals (S - model)/yerr,
            robust (Cauchy) threshold of `robust_scale` sigma, and absolute
            parameter errors from the covariance at the best fit. These
            stay defined when A_ff (or A_sy) sits on its lower bound of 0,
            where the error is the scale of an upper limit.
        'log': the previous behaviour; the residual pixel value is the noise
            term, used as mild log-weights; errors are rescaled by the
            scatter about the fit and are missing (NaN) at the bounds.
    robust_scale : float
        'chi2' only: points beyond ~robust_scale sigma are down-weighted.
    noise_rms : None, array (nfreq,) or array (ny, nx, nfreq)
        'chi2' only: noise per band. None: mad_std of each residual image
        outside the mask.
    alpha_nt_fixed : float
        Value of alpha_nt when fix_alpha_nt=True.
    alpha_nt_bounds : (min, max)
        Bounds of alpha_nt when it is free. The maximum is kept below the
        free-free index (-0.1): above it the two components swap roles.
    use_model_images, modellist, target_beam, add_residual,
    smooth_residual_to_target, residual_beams, model_out_dir :
        Build the cube from the clean-component models restored with a
        common beam; see `specidx_map`.
    out_dir : str, optional
        Directory for the FITS outputs. Default (None): one directory above
        the input images.
    save_fits : bool
        False: write no FITS files.

    Returns:
    --------
    dict : Comprehensive results dictionary containing:
        - parameter_maps: alpha_nt, A_sy, A_ff, S_nu0 (total at nu0) and
          their '_err'; cov_sy_ff, cov_sy_alpha, cov_ff_alpha; chi2_red;
          n_bands; fit_flag. Amplitudes are in mJy/beam.
          fit_flag marks fits that ended on a bound: 1 = A_sy at 0,
          2 = A_ff at 0 (no free-free detected; A_ff_err is then the scale
          of an upper limit), 4 = alpha_nt at one of alpha_nt_bounds (the
          decomposition is then not constrained by the data and the errors
          are not meaningful); the values add up.
        - frequency_maps: synchrotron, free_free, total, thermal_fraction
          (and '_err') at each of output_frequencies
        - data_products: conv_cube, masked_cube, masked_cube_res, mask,
          bin_map (None without binning), noise_rms (per band, 'chi2')
        - metadata: processing information, including the files written

    FITS outputs (file inputs only) get a tag for the binning, the weighting
    and a fixed alpha_nt ('-alphafix'), e.g.
    '<name>-fth_33.0GHz-contourbin-chi2-alphafix.fits'.
    """
    if sed_model != 'S2':
        raise NotImplementedError(f"sed_model={sed_model!r}: only 'S2' "
                                  "(synchrotron + free-free) is implemented")

    # 1-4. Cubes, mask, convolution / model restoring, frequency averaging
    prep = _prepare_spec_cube(
        imagelist, residuallist, ref_image_conv=ref_image_conv, freqs=freqs,
        ref_image_mask=ref_image_mask, mask=mask,
        sigma_global_mask=sigma_global_mask, iterations=iterations,
        dilation_size=dilation_size, needs_convolution=needs_convolution,
        conv_task=conv_task, do_average_cube=do_average_cube,
        bin_size=bin_size, beam_area_px=beam_area_px, verbose=verbose,
        use_model_images=use_model_images, modellist=modellist,
        target_beam=target_beam, add_residual=add_residual,
        smooth_residual_to_target=smooth_residual_to_target,
        residual_beams=residual_beams, model_out_dir=model_out_dir)
    conv_cube = prep['conv_cube']
    masked_cube = prep['masked_cube']
    masked_cube_res = prep['masked_cube_res']
    mask = prep['mask']
    ref_image_conv = prep['ref_image_conv']
    target_beam = prep['target_beam']
    if prep['freqs'] is None:
        raise ValueError("Frequencies must be provided when using pre-loaded cubes")
    frequencies_ghz = prep['freqs'] / 1e9
    file_inputs = isinstance(imagelist[0], str)
    if ref_image_conv is None:
        ref_image_conv = (sort_list_by_beam_size(imagelist=imagelist,
                                                 return_df=False)[0][::-1][0]
                          if file_inputs else "ref_image")

    # 5. What is fitted: pixels or bins, each with its noise term
    (bin_labels, fit_idx, fit_cube, fit_cube_res,
     beam_area_px, noise_term) = _prepare_fit_spectra(
        imagelist, prep['conv_cube_res'], masked_cube, masked_cube_res, mask,
        weighting=weighting, noise_rms=noise_rms,
        do_average_cube=do_average_cube, bin_size=bin_size,
        binning=binning, binning_kwargs=binning_kwargs,
        beam_area_px=prep['beam_area_px'], verbose=verbose)
    nspec = len(fit_idx)
    n_free = 2 if fix_alpha_nt else 3
    if len(frequencies_ghz) < n_free:
        raise ValueError(f"{len(frequencies_ghz)} frequencies cannot constrain "
                         f"{n_free} free parameters")
    if verbose > 0 and len(frequencies_ghz) == n_free:
        print(f" > WARNING: {n_free} free parameters and {n_free} frequencies: "
              "the fit has no degrees of freedom.")

    # 6. Fit
    chunk_results = _fit_pixels_in_chunks(
        fit_spectra_SY_FF_cov_chunk, fit_idx, frequencies_ghz,
        fit_cube, fit_cube_res,
        (nu0, flux_sys_error_frac, fix_alpha_nt, weighting, robust_scale,
         alpha_nt_fixed, tuple(alpha_nt_bounds)),
        n_jobs=n_jobs, chunk_size=chunk_size, verbose=verbose,
        label=('t_specidx_map_SY_FF' if binning is None
               else f't_specidx_map_SY_FF [{binning} binning]'),
        unit='pixels' if binning is None else 'bins')
    fit_values = (np.concatenate(chunk_results, axis=0) if nspec
                  else np.empty((0, len(SY_FF_KEYS))))

    # 7. Parameter maps: every pixel of a bin gets the fit of the bin
    image_shape = conv_cube.shape[:2]

    def _to_map(values):
        out_map = np.full(image_shape, np.nan)
        if bin_labels is None:
            out_map[fit_idx[:, 0], fit_idx[:, 1]] = values
        else:
            in_bin = bin_labels >= 0
            out_map[in_bin] = values[bin_labels[in_bin]]
        return out_map

    output_maps = {key: _to_map(fit_values[:, k])
                   for k, key in enumerate(SY_FF_KEYS)}
    # total flux density at nu0 (A_sy + A_ff) with its covariance
    output_maps['S_nu0'] = output_maps['A_sy'] + output_maps['A_ff']
    with np.errstate(invalid='ignore'):
        output_maps['S_nu0_err'] = np.sqrt(np.maximum(
            output_maps['A_sy_err'] ** 2 + output_maps['A_ff_err'] ** 2
            + 2.0 * np.nan_to_num(output_maps['cov_sy_ff']), 0.0))

    # Fits that ended on a bound: 1 = A_sy at 0, 2 = A_ff at 0 (no free-free
    # detected), 4 = alpha_nt at one of alpha_nt_bounds (values add up). A
    # first-order error is not reliable there; with alpha_nt on a bound the
    # decomposition itself is not constrained by the data.
    fitted = np.isfinite(fit_values[:, 2])
    with np.errstate(invalid='ignore'):
        amp_tol = np.maximum(1e-4 * (np.abs(fit_values[:, 2]) + np.abs(fit_values[:, 4])), 2e-6)
        fit_flag = (1 * (fit_values[:, 2] <= amp_tol)
                    + 2 * (fit_values[:, 4] <= amp_tol)).astype(float)
        if not fix_alpha_nt:
            fit_flag += 4 * ((fit_values[:, 0] <= alpha_nt_bounds[0] + 1e-3)
                             | (fit_values[:, 0] >= alpha_nt_bounds[1] - 1e-3))
    fit_flag[~fitted] = np.nan
    output_maps['fit_flag'] = _to_map(fit_flag)

    if verbose > 0 and nspec:
        unit = 'pixels' if binning is None else 'bins'
        flags = fit_flag[fitted].astype(int)
        print(f" > t_specidx_map_SY_FF: {int(fitted.sum())}/{nspec} {unit} "
              f"fitted; median chi2_red = {np.nanmedian(fit_values[:, 9]):.2f}; "
              f"on a bound: A_ff=0 in {int(np.sum(flags & 2 > 0))}, "
              f"A_sy=0 in {int(np.sum(flags & 1 > 0))}"
              + ('' if fix_alpha_nt else
                 f", alpha_nt in {int(np.sum(flags & 4 > 0))}") + f" {unit}.")

        # ten example fits, re-fitted here (deterministic, same result) to
        # get the report
        for p in np.unique(np.linspace(0, nspec - 1, min(10, nspec)).astype(int)):
            i, j = fit_idx[p]
            alpha_p, A_sy_p, A_ff_p = fit_values[p, 0], fit_values[p, 2], fit_values[p, 4]
            if not np.isfinite(A_sy_p):
                continue
            y, yerr = spectral_fit.spectrum_yerr(fit_cube[i, j, :],
                                                 fit_cube_res[i, j, :],
                                                 flux_sys_error_frac)
            results_fit = do_fit_spec_SY_FF_map(
                frequencies_ghz, y, yerr, nu0, fix_alpha_nt,
                weighting=weighting, robust_scale=robust_scale,
                alpha_nt_fixed=alpha_nt_fixed,
                alpha_nt_bounds=tuple(alpha_nt_bounds))
            nu_plot = np.geomspace(0.8 * np.nanmin(frequencies_ghz),
                                   1.25 * np.nanmax(frequencies_ghz), 100)
            plt.figure(figsize=(4, 5))
            plt.errorbar(frequencies_ghz, y, yerr=abs(yerr), fmt='o',
                        label='Data', color='k', ecolor='gray', alpha=0.5)
            plt.plot(nu_plot, RC_function_SY_FF(nu_plot, A_sy_p, A_ff_p, alpha_p, nu0),
                     'r-.', label='Best-fit model')
            plt.plot(nu_plot, RC_function_SY_FF(nu_plot, A_sy_p, 0.0, alpha_p, nu0),
                     'b-.', label='Synchrotron')
            plt.plot(nu_plot, RC_function_SY_FF(nu_plot, 0.0, A_ff_p, 0.0, nu0),
                     'orange', linestyle='-.', label='Free-free')
            y_pos = y[np.isfinite(y) & (y > 0)]
            if y_pos.size:
                plt.ylim(y_pos.min() * 0.1, y_pos.max() * 10)
            plt.semilogy()
            plt.semilogx()
            plt.xlabel('Frequency [GHz]')
            plt.ylabel('Pixel Flux Density [mJy/Beam]' if bin_labels is None
                       else 'Bin mean Flux Density [mJy/Beam]')
            plt.legend()
            plt.title(f'Fit at pixel ({i}, {j})' if bin_labels is None
                      else f'Fit of bin {i}')
            plt.show()
            print(lmfit.fit_report(results_fit.params))

    # 8. Create frequency-dependent maps
    output_frequencies = list(output_frequencies)
    frequency_maps = create_frequency_maps(
        output_maps['A_sy'], output_maps['A_sy_err'],
        output_maps['A_ff'], output_maps['A_ff_err'],
        output_maps['alpha_nt'], output_maps['alpha_nt_err'],
        output_frequencies, nu0,
        cov_sy_ff=np.nan_to_num(output_maps['cov_sy_ff']),
        cov_sy_alpha=np.nan_to_num(output_maps['cov_sy_alpha']),
        cov_ff_alpha=np.nan_to_num(output_maps['cov_ff_alpha'])
    )

    # 9. Save maps if input was file paths
    tag = (_spec_output_tag(binning, weighting, use_model_images, add_residual,
                            smooth_residual_to_target)
           + ('-alphafix' if fix_alpha_nt else ''))
    files_written = []
    if save_fits:
        files_written = save_parameter_maps(
            output_maps, ref_image_conv, imagelist, tag=tag,
            frequency_maps=frequency_maps, bin_labels=bin_labels,
            target_beam=target_beam if use_model_images else None,
            out_dir=out_dir)
    if verbose > 0:
        _report_files_written(files_written, save_fits, file_inputs)

    # 10. Compile comprehensive results
    # noise rms per band used by 'chi2' (the full cube if it varies spatially)
    noise_per_band = None
    if weighting == 'chi2' and mask.any():
        noise_in_mask = np.asarray(noise_term[mask], dtype=float)
        noise_per_band = noise_in_mask[0]
        if not np.allclose(noise_in_mask, noise_per_band, equal_nan=True):
            noise_per_band = np.array(noise_term, dtype=float)
    results_dict = {
        'parameter_maps': output_maps,
        'frequency_maps': frequency_maps,
        'data_products': {
            'conv_cube': conv_cube,
            'masked_cube': masked_cube,
            'masked_cube_res': masked_cube_res,
            'mask': mask,
            'bin_map': bin_labels,
            'noise_rms': noise_per_band
        },
        'metadata': {
            'frequencies_ghz': frequencies_ghz,
            'output_frequencies': output_frequencies,
            'nu0': nu0,
            'weighting': weighting,
            'binning': binning,
            'beam_area_px': beam_area_px,
            'fix_alpha_nt': fix_alpha_nt,
            'alpha_nt_fixed': alpha_nt_fixed if fix_alpha_nt else None,
            'flux_sys_error_frac': flux_sys_error_frac,
            'units': 'A_sy, A_ff, S_nu0 and the frequency maps are in mJy/beam',
            'files_written': files_written,
            'processing_stats': {
                'total_pixels_fitted': int(np.count_nonzero(np.isfinite(output_maps['A_sy']))),
                'n_fitted': nspec,
                'cube_shape': conv_cube.shape,
                'averaging_applied': do_average_cube
            }
        }
    }

    return results_dict


# ---------------------------------------------------------------------------
# Resolved synchrotron + free-free maps with free-free absorption (FFA)
# ---------------------------------------------------------------------------

def emission_measure_from_turnover(nu_t_GHz, T_e=1e4):
    """
    Emission measure [pc cm^-6] of the gas that has free-free optical depth
    1 at nu_t (Mezger & Henderson 1967):
        EM_t = nu_t**2.1 (T_e/1e4 K)**1.35 / 3.28e-7  ~ 3.05e6 nu_t**2.1.
    Use the fitted turnover frequency, never the frequency of the spectral
    peak (which is ~1.5 nu_t behind a screen and would give EM ~2.5x too high).
    """
    return np.asarray(nu_t_GHz, dtype=float) ** 2.1 * (T_e / 1e4) ** 1.35 / FF_TAU_COEFF


def brightness_temperature_beam(S_mJy_beam, nu_GHz, bmaj_arcsec, bmin_arcsec, z=0.0):
    """
    Rayleigh-Jeans brightness temperature [K] of a surface brightness in
    mJy/beam, for a Gaussian beam with FWHM axes in arcsec:
        T_b = (2 ln2 / (pi k_B)) c^2 (1+z) S / (nu^2 bmaj bmin)
            = 1222.12 K (1+z) (S/mJy) (nu/GHz)**-2 (bmaj*bmin/arcsec**2)**-1
    (same constant as `Tb_source` / `brightness_temperature_props`; Lucatelli
    et al. 2024). z = 0 (default) gives the observed-frame value.
    """
    return (1222.12 * (1.0 + z) * np.asarray(S_mJy_beam, dtype=float)
            / (nu_GHz ** 2 * bmaj_arcsec * bmin_arcsec))


def emission_measure_from_thin_ff(S_ff_mJy_beam, nu_GHz, bmaj_arcsec,
                                  bmin_arcsec, T_e=1e4):
    """
    Beam-averaged emission measure [pc cm^-6] implied by an optically thin
    free-free surface brightness: T_b = T_e tau_ff, with tau_ff from
    Mezger & Henderson (1967). Returns (EM_thin, T_b_ff).
    """
    Tb = brightness_temperature_beam(S_ff_mJy_beam, nu_GHz, bmaj_arcsec, bmin_arcsec)
    em = Tb / (T_e * FF_TAU_COEFF * (T_e / 1e4) ** (-1.35) * nu_GHz ** (-2.1))
    return em, Tb


def estimate_subband_scatter(imagelist, residuallist, band_groups, freqs=None,
                             mask=None, min_snr=50.0):
    """
    Measure how much the sub-band images of each band scatter beyond
    thermal noise, to set `flux_ind_error_frac` of `t_specidx_map_SY_FF_FFA`.

    In every bright pixel a power law is fitted to the sub-bands of one band
    (in log space); the fractional rms about it, corrected for the two
    fitted parameters, is compared with the fractional thermal noise. What
    is left is the image-to-image error that a common flux-scale factor of
    the band cannot describe (sub-band imaging, primary beam, in-band
    calibration).

    Parameters
    ----------
    imagelist, residuallist : lists of FITS names, or cubes (ny, nx, nfreq)
    band_groups : array (nfreq,) of group ids (same id = same band)
    freqs : array (nfreq,) in Hz; read from the headers if None
    mask : bool (ny, nx), optional; pixels to use (default: all)
    min_snr : float
        Only pixels above min_snr in every sub-band of the band are used.

    Returns
    -------
    dict of arrays with one entry per group (NaN for groups with fewer than
    three images or without bright pixels):
        group : the group ids;
        scatter : median fractional scatter about the in-band power law;
        noise : median fractional thermal noise of the same pixels;
        excess : sqrt(scatter**2 - noise**2), the suggested
            flux_ind_error_frac of the images of that group;
        alpha_inband : median in-band spectral index;
        n_pix : pixels used;
    and 'per_image' (nfreq,): `excess` of the group of each image.
    """
    if isinstance(imagelist[0], str):
        cube = makecube(imagelist)
        cube_res = makecube(residuallist)
        if freqs is None:
            freqs = getfreqs(imagelist)
    else:
        cube, cube_res = np.asarray(imagelist), np.asarray(residuallist)
    if freqs is None:
        raise ValueError("freqs must be given when the inputs are arrays")
    freqs = np.asarray(freqs, dtype=float)
    band_groups = np.asarray(band_groups)
    sel = (np.ones(cube.shape[:2], dtype=bool) if mask is None
           else np.asarray(mask, dtype=bool))
    rms = np.array([mad_std(cube_res[:, :, k][np.isfinite(cube_res[:, :, k])])
                    for k in range(cube.shape[2])])
    groups = np.unique(band_groups)
    out = {key: np.full(groups.size, np.nan)
           for key in ('scatter', 'noise', 'excess', 'alpha_inband')}
    out['group'] = groups
    out['n_pix'] = np.zeros(groups.size, dtype=int)
    for g_idx, g in enumerate(groups):
        k = np.where(band_groups == g)[0]
        if k.size < 3:
            continue
        Y = np.asarray(cube[sel][:, k], dtype=float)
        ok = np.all(np.isfinite(Y) & (Y > min_snr * rms[k]), axis=1)
        if ok.sum() < 10:
            continue
        Y = Y[ok]
        A = np.vstack([np.log(freqs[k]), np.ones(k.size)]).T
        coef, *_ = np.linalg.lstsq(A, np.log(Y).T, rcond=None)
        resid = np.log(Y) - (A @ coef).T
        scatter = (np.sqrt(np.mean(resid ** 2, axis=1))
                   * np.sqrt(k.size / (k.size - 2.0)))
        noise = np.sqrt(np.mean((rms[k] / Y) ** 2, axis=1))
        out['scatter'][g_idx] = np.median(scatter)
        out['noise'][g_idx] = np.median(noise)
        out['excess'][g_idx] = np.sqrt(max(np.median(scatter) ** 2
                                           - np.median(noise) ** 2, 0.0))
        out['alpha_inband'][g_idx] = np.median(coef[0])
        out['n_pix'][g_idx] = int(ok.sum())
    out['per_image'] = out['excess'][np.searchsorted(groups, band_groups)]
    return out


def _propagate_errors(func, P, COV, rel_step=1e-4):
    """
    First-order error of func(P) for parameter vectors P (n, npar) with
    covariances COV (n, npar, npar): sqrt(g^T C g), g the numerical gradient
    (central differences). func must be vectorised over the rows of P.
    """
    n, npar = P.shape
    grad = np.zeros((n, npar))
    for a in range(npar):
        h = rel_step * np.maximum(np.abs(P[:, a]), 1e-3)
        Pp, Pm = P.copy(), P.copy()
        Pp[:, a] += h
        Pm[:, a] -= h
        with np.errstate(all='ignore'):
            grad[:, a] = (func(Pp) - func(Pm)) / (2.0 * h)
    grad = np.where(np.isfinite(grad), grad, 0.0)
    var = np.einsum('na,nab,nb->n', grad, COV, grad)
    with np.errstate(invalid='ignore'):
        return np.sqrt(np.maximum(var, 0.0))


def plot_SY_FF_FFA_fit(x, y, sigma, fit, nu0=10.0, geometry='thesis',
                       f_cov=1.0, flux_sys_error_frac=0.1, ax=None, title=None,
                       ylabel='Flux Density [mJy/Beam]',
                       flux_ind_error_frac=0.0):
    """
    Plot one spectrum with its FFA fit (`fit` is the output of
    `fit_SY_FF_FFA_spectrum`): data, absorbed model and its two terms, the
    intrinsic (unabsorbed) spectrum and the best unabsorbed fit.

    The error bars are sqrt(sigma**2 + f**2 |data| model), with f the
    quadrature sum of flux_sys_error_frac and flux_ind_error_frac (the
    fractional errors refer to the geometric mean of data and model, as in
    the fit).
    """
    x = np.asarray(x, dtype=float)
    if ax is None:
        _, ax = plt.subplots(figsize=(5.5, 4.5))
    fitted_fcov = isinstance(f_cov, str) or fit.get('f_cov_err', 0.0) != 0.0
    if isinstance(f_cov, str):      # 'fit': use the fitted covering factor
        f_cov = fit['f_cov']
    nu = np.geomspace(0.7 * np.nanmin(x), 1.4 * np.nanmax(x), 300)
    A_sy, A_ff, al, nut = fit['A_sy'], fit['A_ff'], fit['alpha_nt'], fit['nu_t']
    model = RC_function_SY_FF_FFA(nu, A_sy, A_ff, al, nut, nu0, geometry, f_cov)
    sy, ff = RC_function_SY_FF_FFA(nu, A_sy, A_ff, al, nut, nu0, geometry,
                                   f_cov, components=True)
    intrinsic = RC_function_SY_FF_FFA(nu, A_sy, A_ff, al, 0.0, nu0, geometry, f_cov)
    null = RC_function_SY_FF_FFA(nu, fit['A_sy_null'], fit['A_ff_null'],
                                 fit['alpha_nt_null'], 0.0, nu0, geometry, f_cov)
    model_x = RC_function_SY_FF_FFA(x, A_sy, A_ff, al, nut, nu0, geometry, f_cov)
    fsys = np.broadcast_to(np.asarray(flux_sys_error_frac, dtype=float), x.shape)
    find = np.broadcast_to(np.asarray(
        0.0 if flux_ind_error_frac is None else flux_ind_error_frac,
        dtype=float), x.shape)
    yerr = np.sqrt(np.asarray(sigma, dtype=float) ** 2
                   + (fsys ** 2 + find ** 2) * np.abs(np.asarray(y, dtype=float) * model_x))
    ax.errorbar(x, y, yerr=yerr, fmt='o', color='k', ecolor='gray', ms=4,
                alpha=0.7, label='Data', zorder=5)
    ax.plot(nu, model, 'r-', lw=1.8, label='SY+FF+FFA')
    ax.plot(nu, sy, 'b--', lw=1, label='Synchrotron (absorbed)')
    ax.plot(nu, ff, color='orange', ls='--', lw=1, label='Free-free (absorbed)')
    ax.plot(nu, intrinsic, color='r', ls=':', lw=1, label='Intrinsic (no absorber)')
    ax.plot(nu, null, color='gray', ls='-.', lw=1, label='Best SY+FF, no FFA')
    if np.isfinite(nut):
        ax.axvline(nut, color='purple', ls=':', lw=0.8)
    y_pos = np.asarray(y, dtype=float)
    y_pos = y_pos[np.isfinite(y_pos) & (y_pos > 0)]
    if y_pos.size:
        ax.set_ylim(0.2 * y_pos.min(), 3.0 * y_pos.max())
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Frequency [GHz]')
    ax.set_ylabel(ylabel)
    dof = max(fit['n_bands'] - fit['n_free'], 1)
    lo, hi = 10.0 ** fit['log_nu_t_lo'], 10.0 ** fit['log_nu_t_hi']
    text = (rf"$\nu_t$ = {nut:.2f} [{lo:.2f}, {hi:.2f}] GHz" "\n"
            rf"$\alpha_{{\rm sy}}$ = {al:.2f}, $\chi^2_\nu$ = {fit['chi2'] / dof:.2f}" "\n"
            rf"$\Delta$BIC = {fit['delta_bic']:.1f}")
    if fitted_fcov or fit.get('f_cov', 1.0) != 1.0:
        text += "\n" + rf"$\Xi$ = {fit.get('f_cov', f_cov):.2f}"
        if fit.get('f_cov_err', 0.0) > 0:
            text += rf" $\pm$ {fit['f_cov_err']:.2f}"
    ax.text(0.97, 0.97, text, transform=ax.transAxes, ha='right', va='top',
            fontsize=8)
    ax.legend(fontsize=6, loc='lower left')
    if title:
        ax.set_title(title, fontsize=9)
    return ax


def t_specidx_map_SY_FF_FFA(imagelist, residuallist,
                            ref_image_conv=None,
                            freqs=None,
                            ref_image_mask=None,
                            flux_sys_error_frac=0.1,
                            nu0=10.0,
                            mask=None,
                            sigma_global_mask=6,
                            iterations=1,
                            dilation_size=2,
                            geometry='thesis',
                            f_cov=1.0,
                            fix_alpha_nt=False,
                            alpha_nt_fixed=-0.85,
                            alpha_nt_bounds=(-2.0, -0.3),
                            alpha_nt_prior=None,
                            nu_t_bounds=None,
                            band_groups=None,
                            flux_ind_error_frac=None,
                            delta_bic_min=6.0,
                            needs_convolution=False,
                            conv_task='fft',
                            do_average_cube=False,
                            bin_size=2,
                            n_jobs=1,
                            chunk_size=None,
                            verbose=0,
                            output_frequencies=(1.4, 6.0, 33.0),
                            binning=None,
                            binning_kwargs=None,
                            beam_area_px=None,
                            noise_rms=None,
                            T_e=1e4,
                            fth_ref_freq=33.0,
                            beam_arcsec=None,
                            use_model_images=False,
                            modellist=None,
                            target_beam=None,
                            add_residual=True,
                            smooth_residual_to_target=False,
                            residual_beams=None,
                            model_out_dir=None,
                            out_dir=None,
                            save_fits=True,
                            fit_options=None,
                            n_examples=10,
                            example_panels='SNR_fit'):
    """
    Resolved SED maps of synchrotron + free-free emission with free-free
    absorption (FFA):

        S(nu) = A_sy (nu/nu0)**alpha_nt g_sy(tau) + A_ff (nu/nu0)**-0.1 g_ff(tau),
        tau = (nu/nu_t)**-2.1,

    fitted to the spectrum of every pixel, or of every spatial bin. A_sy and
    A_ff are the UNABSORBED surface brightnesses at nu0 and nu_t is the
    turnover frequency (tau = 1), not the frequency of the spectral peak.

    Geometry (`geometry`, `f_cov`; see `spectral_fit.ffa_factors`)
        'thesis' (default): synchrotron behind a foreground screen,
            g_sy = 1 - f_cov (1 - e^-tau) (= e^-tau for the homogeneous
            screen, f_cov = 1); free-free emission absorbed by the gas that
            emits it, g_ff = (1 - e^-tau)/tau.
        'screen': both terms behind the screen.  'mixed': both terms mixed
            with the absorber, (1 - e^-tau)/tau.
        f_cov is a number (fixed, default 1), or 'fit' (geometry 'thesis'
        only): the covering factor is then fitted for every pixel or bin.
        The screen term splits into a covered part, f_cov e^-tau, and an
        uncovered part, 1 - f_cov, which are two non-negative amplitudes of
        the synchrotron spectrum; no extra non-linear parameter is searched.
        With f_cov < 1 the spectrum does not cut off below nu_t but levels
        to (1 - f_cov) of the unabsorbed synchrotron, so a point well below
        nu_t measures 1 - f_cov and the bands around nu_t measure nu_t. Use
        it with enough data: it adds one free parameter (see the warning
        printed when fewer than three degrees of freedom remain).

    Method (`spectral_fit.fit_SY_FF_FFA_spectrum`)
        Separable least squares: for each trial (alpha_nt, nu_t) the
        amplitudes are solved in closed form with A >= 0; the two non-linear
        parameters are searched on a grid over their whole range and then
        refined, so the result does not depend on a starting point. The
        error model depends on the model flux, so each solution is judged
        with its own error model (-2 ln L = chi2 + ln det C); the unabsorbed
        model is fitted the same way and gives delta_bic, the evidence for
        the absorption.

    The cube preparation, the spatial binning and the noise model are those
    of `specidx_map` / `t_specidx_map_SY_FF` (noise rms per band, i.e.
    their weighting='chi2'); arguments with the same name mean the same.

    Parameters specific to this function
    ------------------------------------
    flux_sys_error_frac : float or array (nfreq,)
        Fractional flux-scale error, one value or one per image (e.g. larger
        for LOFAR). It refers to the geometric mean of the observed and the
        model flux density, so over- and under-prediction cost the same.
    band_groups : array (nfreq,) of group ids, optional
        Images with the same id share ONE flux-scale error (the sub-band
        images of one observation). None: every image has its own, which
        over-weights the in-band slopes when several sub-bands per band are
        fitted.
    flux_ind_error_frac : float or array (nfreq,), optional
        Fractional error that stays independent from image to image, added
        on top of the shared flux-scale error. Default (None): 0.03 when
        band_groups is given, 0 otherwise. Sub-band images of one band
        differ by more than thermal noise (imaging, primary beam, in-band
        calibration); without this term their scatter is judged against the
        noise alone and chi2_red becomes very large.
    fix_alpha_nt, alpha_nt_fixed, alpha_nt_bounds :
        Synchrotron index: free within the bounds (default) or fixed.
    alpha_nt_prior : (mean, sigma), optional
        Gaussian prior on a free alpha_nt, e.g. (-0.85, 0.15).
    nu_t_bounds : (min, max) in GHz, optional
        Range searched for nu_t. Default: (lowest frequency / 10,
        min(highest frequency, 20 GHz)).
    delta_bic_min : float
        Absorption counts as detected where delta_bic >= delta_bic_min
        (6: strong, 10: very strong evidence).
    T_e : float
        Electron temperature [K] for the emission measures.
    fth_ref_freq : float
        Frequency [GHz] at which the free-free brightness is compared with
        the emission measure of the absorber (see physical_maps).
    beam_arcsec : (bmaj, bmin), optional
        Beam FWHM [arcsec] for the brightness temperature. Default: the
        target beam (use_model_images) or the header of the reference image.
    out_dir, save_fits :
        Where the FITS maps go (default: one directory above the images),
        or no files at all.
    fit_options : dict, optional
        Further entries of `spectral_fit.FFA_FIT_DEFAULTS` (n_alpha, n_nu_t,
        n_iter).
    n_examples : int
        With verbose > 0, three figures of n_examples fits each are shown
        (see example_panels). 0 shows none.
    example_panels : 'SNR_fit' (default) or 'SNR_pix'
        How the example fits are chosen.
        'SNR_fit': by the quality of the fit: the best (lowest chi2_red),
            intermediate and worst (highest chi2_red) fits.
        'SNR_pix': by the signal-to-noise of the data: the highest,
            intermediate and lowest S/N spectra (S/N = median over the
            bands of flux density / noise).
        Every panel title gives both the S/N of the spectrum and, in the
        panel text, its chi2_red and delta_bic.

    Returns
    -------
    dict with
    parameter_maps
        A_sy, A_ff [mJy/beam, unabsorbed, at nu0], alpha_nt, nu_t [GHz] and
        their '_err' (first-order, from the covariance; for nu_t the
        half-width of the profile interval where that is closed);
        nu_t_lo, nu_t_hi : 68% interval of nu_t from the chi2 profile (the
            better error when nu_t is poorly constrained; lo = 0 where the
            data are compatible with no absorption, hi = inf where not
            bounded from above);
        nu_t_upper95 : 95% upper limit of nu_t;
        nu_t_detected : nu_t where delta_bic >= delta_bic_min, NaN elsewhere;
        f_cov, f_cov_err : covering factor (the fixed input with error 0,
            or the fitted value; a fitted value on its limit of 1 has no
            first-order error and is flagged 128);
        delta_bic (the absorbed model has one parameter more than the
            unabsorbed one, two with f_cov='fit'), chi2_red, chi2_red_null,
            alpha_nt_null, n_bands,
        the covariances ('cov_*', see spectral_fit.FFA_COV_KEYS) and fit_flag.
    fit_flag (values add up)
        1 A_sy = 0; 2 A_ff = 0; 4 alpha_nt on a bound; 8 nu_t at its lower
        bound (no absorption needed); 16 nu_t at its upper bound;
        32 absorption not significant (delta_bic < delta_bic_min);
        64 Xi_max < f_cov where absorption is detected (see physical_maps);
        128 fitted f_cov = 1 (nothing uncovered).
    frequency_maps['<freq>GHz']
        synchrotron, free_free, total_intrinsic, thermal_fraction: the
            UNABSORBED quantities, with '_err';
        total (+ '_err'), synchrotron_abs, free_free_abs,
            thermal_fraction_abs: as observed, after absorption;
        tau: free-free optical depth at that frequency.
    physical_maps  (the "absorber must emit" test)
        EM_t, EM_t_lo, EM_t_hi : emission measure [pc cm^-6] of the absorber
            from nu_t;
        Tb_ff : brightness temperature [K] of the fitted free-free emission
            at fth_ref_freq; it must stay well below T_e for the -0.1 slope
            to hold;
        EM_thin : beam-averaged emission measure that this free-free
            emission implies;
        Xi_max = EM_thin / EM_t : the largest covering factor of a uniform
            absorber that would not emit more free-free than is fitted.
            Xi_max < f_cov means the absorber cannot cover that fraction of
            the beam uniformly (clumpy or cooler gas, or the curvature is
            not free-free absorption).
    data_products : conv_cube, masked_cube, masked_cube_res, mask, bin_map,
        noise_rms.
    metadata : options used, frequencies, files written.
    """
    # 1. Cubes, mask, convolution / model restoring, frequency averaging
    prep = _prepare_spec_cube(
        imagelist, residuallist, ref_image_conv=ref_image_conv, freqs=freqs,
        ref_image_mask=ref_image_mask, mask=mask,
        sigma_global_mask=sigma_global_mask, iterations=iterations,
        dilation_size=dilation_size, needs_convolution=needs_convolution,
        conv_task=conv_task, do_average_cube=do_average_cube,
        bin_size=bin_size, beam_area_px=beam_area_px, verbose=verbose,
        use_model_images=use_model_images, modellist=modellist,
        target_beam=target_beam, add_residual=add_residual,
        smooth_residual_to_target=smooth_residual_to_target,
        residual_beams=residual_beams, model_out_dir=model_out_dir)
    conv_cube = prep['conv_cube']
    masked_cube = prep['masked_cube']
    masked_cube_res = prep['masked_cube_res']
    mask = prep['mask']
    ref_image_conv = prep['ref_image_conv']
    target_beam = prep['target_beam']
    if prep['freqs'] is None:
        raise ValueError("Frequencies must be provided when using pre-loaded cubes")
    frequencies_ghz = prep['freqs'] / 1e9
    nfreq = len(frequencies_ghz)
    file_inputs = isinstance(imagelist[0], str)
    if ref_image_conv is None and file_inputs:
        ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                                return_df=False)[0][::-1][0]

    # per-image inputs follow the frequency averaging
    def _per_image(values, name):
        if values is None or np.ndim(values) == 0:
            return values
        values = np.asarray(values)
        if do_average_cube and values.size == nfreq * bin_size:
            values = values[::bin_size]
        if values.size != nfreq:
            raise ValueError(f"{name} must have one entry per image "
                             f"({nfreq}), not {values.size}")
        return values
    flux_sys_error_frac = _per_image(flux_sys_error_frac, 'flux_sys_error_frac')
    band_groups = _per_image(band_groups, 'band_groups')
    flux_ind_error_frac = _per_image(flux_ind_error_frac, 'flux_ind_error_frac')

    options = dict(fit_options or {})
    options.update(geometry=geometry, f_cov=f_cov, fix_alpha_nt=fix_alpha_nt,
                   alpha_nt_fixed=alpha_nt_fixed,
                   alpha_nt_bounds=tuple(alpha_nt_bounds),
                   alpha_nt_prior=alpha_nt_prior, nu_t_bounds=nu_t_bounds,
                   band_groups=band_groups,
                   flux_ind_error_frac=flux_ind_error_frac)
    unknown = set(options) - set(FFA_FIT_DEFAULTS)
    if unknown:
        raise TypeError(f"unknown fit_options: {sorted(unknown)}")
    fit_fcov = isinstance(f_cov, str)
    if fit_fcov and (f_cov != 'fit' or geometry != 'thesis'):
        raise ValueError("f_cov must be a number, or 'fit' with "
                         "geometry='thesis'")
    n_free = (3 if fix_alpha_nt else 4) + (1 if fit_fcov else 0)
    if nfreq < n_free:
        raise ValueError(f"{nfreq} frequencies cannot constrain {n_free} "
                         "free parameters")
    if nu_t_bounds is None:
        nu_t_bounds = (frequencies_ghz.min() / 10.0,
                       min(frequencies_ghz.max(), 20.0))
    if verbose > 0:
        if nfreq <= n_free:
            print(f" > WARNING: {n_free} free parameters and {nfreq} "
                  "frequencies: the fit has no degrees of freedom.")
        if fit_fcov and nfreq - n_free < 3:
            print(f" > WARNING: f_cov='fit' with {nfreq} frequencies leaves "
                  f"{nfreq - n_free} degree(s) of freedom; fix alpha_nt, use "
                  "alpha_nt_prior, or add sub-band images.")
        n_below = int(np.sum(frequencies_ghz < 1.0))
        print(f" > t_specidx_map_SY_FF_FFA: geometry '{geometry}', f_cov = "
              f"{f_cov}; nu_t searched in [{nu_t_bounds[0]:.3g}, "
              f"{nu_t_bounds[1]:.3g}] GHz; {nfreq} frequencies "
              f"({frequencies_ghz.min():.3g}-{frequencies_ghz.max():.3g} GHz, "
              f"{n_below} below 1 GHz).")

    # 2. What is fitted: pixels or bins, each with its noise rms
    (bin_labels, fit_idx, fit_cube, fit_cube_res,
     beam_area_px, noise_term) = _prepare_fit_spectra(
        imagelist, prep['conv_cube_res'], masked_cube, masked_cube_res, mask,
        weighting='chi2', noise_rms=noise_rms,
        do_average_cube=do_average_cube, bin_size=bin_size,
        binning=binning, binning_kwargs=binning_kwargs,
        beam_area_px=prep['beam_area_px'], verbose=verbose)
    nspec = len(fit_idx)
    unit = 'pixels' if binning is None else 'bins'

    # 3. Fit
    chunk_results = _fit_pixels_in_chunks(
        fit_spectra_SY_FF_FFA_chunk, fit_idx, frequencies_ghz,
        fit_cube, fit_cube_res,
        (nu0, flux_sys_error_frac, options),
        n_jobs=n_jobs, chunk_size=chunk_size, verbose=verbose,
        label=('t_specidx_map_SY_FF_FFA' if binning is None
               else f't_specidx_map_SY_FF_FFA [{binning} binning]'),
        unit=unit)
    V = (np.concatenate(chunk_results, axis=0) if nspec
         else np.empty((0, len(FFA_KEYS))))
    col = {key: V[:, k] for k, key in enumerate(FFA_KEYS)}
    fitted = np.isfinite(col['A_sy'])
    # covering factor of each spectrum: the fixed input or the fitted value
    f_cov_arr = np.where(fitted, col['f_cov'], np.nan)

    # parameter vector and covariance of (A_sy, A_ff, alpha_nt, log10 nu_t)
    P = V[:, :4]
    COV = np.zeros((nspec, 4, 4))
    iu = np.triu_indices(4)
    COV[:, iu[0], iu[1]] = V[:, 4:14]
    COV[:, iu[1], iu[0]] = V[:, 4:14]
    no_cov = ~np.isfinite(COV).all(axis=(1, 2))
    COV_safe = np.where(no_cov[:, None, None], 0.0, COV)

    def _err(func):
        e = _propagate_errors(func, np.where(fitted[:, None], P, 0.0), COV_safe)
        e[no_cov | ~fitted] = np.nan
        return e

    image_shape = conv_cube.shape[:2]

    def _to_map(values):
        out_map = np.full(image_shape, np.nan)
        if bin_labels is None:
            out_map[fit_idx[:, 0], fit_idx[:, 1]] = values
        else:
            in_bin = bin_labels >= 0
            out_map[in_bin] = values[bin_labels[in_bin]]
        return out_map

    # 4. Parameters
    with np.errstate(all='ignore'):
        nu_t = 10.0 ** col['log_nu_t']
        sig_log = np.sqrt(COV[:, 3, 3])
        nu_t_lo = np.where(np.isneginf(col['log_nu_t_lo']), 0.0,
                           10.0 ** col['log_nu_t_lo'])
        nu_t_hi = 10.0 ** col['log_nu_t_hi']
        nu_t_up = 10.0 ** col['log_nu_t_up95']
        # error of nu_t: half-width of the profile interval where it is
        # closed on both sides, first-order otherwise
        closed = np.isfinite(nu_t_hi) & (nu_t_lo > 0)
        nu_t_err = np.where(closed, 0.5 * (nu_t_hi - nu_t_lo),
                            nu_t * np.log(10.0) * sig_log)
        dof = col['n_bands'] - n_free
        chi2_red = np.where(dof > 0, col['chi2'] / dof, np.nan)
        dof0 = col['n_bands'] - (n_free - 1)
        chi2_red_null = np.where(dof0 > 0, col['chi2_null'] / dof0, np.nan)
        detected = fitted & (col['delta_bic'] >= delta_bic_min)
    values = {
        'A_sy': col['A_sy'], 'A_sy_err': np.sqrt(COV[:, 0, 0]),
        'A_ff': col['A_ff'], 'A_ff_err': np.sqrt(COV[:, 1, 1]),
        'alpha_nt': col['alpha_nt'], 'alpha_nt_err': np.sqrt(COV[:, 2, 2]),
        'nu_t': nu_t, 'nu_t_err': nu_t_err,
        'nu_t_lo': nu_t_lo, 'nu_t_hi': nu_t_hi, 'nu_t_upper95': nu_t_up,
        'nu_t_detected': np.where(detected, nu_t, np.nan),
        'f_cov': f_cov_arr, 'f_cov_err': col['f_cov_err'],
        'delta_bic': col['delta_bic'], 'chi2_red': chi2_red,
        'chi2_red_null': chi2_red_null, 'alpha_nt_null': col['alpha_nt_null'],
        'n_bands': col['n_bands'],
    }
    for k, key in enumerate(spectral_fit.FFA_COV_KEYS):
        values[key] = V[:, 4 + k]

    # 5. Maps at the output frequencies (errors by propagation)
    def _sy(Pv, r):
        return Pv[:, 0] * r ** Pv[:, 2]

    def _ff(Pv, r):
        return Pv[:, 1] * r ** (-0.1)

    def _absorbed(Pv, nu):
        return RC_function_SY_FF_FFA(nu, Pv[:, 0], Pv[:, 1], Pv[:, 2],
                                     10.0 ** Pv[:, 3], nu0, geometry,
                                     np.nan_to_num(f_cov_arr, nan=1.0))

    output_frequencies = list(output_frequencies)
    freq_values = {}
    for freq in output_frequencies:
        r = freq / nu0
        with np.errstate(all='ignore'):
            sy, ff = _sy(P, r), _ff(P, r)
            sy_abs, ff_abs = RC_function_SY_FF_FFA(
                freq, P[:, 0], P[:, 1], P[:, 2], nu_t, nu0, geometry,
                f_cov_arr, components=True)
            # same keys as t_specidx_map_SY_FF ('33.0GHz'); frequencies that
            # one decimal cannot represent keep their digits ('0.144GHz')
            freq_key = (f'{freq:.1f}GHz' if abs(freq - round(freq, 1)) < 1e-9
                        else f'{freq:g}GHz')
            freq_values[freq_key] = {
                'synchrotron': sy,
                'synchrotron_err': _err(lambda Pv: _sy(Pv, r)),
                'free_free': ff,
                'free_free_err': _err(lambda Pv: _ff(Pv, r)),
                'total_intrinsic': sy + ff,
                'total_intrinsic_err': _err(lambda Pv: _sy(Pv, r) + _ff(Pv, r)),
                'thermal_fraction': ff / (sy + ff),
                'thermal_fraction_err': _err(
                    lambda Pv: _ff(Pv, r) / (_sy(Pv, r) + _ff(Pv, r))),
                'total': sy_abs + ff_abs,
                'total_err': _err(lambda Pv: _absorbed(Pv, freq)),
                'synchrotron_abs': sy_abs,
                'free_free_abs': ff_abs,
                'thermal_fraction_abs': ff_abs / (sy_abs + ff_abs),
                'tau': ffa_tau(freq, nu_t),
            }

    # 6. Physics of the absorber: it must also emit
    if beam_arcsec is None:
        if use_model_images and target_beam is not None:
            beam_arcsec = (target_beam[0], target_beam[1])
        elif isinstance(ref_image_conv, str):
            try:
                beam_hdr, _, _ = _fits_beam_cell(ref_image_conv)
                beam_arcsec = (beam_hdr[0], beam_hdr[1])
            except Exception:
                beam_arcsec = None
    with np.errstate(all='ignore'):
        EM_t = emission_measure_from_turnover(nu_t, T_e)
        phys = {'EM_t': EM_t,
                'EM_t_lo': emission_measure_from_turnover(nu_t_lo, T_e),
                'EM_t_hi': emission_measure_from_turnover(nu_t_hi, T_e)}
        if beam_arcsec is not None:
            S_ff_ref = _ff(P, fth_ref_freq / nu0)
            EM_thin, Tb_ff = emission_measure_from_thin_ff(
                S_ff_ref, fth_ref_freq, beam_arcsec[0], beam_arcsec[1], T_e)
            phys.update(EM_thin=EM_thin, Tb_ff=Tb_ff, Xi_max=EM_thin / EM_t)
        else:
            nan = np.full(nspec, np.nan)
            phys.update(EM_thin=nan, Tb_ff=nan, Xi_max=nan)
            if verbose > 0:
                print(" > no beam available (give beam_arcsec): EM_thin, "
                      "Tb_ff and Xi_max are not computed.")
        fit_flag = np.where(fitted, col['at_bound'], np.nan)
        fit_flag = fit_flag + 32 * (fitted & ~detected)
        fit_flag = fit_flag + 64 * (detected & (phys['Xi_max'] < f_cov_arr))
    values['fit_flag'] = fit_flag

    output_maps = {key: _to_map(v) for key, v in values.items()}
    frequency_maps = {fk: {key: _to_map(v) for key, v in d.items()}
                      for fk, d in freq_values.items()}
    physical_maps = {key: _to_map(v) for key, v in phys.items()}

    # 7. Summary and example fits
    if verbose > 0 and nspec:
        nfit = int(fitted.sum())
        flags = np.nan_to_num(fit_flag[fitted]).astype(int)
        with np.errstate(all='ignore'):
            med_nut = np.nanmedian(nu_t[detected]) if detected.any() else np.nan
        print(f" > t_specidx_map_SY_FF_FFA: {nfit}/{nspec} {unit} fitted; "
              f"median chi2_red = {np.nanmedian(chi2_red):.2f} "
              f"(no absorber: {np.nanmedian(chi2_red_null):.2f}).")
        print(f"   absorption detected (delta_bic >= {delta_bic_min:g}) in "
              f"{int(detected.sum())} {unit}, median nu_t = {med_nut:.3g} GHz; "
              f"of these, Xi_max < f_cov in {int(np.sum(flags & 64 > 0))}.")
        if fit_fcov and detected.any():
            fc_det = f_cov_arr[detected]
            print(f"   fitted covering factor where detected: median "
                  f"{np.nanmedian(fc_det):.2f} (p16 {np.nanpercentile(fc_det, 16):.2f}, "
                  f"p84 {np.nanpercentile(fc_det, 84):.2f}); f_cov = 1 in "
                  f"{int(np.sum(flags & 128 > 0))} {unit}.")
        print(f"   on a bound: A_ff=0 in {int(np.sum(flags & 2 > 0))}, "
              f"A_sy=0 in {int(np.sum(flags & 1 > 0))}, alpha_nt in "
              f"{int(np.sum(flags & 4 > 0))}, nu_t upper in "
              f"{int(np.sum(flags & 16 > 0))} {unit}.")

        # example fits in three figures: the highest, intermediate and
        # lowest signal-to-noise spectra (S/N = median over the bands of
        # flux density / noise)
        with np.errstate(all='ignore'):
            snr = np.nanmedian(
                np.asarray(fit_cube[fit_idx[:, 0], fit_idx[:, 1], :], dtype=float)
                / np.asarray(fit_cube_res[fit_idx[:, 0], fit_idx[:, 1], :], dtype=float),
                axis=1)
        if example_panels == 'SNR_pix':
            ranked = np.nonzero(fitted & np.isfinite(snr))[0]
            ranked = ranked[np.argsort(-snr[ranked])]       # highest S/N first
            panel_names = ('Highest S/N', 'Intermediate S/N', 'Lowest S/N')
        elif example_panels == 'SNR_fit':
            ranked = np.nonzero(fitted & np.isfinite(chi2_red))[0]
            ranked = ranked[np.argsort(chi2_red[ranked])]   # best fit first
            panel_names = ('Best fits (lowest chi2_red)', 'Intermediate fits',
                           'Worst fits (highest chi2_red)')
        else:
            raise ValueError("example_panels must be 'SNR_fit' or 'SNR_pix', "
                             f"not {example_panels!r}")
        n_ex = int(min(max(n_examples, 0), len(ranked)))
        groups_ex = []
        if n_ex > 0:
            mid0 = max((len(ranked) - n_ex) // 2, 0)
            used = set()
            for label, sel in ((panel_names[0], ranked[:n_ex]),
                               (panel_names[1], ranked[mid0:mid0 + n_ex]),
                               (panel_names[2], ranked[len(ranked) - n_ex:])):
                sel = [p for p in sel if p not in used]     # few spectra: no repeats
                used.update(sel)
                if sel:
                    groups_ex.append((label, sel))
        f_ind_plot = flux_ind_error_frac
        if f_ind_plot is None:
            f_ind_plot = 0.03 if band_groups is not None else 0.0
        for label, picks in groups_ex:
            ncol = min(5, len(picks))
            nrow = int(np.ceil(len(picks) / ncol))
            fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 4.0 * nrow),
                                     squeeze=False)
            for ax, p_idx in zip(axes.ravel(), picks):
                i, j = fit_idx[p_idx]
                y = np.asarray(fit_cube[i, j, :], dtype=float) * 1000
                sig = np.asarray(fit_cube_res[i, j, :], dtype=float) * 1000
                fit = fit_SY_FF_FFA_spectrum(frequencies_ghz, y, sig, nu0,
                                             flux_sys_error_frac, **options)
                where = f'pixel ({i}, {j})' if bin_labels is None else f'bin {i}'
                plot_SY_FF_FFA_fit(frequencies_ghz, y, sig, fit, nu0, geometry,
                                   fit['f_cov'], flux_sys_error_frac, ax=ax,
                                   flux_ind_error_frac=f_ind_plot,
                                   title=f'{where}, S/N = {snr[p_idx]:.0f}')
            for ax in axes.ravel()[len(picks):]:
                ax.set_visible(False)
            if example_panels == 'SNR_pix':
                span = f'S/N {snr[picks[-1]]:.0f} to {snr[picks[0]]:.0f}'
            else:
                span = (f'chi2_red {chi2_red[picks[0]]:.2f} to '
                        f'{chi2_red[picks[-1]]:.2f}')
            fig.suptitle(f'{label}: {unit} ({span})', fontsize=13)
            fig.tight_layout()
            plt.show()

    # 8. FITS
    tag = (_spec_output_tag(binning, 'chi2', use_model_images, add_residual,
                            smooth_residual_to_target)
           + f'-ffa{geometry}' + ('-fcovfit' if fit_fcov else '')
           + ('-alphafix' if fix_alpha_nt else ''))
    files_written = []
    if save_fits and file_inputs and isinstance(ref_image_conv, str):
        to_write = [
            ('alpha_nt', output_maps['alpha_nt'], ''),
            ('alpha_nt_error', output_maps['alpha_nt_err'], ''),
            ('A_sy', output_maps['A_sy'], 'mJy/beam'),
            ('A_sy_error', output_maps['A_sy_err'], 'mJy/beam'),
            ('A_ff', output_maps['A_ff'], 'mJy/beam'),
            ('A_ff_error', output_maps['A_ff_err'], 'mJy/beam'),
            ('nu_t', output_maps['nu_t'], 'GHz'),
            ('nu_t_lo', output_maps['nu_t_lo'], 'GHz'),
            ('nu_t_hi', output_maps['nu_t_hi'], 'GHz'),
            ('nu_t_upper95', output_maps['nu_t_upper95'], 'GHz'),
            ('nu_t_detected', output_maps['nu_t_detected'], 'GHz'),
            ('f_cov', output_maps['f_cov'], ''),
            ('f_cov_error', output_maps['f_cov_err'], ''),
            ('delta_bic', output_maps['delta_bic'], ''),
            ('chi2red', output_maps['chi2_red'], ''),
            ('chi2red_noffa', output_maps['chi2_red_null'], ''),
            ('fitflag', output_maps['fit_flag'], ''),
            ('EM_t', physical_maps['EM_t'], 'pc cm-6'),
            ('EM_thin', physical_maps['EM_thin'], 'pc cm-6'),
            ('Xi_max', physical_maps['Xi_max'], ''),
        ]
        for freq_key, maps in frequency_maps.items():
            to_write += [
                # intrinsic (unabsorbed) components and their sum
                (f'S_sy_{freq_key}', maps['synchrotron'], 'mJy/beam'),
                (f'S_sy_error_{freq_key}', maps['synchrotron_err'], 'mJy/beam'),
                (f'S_ff_{freq_key}', maps['free_free'], 'mJy/beam'),
                (f'S_ff_error_{freq_key}', maps['free_free_err'], 'mJy/beam'),
                (f'S_total_{freq_key}', maps['total_intrinsic'], 'mJy/beam'),
                (f'S_total_error_{freq_key}', maps['total_intrinsic_err'], 'mJy/beam'),
                # as observed, after absorption
                (f'S_sy_abs_{freq_key}', maps['synchrotron_abs'], 'mJy/beam'),
                (f'S_ff_abs_{freq_key}', maps['free_free_abs'], 'mJy/beam'),
                (f'S_total_abs_{freq_key}', maps['total'], 'mJy/beam'),
                (f'S_total_abs_error_{freq_key}', maps['total_err'], 'mJy/beam'),
                (f'fth_{freq_key}', maps['thermal_fraction'], ''),
                (f'fth_error_{freq_key}', maps['thermal_fraction_err'], ''),
                (f'fth_abs_{freq_key}', maps['thermal_fraction_abs'], ''),
                (f'tau_{freq_key}', maps['tau'], ''),
            ]
        if bin_labels is not None:
            to_write.append(('bins', bin_labels.astype(np.int32), ''))
        files_written = _save_spec_maps(
            to_write, ref_image_conv, tag=tag,
            target_beam=target_beam if use_model_images else None,
            out_dir=out_dir)
    if verbose > 0:
        _report_files_written(files_written, save_fits, file_inputs)

    # 9. Results
    noise_per_band = None
    if mask.any():
        noise_in_mask = np.asarray(noise_term[mask], dtype=float)
        noise_per_band = noise_in_mask[0]
        if not np.allclose(noise_in_mask, noise_per_band, equal_nan=True):
            noise_per_band = np.array(noise_term, dtype=float)
    return {
        'parameter_maps': output_maps,
        'frequency_maps': frequency_maps,
        'physical_maps': physical_maps,
        'data_products': {
            'conv_cube': conv_cube,
            'masked_cube': masked_cube,
            'masked_cube_res': masked_cube_res,
            'mask': mask,
            'bin_map': bin_labels,
            'noise_rms': noise_per_band,
        },
        'metadata': {
            'frequencies_ghz': frequencies_ghz,
            'output_frequencies': output_frequencies,
            'nu0': nu0,
            'geometry': geometry,
            'f_cov': f_cov,
            'fix_alpha_nt': fix_alpha_nt,
            'alpha_nt_fixed': alpha_nt_fixed if fix_alpha_nt else None,
            'alpha_nt_prior': alpha_nt_prior,
            'nu_t_bounds': tuple(nu_t_bounds),
            'band_groups': band_groups,
            'flux_sys_error_frac': flux_sys_error_frac,
            'flux_ind_error_frac': flux_ind_error_frac,
            'delta_bic_min': delta_bic_min,
            'T_e': T_e,
            'fth_ref_freq': fth_ref_freq,
            'beam_arcsec': beam_arcsec,
            'binning': binning,
            'beam_area_px': beam_area_px,
            'fit_options': options,
            'units': 'amplitudes and frequency maps in mJy/beam, nu_t in GHz, '
                     'EM in pc cm^-6',
            'files_written': files_written,
            'processing_stats': {
                'n_fitted': int(fitted.sum()),
                'n_spectra': nspec,
                'n_detected': int(detected.sum()),
                'cube_shape': conv_cube.shape,
                'averaging_applied': do_average_cube,
            },
        },
    }
