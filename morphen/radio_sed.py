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


def RC_function_SY_FF(nu, A_sy, A_ff, alpha_nt, nu0):
    # return A_ff * (nu0 ** (-0.1)) * ((nu / nu0) ** (-0.1)) + A_sy * (nu0 ** (alpha_nt)) * ((nu / nu0) ** (alpha_nt))
    return A_ff * ((nu / nu0) ** (-0.1)) + A_sy * ((nu / nu0) ** (alpha_nt))


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
        res = (y - model) / (np.log(yerr))
        return res.copy()


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
    text = (r"$\alpha_{\rm nth}"+f"= {(result.params['alpha'].value):.2f}\pm "
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
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, figsize=(6, 3),
                                   gridspec_kw={'height_ratios': [3, 1]})

    # fig, ax = plt.subplots()
    x_resample = np.linspace(np.min(x)*0.9, np.max(x)*1.1, 500)
    ax1.errorbar(x, y, yerr=yerr, 
                 fmt='o',label='Data', color='k', ecolor='gray',alpha=0.5)

    # nwalkers = int(len(y) * 2 * 20)
    
    
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
                             model_mean - sigma_shade*model_std,
                             model_mean + sigma_shade*model_std, color='lightgray',
                            alpha=0.7)
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

            median_prediction = np.median(model_predictions, axis=0)
            std_prediction = np.std(model_predictions, axis=0)

            ax1.fill_between(x_resample, median_prediction - sigma_shade*std_prediction,
                             median_prediction + sigma_shade*std_prediction,
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
    text = (r"$\alpha"+f"= {(result.params['alpha'].value):.2f}\pm "
            rf"{(result.params['alpha'].stderr):.2f}$")
    
    text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
    text_bbox = plt.text(text_x, text_y, text,
                        # ha='center', va='center',
                        fontsize=12, color='black',
                        bbox=text_bbox_props, transform=fig.transFigure)
    if do_mcmc_fit == True:
        if mcmc_version == 'general':
            text_x, text_y = 0.65, 0.77
            text = (r"$\alpha^{\rm MCMC}"+f"= {(param_dict['alpha']['best']):.2f}"
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
    plt.show()
    plt.clf()
    plt.close()
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
            
            plt.show()
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
        res = (y - model) / (y*(np.log(yerr)))  # okay
        return res.copy()

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
                   burn_in = 3000,
                   nsteps = 10000,
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
            ax1.legend(loc=(0.05, 0.05),frameon=True,prop={'size': 11},ncol=2)
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
            text = (r"$\alpha_{\rm sy}"+f"= {(result.params['alpha_nt'].value):.2f}\pm "
                    rf"{(result.params['alpha_nt'].stderr):.2f}$")
        text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
        text_bbox = plt.text(text_x, text_y, text,
                            # ha='center', va='center',
                            fontsize=12, color='black',
                            bbox=text_bbox_props, transform=fig.transFigure)
        
        text_x, text_y = 0.40, 0.84
        text =r"$f_{\rm th}$"f"$({nu_th}) = {(thermal_fraction_nu_th):.2f}\pm{(thermal_fraction_nu_th_err):.2f}$"
        text_bbox_props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5)
        text_bbox = plt.text(text_x, text_y, text,
                            # ha='center', va='center',
                            fontsize=12, color='black',
                            # bbox=text_bbox_props, 
                            transform=fig.transFigure)
        
        
        if do_mcmc_fit == True:
            if mcmc_version == 'general':
                text_x, text_y = 0.62, 0.80
                text = (r"$\alpha_{\rm sy}^{\rm MCMC}"+f"= {(param_dict['alpha_nt']['best']):.2f}"
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
    if do_mcmc_fit:
        return(mini,result,thermal_fraction,samples_emcee, param_dict)
    else:
        return(mini,result,thermal_fraction)



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
    fth0_init = 0.1
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
        text = (r"$\alpha_{\rm nt}"+f"= {(result.params['alpha_nt'].value):.2f}\pm "
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
        text = (r"$\alpha_{\rm nt}"+f"= {(result.params['alpha_nt'].value):.2f}\pm "
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
        res = (y - model) / np.log(y+yerr)
        # res = (y - model) / (np.log(abs(yerr)+1.0)) # okay 9
        # res = data - RC_function_S2(nu, A1l, alpha_nt)
        return res.copy()

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
        res = (y - model) / np.log(y+yerr)
        # res = (y - model) / (np.log(abs(yerr)+1.0)) # okay 9
        # res = data - RC_function_S2(nu, A1l, alpha_nt)
        return res.copy()

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

def linear_function(x, alpha, b,x0=1.0):
    return b*((x/x0)**alpha) #+ c*x*x


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

def do_fit_spec_map(freqs,fluxes,fluxes_err,nu0=1.0,verbose=0):

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
    
    def min_func(params):
        """
        More robust implementation for uncertainties of y (yerr).
        """
        alpha = params['alpha']
        b = params['b']
        model = linear_function(x, alpha, b, nu0)
        
        # Log-based weighting
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights
    
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
    
    
    fit_params = lmfit.Parameters()
    fit_params.add("alpha", value=-0.5, min=-3, max=5)
    fit_params.add("b", value=5.0, min=-1, max=500)
    
    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                        nan_policy='omit', reduce_fcn='neglogcauchy')
    
    result_1 = mini.minimize(method='least_squares',
                           max_nfev=15000, #f_scale = 1.0,
                        #    loss="huber", 
                           loss="cauchy",
                           tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12, 
                           verbose=verbose
                           )
    
    # mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
    #                     nan_policy='omit')

    # result_1 = mini.minimize(method='least_squares',
    #                     loss='soft_l1',  # Alternative to Cauchy
    #                     f_scale=0.1,     # Tune this parameter
    #                     tr_solver='exact',
    #                     ftol=1e-12, xtol=1e-12, gtol=1e-12)
    
    second_run_params = result_1.params
    
    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=15000, #f_scale = 1.0,
                        #    loss="huber", 
                           loss="cauchy", 
                           tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12, 
                           verbose=0
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


def do_fit_spec_SY_FF_map(freqs,fluxes,fluxes_err,nu0=None,
                          fix_alpha_nt=False,
                          verbose=0):
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
        # Log-based weighting
        log_weights = 1.0 / (1.0 + np.log1p(yerr / np.median(yerr)))
        return (y - model) * log_weights

    fit_params = lmfit.Parameters()
    fit_params.add("A_sy", value=1.0, min=1.0e-6, max=1000)
    fit_params.add("A_ff", value=0.1, min=1.0e-6, max=100)
    # fit_params.add("A1", value=0.5, min=-10, max=500)
    # fit_params.add("A2", value=0.5, min=-10, max=5000)
    if fix_alpha_nt == True:
        fit_params.add("alpha_nt", value=-0.85, min=-2.0, max=0.0, vary=False)
    else:
        fit_params.add("alpha_nt", value=-0.85, min=-3.0, max=3.0)
    # fit_params.add("alpha_nt", value=-0.85, min=-0.9, max=-0.8)

    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=15000,
                           nan_policy='omit', reduce_fcn='neglogcauchy')

    result_1 = mini.minimize(method='least_squares',
                             max_nfev=200000,  # f_scale = 1.0,
                            #  loss="huber", 
                             loss="cauchy",
                             tr_solver="exact",
                             ftol=1e-12, xtol=1e-12, gtol=1e-12,
                             verbose=verbose
                             )
    second_run_params = result_1.params

    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000,  # f_scale = 1.0,
                        #    loss="huber", 
                           loss="cauchy",
                           tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12,
                           verbose=verbose
                           )
    
    return result

def specidx_map(imagelist,residuallist,
                ref_image_conv=None,
                freqs=None,
                ref_image_mask = None,
                flux_sys_error_frac = 0.05,
                nu0=10.0,
                mask=None,sigma_global_mask=6,
                iterations=1,
                dilation_size=2,
                needs_convolution=False,conv_task='fft',
                return_only_cube=False,
                do_average_cube=False,bin_size=int(2),
                n_jobs=1,
                verbose=0):
    
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
    verbose (int): Verbosity level to use.
    
    """
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
    
    if needs_convolution:
        if ref_image_conv is None:
            ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                                return_df = False)[0][::-1][0]
        else:
            ref_image_conv = ref_image_conv
    
    
    if isinstance(ref_image_conv, str):
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
    
    if return_only_cube:
        return(None,None,conv_cube,masked_cube_res,masked_cube)
    else:
        # cube_image * mask_3d
        idx = np.column_stack(np.where(mask==True))
        if freqs is None:
            freqs = getfreqs(imagelist)
        alphaimage = np.empty_like(conv_cube[:,:,0])
        alphaimage_error = np.empty_like(conv_cube[:,:,0])
        
        alphaimage[:] = np.nan
        alphaimage_error[:] = np.nan
        
        # x = np.log10(freqs/1e9)
        x = freqs/1e9
        
        nspec = len(idx)
        tenpc = int(nspec / 10.0)
        count = 0
        pcount = 0
        
        
        def compute_pixel_spectral_index(i, j, x, masked_cube, masked_cube_res, nu0):
            y = masked_cube[i, j, :] * 1000
            yerr = np.sqrt((masked_cube_res[i, j, :])**2.0 + (flux_sys_error_frac * masked_cube[i, j, :])**2.0) * 1000
            results_fit = do_fit_spec_map(x, y, yerr, nu0)
            return (i, j, 
                    results_fit.params['alpha'].value, 
                    results_fit.params['alpha'].stderr,
                    results_fit.params['b'].value,
                    y,yerr,results_fit)
        
        
        pixel_indices = [(i, j) for i, j in idx]
        results = Parallel(n_jobs=n_jobs)(
        delayed(compute_pixel_spectral_index)(i, j, x, masked_cube, masked_cube_res, nu0) 
        for i, j in tqdm(pixel_indices, total=len(pixel_indices))
        )
        
        
        
        # for i, j in tqdm(idx):
        # # for i, j in idx:
        #     # if count == 0:
        #     #     print(str(pcount) + '%...')
            
            
            
        #     # yerr = (masked_cube_res[i,j,:]*1000)/(masked_cube[i,j,:]*1000 * np.log(10))
        #     # y = np.log10(masked_cube[i,j,:]*1000)
        #     y = masked_cube[i,j,:]*1000
        #     yerr = np.sqrt((masked_cube_res[i,j,:])**2.0 + (flux_sys_error_frac*masked_cube[i,j,:])**2.0)*1000
            
        
        #     results_fit = do_fit_spec_map(x,y,yerr,nu0)
        for (i, j, 
            alpha_value, alpha_err,
            b_value,y,yerr,results_fit) in results:
            alphaimage[i, j] = alpha_value
            alphaimage_error[i, j] = alpha_err
            
            # alphaimage[i, j] = results_fit.params['alpha'].value
            # alphaimage_error[i, j] = results_fit.params['alpha'].stderr
            count += 1
            if count == tenpc:
                count = 0
                pcount += 10
                if verbose>0:
                    model_best = linear_function(x, 
                                                alpha_value, 
                                                b_value,
                                                nu0)
                    plt.figure()

                    plt.errorbar(x, 
                                y, 
                                yerr=abs(yerr), 
                                fmt='o',
                                label='Data', color='k', ecolor='gray',
                                alpha=0.5)

                    plt.plot(x, model_best,
                            color='red', ls='-.', label='Best-fit model')
                    
                    plt.ylim(abs(np.nanmin(y)*0.1),np.nanmax(y)*5)
                    plt.semilogy()
                    plt.semilogx()
                    
                    plt.legend()
                    plt.show()
                    print(lmfit.fit_report(results_fit.params))

        
        if isinstance(imagelist[0], str):
            """
            Only save the alpha image if the input is a list of strings.
            """
            # alphaimage_name = imagelist[0].replace('-0000-image.fits','-alpha.fits')
            # alphaimage_error_name = imagelist[0].replace('-0000-image.fits','-alphaerror.fits')
            if ref_image_conv is None:
                ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                                    return_df = False)[0][::-1][0]
                
            _alphaimage_name = (ref_image_conv.replace('-image.','-alpha.').
                            replace('-image.cutout','-alpha.cutout').
                            replace('-image-pb','-alpha-pb'))
            
            _alphaimage_error_name = (ref_image_conv.replace('-image.','-alphaerror.').
                                    replace('-image.cutout','-alphaerror.cutout').
                                    replace('-image-pb','-alphaerror-pb'))
            
            alphaimage_name = os.path.dirname(os.path.dirname(_alphaimage_name))+'/'+os.path.basename(_alphaimage_name)
            alphaimage_error_name = os.path.dirname(os.path.dirname(_alphaimage_error_name))+'/'+os.path.basename(_alphaimage_error_name)
            
            
            pf.writeto(alphaimage_name,alphaimage,overwrite=True)
            pf.writeto(alphaimage_error_name,alphaimage_error,overwrite=True)
            copy_header(ref_image_conv,alphaimage_name,alphaimage_name)
            copy_header(ref_image_conv,alphaimage_error_name,alphaimage_error_name)
        
        return(alphaimage, alphaimage_error,conv_cube,masked_cube_res,masked_cube)


def specidx_map_SY_FF(imagelist, residuallist,
                      ref_image_conv=None,
                      freqs=None,
                      ref_image_mask=None,
                      flux_sys_error_frac=0.1,
                      nu0=10.0,
                      mask=None, sigma_global_mask=6,
                      iterations=1,
                      dilation_size=2,
                      sed_model='S2',
                      fix_alpha_nt=False,
                      needs_convolution=False, conv_task='fft',
                      do_average_cube=False, bin_size=int(2),
                      n_jobs=1,
                      verbose=0):
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
        _, mask = mask_dilation(ref_image_mask,
                                rms=mad_std(load_fits_data(residuallist[-1])),
                                show_figure=True,
                                PLOT=True,
                                iterations=iterations,
                                dilation_size=dilation_size,
                                sigma=sigma_global_mask)
    mask_3d = mask[:, :, np.newaxis]
    inv_mask = ~mask
    inv_mask_3d = inv_mask[:, :, np.newaxis]

    if needs_convolution == False:
        if ref_image_conv is None:
            ref_image_conv = sort_list_by_beam_size(imagelist=imagelist,
                                                    return_df=False)[0][::-1][0]
        else:
            ref_image_conv = imagelist[0]

    if isinstance(ref_image_conv, str):
        psf_zise = int(get_beam_size_px(ref_image_conv)[0])
        psf_image_size = int(load_fits_data(ref_image_conv).shape[0])
        print(f"PSF BEAM SIZE is --> {psf_zise} px")
        print(f"PSF IMAGE SIZE is --> {psf_image_size} px")
        larger_beam_image_data = load_fits_data(ref_image_conv)

        psf_name = tcreate_beam_psf(ref_image_conv,
                                    size=(psf_image_size, psf_image_size),
                                    aspect='elliptical')
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
                factor_conv_i = larger_beam_area / beam_area_i

                print(f"Factor Convolution is --> {factor_conv_i}")
                # factor_conv = 1.0

                conv_cube[:, :, i] = conv_image_i_uncor * factor_conv_i
                conv_image_i = conv_cube[:, :, i]

                conv_cube_res[:, :, i] = conv_residual_i_uncor * factor_conv_i
                conv_res_i = conv_cube_res[:, :, i]

                conv_name = imagelist[i].replace('-image.', '-image-conv.').replace('-image.cutout',
                                                                                    '-image-conv.cutout')
                pf.writeto(conv_name,
                           conv_image_i,
                           overwrite=True)
                copy_header(ref_image_conv, conv_name, conv_name)

                conv_res_name = residuallist[i].replace('-residual.', '-residual-conv.').replace(
                    '-residual.cutout', '-residual-conv.cutout')
                pf.writeto(conv_res_name,
                           conv_res_i,
                           overwrite=True)
                copy_header(ref_image_conv, conv_res_name, conv_res_name)
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
        averaged_cube_res = np.sqrt(np.sum(reshaped_residual ** 2, axis=3)) / bin_size

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
    idx = np.column_stack(np.where(mask == True))
    if freqs is None:
        freqs = getfreqs(imagelist)
    alphaimage = np.empty_like(conv_cube[:, :, 0])
    alphaimage_error = np.empty_like(conv_cube[:, :, 0])
    A_sy_map = np.empty_like(conv_cube[:, :, 0])
    A_sy_map_err = np.empty_like(conv_cube[:, :, 0])
    A_ff_map = np.empty_like(conv_cube[:, :, 0])
    A_ff_map_err = np.empty_like(conv_cube[:, :, 0])
    Snu0 = np.empty_like(conv_cube[:, :, 0])
    Snu0_err = np.empty_like(conv_cube[:, :, 0])
    f_th_33 = np.empty_like(conv_cube[:, :, 0])
    f_th_33_err = np.empty_like(conv_cube[:, :, 0])
    S_tot_33 = np.empty_like(conv_cube[:, :, 0])
    S_tot_33_err = np.empty_like(conv_cube[:, :, 0])

    sy_map_33 = np.empty_like(conv_cube[:, :, 0])
    sy_map_33_err = np.empty_like(conv_cube[:, :, 0])
    ff_map_33 = np.empty_like(conv_cube[:, :, 0])
    ff_map_33_err = np.empty_like(conv_cube[:, :, 0])

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

    #     x = np.log10(freqs)
    x = freqs.copy() / 1e9

    nspec = len(idx)
    tenpc = int(nspec / 10.0)
    count = 0
    pcount = 0

    masked_cube = masked_cube.astype(np.float32)
    masked_cube_res = masked_cube_res.astype(np.float32)

    def compute_pixel_nth_spectral_index(i, j, x, masked_cube, masked_cube_res, nu0):
        y = masked_cube[i, j, :] * 1000
        yerr = np.sqrt((1 * masked_cube_res[i, j, :]) ** 2.0 + (
                    flux_sys_error_frac * masked_cube[i, j, :]) ** 2.0) * 1000
        results_fit = do_fit_spec_SY_FF_map(x, y, yerr, nu0, fix_alpha_nt)
        return (i, j,
                results_fit.params['alpha_nt'].value,
                results_fit.params['alpha_nt'].stderr,
                results_fit.params['A_sy'].value,
                results_fit.params['A_sy'].stderr,
                results_fit.params['A_ff'].value,
                results_fit.params['A_ff'].stderr,
                y, yerr, results_fit
                )

    # pixel_indices = [(i, j) for i, j in idx]
    # results = Parallel(n_jobs=n_jobs)(
    # delayed(compute_pixel_nth_spectral_index)(i, j, x, masked_cube, masked_cube_res, nu0)
    # for i, j in tqdm(pixel_indices, total=len(pixel_indices))
    # )

    pixel_indices = [(i, j) for i, j in idx]
    with Parallel(n_jobs=n_jobs) as parallel:
        results = parallel(
            delayed(compute_pixel_nth_spectral_index)(
                i, j, x, masked_cube, masked_cube_res, nu0
            ) for i, j in tqdm(pixel_indices, total=len(pixel_indices))
        )

    if sed_model == 'S2':
        # for i, j in tqdm(idx):
        # for i, j in idx:
        for (i, j, alpha_nt_value, alpha_nt_err,
             A_sy_value, A_sy_err,
             A_ff_value, A_ff_err,
             y, yerr, results_fit) in results:
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

            sy_map_33[i, j] = A_sy_value * (33 / nu0) ** alpha_nt_value
            ff_map_33[i, j] = A_ff_value * (33 / nu0) ** (-0.1)
            sy_map_33_err[i, j] = sy_map_33[i, j] * np.sqrt(
                (A_sy_err / A_sy_value) ** 2.0 + (alpha_nt_err * np.log(33 / nu0)) ** 2.0)
            ff_map_33_err[i, j] = ff_map_33[i, j] * np.sqrt(
                (A_ff_err / A_ff_value) ** 2.0 + (0.1 * np.log(33 / nu0)) ** 2)

            Snu0[i, j] = RC_function_SY_FF(nu0,
                                           A_sy_value,
                                           A_ff_value,
                                           alpha_nt_value,
                                           nu0)
            S_tot_33[i, j] = RC_function_SY_FF(33,
                                               A_sy_value,
                                               A_ff_value,
                                               alpha_nt_value,
                                               nu0)

            S_tot_33_err = RC_function_SY_FF(33,
                                             A_sy_err,
                                             A_ff_err,
                                             alpha_nt_value,
                                             nu0)

            S_sy_33 = RC_function_SY_FF(33,
                                        A_sy_value,
                                        0.0,
                                        alpha_nt_value,
                                        nu0)
            S_sy_33_err = RC_function_SY_FF(33,
                                            A_sy_err,
                                            0.0,
                                            alpha_nt_value,
                                            nu0)

            S_ff_33 = RC_function_SY_FF(33,
                                        0.0,
                                        A_ff_value,
                                        0.0,
                                        nu0)
            S_ff_33_err = RC_function_SY_FF(33,
                                            0.0,
                                            A_ff_err,
                                            0.0,
                                            nu0)

            sed_sy_ff = {'S_sy_33': S_sy_33,
                         'S_sy_33_err': S_sy_33_err,
                         'S_ff_33': S_ff_33,
                         'S_ff_33_err': S_ff_33_err,
                         'sy_map_33': sy_map_33,
                         'sy_map_33_err': sy_map_33_err,
                         'ff_map_33': ff_map_33,
                         'ff_map_33_err': ff_map_33_err,
                         'S_tot_33': S_tot_33}

            f_th_33[i, j] = S_ff_33 / S_tot_33[i, j]

            f_th_33_err[i, j] = f_th_33[i, j] * np.sqrt(
                (S_ff_33_err / S_ff_33) ** 2.0 + (S_tot_33_err / (S_ff_33 + S_sy_33)) ** 2.0)

            try:
                Snu0_err[i, j] = RC_function_SY_FF(nu0,
                                                   A_sy_err,
                                                   A_ff_err,
                                                   alpha_nt_value,
                                                   nu0)
            except:
                Snu0_err[i, j] = np.nan

            count += 1
            if count == tenpc:
                count = 0
                pcount += 10
                if verbose > 0:
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

                    plt.ylim(abs(np.nanmin(y) * 0.1), np.nanmax(y) * 100)
                    plt.semilogy()
                    plt.semilogx()
                    plt.xlabel('Frequency [GHz]')
                    plt.ylabel('Pixel Flux Density [mJy/Beam]')
                    plt.legend()
                    plt.show()
                    print(lmfit.fit_report(results_fit.params))

    if isinstance(imagelist[0], str):
        """
        Only save the alpha image if the input is a list of strings.
        """
        # alphaimage_name = imagelist[0].replace('-0000-image.fits','-alpha.fits')
        # alphaimage_error_name = imagelist[0].replace('-0000-image.fits','-alphaerror.fits')
        _alphaimage_name = (ref_image_conv.replace('-image.', '-alpha_nt.').
                            replace('-image.cutout', '-alpha_nt.cutout').
                            replace('-image-pb', '-alpha_nt-pb'))

        _alphaimage_error_name = (ref_image_conv.replace('-image.', '-alpha_nt_error.').
                                  replace('-image.cutout', '-alpha_nt_error.cutout').
                                  replace('-image-pb', '-alpha_nt_error-pb'))

        alphaimage_name = os.path.dirname(
            os.path.dirname(_alphaimage_name)) + '/' + os.path.basename(_alphaimage_name)
        alphaimage_error_name = os.path.dirname(
            os.path.dirname(_alphaimage_error_name)) + '/' + os.path.basename(
            _alphaimage_error_name)

        pf.writeto(alphaimage_name, alphaimage, overwrite=True)
        pf.writeto(alphaimage_error_name, alphaimage_error, overwrite=True)
        copy_header(ref_image_conv, alphaimage_name, alphaimage_name)
        copy_header(ref_image_conv, alphaimage_error_name, alphaimage_error_name)

    return (alphaimage, alphaimage_error, Snu0, Snu0_err,
            A_ff_map, A_ff_map_err,
            A_sy_map, A_sy_map_err,
            f_th_33, f_th_33_err,
            conv_cube, masked_cube_res, masked_cube, sed_sy_ff)