"""
 __  __       _   _
|  \/  | __ _| |_| |__
| |\/| |/ _` | __| '_ \
| |  | | (_| | |_| | | |
|_|  |_|\__,_|\__|_| |_|

 _____                 _   _
|  ___|   _ _ __   ___| |_(_) ___  _ __  ___
| |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
|  _|| |_| | | | | (__| |_| | (_) | | | \__ \
|_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/

"""
def gaussian2D(x0, y0, a, fwhm, q, c, PA, size):
    """
    Creates a 2D gaussian model.

    Parameters
    ----------
    x0,y0 : float float
        center position in pixels
    a : float
        amplitude of the gaussian function, arbitrary units
        [0, inf]
    fwhm : float
        full width at half maximum of the gaussian, in pixels
        [0, inf]
    q : float
        axis ratio, q = b/a; e = 1 -q
        q in [0,1]
    c : float
        geometric parameter that controls how boxy the ellipse is
        c in [-2, 2]
    PA : float
        position angle in degrees of the meshgrid
        [-180, +180]
    size : tuple float
        size of the 2D image data array

    Returns
    -------
    numpt.ndarray 2D
        2D gaussian function image
    """
    # print(size)
    x, y = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
    x, y = rotation(PA, x0, y0, x, y)
    r = (abs(x) ** (c + 2.0) + ((abs(y)) / (q)) ** (c + 2.0)) ** (1.0 / (c + 2.0))
    # mask = 1./np.sqrt(2.*np.pi*sigma**2.) * np.exp(-r2/(2.*sigma**2.))
    gaussian_2D_model = a * np.exp(-4 * (np.log(2)) * (r) / (fwhm ** 2.0))
    return (gaussian_2D_model)


def mf_gaussian2D(x0, y0, sigma, M, N):
    x, y = np.meshgrid(np.arange(N) - x0, np.arange(M) - y0)
    r2 = (x) ** 2 + (y) ** 2
    mask = 1.0 / np.sqrt(2.0 * np.pi * sigma**2.0) * np.exp(-r2 / (2.0 * sigma**2.0))
    return mask


def rotation(PA, x0, y0, x, y):
    """
    Rotate an input image array. It can be used to modify
    the position angle (PA).

    Params
    ------
        x0,y0: center position
        PA: position angle of the meshgrid
        x,y: meshgrid arrays
    Returns
    -------
        tuple float
            rotated meshgrid arrays
    """
    # gal_center = (x0+0.01,y0+0.01)
    x0 = x0 + 0.25
    y0 = y0 + 0.25
    # convert to radians
    t = (PA * np.pi) / 180.0
    return ((x - x0) * np.cos(t) + (y - y0) * np.sin(t),
            -(x - x0) * np.sin(t) + (y - y0) * np.cos(t))


def Rn_R0(R0,n):
    """
    Determine the effective radius from the scale radius and the sersic index.
    """
    return(R0 * (bn_cpu(n)**n))
    

def bn_cpu(n):
    """
    bn function from Cioti .... (1997);
    Used to define the relation between Rn (half-light radii) and total
    luminosity

    Parameters:
        n: sersic index
    """
    return 2. * n - 1. / 3. + 4/(405*n) + 46. / (25515. * n**2.0)

def sersic2D(xy, x0, y0, PA, ell, n, In, Rn,cg=0.0, 
            #  Rtrunc=10000, 
            #  delta_r=1.0
             ):
    """
    Parameters
    ----------
    xy : tuple float
        meshgrid arrays
    x0,y0 : float float
        center position in pixels
    PA : float
        position angle in degrees of the meshgrid
        [-180, +180]
    ell : float
        ellipticity, e = 1 - q
        ell in [0,1]
    n : float
        sersic index
        n in [0, inf]
   Rn : float
        half-light radius
        Rn in [0, inf]
    In : float
        intensity at Rn
        In in [0, inf]
    cg : float
        geometric parameter that controls how boxy the ellipse is
        c in [-2, 2]
    Returns
    -------
    model : 2D array
        2D sersic function image
    """
    q = 1 - ell
    x, y = xy
    # x,y   = np.meshgrid(np.arange((size[1])),np.arange((size[0])))
    xx, yy = rotation(PA, x0, y0, x, y)
    # r     = (abs(xx)**(c+2.0)+((abs(yy))/(q))**(c+2.0))**(1.0/(c+2.0))
    r = np.sqrt((abs(xx) ** (cg+2.0) + ((abs(yy)) / (q)) ** (cg+2.0)))
    model = In * np.exp(-bn_cpu(n) * ((r / (Rn)) ** (1.0 / n) - 1.))
    # delta_r = Rtrunc * 0.5
    # truncation_model = 1.0 / (1.0 + np.exp((r - Rtrunc) / delta_r))
    # truncation_model = jnp.exp(-(r/Rtrunc)**delta_r)
    truncation_model = 1.0
    return (model * truncation_model)

def FlatSky_cpu(data_level, a):
    """
    Parameters
    ----------
    data_level : float
        data level, usually the std of the image.
        data_level in [0, inf]
    a : float
        flat sky level factor, to multiply the data_level.
        a in [0, inf]
    Returns
    -------
    float
        flat sky level

    """
    return (a * data_level)


# def break_tangent(r_brk,rs):
#     B = 2.65-4.98*(r_brk/(r_brk-rs))
#     return(B)
# def profile_break_tangent(r,r_brk,rs):
#     func = 0.5*(np.tanh( (2.0-break_tangent(r_brk,rs))*(r/r_brk)+break_tangent(r_brk,rs)
#                       ) +1 )
#     return(func)


def deconvolve_fft(image, psf):
    """
    Simple deconvolution in input array image. 
    
    CAUTION: This is just indented to simulate how a convolved residual map 
    would look like as deconvolved. It is not a real deconvolution.
    
    This was designed to provide a residual map to be used as input for the 
    Sersic fitting. Instead of providing the convolved residual map, it is 
    more correct to provide a deconvolved residual map.

    Parameters
    ----------
    image : 2D array
        Input image array.
    psf : 2D array
        Input psf array.
    Returns
    -------
    deconvolved_scaled : 2D array
        Deconvolved image array.
    deconvolved_norm : 2D array
        Deconvolved image array, normalised.

    
    """
    padded_shape = (image.shape[0] + psf.shape[0] - 1,
                    image.shape[1] + psf.shape[1] - 1)

    # Pad both image and psf to the new shape
    pad_shape = [(0, ts - s) for s, ts in zip(image.shape, padded_shape)]
    image_padded = np.pad(image, pad_shape, mode='constant')
    pad_shape = [(0, ts - s) for s, ts in zip(psf.shape, padded_shape)]
    psf_padded = np.pad(psf, pad_shape, mode='constant')
    
    
    image_fft = scipy.fftpack.fft2(image_padded)
    psf_fft = scipy.fftpack.fft2(psf_padded)
    deconvolved_fft_full = image_fft / psf_fft
    
    deconvolved_fft = deconvolved_fft_full[psf.shape[0] // 2:image.shape[0] + psf.shape[0] // 2,
            psf.shape[1] // 2:image.shape[1] + psf.shape[1] // 2]
    deconvolved = np.abs(scipy.fftpack.ifft2(deconvolved_fft))
    deconvolved_norm = deconvolved/np.sum(deconvolved)
    deconvolved_scaled = (image/np.mean(image)) * deconvolved_norm
    return deconvolved_scaled, deconvolved_norm

"""
 __  __       _   _
|  \/  | __ _| |_| |__
| |\/| |/ _` | __| '_ \
| |  | | (_| | |_| | | |
|_|  |_|\__,_|\__|_| |_|

 _____                 _   _
|  ___|   _ _ __   ___| |_(_) ___  _ __  ___
| |_ | | | | '_ \ / __| __| |/ _ \| '_ \/ __|
|  _|| |_| | | | | (__| |_| | (_) | | | \__ \
|_|   \__,_|_| |_|\___|\__|_|\___/|_| |_|___/


  ____ ____  _   _        _____             _     _          _ 
 / ___|  _ \| | | |      | ____|_ __   __ _| |__ | | ___  __| |
| |  _| |_) | | | |      |  _| | '_ \ / _` | '_ \| |/ _ \/ _` |
| |_| |  __/| |_| |      | |___| | | | (_| | |_) | |  __/ (_| |
 \____|_|    \___/       |_____|_| |_|\__,_|_.__/|_|\___|\__,_|

"""
try:
    @jit
    def bn(n):
        """
        bn function from Cioti .... (1997);
        Used to define the relation between Rn (half-light radii) and total
        luminosity

        Parameters:
            n: sersic index
        """
        return 2. * n - 1. / 3. + 4/(405*n) + 46. / (25515. * n**2.0)
except:
    def bn(n):
        """
        bn function from Cioti .... (1997);
        Used to define the relation between Rn (half-light radii) and total
        luminosity

        Parameters:
            n: sersic index
        """
        return 2. * n - 1. / 3. + 4/(405*n) + 46. / (25515. * n**2.0)

try:
    @jit
    def sersic2D_GPU(xy, x0=256, y0=256, PA=10, ell=0.9,
                    n=1.0, In=0.1, Rn=10.0, cg=0.0, 
                    # Rtrunc=10000, 
                    # delta_r=1.0
                    ):
        """
        Using Jax >> 10x to 100x faster.

        Parameters
        ----------
        xy : tuple float
            meshgrid arrays
        x0,y0 : float float
            center position in pixels
        PA : float
            position angle in degrees of the meshgrid
            [-180, +180]
        ell : float
            ellipticity, e = 1 - q
            ell in [0,1]
        n : float
            sersic index
            n in [0, inf]
        Rn : float
            half-light radius
            Rn in [0, inf]
        In : float
            intensity at Rn
            In in [0, inf]
        cg : float
            geometric parameter that controls how boxy the ellipse is
            c in [-2, 2]
        Returns
        -------
        model : 2D Jax array
            2D sersic function image
        """

        q = 1 - ell
        x, y = xy

        xx, yy = rotation_GPU(PA, x0, y0, x, y)
        # r     = (abs(xx)**(c+2.0)+((abs(yy))/(q))**(c+2.0))**(1.0/(c+2.0))
        r = jnp.sqrt((abs(xx) ** (cg + 2.0) + ((abs(yy)) / (q)) ** (cg + 2.0)))
        model = In * jnp.exp(-bn(n) * ((r / (Rn)) ** (1.0 / n) - 1.))
        # delta_r = Rtrunc * 0.5
        # delta_r = 1.0
        # truncation_model = 1.0 / (1.0 + jnp.exp((r - Rtrunc) / delta_r))
        # truncation_model = jnp.where(r < Rtrunc, 
        #                    1.0,
        #                    jnp.exp(-((r/Rtrunc)**delta_r - 1)))
        truncation_model = 1.0
        return (model * truncation_model)
except:
    def sersic2D_GPU(xy, x0=256, y0=256, PA=10, ell=0.9,
                    n=1.0, In=0.1, Rn=10.0, cg=0.0, 
                    # Rtrunc=10000, 
                    # delta_r=1.0
                    ):
        """
        Using Jax >> 10x to 100x faster.

        Parameters
        ----------
        xy : tuple float
            meshgrid arrays
        x0,y0 : float float
            center position in pixels
        PA : float
            position angle in degrees of the meshgrid
            [-180, +180]
        ell : float
            ellipticity, e = 1 - q
            ell in [0,1]
        n : float
            sersic index
            n in [0, inf]
        Rn : float
            half-light radius
            Rn in [0, inf]
        In : float
            intensity at Rn
            In in [0, inf]
        cg : float
            geometric parameter that controls how boxy the ellipse is
            c in [-2, 2]
        Returns
        -------
        model : 2D Jax array
            2D sersic function image
        """

        q = 1 - ell
        x, y = xy

        xx, yy = rotation_GPU(PA, x0, y0, x, y)
        # r     = (abs(xx)**(c+2.0)+((abs(yy))/(q))**(c+2.0))**(1.0/(c+2.0))
        r = jnp.sqrt((abs(xx) ** (cg + 2.0) + ((abs(yy)) / (q)) ** (cg + 2.0)))
        model = In * jnp.exp(-bn(n) * ((r / (Rn)) ** (1.0 / n) - 1.))
        # delta_r = Rtrunc * 0.5
        # truncation_model = 1.0 / (1.0 + jnp.exp( (r - Rtrunc) / delta_r) )
        # truncation_model = jnp.exp(-(r/Rtrunc)**delta_r)
        truncation_model = 1.0
        return (model * truncation_model)

def sersic2D_GPU_new(xy, params):
    """
    Using Jax >> 10x to 100x faster.

    Parameters
    ----------
    xy : tuple float
        meshgrid arrays
    x0,y0 : float float
        center position in pixels
    PA : float
        position angle in degrees of the meshgrid
        [-180, +180]
    ell : float
        ellipticity, e = 1 - q
        ell in [0,1]
    n : float
        sersic index
        n in [0, inf]
    Rn : float
        half-light radius
        Rn in [0, inf]
    In : float
        intensity at Rn
        In in [0, inf]
    cg : float
        geometric parameter that controls how boxy the ellipse is
        c in [-2, 2]
    Returns
    -------
    model : 2D Jax array
        2D sersic function image
    """
    print(params)
    # print(params.shape)
    x0, y0, PA, ell, n, In, Rn, cg = params
    q = 1 - ell
    x, y = xy

    xx, yy = rotation_GPU(PA, x0, y0, x, y)
    # r     = (abs(xx)**(c+2.0)+((abs(yy))/(q))**(c+2.0))**(1.0/(c+2.0))
    r = jnp.sqrt((abs(xx) ** (cg + 2.0) + ((abs(yy)) / (q)) ** (cg + 2.0)))
    model = In * jnp.exp(-bn(n) * ((r / (Rn)) ** (1.0 / n) - 1.))
    return (model)

try:
    @jit
    def rotation_GPU(PA, x0, y0, x, y):
        """
        Rotate an input image array. It can be used to modify
        the position angle (PA).

        Using Jax >> 10-100x faster.

        Params:
            x0,y0: center position
            PA: position angle of the meshgrid
            x,y: meshgrid arrays
        """
        # gal_center = (x0+0.01,y0+0.01)
        x0 = x0 + 0.25
        y0 = y0 + 0.25
        # convert to radians
        t = (PA * jnp.pi) / 180.0
        return ((x - x0) * jnp.cos(t) + (y - y0) * jnp.sin(t),
                -(x - x0) * jnp.sin(t) + (y - y0) * jnp.cos(t))

    @jit
    def FlatSky(background_data, a):
        """
        A simple model for the background.

        Parameters
        ----------
        background_data : 2D array
            Input background array.
        a : float
            flat sky level factor, to multiply the background_data.
        """
        return (a * background_data)

    @jit
    def _fftconvolve_jax(image, psf):
        """
        2D Image convolution using the analogue of scipy.signal.fftconvolve,
        but with Jax. This function is decorated to speed up things.
        """
        return jax.scipy.signal.fftconvolve(image, psf, mode='same')

except:
    def rotation_GPU(PA, x0, y0, x, y):
        """
        Rotate an input image array. It can be used to modify
        the position angle (PA).

        Using Jax >> 10-100x faster.

        Params:
            x0,y0: center position
            PA: position angle of the meshgrid
            x,y: meshgrid arrays
        """
        # gal_center = (x0+0.01,y0+0.01)
        x0 = x0 + 0.25
        y0 = y0 + 0.25
        # convert to radians
        t = (PA * jnp.pi) / 180.0
        return ((x - x0) * jnp.cos(t) + (y - y0) * jnp.sin(t),
                -(x - x0) * jnp.sin(t) + (y - y0) * jnp.cos(t))

    def FlatSky(background_data, a):
        """
        A simple model for the background.

        Parameters
        ----------
        background_data : 2D array
            Input background array.
        a : float
            flat sky level factor, to multiply the background_data.
        """
        return (a * background_data)

    def _fftconvolve_jax(image, psf):
        """
        2D Image convolution using the analogue of scipy.signal.fftconvolve,
        but with Jax. This function is decorated to speed up things.
        """
        return jax.scipy.signal.fftconvolve(image, psf, mode='same')


def read_imfit_params(fileParams,return_names=False):
    dlines = [ line for line in open(fileParams) if len(line.strip()) > 0 and line[0] != "#" ]
    values=[]
    temp=[]
    param_names = []
    for line in dlines:
#         print(line)
        if line.split()[0]=='FUNCTION' or line.split()[0]=='GAIN' or line.split()[0]=='READNOISE':
            pass
        else:
#             print(float(line.split()[1]))
            temp.append(float(line.split()[1]))
            param_names.append(line.split()[0])
        if line.split()[0]=='R_e' or line.split()[0]=='r_e':
    #         values['c1'] = {}
            values.append(np.asarray(temp))
            temp = []

    if dlines[-2].split()[1]=='FlatSky':
        values.append(np.asarray(float(dlines[-1].split()[1])))
    if return_names == True:
        return(values,param_names)
    else:
        return(values)



"""

                        ___
                       |_ _|_ __ ___   __ _  __ _  ___
                        | || '_ ` _ \ / _` |/ _` |/ _ \
                        | || | | | | | (_| | (_| |  __/
                       |___|_| |_| |_|\__,_|\__, |\___|
                                            |___/
        ____                                           _ _   _
       |  _ \  ___  ___ ___  _ __ ___  _ __   ___  ___(_) |_(_) ___  _ __
       | | | |/ _ \/ __/ _ \| '_ ` _ \| '_ \ / _ \/ __| | __| |/ _ \| '_ \
       | |_| |  __/ (_| (_) | | | | | | |_) | (_) \__ \ | |_| | (_) | | | |
       |____/ \___|\___\___/|_| |_| |_| .__/ \___/|___/_|\__|_|\___/|_| |_|
                                      |_|

"""

def setup_model_components(n_components=2):
    """
        Set up a single sersic component or a composition of sersic components.

        Uses the LMFIT objects to easily create model components.

        fi_ is just a prefix to distinguish the set of parameters for each component.

    """
    if n_components == 1:
        smodel2D = Model(sersic2D, prefix='f1_') + Model(FlatSky, prefix='s_')
    if n_components > 1:
        smodel2D = Model(sersic2D, prefix='f1_')
        for i in range(2, n_components + 1):
            smodel2D = smodel2D + Model(sersic2D, prefix='f' + str(i) + '_')
        smodel2D = smodel2D + Model(FlatSky, prefix='s_')
    return (smodel2D)


def construct_model_parameters(n_components, params_values_init_IMFIT=None,
                               init_constraints=None,observation_type='radio',
                               constrained=True, fix_n=False, 
                               fix_value_n=False,
                               fix_max_value_n=False,
                               fix_min_value_n=False,
                               fix_x0_y0=False,dr_fix = None,fix_geometry=None,
                               force_circular=None,
                               fix_max_value_Rn = False,
                               fix_min_value_Rn = False,
                               trunc=False, verbose=0,
                               init_params=0.25, final_params=4.0):
    """
    This function creates a single or multi-component Sersic model to be fitted
    onto an astronomical image.

    It uses the function setup_model_components to create the model components and specify/constrain
    the parameters space in which each parameter will vary during the fit.

    DEV NOTES:

        Note that this function handles parameter/model generation in four different ways:
            -- free parameters (params_values_init_IMFIT=None, init_constraints=None,
            constrained=False)
            -- constrained parameters from IMFIT (params_values_init_IMFIT=np.array of IMFIT
            parameters, init_constraints=None, constrained=True)
            -- initial parameter from a source extraction object and no constraints
            (params_values_init_IMFIT=None, init_constraints=SE.object,
            constrained=False)
            -- initial and constrained parameters from a source extraction object
            (params_values_init_IMFIT=None, init_constraints=SE.object, constrained=True)

        These are the four possible combinations of parameters and constraints that can be used.
        However, only the last one was tested extensively and is currently being used as default.
        It showed to be the most robust and reliable way to fit the model to the data.
        The other methods need some more testing and improvements.


    Note:

    Parameters
    ----------
    n_components : int, optional
        Number of components to be fitted. The default is None.
    params_values_init_IMFIT : list, optional
        List of initial parameters from a IMFIT config file to be used as initial guess for the fit.
        The default is None.
    init_constraints : dict, optional
        Dictionary containing initial constraints to be used as initial guess
        for the fit. The default is None.
    constrained : bool, optional
        If True, then the fit will be constrained. The default is True.
    fix_n : bool, optional
        If True, then the Sersic index will be fixed to 0.5. The default is False.
    fix_value_n : float, optional
        If True, then the Sersic index will be fixed to this value. The default is False.
    fix_x0_y0 : bool, optional
        If True, then the centre position will be fixed to the initial guess
        value. The default is False.
    dr_fix : float, optional
        If True, then the centre position will be fixed to the initial guess
        value. The default is False.
    fix_geometry : bool, optional
        If True, then the geometry of the components will be fixed to the
        initial guess value. The default is True.
    force_circular : bool or float or dict or list, optional
        Per-component circularity. `True` pins `ell` to 0 and freezes `PA` and
        `cg`; a float caps `ell` at that value instead. Normally set on the
        component itself (`ids_types` / `ids_to_add`) and read back off
        `init_constraints`; this argument is the by-final-index escape hatch for
        a direct call. The default None means "whatever the component declared".

    ----------------------------
    These will be removed in a future version.
    init_params : float, optional
        Initial parameter value. The default is 0.25.
    final_params : float, optional
        Final parameter value. The default is 4.0.
    """
    
    if n_components is None:
        n_components = len(params_values_init_IMFIT) - 1

    smodel2D = setup_model_components(n_components=n_components)
    # print(smodel2D)
    model_temp = Model(sersic2D)
    dr = 10.0  # +/- pixels bounds for x0,y0
    sky_init_bound = 0.99
    # sky_min_bound = -10.0
    sky_min_bound = -1e-3
    sky_max_bound = +1.0
    max_n_sersic = 15.0
    # sky_init_bound = 0.998
    # sky_min_bound = 0.99
    # sky_max_bound = 1.0


    # params_values_init_IMFIT = [] #grid of parameter values, each row is the
    # parameter values of a individual component

    if params_values_init_IMFIT is not None:
        """This takes the values from an IMFIT config file as init
        params and set number of components. This is useful to use results 
        from IMFIT, for example. 
        
        WARNING: This portion of the code was not revised and tested properly 
        since it was implemented. It will remain here for practical reasons 
        and for future improvements and experiments. 
        """
        for i in range(0, n_components):
            # x0, y0, PA, ell, n, In, Rn = params_values_init_IMFIT[i]
            x0, y0, PA, ell, n, In, Rn = params_values_init_IMFIT[i]
            if fix_x0_y0 is not False:
                fix_x0_y0_i = fix_x0_y0[i]
                dr_fix_i = dr_fix[i]
            else:
                fix_x0_y0_i = False
                dr_fix_i = False

            if fix_n is not False:
                fix_n_i = fix_n[i]
            else:
                fix_n_i = False

            ii = str(i + 1)
            if constrained == True:
                for param in model_temp.param_names:
                    # apply bounds to each parameter.
                    smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                            value=eval(param),
                                            min=init_params * eval(param),
                                            max=final_params * eval(param))

                    # still, some of them must be treated in particular.
                    if param == 'n':
                        if fix_n_i == True:
                            print('++==>> Fixing sersic index of component',i+1,' to 0.5')
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=2.0, min=0.49, max=0.51)
                        else:
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=eval(param), min=0.3,
                                max=8.0)
                    if param == 'x0':
                        if fix_x0_y0_i is not False:
                            """
                            Fix centre position by no more than dr_fix.
                            """
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=eval(param),
                                min=eval(param) - dr_fix_i,
                                max=eval(param) + dr_fix_i)
                        else:
                            if (init_constraints is not None) and (
                                    init_constraints['ncomps'] == n_components):
                                """
                                If initial constraints using Petro analysis are
                                provided, then use!
                                """
                                ddxx = 3  # the offset on x direction from Petro centre.
                                x0 = init_constraints['c' + ii + '_x0c']
                                x0_max = x0 + ddxx
                                x0_min = x0 - ddxx
                                
                                if verbose:
                                    print('Limiting ', param)
                                smodel2D.set_param_hint(
                                    'f' + str(i + 1) + '_' + param,
                                    value=x0,
                                    min=x0_min,
                                    max=x0_max)
                            else:
                                """
                                Then, consider that input File is good, then
                                give some bound
                                around those values.
                                """
                                if verbose:
                                    print('Limiting ', param)
                                smodel2D.set_param_hint(
                                    'f' + str(i + 1) + '_' + param,
                                    value=eval(param),
                                    min=eval(param) - dr,
                                    max=eval(param) + dr)
                    if param == 'y0':
                        if fix_x0_y0_i is not False:
                            """
                            Fix centre position by no more than dr_fix_i.
                            """
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=eval(param),
                                min=eval(param) - dr_fix_i,
                                max=eval(param) + dr_fix_i)
                        else:
                            if (init_constraints is not None) and (
                                    init_constraints['ncomps'] == n_components):
                                """
                                If initial constraints is using Petro analysis
                                are provided, then use!
                                """
                                ddyy = 3  # the offset on x direction from Petro centre.
                                y0 = init_constraints['c' + ii + '_y0c']
                                y0_max = y0 + ddyy
                                y0_min = y0 - ddyy
                                if verbose:
                                    print('Limiting ', param)
                                smodel2D.set_param_hint(
                                    'f' + str(i + 1) + '_' + param,
                                    value=y0,
                                    min=y0_min,
                                    max=y0_max)
                            else:
                                """
                                Then, consider that input File is good, then give
                                some bound around those values.
                                """
                                if verbose:
                                    print('Limiting ', param)
                                smodel2D.set_param_hint(
                                    'f' + str(i + 1) + '_' + param,
                                    value=eval(param),
                                    min=eval(param) - dr,
                                    max=eval(param) + dr)
                    if param == 'ell':
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=eval(param), min=0.001,
                                                max=0.8)
                    if param == 'PA':
                        if (init_constraints is not None) and (
                                init_constraints['ncomps'] == n_components):
                            _PA = init_constraints['c' + ii + '_PA']
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=_PA, min=_PA - 90,
                                max=_PA + 90)
                        else:
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=eval(param), min=0.0,
                                max=180.0)
                    if param == 'In':
                        if (init_constraints is not None) and (
                                init_constraints['ncomps'] == n_components):
                            I50 = init_constraints['c' + ii + '_I50']
                            I50_max = I50 * 10
                            I50_min = I50 * 0.1
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=I50_max, min=I50_min, max=I50_max)
                        else:
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=eval(param),
                                min=init_params * eval(param),
                                max=10 * final_params * eval(param))
            if constrained == False:
                for param in model_temp.param_names:
                    smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                            value=eval(param), min=0.000001)
                    if param == 'n':
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=0.5, min=0.3, max=8)
                    if param == 'PA':
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=45, min=-50.0, max=190)
                    if param == 'ell':
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=eval(param), min=0.001,
                                                max=0.99)
                    if param == 'In':
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=eval(param), min=0.0000001,
                                                max=10.0)
                    if param == 'Rn':
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=eval(param), min=0.5,
                                                max=300.0)
                    if param == 'x0':
                        if verbose:
                            print('Limiting ', param)
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=eval(param),
                                                min=eval(param) - dr * 5,
                                                max=eval(param) + dr * 5)
                    if param == 'y0':
                        # print('Limiting ',param)
                        smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                                                value=eval(param),
                                                min=eval(param) - dr * 5,
                                                max=eval(param) + dr * 5)

        # smodel2D.set_param_hint('s_a', value=1, min=0.99, max=1.01)
        smodel2D.set_param_hint('s_a', value=sky_init_bound, min=sky_min_bound, max=sky_max_bound)
        # smodel2D.set_param_hint('s_a', value=1, min=0.3, max=6.0)
    else:
        if init_constraints is not None:
            """
            This is the default option to use, and the more robust.
            """
            if constrained == True:
                """
                This is the default option to use, and the more robust.
                """
                ncomps_ = init_constraints['ncomps']

                # Every per-component argument is normalised once, up front, into a
                # dict keyed 1..ncomps. This accepts all the shapes these arguments
                # have ever taken -- a positional list (the historical form), a
                # dict keyed by component number (the canonical form now), or a
                # bare scalar. The scalar case is what makes `fix_geometry=True`
                # work: both fit drivers default to it, and the old
                # `fix_geometry[j]` raised `TypeError: 'bool' object is not
                # subscriptable` on it.
                _maps = {
                    'fix_n': as_component_map(
                        False if fix_n is False else fix_n, ncomps_,
                        default=False, name='fix_n'),
                    'fix_value_n': as_component_map(
                        None if fix_value_n is False else fix_value_n, ncomps_,
                        default=1.0, name='fix_value_n'),
                    'fix_max_value_n': as_component_map(
                        False if fix_max_value_n is False else fix_max_value_n,
                        ncomps_, default=False, name='fix_max_value_n'),
                    'fix_min_value_n': as_component_map(
                        False if fix_min_value_n is False else fix_min_value_n,
                        ncomps_, default=False, name='fix_min_value_n'),
                    'fix_max_value_Rn': as_component_map(
                        False if fix_max_value_Rn is False else fix_max_value_Rn,
                        ncomps_, default=False, name='fix_max_value_Rn'),
                    # Previously read as a single global scalar, so one component
                    # could not be pinned while another stayed capped-but-free.
                    'fix_min_value_Rn': as_component_map(
                        fix_min_value_Rn, ncomps_, default=False,
                        name='fix_min_value_Rn'),
                    'fix_x0_y0': as_component_map(
                        False if fix_x0_y0 is False else fix_x0_y0, ncomps_,
                        default=False, name='fix_x0_y0'),
                    'dr_fix': as_component_map(
                        False if dr_fix is None else dr_fix, ncomps_,
                        default=False, name='dr_fix'),
                    'fix_geometry': as_component_map(
                        False if fix_geometry is False else fix_geometry, ncomps_,
                        default=False, name='fix_geometry'),
                    'trunc': as_component_map(
                        False if trunc is False else trunc, ncomps_,
                        default=False, name='trunc'),
                    # `None`, not `False`, is the unset marker here: it is what
                    # lets the legacy `c{i}_force_circular` written by
                    # `prepare_fit` act as the fallback, below.
                    'force_circular': as_component_map(
                        force_circular, ncomps_, default=None,
                        name='force_circular'),
                }

                # Bound widths, and any fix/free setting the caller left out, come
                # from each component's preset. `init_constraints` carries
                # `c{i}_type` only if it was built by the current `prepare_fit`;
                # anything older gets no preset and falls back to the historical
                # hardcoded widths below, so existing results are reproducible.
                #
                # The fit drivers normally materialise the fix/free arguments
                # themselves (`build_fit_control_maps`); filling them in here too
                # means a direct `do_fit2D` call on a typed source cannot end up
                # with a component's bounds honouring its preset while its locked
                # Sersic index quietly does not.
                _bounds = {}
                _preset_fit = {}
                for _i in range(1, ncomps_ + 1):
                    _ctype = init_constraints.get(f'c{_i}_type')
                    if _ctype is None:
                        _bounds[_i] = {}
                        _preset_fit[_i] = {}
                        continue
                    # `observation_type` matters here: the In-bound factors are
                    # written per observation type, and leaving them as a nested
                    # dict would make `I50 * bounds_j['In_min_factor']` a TypeError.
                    _preset = resolve_component_preset(
                        _ctype, warn=False, component_id=_i,
                        observation_type=observation_type,
                        psf_fwhm_px=init_constraints.get('psf_fwhm_px'),
                        overrides=(init_constraints.get('component_overrides')
                                   or {}).get(_i, {}))
                    _bounds[_i] = {k: _preset[k] for k in PRESET_BOUND_KEYS
                                   if k in _preset}
                    _preset_fit[_i] = {k: _preset[k] for k in PRESET_FIT_KEYS
                                       if k in _preset}

                _given = {'fix_n': fix_n, 'fix_value_n': fix_value_n,
                          'fix_max_value_n': fix_max_value_n,
                          'fix_min_value_n': fix_min_value_n,
                          'fix_max_value_Rn': fix_max_value_Rn,
                          'fix_min_value_Rn': fix_min_value_Rn,
                          'dr_fix': dr_fix, 'fix_x0_y0': fix_x0_y0,
                          'fix_geometry': fix_geometry, 'trunc': trunc,
                          'force_circular': force_circular}
                for _name, _raw in _given.items():
                    if _raw is not None and _raw is not False:
                        continue  # caller supplied it; their value wins
                    for _i in range(1, ncomps_ + 1):
                        if _name in _preset_fit[_i]:
                            _maps[_name][_i] = _preset_fit[_i][_name]

                # Detected components occupy 1..nIDs; anything beyond that was
                # appended by `ids_to_add`. The two are bounded differently -- see
                # the `Rn` branch below.
                _n_IDs = init_constraints.get('nIDs', ncomps_)

                for j in range(ncomps_):
                    jj_id = j + 1
                    bounds_j = _bounds[jj_id]
                    _is_detected = jj_id <= _n_IDs
                    _has_preset = bool(_preset_fit[jj_id])

                    # Nesting context. `add_extra_component` records which
                    # detected region a component was appended to (`c{i}_parent`);
                    # `resolve_component_types` guarantees that ID is a detected
                    # region, so a parent is never itself an added component.
                    # Detected components carry `c{i}_parent == i` and are
                    # excluded by the identity test.
                    # `jj` itself is only bound further down, so key off
                    # `jj_id`, which is the same number.
                    _parent_j = init_constraints.get(
                        'c' + str(jj_id) + '_parent')
                    # `In` precedes `Rn` in `param_names`, so the seed radius has
                    # to be read here rather than borrowed from the `Rn` branch.
                    _R50_j = init_constraints.get('c' + str(jj_id) + '_R50')
                    _R50_parent_j = None
                    _I50_parent_j = None
                    if (NEST_ADDED_COMPONENT_BOUNDS and not _is_detected
                            and _parent_j is not None
                            and int(_parent_j) != jj_id):
                        _R50_parent_j = init_constraints.get(
                            'c' + str(int(_parent_j)) + '_R50')
                        _I50_parent_j = init_constraints.get(
                            'c' + str(int(_parent_j)) + '_I50')

                    fix_n_j = _maps['fix_n'][jj_id]
                    fix_value_n_j = _maps['fix_value_n'][jj_id]
                    if fix_value_n_j is None:
                        fix_value_n_j = 1.0  # will be skipped.
                    fix_max_value_n_j = _maps['fix_max_value_n'][jj_id]
                    fix_min_value_n_j = _maps['fix_min_value_n'][jj_id]
                    fix_max_value_Rn_j = _maps['fix_max_value_Rn'][jj_id]
                    fix_min_value_Rn_j = _maps['fix_min_value_Rn'][jj_id]

                    # `False <= 3.0` is True in Python, so the no-cap case used to
                    # fall into the delta-function intensity multipliers by
                    # accident. Test for a cap first, then for its size.
                    _has_Rn_cap = fix_max_value_Rn_j is not False \
                        and fix_max_value_Rn_j is not None
                    _is_delta_cap = _has_Rn_cap and fix_max_value_Rn_j <= 3.0
                    if _is_delta_cap or not _has_Rn_cap:
                        intensity_multiplier_max_factor = 10000
                        intensity_multiplier_min_factor = 1.0
                    else:
                        intensity_multiplier_max_factor = 100
                        intensity_multiplier_min_factor = 0.01
                    intensity_multiplier_min_factor = bounds_j.get(
                        'In_min_factor', intensity_multiplier_min_factor)
                    intensity_multiplier_max_factor = bounds_j.get(
                        'In_max_factor', intensity_multiplier_max_factor)

                    fix_x0_y0_j = _maps['fix_x0_y0'][jj_id]
                    dr_fix_j = _maps['dr_fix'][jj_id]
                    fix_geometry_j = _maps['fix_geometry'][jj_id]
                    trunc_j = _maps['trunc'][jj_id]

                    jj = str(j + 1)

                    # Circularity is a statement about the SHAPE of the component,
                    # so it has to reach the parameter bounds. Setting only the
                    # initial `q` (which is all `prepare_fit` used to do) left the
                    # preset's own `ell_max` in charge, and the fit walked straight
                    # back to an ellipse.
                    #
                    # What the component declared wins; the coarse
                    # `force_circular` / `force_circular_all` flags of
                    # `source_extraction`, which `prepare_fit` records as
                    # `c{i}_force_circular`, are only the fallback.
                    _fc_j = _maps['force_circular'][jj_id]
                    if _fc_j is None:
                        _fc_j = init_constraints.get(
                            'c' + jj + '_force_circular', False)
                    if _fc_j is None:
                        _fc_ell_max_j = None
                    elif isinstance(_fc_j, (bool, np.bool_)):
                        # A flag defers to the module-wide ceiling; a number is
                        # this one component's own ceiling.
                        _fc_ell_max_j = (float(FORCE_CIRCULAR_ELL_MAX)
                                         if _fc_j else None)
                    else:
                        _fc_ell_max_j = float(_fc_j)
                    _force_circular_j = _fc_ell_max_j is not None
                    _circular_frozen_j = (_force_circular_j
                                          and _fc_ell_max_j <= 0.0)

                    for param in model_temp.param_names:
                        #                         smodel2D.set_param_hint('f' + str(i + 1) + '_' + param,
                        #                                                 value=eval(param), min=0.000001)
                        if (param == 'n'):
                            if (fix_n_j == True):
                                print(f'++==>> Fixing sersic index of component {j+1} to {fix_value_n_j}.')
                                dn = 0.01
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=fix_value_n_j,
                                    min=fix_value_n_j-dn, max=fix_value_n_j+dn)
                            else:
                                # A floor on a free index, so a preset can say
                                # "disk-like but not fixed" (`disk-f`: n in
                                # [0.3, 1.5]). Historically hardcoded to 0.2.
                                n_min = 0.2
                                if (fix_min_value_n_j is not False
                                        and fix_min_value_n_j is not None):
                                    n_min = fix_min_value_n_j
                                if fix_max_value_n_j is not False:
                                    if verbose:
                                        print(f" ++==>> Limiting {param}_{jj} to {fix_max_value_n_j}")
                                    if _has_preset:
                                        # Start from the preset's preferred index
                                        # rather than sitting on the ceiling, the
                                        # same reasoning as the Rn cap above.
                                        n_init = float(np.clip(fix_value_n_j, n_min,
                                                               fix_max_value_n_j))
                                    else:
                                        n_init = fix_max_value_n_j * 1.00
                                    smodel2D.set_param_hint(
                                        'f' + str(j + 1) + '_' + param,
                                        value=n_init,
                                        min=n_min, max=fix_max_value_n_j)
                                else:
                                    smodel2D.set_param_hint(
                                        'f' + str(j + 1) + '_' + param,
                                        value=max(0.5, n_min),
                                        min=n_min, max=max_n_sersic)

                        if param == 'Rn':
                            R50 = init_constraints['c' + jj + '_R50']
                            if _has_Rn_cap:
                                """
                                Fix the value of Rn to a maximum value.
                                This is useful if one wants to fit a delta function.
                                That is, Rn=1.0 ~ delta functio, 1 pixel component.
                                """
                                if verbose:
                                    print(f" ++==>> Limiting {param}_{jj} to {fix_max_value_Rn_j}")
                                # The two cap regimes (delta-function vs a moderate
                                # ceiling) live in `effective_initial_Rn`, so the
                                # printed component summary starts from exactly the
                                # same numbers the fit does.
                                Rn_init, R50_min = effective_initial_Rn(
                                    R50, fix_max_value_Rn=fix_max_value_Rn_j,
                                    fix_min_value_Rn=fix_min_value_Rn_j)
                                R50_max = fix_max_value_Rn_j

                            else:
                                # R50_max = R50 * 4.0
                                # R50_max = init_constraints['c' + jj + '_Rp']
                                if ('Rn_min_factor' in bounds_j
                                        or 'Rn_max_factor' in bounds_j):
                                    """
                                    Preset-driven widths. These replace the
                                    `j == 0` special case below, which hardwired
                                    "component 1 is the compact one" regardless of
                                    what the component actually is -- with a preset
                                    the component says so itself. Untyped
                                    components keep the old index-based rule.
                                    """
                                    R50_max = R50 * bounds_j.get(
                                        'Rn_max_factor', 5.0)
                                    if _is_detected:
                                        """
                                        A *detected* component's R50 comes from the
                                        source extraction, which measures it on the
                                        convolved image -- so it is systematically
                                        oversized as a prior for the deconvolved
                                        model. Anchoring the lower bound to a
                                        fraction of it puts a floor under the fit
                                        that the component cannot get below, and in
                                        practice components were converging to
                                        exactly that floor. One pixel is the only
                                        physically meaningful lower limit here.

                                        Added components keep the R50-relative
                                        floor: their R50 is a deliberate scaling of
                                        the parent, not a measurement, so it is a
                                        legitimate thing to bound against.
                                        """
                                        R50_min = min(1.0, R50_max * 0.5)
                                    else:
                                        R50_min = max(R50 * bounds_j.get(
                                            'Rn_min_factor', 0.1), 0.5)
                                elif observation_type == 'radio':
                                    if j == 0:
                                        # R50_min = R50 * 0.01
                                        R50_min = 0.5
                                    else:
                                        R50_min = R50 * 0.1 #should be small.
                                    R50_max = R50 * 5.0
                                    # if R50_max < 1.0:
                                    #     intensity_multiplier_max_factor = 5000
                                else:
                                    if j == 0:
                                        R50_min = 1.0
                                        R50_max = R50 * 5.0
                                    else:
                                        R50_min = R50 * 0.5
                                        R50_max = R50 * 10.0
                                Rn_init = R50

                            """
                            Nest an added component against the one it was added
                            to. Seeded larger, it is the OUTER component -- it
                            exists to soak up emission the deblended core cannot
                            reach -- so it must not be allowed to collapse back
                            inside that core. Seeded smaller (`compact`), it is the
                            inner one and must not swell past its parent. Without
                            this the two windows overlap and the minimiser is free
                            to swap the components, which fits the data equally
                            well and means nothing physically.
                            """
                            if _R50_parent_j is not None:
                                if R50 > _R50_parent_j:
                                    R50_min = min(max(R50_min, _R50_parent_j),
                                                  R50_max)
                                elif R50 < _R50_parent_j:
                                    R50_max = max(min(R50_max, _R50_parent_j),
                                                  R50_min)
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=float(np.clip(Rn_init, R50_min, R50_max)),
                                min=R50_min, max=R50_max)


                        if param == 'In':
                            I50 = init_constraints['c' + jj + '_I50']
                            """
                            A high value of I50 is required because the
                            deconvolved model has a higher peak intensity
                            (and therefore the same for the I50 region) than the
                            convolved model. The PSF convolution attenuates significantly
                            the peak intensity.
                            
                            To-do:
                            1. Use a more robust approach, from the theoretical prediction of a 
                            deconvolved signal from a convolved signal with a Gaussian 
                            kernel.
                            """
                            if observation_type == 'radio':
                                I50_max = I50 * intensity_multiplier_max_factor
                                I50_min = I50 * intensity_multiplier_min_factor
                            else:
                                # A preset speaks for both observation types; only
                                # fall back to the optical hardcoded window when it
                                # is silent.
                                I50_max = I50 * bounds_j.get('In_max_factor', 100)
                                I50_min = I50 * bounds_j.get('In_min_factor', 0.01)

                            """
                            The intensity half of the nesting above. The more
                            extended component of a pair carries the LOWER central
                            surface brightness, so an added outer component is
                            capped at its parent's I50 and an added inner one
                            (`compact`) is floored there. The large `In_max_factor`
                            values exist to give a delta function headroom against
                            PSF attenuation (`point-like`: 1.0 / 1e4); an extended
                            outer component does not need it, and leaving it there
                            is precisely what let the outer component take the
                            core's intensity and the two swap places.

                            The `_I50_parent_j` vs `I50` tests keep the clamp
                            strictly outside this component's own seed, so an
                            inline override that inverts the intensity scaling can
                            never produce an empty window.
                            """
                            if _I50_parent_j is not None and _R50_j is not None:
                                if (_R50_j > _R50_parent_j
                                        and _I50_parent_j > I50):
                                    I50_max = min(I50_max, _I50_parent_j)
                                elif (_R50_j < _R50_parent_j
                                        and _I50_parent_j < I50):
                                    I50_min = max(I50_min, _I50_parent_j)
                                # A clamp must never invert the window, even if an
                                # inline `In_min_factor` above 1.0 put the floor
                                # over the seed to begin with.
                                I50_max = max(I50_max, I50_min)
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=I50, min=I50_min, max=I50_max)

                        """
                        Constraining PA and q from the pre-analysis of the image
                        (e.g. petro analysys) is not robust, since that image is
                        already convolved with the restoring beam, which can be
                        rotated. So the PA and q of a DECONVOLVED_MODEL
                        (the actual minimization problem here) can be different
                        from the PA and q of a CONVOLVED_MODEL
                        (as well Rn_conv > Rn_decon; In_conv< In_deconv).
                        So, at least we give some large bound.
                        """
                        if param == 'PA':
                            dO = bounds_j.get('PA_window', 110)
                            _PA = init_constraints['c' + jj + '_PA'] + 30
                            PA_max = _PA + dO
                            PA_min = _PA - dO
                            # A circle has no position angle, so leaving PA free
                            # only adds a fully degenerate direction to the fit
                            # (that is where the ~100% errors on PA came from).
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=_PA, min=PA_min, max=PA_max,
                                vary=not _circular_frozen_j)

                        if param == 'ell':
                            ell = 1 - init_constraints['c' + jj + '_q']
                            # ell_min = ell * 0.2
                            ell_min = 0.0
                            #                         if ell + dell <= 1.0:
                            # if ell * 1.0 <= 0.5:
                            if _has_Rn_cap:
                                if _is_delta_cap:
                                    ell_max = 0.03
                                    ell = 0.01
                                else:
                                    ell_max = 0.6
                                    ell = 0.01
                            else:
                                if ell <= 0.15:
                                    ell_max = 0.15
                                elif ell <= 0.02:
                                    ell_max = 0.02
                                else:
                                    if observation_type == 'radio':
                                        ell_max = 0.6
                                    else:
                                        ell_max = 0.6

                            # A preset states how round its component should be
                            # directly, rather than having it inferred from whether
                            # an Rn cap happens to be set.
                            if 'ell_max' in bounds_j:
                                ell_max = bounds_j['ell_max']
                                ell = min(ell, ell_max)

                            # Forced-circular wins over both the preset and the
                            # ladder above -- it is the caller's explicit request.
                            if _force_circular_j:
                                ell_min = 0.0
                                ell_max = _fc_ell_max_j
                                # Start inside the window rather than on its edge.
                                ell = 0.5 * ell_max

                            if _circular_frozen_j:
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=0.0, min=0.0, max=1.0, vary=False)
                            else:
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=ell, min=ell_min, max=ell_max)
                            # print(f" -- Init value for ell: {ell}")
                            # print(f" -- min value for ell: {ell_min}")
                            # print(f" -- max value for ell: {ell_max}")
                            # smodel2D.set_param_hint(
                            #     'f' + str(j + 1) + '_' + param,
                            #     value=ell, min=0.01, max=0.9)


                        if param == 'cg':
                            if _circular_frozen_j:
                                # `cg` bends the isophote away from an ellipse, so
                                # any non-zero value is not a circle either, however
                                # small the window `fix_geometry` would allow.
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=0.0, min=-0.01, max=0.01, vary=False)
                            elif fix_geometry_j == True:
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=0.0, min=-0.01, max=0.01)
                            else:
                                print('Using general elliptical geometry during '
                                      'fitting... may take longer.')
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=0.0, min=-2.0, max=2.0)
                                
                        # if param == 'Rtrunc':
                        #     if trunc_j == True:
                        #         smodel2D.set_param_hint(
                        #             'f' + str(j + 1) + '_' + param,
                        #             value=100.0, min=1.0, max=5000.0,
                        #             # expr = f"f{j+1}_Rn * 1"
                        #             )
                        #     else:
                        #         smodel2D.set_param_hint(
                        #             'f' + str(j + 1) + '_' + param,
                        #             value=10000.0, min=9999.99, max=10000.0001)
                        
                        # # if param == 'delta_r':
                        # #     if trunc_j == True:
                        # #         smodel2D.set_param_hint(
                        # #             'f' + str(j + 1) + '_' + param,
                        # #             value=1.0, min=0.9999, max=1.0001)
                        # #     else:
                        # #         smodel2D.set_param_hint(
                        # #             'f' + str(j + 1) + '_' + param,
                        # #             value=1.0, min=0.9999, max=1.0001)

                        if param == 'x0':
                            if fix_x0_y0_j is not False:
                                """
                                Fix centre position by no more than dr_fix.
                                """
                                x0c = init_constraints['c' + jj + '_x0c']
                                x0_max = x0c + dr_fix_j
                                x0_min = x0c - dr_fix_j
                                
                                if verbose:
                                    print(f" ++==>> Limiting {param}={x0c}+/-{dr_fix_j}")
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=x0c,
                                    min=x0_min,
                                    max=x0_max)
                            else:
                                ddxx = 10
                                x0c = init_constraints['c' + jj + '_x0c']
                                x0_max = x0c + ddxx
                                x0_min = x0c - ddxx
                                if verbose:
                                    print(f" ++==>> Limiting {param}={x0c}+/-{ddxx}")
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=x0c,
                                    min=x0_min,
                                    max=x0_max)
                        if param == 'y0':
                            if fix_x0_y0_j is not False:
                                """
                                Fix centre position by no more than dr_fix.
                                """
                                y0c = init_constraints['c' + jj + '_y0c']
                                y0_max = y0c + dr_fix_j
                                y0_min = y0c - dr_fix_j
                                if verbose:
                                    print(f" ++==>> Limiting {param}={y0c}+/-{dr_fix_j}")
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=y0c,
                                    min=y0_min,
                                    max=y0_max)
                            else:
                                ddyy = 10
                                y0c = init_constraints['c' + jj + '_y0c']
                                y0_max = y0c + ddyy
                                y0_min = y0c - ddyy
                                if verbose:
                                    print(f" ++==>> Limiting {param}={y0c}+/-{ddyy}")
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=y0c,
                                    min=y0_min,
                                    max=y0_max)
            if constrained == False:
                for j in range(init_constraints['ncomps']):
                    jj = str(j + 1)
                    for param in model_temp.param_names:
                        smodel2D.set_param_hint('f' + str(j + 1) + '_' + param,
                                                value=eval(param), min=0.000001,
                                                max=0.5)
                        if param == 'n':
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=0.5, min=0.3, max=8)
                        if param == 'PA':
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=45, min=-50.0, max=190)
                        if param == 'ell':
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=0.2, min=0.001, max=0.9)
                        if param == 'In':
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=0.1, min=0.0000001, max=10.0)
                        if param == 'Rn':
                            Rp = init_constraints['c' + jj + '_x0c']
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=10, min=2.0, max=2 * Rp)

                        """
                        This is not contrained, but at least is a good idea to
                        give some hints to the centre (x0,y0).
                        """
                        if param == 'x0':
                            ddxx = 20
                            x0c = init_constraints['c' + jj + '_x0c']
                            x0_max = x0c + ddxx
                            x0_min = x0c - ddxx
                            if verbose:
                                print('Limiting ', param)
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=x0c,
                                min=x0_min,
                                max=x0_max)
                        if param == 'y0':
                            ddyy = 20
                            y0c = init_constraints['c' + jj + '_y0c']
                            y0_max = y0c + ddyy
                            y0_min = y0c - ddyy
                            if verbose:
                                print('Limiting ', param)
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=y0c,
                                min=y0_min,
                                max=y0_max)

            # smodel2D.set_param_hint('s_a', value=1, min=0.99, max=1.01)
            smodel2D.set_param_hint('s_a', value=sky_init_bound, min=sky_min_bound, max=sky_max_bound)
            # smodel2D.set_param_hint('s_a', value=1, min=0.3, max=6.0)
        else:
            '''
            Run a complete free-optimization.
            '''
            try:
                for j in range(n_components):
                    jj = str(j + 1)
                    for param in model_temp.param_names:
                        smodel2D.set_param_hint('f' + str(j + 1) + '_' + param,
                                                value=0.5, min=0.000001)
                        if param == 'n':
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=0.5, min=0.3, max=6)
                # smodel2D.set_param_hint('s_a', value=1, min=0.99, max=1.01)
                smodel2D.set_param_hint('s_a', value=sky_init_bound, min=sky_min_bound, max=sky_max_bound)
                # smodel2D.set_param_hint('s_a', value=1, min=0.3, max=6.0)
            except:
                print('Please, if not providing initial parameters file,')
                print('provide basic information for the source.')
                return (ValueError)

    params = smodel2D.make_params()
    # print(smodel2D.param_hints)
    
    if verbose:
        # Create a PrettyTable object
        table = PrettyTable()
        table.field_names = ["Parameter", "Value", "Min", "Max"]

        # Add rows to the table with formatted values
        for key, val in smodel2D.param_hints.items():
            table.add_row([
                key, 
                f"{val['value']:.3f}", 
                f"{val['min']:.3f}", 
                f"{val['max']:.3f}"
            ])

        # Print the table
        print(table)
    
    return (smodel2D, params)



def constrain_nelder_mead_params(params,
                                 max_factor = 1.03,
                                 min_factor = 0.97):
    """
    Constrain Nelder-Mead optimised parameters.
    Since Nelder-Mead is robust, we can feed these values
    into Least-Squares.

    This is a workaround since nelder-mead does not provide statistical errors.
    SO, this is an attempt to produce a set of parameter distributions around the best fit ones.
    """
    params_copy = params.copy()
    for name, param in params_copy.items():
        value = param.value
        if param.value > 0:
            param.max = value * max_factor
        if param.value < 0:
            param.max = value * min_factor
        if param.value > 0:
            param.min = value * min_factor
        if param.value < 0:
            param.min = value * max_factor
        if param.value == 0:
            param.max = 0.01
            param.min = -0.01
    return(params_copy)


def generate_random_params_uniform_old(params, param_errors):
    # Generate a set of random numbers from a normal distribution with mean 0 and standard deviation 1
    # try:
    #     # Scale the random numbers by the standard errors of the parameters
    param_errors_corr = param_errors.copy()
    random_nums = np.random.uniform(-5, 5, size=len(params))
    scaled_random_nums = random_nums * param_errors
    random_params = params + scaled_random_nums
    return random_params

def generate_random_params_uniform(params, param_errors, sigma_errors=5.0):
    # Generate a set of random numbers from a normal distribution with mean 0 and standard deviation 1

    random_params = np.random.uniform(params-sigma_errors*param_errors, 
                                    params+sigma_errors*param_errors, 
                                    size=len(params))
    return random_params


def generate_random_params_normal(params, param_errors, sigma_errors=3.0):
    # Generate a set of random numbers from a normal distribution with mean 0 and standard deviation 1
    param_errors_corr = param_errors.copy()

    #     ndim_params = int((len(params)-1)/(len(params[0:8])))
    #     weights = np.asarray([3,3,5,0.05,0.1,0.00001,5,0.01])
    #     weights_m = np.tile(weights, (ndim_params, 1))
    #     weights_f = weights_m.flatten()
    #     weights_f = np.append(weights_f,np.asarray([0.1]))
    # #     np.random.seed(123)

    #     # Generate a random distribution of values between -1 and 1
    #     random_noise = np.random.uniform(low=-1, high=1, size=len(weights_f))
    # #     random_noise = np.random.random(len(weights_f)) * weights_f

    random_params = np.random.normal(params, 
                                    #  abs(sigma_errors*params*0.3), 
                                     sigma_errors*param_errors, 
                                     size=len(params))
    
    # scaled_random_nums = random_nums * param_errors  # + random_noise
    #     random_nums = np.random.normal(0.0, 0.1, size=len(params))
    # scaled_random_nums = random_nums * params
    # random_params = scaled_random_nums
    # random_params = params + scaled_random_nums
    return random_params


def generate_random_params_tukeylambda(params, param_errors):
    from scipy.stats import tukeylambda
    # Generate a set of random numbers from a tukeylambda distribution.
    param_errors_corr = param_errors.copy()
    random_nums = tukeylambda.rvs(0.5, size=len(params)) * 10
    scaled_random_nums = random_nums * param_errors
    #     random_nums =  tukeylambda.rvs(0.5, size=len(params)) *0.1
    #     scaled_random_nums = random_nums * params
    random_params = params + scaled_random_nums
    return random_params


# =============================================================================
# Component type presets
# =============================================================================
"""
A "component type" is a named PRESET, not a new model function.

Every component is still a `sersic2D` (`sersic2D_GPU` on the GPU path): a Gaussian
is simply a Sersic with n=0.5, so nothing new needs to be minimised. What a preset
does is bundle values for knobs that already exist, and feed them to the two
subsystems that shape a component:

  (a) `add_extra_component`        -- the initial guess (R50/I50 scale factors);
  (b) `construct_model_parameters` -- the fit bounds and which parameters are fixed.

The two are not independent: every Rn and In bound is `seed * factor`, so the scale
factors in (a) set the windows in (b) as well. A growth preset therefore states an
invariant, not just a starting point -- the component is seeded outside its parent
AND bounded to stay there (see `NEST_ADDED_COMPONENT_BOUNDS`). A preset that scales
by 1.0/1.0 does neither, and leaves the fit free to swap the two components.

Each preset classifies every parameter it touches into one of three tiers:

  locked     the preset owns it. A user override is warned about and ignored, which
             is what makes the label a guarantee: `disk` means n=1, always. To take
             control of a locked parameter, switch the component to `sersic`.
  defaulted  the preset supplies a starting point; an explicit user value wins
             silently.
  free       the preset says nothing; the user (or the global default) decides.

`sersic` is by construction the preset whose `locked` set is empty -- the escape
hatch.

Extensibility. The `model` field exists so that a future genuinely-different
profile (a real `gaussian2D`, a non-Sersic form) can be selected per component.
Adding a new *type* today is one registry entry and nothing else. Adding a new
*model function* is a larger job than it looks: the profile is spelled out in four
places -- `setup_model_components`, the CPU and GPU residual closures in
`_create_minimizer_functions`, and the per-component rebuild in
`_build_model_dict`. `_check_supported_models` below fails loudly rather than
silently fitting the wrong profile if a preset ever names something else.
"""

DEFAULT_COMPONENT_TYPE = 'sersic'

#: keys consumed by `add_extra_component` (initial guess)
PRESET_GUESS_KEYS = ('radius_scale_factor', 'intensity_scale_factor')

#: keys consumed by `construct_model_parameters` as per-component fix/bound args
PRESET_FIT_KEYS = ('fix_n', 'fix_value_n', 'fix_max_value_n', 'fix_min_value_n',
                   'fix_max_value_Rn', 'fix_min_value_Rn',
                   'dr_fix', 'fix_x0_y0', 'fix_geometry', 'trunc',
                   'force_circular')

#: keys consumed by `construct_model_parameters` to size the parameter bounds.
#: When a preset is silent about one of these, the historical hardcoded value is
#: used, so behaviour is unchanged for anything that does not carry a preset.
PRESET_BOUND_KEYS = ('Rn_min_factor', 'Rn_max_factor',
                     'In_min_factor', 'In_max_factor',
                     'ell_max', 'PA_window')

#: Ellipticity a component is still allowed when it is marked round but does not
#: state its own ceiling (`force_circular=True` rather than `force_circular=0.02`).
#: At the default 0.0 the component is genuinely circular: `ell` is pinned to 0 and
#: `PA` and `cg` are frozen with it, since a circle has no orientation and a boxy
#: isophote (`cg` != 0) is not a circle either. Set it to a small positive number to
#: cap `ell` at that value and leave the three parameters free instead.
FORCE_CIRCULAR_ELL_MAX = 0.1

"""
Circularity is a per-component fit control, declared with the component::

    ids_types  = {'1': ('point-like', {'force_circular': True})}
    ids_to_add = {'1': [('sersic', {'force_circular': 0.05})]}

`True` means exactly round (see `FORCE_CIRCULAR_ELL_MAX`); a number is this one
component's ellipticity ceiling, so a compact core can be pinned perfectly round
while a disk is merely kept from going too flat. `False` -- the default of every
preset, since no preset declares it -- leaves `ell`, `PA` and `cg` alone.

`source_extraction`'s `force_circular` / `force_circular_all` are the older,
coarser way of saying the same thing (all detected components, or all components).
They still work and are still honoured, but only where the component itself said
nothing, and they are due for removal: prefer declaring it on the component.
"""

#: Whether an added component's radius and intensity windows are nested against
#: the component it was added to. A component appended with `ids_to_add` is a
#: statement about WHERE it sits: a growth preset (`sersic`, `gaussian`, `disk`,
#: ...) seeds it larger and fainter than its parent, so it is the outer component
#: and must not be allowed to collapse back inside the core; `compact` seeds it
#: smaller and brighter, so it is the inner one and must not swell past the
#: parent. Seeds alone only bias the minimiser -- the windows still overlap, and
#: an ill-conditioned source can converge with the two components swapped, which
#: fits the data just as well and means nothing physically.
#:
#: Set `mlibs.NEST_ADDED_COMPONENT_BOUNDS = False` to restore the unnested
#: windows (both components bounded only by their own seed) for A/B comparison.
NEST_ADDED_COMPONENT_BOUNDS = True

COMPONENT_TYPE_PRESETS = {
    'point-like': {
        'model': 'sersic2D',
        'doc': 'Unresolved point source: Rn pinned to ~1 pixel, n locked to 0.5.',
        'locked': {
            'fix_n': True, 'fix_value_n': 0.5,
            # 1 px is scale-free, hence locked: it *is* the definition.
            'fix_max_value_Rn': 1.0, 'fix_min_value_Rn': True,
        },
        'defaulted': {
            'fix_max_value_n': False,
            'fix_min_value_n': False, 'fix_x0_y0': True, 'fix_geometry': True,
            'trunc': False, 'dr_fix': 5,
            'radius_scale_factor': 1.0, 'intensity_scale_factor': 1.0,
            # A delta function loses the most peak intensity to the PSF, so it
            # needs the widest headroom on In of any preset.
            'In_min_factor': {'radio': 1.0, 'other': 0.01},
            'In_max_factor': {'radio': 1e4, 'other': 100},
            'ell_max': 0.1, 'PA_window': 110,
        },
    },
    'gaussian': {
        'model': 'sersic2D',
        'doc': 'Resolved Gaussian blob: n locked to 0.5, radius free. The default '
               'added-component type for radio fitting.',
        'locked': {'fix_n': True, 'fix_value_n': 0.5},
        'defaulted': {
            'fix_max_value_n': False,
            'fix_min_value_n': False, 'fix_max_value_Rn': False,
            'fix_min_value_Rn': False, 'fix_x0_y0': True, 'fix_geometry': True,
            'trunc': False, 'dr_fix': 10,
            'radius_scale_factor': 1.5, 'intensity_scale_factor': 0.5,
            'Rn_min_factor': 0.3, 'Rn_max_factor': 3.0,
            'In_min_factor': {'radio': 0.01, 'other': 0.01},
            'In_max_factor': {'radio': 1e3, 'other': 100},
            'ell_max': 0.75, 'PA_window': 110,
        },
    },
    'compact': {
        'model': 'sersic2D',
        'doc': 'Compact/bulge component: the Sersic index is free (that is what '
               'separates it from `gaussian`), but the radius is capped at a few '
               'PSF widths in absolute pixels so it cannot grow into the extended '
               'emission. Effectively `sersic` plus a small Rn ceiling.',
        'locked': {},
        'defaulted': {
            'fix_n': False, 'fix_value_n': 1.0, 'fix_max_value_n': False,
            'fix_min_value_n': False,
            # The ceiling is "a few PSF widths": `Rn_max_psf_factor` is the real
            # statement, and `fix_max_value_Rn` is only the fallback used when the
            # PSF/beam size cannot be determined. Both are `defaulted` rather than
            # `locked` -- how compact "compact" should be is data-dependent, so it
            # must stay retunable without abandoning the preset.
            'Rn_max_psf_factor': 1.0,
            'fix_max_value_Rn': 10.0, 'fix_min_value_Rn': False,
            'fix_x0_y0': True, 'fix_geometry': True, 'trunc': False, 'dr_fix': 5,
            'radius_scale_factor': 0.5, 'intensity_scale_factor': 2.0,
            'In_min_factor': {'radio': 0.05, 'other': 0.01},
            'In_max_factor': {'radio': 1e2, 'other': 100},
            'ell_max': 0.50, 'PA_window': 110,
        },
    },
    'disk': {
        'model': 'sersic2D',
        'doc': 'Extended exponential disk: n locked to 1.0, larger radius, lower '
               'surface brightness.',
        'locked': {'fix_n': True, 'fix_value_n': 1.0},
        'defaulted': {
            'fix_max_value_n': False,
            'fix_min_value_n': False, 'fix_max_value_Rn': False,
            'fix_min_value_Rn': False, 'fix_x0_y0': True, 'fix_geometry': True,
            'trunc': False, 'dr_fix': 20,
            'radius_scale_factor': 2.5, 'intensity_scale_factor': 0.3,
            'Rn_min_factor': 0.5, 'Rn_max_factor': 3.0,
            'In_min_factor': {'radio': 0.01, 'other': 0.01},
            'In_max_factor': {'radio': 1e2, 'other': 100},
            'ell_max': 0.90, 'PA_window': 110,
        },
    },
    'disk-f': {
        'model': 'sersic2D',
        'doc': 'Flexible disk: same size and brightness scaling as `disk`, but the '
               'Sersic index is fitted within a disk-like range (0.3-1.5) instead '
               'of being locked to exactly 1.0. Use it when the outer component is '
               'disk-ish but not a clean exponential -- a thick disk, a disk with a '
               'bar or truncation, a lenticular. Nothing is locked, so the range '
               'itself can be retuned with fix_min_value_n / fix_max_value_n.',
        'locked': {},
        'defaulted': {
            # n free, but confined to the disk-like range. `fix_value_n` is the
            # starting index here rather than a fixed one, since `fix_n` is False.
            'fix_n': False, 'fix_value_n': 1.0,
            'fix_min_value_n': 0.3, 'fix_max_value_n': 1.5,
            'fix_max_value_Rn': False, 'fix_min_value_Rn': False,
            'fix_x0_y0': True, 'fix_geometry': True, 'trunc': False, 'dr_fix': 20,
            'radius_scale_factor': 2.5, 'intensity_scale_factor': 0.3,
            'Rn_min_factor': 0.5, 'Rn_max_factor': 3.0,
            'In_min_factor': {'radio': 0.01, 'other': 0.01},
            'In_max_factor': {'radio': 1e2, 'other': 100},
            'ell_max': 0.90, 'PA_window': 110,
        },
    },
    'sersic': {
        'model': 'sersic2D',
        'doc': 'General Sersic component: nothing is locked. The user-controlled '
               'escape hatch, and the default added-component type for optical '
               'fitting.',
        'locked': {},
        'defaulted': {
            'fix_n': False, 'fix_value_n': 1.0, 'fix_max_value_n': False,
            'fix_min_value_n': False,
            'fix_max_value_Rn': False, 'fix_min_value_Rn': False,
            'fix_x0_y0': True, 'fix_geometry': True, 'trunc': False, 'dr_fix': 10,
            # Added *outside* its parent, like every other growth preset --
            # `gaussian` 1.5/0.5, `disk`/`disk-f` 2.5/0.3, `envelope` 4.0/0.10.
            # These were 1.0/1.0, which copied the parent verbatim and, because
            # the Rn/In windows are `seed * factor`, gave the two components
            # identical bounds as well; the fit was then free to swap them.
            'radius_scale_factor': 2.0, 'intensity_scale_factor': 0.4,
            'Rn_min_factor': 0.2, 'Rn_max_factor': 5.0,
            'In_min_factor': {'radio': 0.01, 'other': 0.01},
            'In_max_factor': {'radio': 1e3, 'other': 100},
            'ell_max': 0.75, 'PA_window': 110,
        },
    },
    # Legacy types. `add_extra_component` has always accepted these and they only
    # ever affected the initial guess, so they keep their historical scale factors
    # and lock nothing.
    'envelope': {
        'model': 'sersic2D',
        'doc': 'Legacy: very extended, faint envelope (initial guess only).',
        'locked': {},
        'defaulted': {
            'fix_n': False, 'fix_value_n': 1.0, 'fix_max_value_n': False,
            'fix_min_value_n': False,
            'fix_max_value_Rn': False, 'fix_min_value_Rn': False,
            'fix_x0_y0': True, 'fix_geometry': True, 'trunc': False, 'dr_fix': 20,
            'radius_scale_factor': 4.0, 'intensity_scale_factor': 0.10,
        },
    },
    'halo': {
        'model': 'sersic2D',
        'doc': 'Legacy: the most extended, faintest component (initial guess only).',
        'locked': {},
        'defaulted': {
            'fix_n': False, 'fix_value_n': 1.0, 'fix_max_value_n': False,
            'fix_min_value_n': False,
            'fix_max_value_Rn': False, 'fix_min_value_Rn': False,
            'fix_x0_y0': True, 'fix_geometry': True, 'trunc': False, 'dr_fix': 20,
            'radius_scale_factor': 7.0, 'intensity_scale_factor': 0.05,
        },
    },
}

#: Model functions a component may name. Extending this is not enough on its own --
#: see `_check_supported_models`.
SUPPORTED_COMPONENT_MODELS = ('sersic2D',)


def normalise_component_type(component_type):
    """
    Canonicalise a component-type name and check it against the registry.

    Case and the hyphen/underscore split are both forgiving, so `'Point_Like'`,
    `'point-like'` and `'pointlike'` all resolve to the same preset.

    Raises
    ------
    ValueError
        If the name is not in `COMPONENT_TYPE_PRESETS`. The message lists the
        valid names, since a typo here is otherwise only visible as a component
        that quietly behaves like a disk.
    """
    if component_type is None:
        return None
    if not isinstance(component_type, str):
        raise ValueError(f"Component type must be a string, got "
                         f"{type(component_type).__name__}: {component_type!r}.")

    key = component_type.strip().lower().replace('_', '-')
    if key in COMPONENT_TYPE_PRESETS:
        return key
    # tolerate 'pointlike' for 'point-like'
    squashed = {k.replace('-', ''): k for k in COMPONENT_TYPE_PRESETS}
    if key.replace('-', '') in squashed:
        return squashed[key.replace('-', '')]

    raise ValueError(
        f"Unknown component type {component_type!r}. "
        f"Valid types are: {', '.join(sorted(COMPONENT_TYPE_PRESETS))}.")


def _resolve_observation_value(value, observation_type):
    """
    Flatten a preset value that is specified per observation type.

    A preset value may be written either as a plain scalar, meaning "the same for
    every observation type", or as a dict keyed by observation type::

        'In_max_factor': {'radio': 1e3, 'other': 100}

    This exists because radio and optical genuinely need different windows on `In`.
    A radio model is fitted against a *beam*-convolved image, and deconvolution
    raises the model peak enormously, so the fit needs orders of magnitude of
    headroom above the measured I50. Optical PSF attenuation is far milder, and the
    optical branch of `construct_model_parameters` has always used a much tighter
    window -- applying the radio numbers there would raise the *floor* tenfold and
    stop a faint component (an outer disk, say) from reaching its true brightness.

    With `observation_type=None` the nested form is left alone; callers that do not
    know the observation type only read fix/free keys, which are never nested.
    """
    if observation_type is None or not isinstance(value, dict):
        return value
    if observation_type in value:
        return value[observation_type]
    # 'other' is the catch-all in morphen's obs_type vocabulary ('radio' / 'other').
    return value.get('other', value.get('radio'))


def resolve_component_preset(component_type, overrides=None, component_id=None,
                             warn=True, observation_type=None, psf_fwhm_px=None):
    """
    Merge a preset with user overrides into one flat dict of component settings.

    Precedence, per the three tiers documented above: start from the preset's
    `defaulted` values, let an explicit user value overwrite them, then reapply the
    preset's `locked` values on top -- so a locked parameter always wins, and the
    attempt to override it is reported rather than silently dropped.

    Parameters
    ----------
    component_type : str
        A key of `COMPONENT_TYPE_PRESETS` (case/underscore insensitive).
    overrides : dict, optional
        User-supplied values for this component. `None` values are treated as
        "not supplied" so callers can pass a sparse dict without special-casing.
    component_id : int, optional
        1-indexed component number, used only to make the warning message useful.
    warn : bool, optional
        Emit a `UserWarning` when an override targets a locked parameter. Set
        False when resolving speculatively.
    observation_type : str, optional
        'radio' or 'other'. When given, preset values written per observation type
        are flattened to the matching scalar (see `_resolve_observation_value`).
        Callers that omit it get the nested form back untouched, which is safe as
        long as they only read fix/free keys -- those are never nested.
    psf_fwhm_px : float, optional
        PSF/beam FWHM in pixels. When a preset carries `Rn_max_psf_factor`, this
        turns it into a concrete `fix_max_value_Rn`; without it the preset's
        absolute pixel default stands.

    Returns
    -------
    dict
        Flat settings dict, always including `'model'` and `'component_type'`.
    """
    key = normalise_component_type(component_type)
    preset = COMPONENT_TYPE_PRESETS[key]
    locked = preset.get('locked', {})

    resolved = dict(preset.get('defaulted', {}))
    resolved['model'] = preset.get('model', 'sersic2D')
    resolved['component_type'] = key

    where = f" for component {component_id}" if component_id is not None else ""
    supplied = set()
    for name, value in (overrides or {}).items():
        if value is None:
            continue
        if name in locked:
            if warn and value != locked[name]:
                warnings.warn(
                    f"'{key}' preset{where} locks {name}={locked[name]!r}; "
                    f"the requested {name}={value!r} is ignored. Use the "
                    f"'sersic' component type to control this parameter.",
                    UserWarning, stacklevel=2)
            continue
        resolved[name] = value
        supplied.add(name)

    # A radius ceiling stated in PSF widths beats the preset's absolute pixel
    # fallback, but still loses to an explicit user value -- it replaces a
    # `defaulted` entry, so it sits at the same tier.
    if (psf_fwhm_px and 'Rn_max_psf_factor' in resolved
            and 'fix_max_value_Rn' not in supplied
            and 'fix_max_value_Rn' not in locked):
        resolved['fix_max_value_Rn'] = float(psf_fwhm_px
                                             * resolved['Rn_max_psf_factor'])

    resolved.update(locked)

    if observation_type is not None:
        resolved = {name: _resolve_observation_value(value, observation_type)
                    for name, value in resolved.items()}
    return resolved


def effective_initial_Rn(R50, fix_max_value_Rn=False, fix_min_value_Rn=False):
    """
    The Rn a component will actually start the fit from, once its preset's radius
    cap is applied.

    `c{i}_R50` is the source-extraction *measurement*; a preset can override it
    outright. `point-like` caps Rn at 1 px, so a component whose measured R50 is
    8 px still starts at 0.99 px -- reporting the measurement alone is misleading.
    The rules here are the ones `construct_model_parameters` applies; this exists
    so the fit and the printed component summary cannot disagree.

    Parameters
    ----------
    R50 : float
        Measured (or seeded) half-light radius of the component, in pixels.
    fix_max_value_Rn : float or False
        The component's radius ceiling in pixels, or False for no cap.
    fix_min_value_Rn : bool
        Only meaningful for a delta-function cap: pin the lower bound just below
        the ceiling instead of leaving the component free down to 0.5 px.

    Returns
    -------
    (Rn_init, Rn_min) : tuple of float
        `Rn_min` is None when there is no cap, since the lower bound is then set
        by the Rn_min_factor / detected-vs-added rules instead.
    """
    _has_cap = fix_max_value_Rn is not False and fix_max_value_Rn is not None
    if not _has_cap:
        return float(R50), None
    cap = float(fix_max_value_Rn)
    if cap <= 3.0:
        """
        A cap this small IS the initial guess: the component is being asked to be
        unresolved.
        """
        Rn_min = cap * 0.98 if fix_min_value_Rn else 0.5
        return cap * 0.99, Rn_min
    """
    A moderate cap (`compact`'s ~10 px bulge ceiling) is a ceiling, not a target.
    Starting at 0.99*cap would put every such component at the cap regardless of
    its measured size.
    """
    return float(np.clip(R50, 0.5, cap)), 0.5


def _check_supported_models(component_types):
    """
    Guard the model-dispatch seam.

    Every preset currently names `sersic2D`, so this never fires today. It exists
    so that adding a preset with a new `model` fails immediately and says what
    still needs wiring, instead of being silently fitted as a Sersic.
    """
    unsupported = sorted({
        COMPONENT_TYPE_PRESETS[normalise_component_type(t)].get('model', 'sersic2D')
        for t in component_types
        if COMPONENT_TYPE_PRESETS[normalise_component_type(t)].get(
            'model', 'sersic2D') not in SUPPORTED_COMPONENT_MODELS})
    if unsupported:
        raise NotImplementedError(
            f"Component model(s) {unsupported} are declared in the preset "
            f"registry but not implemented. Wiring a new profile in requires all "
            f"four sites that spell out sersic2D: setup_model_components, the CPU "
            f"and GPU residual closures in _create_minimizer_functions, and "
            f"_build_model_dict.")


def as_component_map(value, ncomps, default=None, name='argument'):
    """
    Normalise a per-component argument to a dict keyed by 1-indexed component.

    Accepts every shape these arguments have historically taken:

    - `None`      -> `default` for every component;
    - a scalar    -> broadcast to every component (this is what makes
                     `fix_geometry=True`, the drivers' own default, work at all --
                     the old code did `fix_geometry[j]` on it and raised
                     `TypeError: 'bool' object is not subscriptable`);
    - a sequence  -> positional, `value[i]` for component `i+1`; a short sequence
                     is padded with `default` and a long one is truncated, matching
                     the previous tolerance for mis-sized lists;
    - a dict      -> keyed by component number (`1` or `'1'`), sparse entries
                     filled with `default`.

    Returns
    -------
    dict
        `{1: ..., 2: ..., ..., ncomps: ...}`.
    """
    if value is None:
        return {i: default for i in range(1, ncomps + 1)}

    if isinstance(value, dict):
        out = {}
        for raw_key, item in value.items():
            try:
                idx = int(raw_key)
            except (TypeError, ValueError):
                raise ValueError(
                    f"{name}: dict keys must be component numbers, got "
                    f"{raw_key!r}.") from None
            if not 1 <= idx <= ncomps:
                raise ValueError(
                    f"{name}: component {idx} is out of range; this fit has "
                    f"{ncomps} component(s), so valid keys are 1..{ncomps}.")
            out[idx] = item
        return {i: out.get(i, default) for i in range(1, ncomps + 1)}

    if isinstance(value, (list, tuple, np.ndarray)):
        seq = list(value)
        return {i: (seq[i - 1] if i - 1 < len(seq) else default)
                for i in range(1, ncomps + 1)}

    # scalar (bool, number, str)
    return {i: value for i in range(1, ncomps + 1)}


#: override names a component may carry inline, next to its type
INLINE_OVERRIDE_KEYS = (PRESET_FIT_KEYS + PRESET_GUESS_KEYS + PRESET_BOUND_KEYS
                        + ('Rn_max_psf_factor',))


def parse_component_spec(spec, source='component'):
    """
    Split a component specification into its type and any inline overrides.

    A component can be named by its type alone, or by its type plus overrides for
    that one component::

        'disk'                                  # type only
        ('disk-f', {'fix_max_value_n': 3.0})    # type + overrides
        {'type': 'disk-f', 'fix_max_value_n': 3.0}

    The point of the second and third forms is locality: without them a component
    is declared by *parent ID* in `ids_to_add` but retuned by *final component
    index* on the fit driver, and the final index is only knowable after source
    extraction has run. Declaring both in one place removes that bookkeeping.

    Returns
    -------
    (str, dict)
        Canonical type name, and a possibly-empty dict of overrides.
    """
    if isinstance(spec, str):
        return normalise_component_type(spec), {}

    if isinstance(spec, dict):
        overrides = dict(spec)
        ctype = overrides.pop('type', None)
        if ctype is None:
            raise ValueError(
                f"{source}: a dict component specification needs a 'type' key, "
                f"e.g. {{'type': 'disk-f', 'fix_max_value_n': 3.0}}. Got "
                f"{spec!r}.")
    elif isinstance(spec, (tuple, list)) and len(spec) == 2 \
            and isinstance(spec[0], str) and isinstance(spec[1], dict):
        ctype, overrides = spec[0], dict(spec[1])
    else:
        raise ValueError(
            f"{source}: cannot read {spec!r} as a component. Use a type name, "
            f"a (type, overrides) pair, or {{'type': ..., **overrides}}.")

    unknown = sorted(set(overrides) - set(INLINE_OVERRIDE_KEYS))
    if unknown:
        raise ValueError(
            f"{source}: unknown inline override(s) {unknown} for component type "
            f"{ctype!r}. Valid names are: {', '.join(sorted(INLINE_OVERRIDE_KEYS))}.")

    return normalise_component_type(ctype), overrides


def resolve_component_types(n_IDs, ids_to_add=None, ids_types=None,
                            default_component_type=DEFAULT_COMPONENT_TYPE):
    """
    Work out the type of every model component, and what has to be cloned.

    Two independent inputs, both keyed by *parent ID* (a detected region, 1..n_IDs):

    `ids_types` re-types a component that was actually detected, without adding
    anything::

        ids_types = {'1': 'point-like'}      # detected region 1 is unresolved

    `ids_to_add` is purely additive -- the detected component stays as detected,
    and each entry in a region's list appends one *extra* component cloned from
    that region::

        ids_to_add = {'1': ['disk'], '2': ['compact', 'sersic']}

    Any component in either argument may also carry overrides for itself, which
    saves having to work out its final index to retune it on the fit driver::

        ids_to_add = {'1': [('disk-f', {'fix_max_value_n': 3.0})]}

    Ordering is the contract the rest of the codebase depends on: detected
    components keep indices `1..n_IDs`, added ones follow in dict-insertion then
    list order.

    Returns
    -------
    component_types : list of str
        Length `n_IDs + total added`; `component_types[i]` is the type of
        component `i+1`.
    additions : list of (int, str)
        `(parent_id, component_type)` in the order they must be appended.
    overrides : dict
        `{component_index: {name: value}}` for the components that declared any.
    """
    default_type = normalise_component_type(default_component_type
                                            or DEFAULT_COMPONENT_TYPE)

    component_types = [default_type] * n_IDs
    overrides_by_index = {}

    def _parent_index(raw_key, source):
        try:
            idx = int(raw_key)
        except (TypeError, ValueError):
            raise KeyError(
                f"{source}: parent ID {raw_key!r} is not a component number. "
                f"Keys must be detected region IDs in 1..{n_IDs}.") from None
        if not 1 <= idx <= n_IDs:
            raise KeyError(
                f"{source}: parent ID {idx} was not detected. This image has "
                f"{n_IDs} detected region(s), so valid IDs are 1..{n_IDs}.")
        return idx

    for raw_key, spec in (ids_types or {}).items():
        idx = _parent_index(raw_key, 'ids_types')
        ctype, inline = parse_component_spec(spec, source=f"ids_types['{raw_key}']")
        component_types[idx - 1] = ctype
        if inline:
            overrides_by_index[idx] = inline

    additions = []
    if ids_to_add:
        if not isinstance(ids_to_add, dict):
            suggestion = {}
            for entry in ids_to_add:
                suggestion.setdefault(str(entry), []).append(default_type)
            raise ValueError(
                "ids_to_add must be a dict mapping a parent ID to the list of "
                "extra components to add to that region, e.g. "
                "{'1': ['point-like', 'disk'], '2': ['sersic']}. The flat-list "
                f"form is no longer supported; {list(ids_to_add)!r} would "
                f"translate to {suggestion!r}.")

        for raw_key, specs in ids_to_add.items():
            idx = _parent_index(raw_key, 'ids_to_add')
            # A bare string, or a single (type, overrides) pair, is one component;
            # a list is several.
            if isinstance(specs, (str, dict)):
                specs = [specs]
            elif (isinstance(specs, tuple) and len(specs) == 2
                  and isinstance(specs[0], str) and isinstance(specs[1], dict)):
                specs = [specs]
            for spec in specs:
                ctype, inline = parse_component_spec(
                    spec, source=f"ids_to_add['{raw_key}']")
                additions.append((idx, ctype))
                component_types.append(ctype)
                if inline:
                    overrides_by_index[len(component_types)] = inline

    _check_supported_models(component_types)
    return component_types, additions, overrides_by_index


def estimate_psf_fwhm_px(imagename=None, psf_data=None, psf_name=None):
    """
    PSF/beam FWHM in pixels, or None when it cannot be determined.

    Presets that constrain a radius in absolute pixels (`compact`) express the
    constraint as a multiple of the PSF, since "a few resolution elements" is the
    physical statement and a pixel count is only a stand-in for it that silently
    depends on the cell size.

    Tried in order: a PSF array, a PSF file, then the restoring beam recorded in
    the image header. Any of these may be unavailable -- optical data often has no
    PSF file, radio images sometimes have no beam keywords -- so the caller must
    handle None (the presets fall back to their absolute pixel default).
    """
    for candidate in (psf_data,
                      load_fits_data(psf_name) if psf_name else None):
        if candidate is None:
            continue
        try:
            fwhm = float(psf_params(np.asarray(candidate)))
            if np.isfinite(fwhm) and fwhm > 0:
                return fwhm
        except Exception:
            pass

    if imagename is not None:
        try:
            omaj, omin, _, _, _ = beam_shape(imagename)
            cell = get_cell_size(imagename)
            fwhm = float(np.sqrt(omaj * omin) / cell)
            if np.isfinite(fwhm) and fwhm > 0:
                return fwhm
        except Exception:
            pass
    return None


def build_fit_control_maps(sources_photometries, ncomps, overrides=None,
                           default_component_type=DEFAULT_COMPONENT_TYPE,
                           warn=True):
    """
    Materialise the full set of per-component fit-control arguments.

    This is what removes the need to hand-write parallel lists like::

        dr_fix       = [5, 5, 5, 99] + [50] * SE.n_IDs + [100, 100, 15, 15]
        fix_value_n  = [0.5] * 7
        fix_max_value_Rn = [1.0, False, False, False, False, False, False]

    Each component's preset supplies its own values, so the lists are always the
    right length and always aligned with the component they describe. A user value
    for a `defaulted` parameter wins; one for a `locked` parameter is warned about
    and dropped.

    Parameters
    ----------
    sources_photometries : dict
        The `prepare_fit` output. `c{i}_type` is read for each component; a
        component with no recorded type falls back to `default_component_type`.
    ncomps : int
        Final number of model components.
    overrides : dict, optional
        `{arg_name: value}` of user-supplied per-component arguments, in any of
        the shapes `as_component_map` accepts. `None` values mean "not supplied".
    default_component_type : str, optional
        Preset for components without a recorded type.
    warn : bool, optional
        Report attempts to override a locked parameter.

    Returns
    -------
    dict
        `{arg_name: {1: value, ..., ncomps: value}}` covering `PRESET_FIT_KEYS`.
    """
    overrides = overrides or {}
    per_component_overrides = {
        name: as_component_map(value, ncomps, default=None, name=name)
        for name, value in overrides.items() if value is not None
    }
    psf_fwhm_px = sources_photometries.get('psf_fwhm_px')
    inline_overrides = sources_photometries.get('component_overrides') or {}

    resolved = {name: {} for name in PRESET_FIT_KEYS}
    for i in range(1, ncomps + 1):
        ctype = sources_photometries.get(f'c{i}_type') or default_component_type
        # Overrides declared inline with the component type sit under anything the
        # caller passed to the driver, which is the more explicit statement.
        merged = dict(inline_overrides.get(i, {}))
        merged.update({name: values[i]
                       for name, values in per_component_overrides.items()
                       if values[i] is not None})
        settings = resolve_component_preset(
            ctype, overrides=merged,
            component_id=i, warn=warn, psf_fwhm_px=psf_fwhm_px)
        for name in PRESET_FIT_KEYS:
            if name in settings:
                resolved[name][i] = settings[name]

    # Drop arguments no preset spoke to, so the caller can tell "unset" from
    # "set to False everywhere" and leave the historical default in place.
    return {name: values for name, values in resolved.items() if values}


def add_extra_component(petro_properties, copy_from_id, image_shape=None,
                        radius_scale_factor=None, intensity_scale_factor=None,
                        component_type='disk', verbose=True):
    """
    Create an additional Sérsic component by copying and scaling parameters
    from an existing component.
    
    This is useful when a single detected source requires multiple Sérsic 
    functions to model its light distribution (e.g., bulge + disk decomposition,
    or a compact core embedded in extended emission).
    
    Parameters
    ----------
    petro_properties : dict
        Dictionary containing Petrosian photometric properties for N detected
        components. Keys follow the pattern 'c{i}_{param}' where i is the 
        component ID (1-indexed) and param is one of:
        PA, q, area, Re, x0c, y0c, label, R50, Snu, Rp, Rpidx, rlast, I50
        Must also contain 'ncomps' (number of components) and optionally 
        'cg_Rp' (global Petrosian radius).
        
    copy_from_id : int
        Component ID (1-indexed) to copy parameters from. Typically the
        compact/central component when adding a more extended component.
        
    image_shape : tuple of int, optional
        Shape of the image array (ny, nx). Used to constrain maximum radius.
        If None, will attempt to estimate from 'image_shape' key in 
        petro_properties.
        
    radius_scale_factor : float, optional
        Factor to scale R50 of the new component relative to the original.
        If None, automatically estimated based on global Petrosian radius
        and the component type.
        
    intensity_scale_factor : float, optional
        Factor to scale I50 (effective surface brightness) of the new 
        component. If None, estimated from Sérsic profile assumptions.
        
    component_type : str, optional
        A key of `COMPONENT_TYPE_PRESETS`: 'point-like', 'gaussian', 'compact',
        'disk', 'sersic', or the legacy 'envelope'/'halo'. The preset supplies the
        default radius and intensity scale factors (and, elsewhere, the fit
        bounds); see the registry for the values and the rationale.

    verbose : bool, optional
        Print diagnostic information.
        
    Returns
    -------
    dict
        Updated petro_properties dictionary with the new component added.
        The new component ID will be ncomps + 1.
        
    Raises
    ------
    ValueError
        If copy_from_id does not exist in petro_properties.
        
    Notes
    -----
    Physical assumptions:
    - The new component shares the same centroid (x0c, y0c) as the original
    - Position angle (PA) and axis ratio (q) are copied but may need adjustment
    - Effective radius is scaled up (for outer components)
    - Effective intensity is scaled down (flux spread over larger area)
    - Petrosian radius (Rp) is set to the global value when available

    The scale factors are not only a starting point. `construct_model_parameters`
    derives every Rn and In bound from the seed (`R50 * Rn_max_factor`,
    `I50 * In_max_factor`, ...) and, with `NEST_ADDED_COMPONENT_BOUNDS`, clamps the
    added component's window against its parent's seed. So a growth preset states
    that this component IS the outer one and keeps it there; `compact`, the one
    preset that scales below 1.0, states the converse. A preset scaling by exactly
    1.0/1.0 copies the parent outright -- identical seed, identical bounds -- and
    the minimiser is then free to converge with the two components swapped.
    
    Examples
    --------
    >>> # Add a disk component based on a detected bulge (component 1)
    >>> props = add_extra_component(props, copy_from_id=1, component_type='disk')
    
    >>> # Add with custom scaling
    >>> props = add_extra_component(props, copy_from_id=1, 
    ...                             radius_scale_factor=3.0,
    ...                             intensity_scale_factor=0.2)
    """
    # Validate input
    if 'ncomps' not in petro_properties:
        raise ValueError("petro_properties must contain 'ncomps' key")
    
    source_key_prefix = f'c{copy_from_id}_'
    if not any(k.startswith(source_key_prefix) for k in petro_properties.keys()):
        raise ValueError(
            f"Component {copy_from_id} not found in petro_properties. "
            f"Available components: 1 to {petro_properties['ncomps']}"
        )
    
    # Extract unique parameter suffixes (e.g., 'PA', 'q', 'R50', etc.)
    param_suffixes = list(dict.fromkeys(
        key.split('_', 1)[1] 
        for key in petro_properties.keys() 
        if key.startswith('c') and '_' in key and key[1:].split('_')[0].isdigit()
    ))
    
    # Create output dictionary
    new_props = petro_properties.copy()
    new_comp_id = petro_properties['ncomps'] + 1
    
    # Determine image constraints
    if image_shape is not None:
        max_radius = 0.5 * np.sqrt(image_shape[0]**2 + image_shape[1]**2)
    elif 'image_shape' in petro_properties:
        shape = petro_properties['image_shape']
        max_radius = 0.5 * np.sqrt(shape[0]**2 + shape[1]**2)
    else:
        # Fallback: use 5× the global Petrosian radius or 500 pixels.
        # `prepare_fit` now always passes `image_shape`, so this is only reached by
        # direct callers. It is a very weak ceiling -- combined with the R50 clip
        # below it works out to 4x the global Petrosian radius, i.e. no ceiling at
        # all -- which is why the caller supplying the real shape matters.
        max_radius = petro_properties.get('cg_Rp', 100) * 5
        if verbose:
            print(f"Warning: image_shape not provided. Using max_radius = {max_radius:.1f} px")

    # Get source component parameters
    R50_source = petro_properties.get(f'c{copy_from_id}_R50', 10.0)
    I50_source = petro_properties.get(f'c{copy_from_id}_I50', 1.0)
    # The *parent region's* Petrosian radius, not the global one. `cg_Rp` describes
    # the whole source, so measuring a sub-component against it made the auto-scale
    # saturate at its ceiling for any compact parent inside a large source.
    Rp_parent = petro_properties.get(f'c{copy_from_id}_Rp')
    Rp_global = petro_properties.get('cg_Rp', R50_source * 3)
    if Rp_parent is None or not np.isfinite(Rp_parent) or Rp_parent <= 0:
        Rp_parent = Rp_global

    # Scale factors come from the component-type preset. Explicit arguments win.
    preset = resolve_component_preset(
        component_type,
        overrides={'radius_scale_factor': radius_scale_factor,
                   'intensity_scale_factor': intensity_scale_factor},
        component_id=new_comp_id)
    component_type = preset['component_type']
    radius_scale_factor = preset.get('radius_scale_factor')
    intensity_scale_factor = preset.get('intensity_scale_factor')

    # Calculate radius scale factor
    if radius_scale_factor is None:
        # No preset value and none supplied: fall back on the parent region's own
        # extent as a guide.
        if Rp_parent > R50_source:
            ratio = Rp_parent / R50_source
            radius_scale_factor = np.clip(ratio, 1.5, 5.0)
        else:
            radius_scale_factor = 2.5

    # Calculate intensity scale factor
    if intensity_scale_factor is None:
        # For Sérsic profiles, if we assume similar n and total flux contribution,
        # I_e scales roughly as (R_e)^{-2} for fixed flux
        # But new component likely has less flux, so we apply additional reduction
        intensity_scale_factor = 0.30 / (radius_scale_factor / 2.5)
        intensity_scale_factor = np.clip(intensity_scale_factor, 0.01, 0.5)

    if verbose:
        print(f"Adding component {new_comp_id} (type: {component_type}) "
              f"from component {copy_from_id}")
        print(f"  Radius scale factor: {radius_scale_factor:.2f}")
        print(f"  Intensity scale factor: {intensity_scale_factor:.3f}")
    
    # A preset that shrinks its component (`compact`, radius_scale_factor 0.5) must
    # not have that undone by a floor written for components that only ever grow.
    # The floor is only meaningful in the growth direction.
    grows = radius_scale_factor > 1.0

    # Copy and adjust parameters
    for param in param_suffixes:
        source_key = f'c{copy_from_id}_{param}'
        new_key = f'c{new_comp_id}_{param}'

        if source_key not in petro_properties:
            continue

        source_value = petro_properties[source_key]

        # Apply parameter-specific transformations
        if param == 'R50':
            new_value = R50_source * radius_scale_factor
            # Constrain to reasonable bounds
            floor = R50_source * 1.2 if grows else 0.5
            new_value = np.clip(new_value, floor, max_radius * 0.8)
            if verbose and new_value >= max_radius * 0.8:
                print(f"  Warning: R50 clipped to {new_value:.1f} px (80% of max_radius)")

        elif param == 'Re_s':
            # Effective radius follows same scaling as R50
            Re_source = petro_properties.get(f'c{copy_from_id}_Re', R50_source)
            new_value = Re_source * radius_scale_factor
            floor = Re_source * 1.2 if grows else 0.5
            new_value = np.clip(new_value, floor, max_radius * 0.8)

        elif param == 'I50':
            new_value = I50_source * intensity_scale_factor
            
        elif param == 'Snu':
            # Total flux: new component contributes fraction of original
            # Rough estimate based on intensity and area scaling
            flux_fraction = intensity_scale_factor * radius_scale_factor**2
            flux_fraction = np.clip(flux_fraction, 0.1, 2.0)
            new_value = source_value * flux_fraction
            
        elif param == 'Rp':
            # Scale the *parent region's* Petrosian radius rather than adopting the
            # global one: a `compact` component added inside a large source is not
            # as extended as the whole source. Only consumed for plot extents
            # (plotting.py:800), never by the fit itself.
            new_value = min(Rp_parent * radius_scale_factor, max_radius * 0.8)

        elif param == 'Rpidx':
            # Update Petrosian radius index
            new_value = int(2 * min(Rp_parent * radius_scale_factor,
                                    max_radius * 0.8))

        elif param == 'area':
            # Scale area with R50²
            new_value = source_value * radius_scale_factor**2
            
        elif param == 'label':
            # Assign new label
            new_value = new_comp_id
            
        else:
            # Copy unchanged: x0c, y0c, PA, q, rlast, etc.
            new_value = source_value
            
        new_props[new_key] = new_value
    
    # Update component count
    new_props['ncomps'] = new_comp_id
    # Record what this component is and where it came from, so the fit drivers can
    # look up its preset later without re-deriving anything.
    new_props[f'c{new_comp_id}_type'] = component_type
    new_props[f'c{new_comp_id}_parent'] = int(copy_from_id)

    if verbose:
        R50_new = new_props.get(f'c{new_comp_id}_R50', 0)
        I50_new = new_props.get(f'c{new_comp_id}_I50', 0)
        print(f"  New component: R50={R50_new:.2f} px, I50={I50_new:.4f}")
    
    return new_props



def sorted_detected_coordinates(reference_x, 
                                reference_y, 
                                detected_x, 
                                detected_y, 
                                reference_coordinate, 
                                tolerance=2):
    """
    This function will sort a new set of detected coordinates 
    by distance to a reference position in the same order as a reference set 
    of coordinates. 
    
    Consider for example, that in a high-resolution image, we have a set of
    detected structures such as ID1, ID2, ID3. 
    
    If in a new image (e.g. a low-resolution image) we detect the same structures
    in addition to other structures, we will have for example: 
        ID1_new, ID2_new, ID3_new, ID4_new, ID5_new.
    However, the order of the detected structures can be any. For example,
    ID1_new may be the same as ID1, but not for the others. 
    
    This function will sort the new coordinates in the same order they appear 
    in the reference coordinates, in addition to the extra detected coordinates.
    
    A typical use case is: 
    - Perform a source detection in a VLA image (reference image), at 33 GHz.
    - Perform a source detection in an e-MERLIN image, at 6 GHz.
    The number of Structures may differ, but some of them are the same.
    This allow us to connect the labels of the structures in both images.
    
    """
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print('+++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    # Reference coordinate
    # Reference coordinate
    x_ref, y_ref = reference_coordinate

    # Calculate distances to the reference coordinate for both sets of coordinates
    reference_distances = np.sqrt((reference_x - x_ref)**2 + (reference_y - y_ref)**2)
    detected_distances = np.sqrt((detected_x - x_ref)**2 + (detected_y - y_ref)**2)

    # Sort the reference and detected indices by distance
    sorted_reference_indices = np.argsort(reference_distances)
    sorted_detected_indices = np.argsort(detected_distances)

    # Initialize lists to store the sorted detected coordinates, distances, and indices
    sorted_detected_coordinates = []
    sorted_detected_distances = []
    sorted_detected_original_indices = []

    # Match reference coordinates to detected coordinates by closest distance
    for ref_idx in sorted_reference_indices:
        if len(sorted_detected_indices) > 0:
            closest_idx = sorted_detected_indices[0]
            sorted_detected_coordinates.append([detected_x[closest_idx], detected_y[closest_idx]])
            sorted_detected_distances.append(detected_distances[closest_idx])
            sorted_detected_original_indices.append(closest_idx)
            # Remove the matched index to avoid duplicate matching
            # sorted_detected_indices = np.delete(sorted_detected_indices, 0)

    # Add any remaining detected coordinates that were not matched
    for det_idx in sorted_detected_indices:
        sorted_detected_coordinates.append([detected_x[det_idx], detected_y[det_idx]])
        sorted_detected_distances.append(detected_distances[det_idx])
        sorted_detected_original_indices.append(det_idx)

    return (np.array(sorted_detected_coordinates), 
            np.array(sorted_detected_distances), 
            np.array(sorted_detected_indices))
    


# def prepare_fit(ref_image, ref_res, z, ids_to_add=[1],
#                 bw=51, bh=51, fw=15, fh=15, sigma=15, ell_size_factor=2.0,
#                 deblend_cont=1e-7, deblend_nthresh=15,
#                 minarea=None,minarea_factor=1.0,npixels=None,
#                 sigma_mask=6,mask=None,dilation_size=None,
#                 show_detection=True,use_extraction_positions=False,
#                 clean_param=0.9,clean=True,sort_by='distance',
#                 apply_mask=False, mask_grow_iterations=3,
#                 obs_type = 'radio',algorithm='SEP',
#                 force_circular=True,force_circular_all=False,
#                 show_petro_plots=False):
#     """
#     Prepare the imaging data to be modeled.

#     This function runs a source extraction, computes basic petrosian properties
#     from the data for each detected source as well as shape morphology
#     (e.g. position angle, axis ration, effective intensity and radii).
#     """
#     crop_image = ref_image
#     crop_residual = ref_res
#     data_2D = load_fits_data(crop_image)
#     if crop_residual is not None:
#         residual_2D = load_fits_data(crop_residual)
#     else:
#         residual_2D = None
#     if minarea is None:
#         try:
#             minarea = int(beam_area2(crop_image))
#         except:
#             minarea = data_2D.shape[0]/30
#     pix_to_pc = pixsize_to_pc(z=z,
#                               cell_size=get_cell_size(crop_image))
#     #     eimshow(crop_image, vmin_factor=5)
#     try:
#         std_res = mad_std(load_fits_data(crop_residual))
#     except:
#         std_res = mad_std(data_2D)

#     if mask is not None:
#         mask_detection = mask
#         apply_mask = False

#     if apply_mask:
#         _, mask = mask_dilation(crop_image, sigma=sigma_mask, dilation_size=None,
#                                 iterations=2, rms=std_res)


#     mask_detection = mask.copy()

#     # if apply_mask == True:
#     #     mask_detection = mask
#     # else:
#     #     mask_detection = None #np.ones(data_2D.shape)
#     # # plt.figure()

#     # _, mask = mask_dilation(crop_image, sigma=6, dilation_size=None,
#     #                         iterations=2)
#     if algorithm == 'SEP':
#         masks, indices, bkg, seg_maps, objects = \
#             sep_source_ext(crop_image, 
#                            residualname=crop_residual,
#                            bw=bw,bh=bh,fw=fw, fh=fh,
#                            minarea=minarea,
#                            minarea_factor=minarea_factor,
#                            segmentation_map=True,
#                            filter_type='matched',
#                            mask = mask_detection,
#                            deblend_nthresh=deblend_nthresh,
#                            deblend_cont=deblend_cont,
#                            clean_param=clean_param,
#                            clean=clean,
#                            sort_by=sort_by,
#                            npixels = npixels,
#                            dilation_size=dilation_size,
#                            iterations = mask_grow_iterations,
#                            sigma=sigma,sigma_mask=sigma_mask,
#                            ell_size_factor=ell_size_factor,
#                            apply_mask=apply_mask,
#                            show_detection=show_detection)
#     if algorithm == 'PF':
#         masks, indices, bkg, seg_maps, objects = \
#             phot_source_ext(crop_image, 
#                            residual=residual_2D,
#                            bw=bw,
#                            bh=bh,
#                            fw=fw, fh=fh,
#                            minarea=minarea,
#                            minarea_factor=minarea_factor,
#                            segmentation_map=True,
#                            filter_type='matched', mask=mask,
#                            deblend_nthresh=deblend_nthresh,
#                            deblend_cont=deblend_cont,
#                            clean_param=clean_param,
#                            clean=clean,
#                            sort_by=sort_by,
#                            sigma=sigma,sigma_mask=sigma_mask,
#                            iterations = mask_grow_iterations,
#                            ell_size_factor=ell_size_factor,
#                            apply_mask=apply_mask,
#                            show_detection=show_detection)



#     sigma_level = 3
#     vmin = 3
#     # i = 0 #to be used in indices[0], e.g. first component
#     sources_photometries = {}  # init dict to store values.
#     # if use_extraction_positions == True:
#     #     for i in range(len(indices)):
#     #         # ii = str(i+1)
#     #         positions = np.array([objects['xc'][i], objects['yc'][i]])
#     #         mask_component = masks[indices[i]]
#     #         data_component = data_2D * mask_component
#     #         sources_photometries = compute_petro_source(data_component,
#     #                                                     mask_component=mask_component,
#     #                                                     sigma_level=1,positions=positions,
#     #                                                     i=i, plot=show_petro_plots,
#     #                                                     source_props=sources_photometries)

#     # else:
    
    
#     sources_photometries = compute_petro_source(data_2D = data_2D,
#                                                 imagename = crop_image,
#                                                 sigma_level=6,
#                                                 nlevels=1, contrast=1,
#                                                 deblend=False, npixels=None,
#                                                 i='g', mask_component = mask_detection,
#                                                 plot=show_petro_plots,
#                                                 source_props=sources_photometries)
    
#     for i in tqdm(range(len(indices))):
#         # ii = str(i+1)
#         mask_component = masks[i]
#         data_component = data_2D * mask_component
#         sources_photometries = compute_petro_source(data_2D = data_component,
#                                                     imagename = crop_image,
#                                                     mask_component=mask_component,
#                                                     obs_type = obs_type,
#                                                     sigma_level=1,nlevels=1, contrast=1,
#                                                     deblend=False, npixels=None,
#                                                     i=i, plot=show_petro_plots,
#                                                     source_props=sources_photometries)


#     if obs_type == 'radio':
#         """
#         PSF image is contained within the header of the original image. 
#         A new psf file will be created (as `psf_name`). 
#         """
#         # omaj, omin, _, _, _ = beam_shape(crop_image)
#         # dilation_size = int(
#         #     np.sqrt(omaj * omin) / (2 * get_cell_size(crop_image)))
#         # psf_image_size = dilation_size*6
#         # psf_image_size = (2 * psf_image_size) // 2 +1
#         psf_image_size = int(data_2D.shape[0])
#         # print('++==>> PSF IMAGE SIZE is', psf_image_size)
#         # creates a psf from the beam shape.
#         psf_name = tcreate_beam_psf(crop_image, size=(
#             psf_image_size, psf_image_size))  # ,app_name='_'+str(psf_image_size)+'x'+str(psf_image_size)+'')
#     if obs_type == 'other':
#         """
#         Provide a psf file.
#         """
#         psf_name = None

#     n_components = len(indices)
#     n_IDs = len(indices)
#     sources_photometries['ncomps'] = n_components
#     sources_photometries['nIDs'] = n_IDs

    
#     print("# of structures (IDs) to be fitted =", n_IDs)
#     # sources_photometies_new = sources_photometies
#     # n_components_new = n_components
#     if ids_to_add is not None:
#         for id_to_add in ids_to_add:
#             sources_photometries = add_extra_component(sources_photometries,
#                                                        copy_from_id=id_to_add)
    
#     if force_circular:
#         print(f"<> Forcing components to be circular.")
#         # for i in range(len(indices)):
#         for i in range(sources_photometries['ncomps']):
#             if force_circular_all:
#                 sources_photometries[f'c{i+1}_q'] = 0.99
#             else:
#                 if i+1 <= n_IDs:
#                     sources_photometries[f'c{i+1}_q'] = 0.99
#                 else:
#                     sources_photometries[f'c{i+1}_q'] = 0.5
#             print(f"<> q_{i+1} {sources_photometries[f'c{i+1}_q']}")
#     else:
#         print(f"<> Leaving components to be elliptical.")
#         for i in range(sources_photometries['ncomps']):
#             # if i+1 > n_IDs:
#             #     sources_photometries[f'c{i+1}_q'] = 0.6
#             sources_photometries[f'c{i+1}_q'] = 0.5
#                 # sources_photometries[f'c{i+1}_R50'] = sources_photometries[f'cg_Rp']
#                 # sources_photometries[f'c{i+1}_Rp'] = sources_photometries[f'cg_Rp']
#             print('<> q = ',sources_photometries[f'c{i+1}_q'])

#     # update variable `n_components`.
#     n_components = sources_photometries['ncomps']
#     print("# of model components (COMPS) to be fitted =", n_components)
#     return (sources_photometries, n_components, n_IDs, masks, indices, objects,
#             psf_name, mask, bkg)

def prepare_fit(ref_image, ref_res, z, ids_to_add=None,
                ids_types=None,
                default_component_type=DEFAULT_COMPONENT_TYPE,
                bw=51, bh=51, fw=15, fh=15, sigma=15,
                psf_data = None,
                ell_size_factor=2.0,
                deblend_cont=1e-7, deblend_nthresh=15,
                minarea=None, minarea_factor=1.0, npixels=None,
                sigma_mask=6, mask=None, dilation_size=None,
                show_detection=True, use_extraction_positions=False,
                clean_param=0.9, clean=True, sort_by='distance',
                apply_mask=False, mask_grow_iterations=3,
                obs_type='radio', algorithm='SEP', SE = None,
                force_circular=True, force_circular_all=False,
                show_petro_plots=False,
                # Robust-specific parameters
                multiscale_levels=3,
                extended_threshold_factor=0.5,
                watershed_connectivity=2,
                merge_threshold=0.3,
                adaptive_background=True,
                preserve_extended=True,
                min_separation=None,
                error_map=None, use_residual_as_error=False):
    """
    Prepare the imaging data to be modeled.

    This function runs a source extraction, computes basic petrosian properties
    from the data for each detected source as well as shape morphology
    (e.g. position angle, axis ration, effective intensity and radii).

    Model components
    ----------------
    Each detected region contributes one model component; `ids_to_add` adds more
    on top of a region, and `ids_types` changes what a detected component *is*.
    Both are keyed by parent ID -- the 1-indexed detected region -- and both take
    names from `COMPONENT_TYPE_PRESETS`::

        ids_types  = {'1': 'point-like'}                        # retype a detection
        ids_to_add = {'1': ['disk'], '2': ['compact', 'sersic']} # add extra components

    With two detections that yields five components: `c1` point-like and `c2`
    `default_component_type` (both detected), then `c3` disk, `c4` compact and `c5`
    sersic. Detected components always keep indices `1..n_IDs`; added ones follow
    in dict-insertion then list order.

    Parameters
    ----------
    ids_to_add : dict, optional
        `{parent_id: [component_type, ...]}`. A region may be given more than one
        extra component. The old flat-list form is no longer accepted; passing one
        raises a `ValueError` that shows the equivalent dict.
    ids_types : dict, optional
        `{parent_id: component_type}` for *detected* components. Untyped
        detections fall back to `default_component_type`.
    force_circular, force_circular_all : bool, optional
        DEPRECATED. Make the detected components (or, with `_all`, every
        component) round. Kept working for now, but they can only speak about
        whole groups; state it on the component instead::

            ids_types  = {'1': ('point-like', {'force_circular': True})}
            ids_to_add = {'1': [('sersic', {'force_circular': 0.05})]}

        A component that declares its own value ignores these flags.
    default_component_type : str, optional
        Preset used for any component not explicitly typed. The fit drivers set
        this per observation kind ('gaussian' for radio, 'sersic' for optical);
        it is deliberately not hardcoded here.
    """
    crop_image = ref_image
    crop_residual = ref_res
    data_2D = load_fits_data(crop_image)
    if crop_residual is not None:
        residual_2D = load_fits_data(crop_residual)
    else:
        residual_2D = None
    _error_for_petro, _ = resolve_flux_error_map(
        data_2D, error_map=error_map,
        residual_map=(residual_2D if use_residual_as_error and error_map is None else None),
        rms=None)
    if minarea is None:
        try:
            minarea = int(beam_area2(crop_image))
        except:
            minarea = data_2D.shape[0]/30
    pix_to_pc = pixsize_to_pc(z=z,
                              cell_size=get_cell_size(crop_image))
    #     eimshow(crop_image, vmin_factor=5)
    try:
        std_res = mad_std(load_fits_data(crop_residual))
    except:
        std_res = mad_std(data_2D)

    if mask is not None:
        mask_detection = mask
        apply_mask = False

    if apply_mask:
        _, mask = mask_dilation(crop_image, sigma=sigma_mask, dilation_size=None,
                                iterations=2, rms=std_res)

    mask_detection = mask.copy()

    if SE is not None:
        print(" ++>> Using provided Source Extractor (SE) object for source extraction.")
        masks = SE.masks
        indices = SE.indices
        bkg = SE.bkg
        objects = SE.objects
    else:
        # Source extraction based on selected algorithm
        if algorithm == 'SEP':
            masks, indices, bkg, seg_maps, objects = \
                sep_source_ext(crop_image, 
                            residualname=crop_residual,
                            bw=bw, bh=bh, fw=fw, fh=fh,
                            minarea=minarea,
                            minarea_factor=minarea_factor,
                            segmentation_map=True,
                            filter_type='matched',
                            mask=mask_detection,
                            deblend_nthresh=deblend_nthresh,
                            deblend_cont=deblend_cont,
                            clean_param=clean_param,
                            clean=clean,
                            sort_by=sort_by,
                            npixels=npixels,
                            dilation_size=dilation_size,
                            iterations=mask_grow_iterations,
                            sigma=sigma,
                            sigma_mask=sigma_mask,
                            ell_size_factor=ell_size_factor,
                            apply_mask=apply_mask,
                            show_detection=show_detection)
        
        if algorithm == 'PF':
            masks, indices, bkg, seg_maps, objects, cat = \
                phot_source_ext(crop_image, 
                            residual=residual_2D,
                            psf_data = psf_data,
                            bw=bw,
                            bh=bh,
                            fw=fw, fh=fh,
                            minarea=minarea,
                            minarea_factor=minarea_factor,
                            segmentation_map=True,
                            filter_type='matched', mask=mask,
                            deblend_nthresh=deblend_nthresh,
                            deblend_cont=deblend_cont,
                            clean_param=clean_param,
                            clean=clean,
                            sort_by=sort_by,
                            sigma=sigma,
                            sigma_mask=sigma_mask,
                            iterations=mask_grow_iterations,
                            ell_size_factor=ell_size_factor,
                            apply_mask=apply_mask,
                            show_detection=show_detection)
                
        if algorithm == 'astphot':
            masks, indices, bkg, seg_maps, objects, cat = \
                astphot_source_ext(crop_image, 
                            residual=residual_2D,
                            psf_data = psf_data,
                            bw=bw,
                            bh=bh,
                            fw=fw, fh=fh,
                            minarea=minarea,
                            minarea_factor=minarea_factor,
                            segmentation_map=True,
                            filter_type='matched', mask=mask,
                            deblend_nthresh=deblend_nthresh,
                            deblend_cont=deblend_cont,
                            clean_param=clean_param,
                            clean=clean,
                            sort_by=sort_by,
                            sigma=sigma,
                            sigma_mask=sigma_mask,
                            iterations=mask_grow_iterations,
                            ell_size_factor=ell_size_factor,
                            apply_mask=apply_mask,
                            show_detection=show_detection)
        
        if algorithm == 'robust':
            masks, indices, bkg, seg_maps, objects = \
                robust_source_ext(crop_image,
                                residualname=crop_residual,
                                bw=bw, bh=bh, fw=fw, fh=fh,
                                minarea=minarea,
                                minarea_factor=minarea_factor,
                                segmentation_map=True,
                                filter_type='matched',
                                mask=mask_detection,
                                deblend_nthresh=deblend_nthresh,
                                deblend_cont=deblend_cont,
                                clean_param=clean_param,
                                clean=clean,
                                sort_by=sort_by,
                                npixels=npixels,
                                dilation_size=dilation_size,
                                iterations=mask_grow_iterations,
                                sigma=sigma,
                                sigma_mask=sigma_mask,
                                ell_size_factor=ell_size_factor,
                                apply_mask=apply_mask,
                                show_detection=show_detection,
                                show_bkg_map=False,
                                # Robust-specific parameters
                                multiscale_levels=multiscale_levels,
                                extended_threshold_factor=extended_threshold_factor,
                                watershed_connectivity=watershed_connectivity,
                                merge_threshold=merge_threshold,
                                adaptive_background=adaptive_background,
                                preserve_extended=preserve_extended,
                                min_separation=min_separation)

    sigma_level = 3
    vmin = 3
    # i = 0 #to be used in indices[0], e.g. first component
    sources_photometries = {}  # init dict to store values.
    
    sources_photometries = compute_petro_source(data_2D=data_2D,
                                                imagename=crop_image,
                                                sigma_level=6,
                                                nlevels=1, contrast=1,
                                                deblend=False, npixels=None,
                                                i='g', #'g' stands for global properties (full image)
                                                mask_component=mask_detection,
                                                plot=show_petro_plots,error=_error_for_petro,
                                                source_props=sources_photometries)
    
    for i in tqdm(range(len(indices))):
        # ii = str(i+1)
        mask_component = masks[i]
        data_component = data_2D * mask_component
        sources_photometries = compute_petro_source(data_2D=data_component,
                                                    imagename=crop_image,
                                                    mask_component=mask_component,
                                                    obs_type=obs_type,
                                                    sigma_level=1, nlevels=2, contrast=1,
                                                    deblend=False, npixels=None,
                                                    i=i, #for each source i in the image, compute the properties.
                                                    plot=show_petro_plots,error=_error_for_petro,
                                                    source_props=sources_photometries)

    if obs_type == 'radio':
        """
        PSF image is contained within the header of the original image. 
        A new psf file will be created (as `psf_name`). 
        """
        psf_image_size = int(data_2D.shape[0])
        # creates a psf from the beam shape.
        psf_name = tcreate_beam_psf(crop_image, size=(
            psf_image_size, psf_image_size))
    if obs_type == 'other':
        """
        Provide a psf file.
        """
        psf_name = None

    n_components = len(indices)
    n_IDs = len(indices)
    sources_photometries['ncomps'] = n_components
    sources_photometries['nIDs'] = n_IDs

    # Recorded so that presets constraining a radius in "a few PSF widths"
    # (`compact`) can turn that into pixels. Carried on the photometry dict because
    # both consumers -- `build_fit_control_maps` in the drivers and
    # `construct_model_parameters` -- already receive it.
    psf_fwhm_px = estimate_psf_fwhm_px(imagename=crop_image, psf_data=psf_data,
                                       psf_name=psf_name)
    sources_photometries['psf_fwhm_px'] = psf_fwhm_px
    if psf_fwhm_px is None:
        print("<> PSF/beam size unavailable; presets that scale a radius with the "
              "PSF will fall back to their absolute pixel defaults.")
    else:
        print(f"<> PSF FWHM = {psf_fwhm_px:.2f} px")

    print("# of structures (IDs) to be fitted =", n_IDs)

    # Work out what every component is before touching the photometry, so that a
    # bad parent ID or a typo'd type fails here rather than half-way through
    # appending components.
    component_types, additions, component_overrides = resolve_component_types(
        n_IDs, ids_to_add=ids_to_add, ids_types=ids_types,
        default_component_type=default_component_type)

    for i, ctype in enumerate(component_types[:n_IDs]):
        sources_photometries[f'c{i+1}_type'] = ctype
        sources_photometries[f'c{i+1}_parent'] = i + 1

    for offset, (parent_id, ctype) in enumerate(additions):
        # Scale factors declared inline apply to this component's initial guess,
        # so they have to reach `add_extra_component` rather than only the fit.
        inline = component_overrides.get(n_IDs + offset + 1, {})
        sources_photometries = add_extra_component(
            sources_photometries,
            copy_from_id=parent_id,
            component_type=ctype,
            image_shape=data_2D.shape,
            radius_scale_factor=inline.get('radius_scale_factor'),
            intensity_scale_factor=inline.get('intensity_scale_factor'))

    sources_photometries['component_types'] = component_types
    sources_photometries['component_overrides'] = component_overrides

    print("<> Component types:")
    for i, ctype in enumerate(component_types):
        origin = ('detected' if i < n_IDs
                  else f'added from ID {additions[i - n_IDs][0]}')
        # Report what the fit will actually start from, not only the measurement:
        # a `point-like` component whose region measures R50 = 8 px still starts
        # at Rn ~ 1 px, because its preset caps the radius. `c{i}_R50` stays the
        # measurement -- `add_extra_component` seeds the child from it.
        _p = resolve_component_preset(
            ctype, warn=False, component_id=i + 1, observation_type=obs_type,
            psf_fwhm_px=psf_fwhm_px,
            overrides=component_overrides.get(i + 1, {}))
        _R50_i = sources_photometries[f'c{i+1}_R50']
        _I50_i = sources_photometries[f'c{i+1}_I50']
        _Rn_i, _ = effective_initial_Rn(
            _R50_i,
            fix_max_value_Rn=_p.get('fix_max_value_Rn', False),
            fix_min_value_Rn=_p.get('fix_min_value_Rn', False))
        _meas = ('' if np.isclose(_Rn_i, _R50_i)
                 else f"  (R50 meas {_R50_i:.2f})")
        print(f"<>   c{i+1}: {ctype:<11s} {origin:<16s} "
              f"Rn={_Rn_i:7.2f} px  I50={_I50_i:.4e}{_meas}")

    # DEPRECATED, and the coarsest of the three ways to say this. Circularity is
    # a per-component control now: declare it with the component
    # (`ids_types = {'1': ('point-like', {'force_circular': True})}`) and it wins
    # over these flags. `c{i}_force_circular` below is the fallback
    # `construct_model_parameters` consults when the component said nothing; the
    # `c{i}_q` values are, as they always were, only the initial guess.
    if force_circular:
        if FORCE_CIRCULAR_ELL_MAX <= 0.0:
            print("<> Forcing components to be circular "
                  "(ell = 0, PA and cg frozen).")
        else:
            print("<> Forcing components to be circular "
                  f"(ell capped at {FORCE_CIRCULAR_ELL_MAX}).")
        print("<> Note: force_circular/force_circular_all are deprecated. "
              "Declare it per component instead, e.g. "
              "ids_types = {'1': ('point-like', {'force_circular': True})}.")
        for i in range(sources_photometries['ncomps']):
            # Without `force_circular_all` only the detected components are made
            # round; anything appended by `ids_to_add` stays free.
            circular_i = bool(force_circular_all or (i + 1 <= n_IDs))
            sources_photometries[f'c{i+1}_q'] = 0.99 if circular_i else 0.5
            sources_photometries[f'c{i+1}_force_circular'] = circular_i
            print(f"<> q_{i+1} {sources_photometries[f'c{i+1}_q']}"
                  f" (circular: {circular_i})")
    else:
        print(f"<> Leaving components to be elliptical.")
        for i in range(sources_photometries['ncomps']):
            sources_photometries[f'c{i+1}_q'] = 0.5
            sources_photometries[f'c{i+1}_force_circular'] = False
            print('<> q = ', sources_photometries[f'c{i+1}_q'])

    # update variable `n_components`.
    n_components = sources_photometries['ncomps']
    print("# of model components (COMPS) to be fitted =", n_components)
    return (sources_photometries, n_components, n_IDs, masks, indices, objects,
            psf_name, mask, bkg)



def mf_gaussian2D(x0, y0, sigma, M, N):
    x, y = np.meshgrid(np.arange(N) - x0, np.arange(M) - y0)
    r2 = (x) ** 2 + (y) ** 2
    mask = 1.0 / np.sqrt(2.0 * np.pi * sigma**2.0) * np.exp(-r2 / (2.0 * sigma**2.0))
    return mask

def PSFfitgaussian(data, param0):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""

    def gaussian(base, height, center_x, center_y, width):
        """
        Returns a 2D gaussian function with the given parameters
        """
        return lambda y, x: base + height * np.exp(
            -(((x - center_x) / width) ** 2.0 + ((y - center_y) / width) ** 2.0) / 2.0
        )

    errorfunction = lambda p: np.ravel(gaussian(*p)(*np.indices(data.shape)) - data)

    p, success = scipy.optimize.leastsq(errorfunction, param0)

    return (p, success)

def psf_params(psf):
    """Calculates the PSF full width half maximum
    from its mean standard deviation across axis.

    Use image moments to evaluate center and then sigma
    (first and second moments)

    <<Morfometryka libs.>>

    """
    mm, nn = psf.shape

    y, x = np.indices(psf.shape)

    # center
    x0 = (x * psf).sum() / psf.sum()
    y0 = (y * psf).sum() / psf.sum()

    # width
    sigmax = np.sqrt(((x - x0) ** 2 * psf).sum() / psf.sum())
    sigmay = np.sqrt(((y - y0) ** 2 * psf).sum() / psf.sum())

    psfsigma = (sigmax + sigmay) / 2.0

    # FIT GAUSSIAN to PSF and make a Synthetic one

    params0 = [0, psf.max(), x0, y0, psfsigma]
    psf_fitparams, fitsucess = PSFfitgaussian(psf, params0)

    if fitsucess > 0:
        psfsigma = np.abs(psf_fitparams[4])
        M = N = int(psfsigma * 11)
        psf = mf_gaussian2D(M / 2.0, N / 2.0, psfsigma, M, N)
        psf = psf / psf.sum()

    psffwhm = 2.35 * psfsigma
    return psffwhm


def _ensure_2D_image(data, name='array', logger=None):
    """
    Drop degenerate leading axes from an image array, returning a 2-D array.

    Interferometric FITS products are routinely written with NAXIS=4, i.e.
    shape (1, 1, ny, nx), the two leading axes being the degenerate
    Stokes/frequency ones. Everything else in morphen reads images through
    `load_fits_data`, which drops those axes, so every mask/map handed to the
    fitting routines is 2-D. `do_fit2D`, however, reads its image and residual
    with a raw `pf.getdata` and takes user arrays as given, so a 4-D input used
    to propagate into the fit: `_setup_coordinate_grid` would build the model
    meshgrid from shape[0] and shape[1] (both 1), producing a 1x1 Sersic model
    broadcast against a 4-D image. With a fitting mask that raised
    `IndexError: boolean index did not match shape of indexed array`; without
    one it silently converged on a 1x1 model. Normalising every array at the
    entry point of `do_fit2D` keeps the whole fit 2-D regardless of how the
    input cutouts were written (`t_cutout_2D_radec` preserves the 4-D shape,
    the older `cutout_2D_radec` collapsed it to 2-D).

    Parameters
    ----------
    data : array, scalar or None
        Array to normalise. None and scalars are returned untouched, so this is
        safe to apply to optional arguments and to flat-sky background levels.
    name : str
        Name used in log/error messages.
    logger : logger, optional
        Logger instance.

    Returns
    -------
    array, scalar or None
        `data` with its degenerate leading axes removed.
    """
    if data is None:
        return None

    ndim = np.ndim(data)
    if ndim <= 2:
        # Scalars (flat-sky levels) and already-2D arrays pass through.
        return data

    shape = np.shape(data)
    if not all(s == 1 for s in shape[:-2]):
        raise ValueError(f"{name}: expected a 2D image, got shape {shape}. "
                         f"Only degenerate (length-1) leading axes can be "
                         f"dropped automatically.")

    if logger is not None:
        logger.debug(f" ==> Dropping degenerate axes of {name}: "
                     f"{shape} -> {shape[-2:]}.")
    # reshape (not squeeze) so that both numpy and jax arrays are handled and
    # a genuine (1, 1) image would not be collapsed away.
    if hasattr(data, 'reshape'):
        return data.reshape(shape[-2:])
    return np.asarray(data).reshape(shape[-2:])


def _initialize_image_data(imagename, data_2D_, convolution_mode, logger=None):
    """
    Load and prepare image data for fitting.

    Parameters
    ----------
    imagename : str
        Path to image file
    data_2D_ : array or None
        Pre-loaded image data
    convolution_mode : str
        'GPU' or 'CPU'
    logger : logger, optional
        Logger instance

    Returns
    -------
    data_2D : array
        Image data in numpy format
    data_2D_gpu : array or None
        Image data in JAX format (if GPU mode)
    """
    if data_2D_ is None:
        data_2D = pf.getdata(imagename)
    else:
        data_2D = data_2D_

    # A 4D (NAXIS=4) image would make the model grid 1x1, see _ensure_2D_image.
    data_2D = _ensure_2D_image(data_2D, name=f'image data ({imagename})',
                               logger=logger)

    # Prepare GPU array if needed
    if convolution_mode == 'GPU':
        data_2D_gpu = jnp.array(data_2D)
    else:
        data_2D_gpu = None

    return data_2D, data_2D_gpu


def _initialize_psf_data(psf_name, convolution_mode, logger=None):
    """
    Load and prepare PSF data for convolution.

    Parameters
    ----------
    psf_name : str or None
        Path to PSF file
    convolution_mode : str
        'GPU' or 'CPU'
    logger : logger, optional
        Logger instance

    Returns
    -------
    PSF_CONV : bool
        Whether PSF convolution is enabled
    PSF_DATA : array or None
        PSF data in appropriate format for computation
    PSF_DATA_raw : array or None
        PSF data in raw numpy format
    """
    if psf_name is not None:
        PSF_CONV = True
        try:
            PSF_DATA_raw = pf.getdata(psf_name)
            if len(PSF_DATA_raw.shape) == 4:
                PSF_DATA_raw = PSF_DATA_raw[0][0]
        except:
            PSF_DATA_raw = load_fits_data(psf_name)

        if convolution_mode == 'GPU':
            if logger is not None:
                logger.debug(f"---------------------------------------")
                logger.debug(f" <<< PERFORMING CONVOLUTION WITH JAX >>> ")
                logger.debug(f"---------------------------------------")
            PSF_DATA = jnp.array(PSF_DATA_raw)
        elif convolution_mode == 'CPU':
            PSF_DATA = PSF_DATA_raw
        else:
            PSF_DATA = PSF_DATA_raw

        # PSF_DATA = pf.getdata(
        #     imagename.replace('-image.cutout.fits', '-beampsf.cutout.fits'))
    else:
        PSF_CONV = False
        PSF_DATA = None
        PSF_DATA_raw = None

    return PSF_CONV, PSF_DATA, PSF_DATA_raw

# def _create_psf_exclusion_mask(PSF_DATA_raw, Npsf, image_shape, convolution_mode, logger=None):
#     """
#     Create a mask that excludes the central region based on PSF FWHM.
    
#     Parameters
#     ----------
#     PSF_DATA_raw : array or None
#         PSF data in raw numpy format
#     Npsf : float
#         Number of PSF FWHM to exclude from center (0 = no exclusion)
#     image_shape : tuple
#         Shape of the image (M, N)
#     convolution_mode : str
#         'GPU' or 'CPU'
#     logger : logger, optional
#         Logger instance
        
#     Returns
#     -------
#     psf_mask : array or None
#         Boolean mask (True = include in fit, False = exclude)
#         None if Npsf <= 0 or no PSF provided
#     """
#     if Npsf <= 0 or PSF_DATA_raw is None:
#         return None
        
#     # Calculate PSF FWHM
#     psf_fwhm = psf_params(PSF_DATA_raw)
    
#     # Calculate exclusion radius in pixels
#     if 0 < Npsf < 1:
#         print("Warning: Npsf < 1, adding 1 pixel to exclusion radius.")
#         # exclusion_radius = np.ceil(Npsf * (psf_fwhm)+1)
#         exclusion_radius = int(np.ceil(Npsf * (psf_fwhm))+1)
#         # exclusion_radius = ((Npsf * (psf_fwhm))+1)
#     else:
#         print("Info: Npsf >= 1, using standard calculation for exclusion radius.")
#         exclusion_radius = np.ceil(Npsf * psf_fwhm)
#         # exclusion_radius = ((Npsf * (psf_fwhm)))
    
#     if logger is not None:
#         logger.info(f" ==> PSF FWHM = {psf_fwhm:.2f} pixels")
#         logger.info(f" ==> Excluding central region with radius = {exclusion_radius:.2f} pixels ({Npsf} x FWHM)")
    
#     # Create coordinate grid centered on image center
#     M, N = image_shape
#     y_center, x_center = M / 2.0, N / 2.0
#     y, x = np.ogrid[0:M, 0:N]
    
#     # Calculate distance from center
#     r = np.sqrt((x - x_center)**2 + (y - y_center)**2)
    
#     # Create mask: True where we WANT to fit (outside exclusion radius)
#     psf_mask = r > exclusion_radius
    
#     # Convert to appropriate format
#     if convolution_mode == 'GPU':
#         psf_mask = jnp.array(psf_mask)

#     # plt.figure(figsize=(6,6))
#     # plt.imshow(psf_mask, cmap='gray')
#     # plt.title('PSF Exclusion Mask')
#     # plt.colorbar(label='Fit Inclusion (True=Fit, False=Exclude)')
#     # plt.show()
    
#     return psf_mask

def psf_exclusion_radius_px(psf_fwhm, Npsf):
    """
    Radius in pixels of the PSF exclusion hole, from `Npsf` resolution elements.

    `Npsf` is the hole **diameter** in PSF FWHM, so the radius is
    `Npsf * FWHM / 2` and `Npsf=1` removes exactly one resolution element across.

    This used to be `Npsf * FWHM` used as a *radius*, which made `Npsf=1` cut a
    hole two FWHM wide -- four times the area the name implies, and enough to
    remove the entire core that constrains `n` and `In` for a compact component.
    There was also a discontinuity: values below 1 had a `+1 px` term the others
    did not, so the hole *shrank* as Npsf crossed 1 (0.99 -> 5.40 px, 1.00 -> 4.44
    px for a 4.44 px FWHM).

    Everything that needs this number calls here, so the fit mask and the
    diagnostic plots cannot drift apart again -- previously they computed it three
    different ways.
    """
    if Npsf is None or Npsf <= 0 or psf_fwhm is None:
        return 0.0
    return float(Npsf) * float(psf_fwhm) / 2.0


def _create_psf_exclusion_mask(PSF_DATA_raw, Npsf, image_shape,
                               centres,
                               convolution_mode, logger=None):
    """
    Mask out a PSF-sized hole around each component centre.

    Parameters
    ----------
    PSF_DATA_raw : array or None
        PSF data in raw numpy format
    Npsf : float
        Hole diameter in PSF FWHM (0 = no exclusion). See
        `psf_exclusion_radius_px`.
    image_shape : tuple
        Shape of the image (M, N)
    centres : sequence of (x, y)
        One centre per model component. A single `(x, y)` pair is accepted for
        convenience. Previously this was a single point -- the global brightest
        pixel -- so in a multi-component fit every component except the brightest
        kept its PSF-dominated core in the fit while the brightest lost it.
    convolution_mode : str
        'GPU' or 'CPU'
    logger : logger, optional
        Logger instance

    Returns
    -------
    psf_mask : array or None
        Boolean mask (True = include in fit, False = exclude)
        None if Npsf <= 0 or no PSF provided
    """
    if Npsf is None or Npsf <= 0 or PSF_DATA_raw is None:
        return None

    psf_fwhm = psf_params(PSF_DATA_raw)
    size = psf_exclusion_radius_px(psf_fwhm, Npsf)
    if size <= 0:
        return None

    # tolerate a bare (x, y)
    centres = np.atleast_2d(np.asarray(centres, dtype=float))

    if logger is not None:
        logger.info(f" ==> PSF FWHM = {psf_fwhm:.2f} pixels")
        logger.info(f" ==> Excluding a hole of radius {size:.2f} px "
                    f"({Npsf} x FWHM across) at {len(centres)} component "
                    f"centre(s).")

    M, N = image_shape
    x, y = np.meshgrid(np.arange(N), np.arange(M))

    psf_mask = np.ones((M, N), dtype=bool)
    for x0, y0 in centres:
        rr = np.sqrt((x - x0) ** 2 + (y - y0) ** 2)
        psf_mask &= (rr >= size)

    # Convert to appropriate format
    if convolution_mode == 'GPU':
        psf_mask = jnp.array(psf_mask)

    return psf_mask


# def _prepare_mask(mask_region, convolution_mode, logger=None):
#     """
#     Prepare mask for constrained fitting.

#     Parameters
#     ----------
#     mask_region : array or None
#         Mask array
#     convolution_mode : str
#         'GPU' or 'CPU'
#     logger : logger, optional
#         Logger instance

#     Returns
#     -------
#     mask_for_fit : array or None
#         Mask in appropriate format
#     """
#     if mask_region is not None:
#         """

#         """
#         if logger is not None:
#             logger.debug(f" ==> Using provided mask region to constrain fit. ")
#             logger.warning(f" !!++==> Fitting with a mask is faster, but experimental!! \n"
#                            f"         Use with caution.")
#         # data_2D = data_2D * mask_region
#         if convolution_mode == 'GPU':
#             mask_for_fit = jnp.array(mask_region)
#         else:
#             mask_for_fit = mask_region
#     else:
#         mask_for_fit = None

#     return mask_for_fit

def _prepare_mask(mask_region, PSF_DATA_raw, Npsf, image_shape,
                  centres,
                  convolution_mode, logger=None):
    """
    Prepare mask for constrained fitting, including PSF-based exclusion.

    Parameters
    ----------
    mask_region : array or None
        User-provided mask array
    PSF_DATA_raw : array or None
        PSF data for calculating exclusion region
    Npsf : float
        Hole diameter in PSF FWHM (see `psf_exclusion_radius_px`)
    image_shape : tuple
        Shape of the image
    centres : sequence of (x, y)
        One centre per model component
    convolution_mode : str
        'GPU' or 'CPU'
    logger : logger, optional
        Logger instance

    Returns
    -------
    mask_for_fit : array or None
        Combined mask in appropriate format (True = fit, False = exclude)
    """
    # Start with user-provided mask if available
    if mask_region is not None:
        if logger is not None:
            logger.debug(f" ==> Using provided mask region to constrain fit.")
            logger.warning(f" !!++==> Fitting with a mask is faster, but experimental!!\n"
                           f"         Use with caution.")
        
        # Convert to appropriate format
        if convolution_mode == 'GPU':
            mask_for_fit = jnp.array(mask_region, dtype=bool)
        else:
            mask_for_fit = np.array(mask_region, dtype=bool)
    else:
        mask_for_fit = None
    
    # Create PSF-based exclusion mask if requested
    psf_exclusion_mask = _create_psf_exclusion_mask(
        PSF_DATA_raw, Npsf, image_shape,
        centres,
        convolution_mode, logger
    )
    
    # Combine masks if both exist
    if psf_exclusion_mask is not None:
        if mask_for_fit is not None:
            # Both masks exist: combine them (AND operation)
            # Only fit where both masks are True
            if convolution_mode == 'GPU':
                mask_for_fit = jnp.logical_and(mask_for_fit, psf_exclusion_mask)
            else:
                mask_for_fit = np.logical_and(mask_for_fit, psf_exclusion_mask)
        else:
            # Only PSF mask exists
            mask_for_fit = psf_exclusion_mask
    
    return mask_for_fit


def _prepare_background_and_residual(residualname, residualdata_2D_, which_residual,
                                     self_bkg, bkg_map, imagename, data_2D,
                                     PSF_CONV, PSF_DATA, PSF_DATA_raw, convolution_mode,
                                     is_bkg_map_conv, logger=None, sky_mode=None):
    """
    Prepare background and residual arrays for fitting.

    This function handles the complex logic of preparing background maps,
    residual images, and RMS estimates based on various input configurations
    optimized for different observation types (radio vs optical).

    Background/Residual Strategy by Observation Type
    -------------------------------------------------

    RADIO IMAGES (residualname provided, which_residual != 'user'):
        - Uses residual map from interferometric cleaning process
        - 'shuffled' mode: Shuffles residual to create noise realization
          This breaks spatial correlations while preserving noise statistics
        - 'natural' mode: Uses residual as-is without shuffling
        - Residual added to model via optimized scaling factor s_a
        - Ensures flux conservation relative to cleaned image noise
        - FlatSky_level computed from residual for initial scaling

    OPTICAL IMAGES (which_residual = 'user'):
        - User provides custom background map via bkg_map parameter
        - Could be photometric error map from pipeline
        - Could be background estimate from source extraction
        - Requires explicit bkg_map parameter (raises ValueError if None)
        - Allows pre-convolved error maps (set is_bkg_map_conv=True)

    SELF-BACKGROUND MODE (self_bkg = True, no residual provided):
        - Estimates background directly from input image using SEP
        - Shuffles background estimate to create noise realization
        - Useful when no external background/residual available
        - Appropriate for well-behaved optical images

    FLAT SKY MODE (default fallback, no inputs):
        - Estimates single RMS value from entire image
        - Uses MAD (Median Absolute Deviation) estimator for robustness
        - Assumes uniform noise across image
        - Simplest case, background is scalar not array

    PSF Convolution of Background/Error Maps
    -----------------------------------------
    The is_bkg_map_conv flag controls convolution behavior:

    - If is_bkg_map_conv = False AND PSF provided:
      Background/error map is convolved with PSF to match model space
      Critical for radio residuals that need beam matching
      Optical error maps from unconvolved data need this

    - If is_bkg_map_conv = True:
      Background/error already convolved, use as-is
      Common for optical pipelines that provide PSF-matched errors
      Avoids double-convolution artifacts

    - If no PSF (PSF_CONV = False):
      No convolution regardless of is_bkg_map_conv setting
      Background/error used directly

    Implementation Notes
    --------------------
    - background_dec stores the unconvolved version for deconvolved model output
    - background stores the (possibly convolved) version for fitting
    - Both are needed to produce conv/deconv model components correctly
    - Scalar backgrounds (FlatSky_level case) are not convolved

    Parameters
    ----------
    residualname : str or None
        Path to residual image (radio: cleaning residual, optical: error map)
    residualdata_2D_ : array or None
        Pre-loaded residual data (avoids re-reading file)
    which_residual : str
        Residual processing mode: 'shuffled', 'natural', or 'user'
    self_bkg : bool
        Whether to use self-background estimation from image
    bkg_map : array or None
        User-provided background map (for 'user' mode)
    imagename : str
        Path to main image (for self-background estimation)
    data_2D : array
        Main image data (for flat sky MAD estimation)
    PSF_CONV : bool
        Whether PSF convolution is enabled
    PSF_DATA : array or None
        PSF data in computation format (GPU or CPU)
    PSF_DATA_raw : array or None
        PSF data in raw numpy format (for CPU convolution)
    convolution_mode : str
        'GPU' or 'CPU'
    is_bkg_map_conv : bool
        Whether the background map is already PSF-convolved
        (True prevents double-convolution)
    logger : logger, optional
        Logger instance

    Returns
    -------
    background : array or scalar
        Background array for fitting (possibly PSF-convolved)
        Array for map-based backgrounds, scalar for flat sky
    background_dec : array or scalar
        Deconvolved/unconvolved background (always unconvolved version)
        Used for generating deconvolved model components
    residual_2D : array
        Residual data array (loaded from file or background fallback)
    FlatSky_level : float or None
        Flat sky RMS estimate (None if using map-based background)
        Used for scaling in some model configurations
    """
    FlatSky_level = None

    if sky_mode == 'flat':
        # A UNIT pedestal, so that FlatSky(1.0, s_a) == s_a and the fitted `s_a`
        # IS the sky level in image units -- the IMFIT/GALFIT convention, where
        # the sky is a data-unit constant rather than a multiple of something.
        #
        # The legacy fallthrough further down instead uses mad_std(data_2D) as
        # the unit, which makes s_a a multiple of the image NOISE and leaves its
        # sign tied to the sign of whatever produced the background. That
        # fallthrough is untouched, so callers that do not ask for a mode keep
        # their numbers.
        #
        # A constant is unchanged by a normalised PSF, so it is returned already
        # "convolved" and the convolution below is skipped for it.
        if logger is not None:
            logger.debug(" ==> sky_mode='flat': fitting a flat pedestal; "
                         "s_a is the sky level in image units.")
        background = 1.0
        background_dec = 1.0
        FlatSky_level = 1.0
        if residualdata_2D_ is not None:
            residual_2D = residualdata_2D_
        else:
            try:
                residual_2D = _ensure_2D_image(pf.getdata(residualname),
                                               name=f'residual ({residualname})',
                                               logger=logger)
            except Exception:
                residual_2D = np.zeros_like(np.asarray(data_2D, dtype=float))
        return background, background_dec, residual_2D, FlatSky_level

    if residualname is not None and which_residual != 'user':
        """
        This is important for radio image fitting.

        It uses the shuffled version of the residual cleaned image
        originated from the interferometric deconvolution.

        This ensures that the best model created here will be on top
        of that rms noise so that flux conservation is maximized.

        However, this residual is not added as model + shuffled_residual
        only, but instead by a multiplication factor,
        e.g. model + const* shuffled_residual, and const will be minimized
        as well during the fitting (here, called `s_a`).
        """
        if residualdata_2D_ is not None:
            residual_2D = residualdata_2D_
        else:
            residual_2D = _ensure_2D_image(pf.getdata(residualname),
                                           name=f'residual ({residualname})',
                                           logger=logger)

        if which_residual == 'shuffled':
            if logger is not None:
                logger.debug(f" ==> Using clean shuffled background for optmization... ")
            residual_2D_to_use = shuffle_2D(residual_2D)
        elif which_residual == 'natural':
            if logger is not None:
                logger.debug(f" ==> Using clean background for optmization... ")
            """            
            if psf_name is not None:
                if logger is not None:
                    logger.debug(f" ====> Deconvolving residual map... ")
                residual_2D_to_use, _ = deconvolve_fft(residual_2D,
                                                            PSF_DATA_raw/PSF_DATA_raw.sum())
            else:
                residual_2D_to_use = residual_2D
            """
            residual_2D_to_use = residual_2D
        else:
            residual_2D_to_use = residual_2D

        FlatSky_level = mad_std(residual_2D_to_use)
        #         background = residual_2D #residual_2D_to_use
        if convolution_mode == 'GPU':
            background = jnp.array(residual_2D_to_use)
        else:
            background = residual_2D_to_use

    else:
        if which_residual == 'user':
            if bkg_map is None:
                print('--==>> A rms map/background mode was selected (user)')
                print('       but no rms/background map was provided.')
                print('       Please, provide a rms/background map.')
                print('||==>> Stopping code now.')
                raise ValueError("bkg_map should not be None when which_residual='user'")
            else:
                if logger is not None:
                    logger.debug(f" ==> Using provided Background map. ")
                background_map = bkg_map
                background = background_map.copy()
        else:
            if self_bkg == True:
                if logger is not None:
                    logger.warning(f" ==> No residual/background provided. Using image bkg map... ")
                background_map = sep_background(imagename)
                background = shuffle_2D(background_map.back())
            else:
                if logger is not None:
                    logger.warning(f" ==> Using only flat sky for rms bkg.")
                FlatSky_level = mad_std(data_2D)
                background = FlatSky_level


    # One hand-off to the GPU covering every branch above. The clean/natural
    # branch already did `jnp.array` on its own; the user-supplied map, the
    # self_bkg map and the flat-sky scalar did not, so a caller's array reached
    # the jitted FlatSky (and _fftconvolve_jax) as raw numpy. jax rejects any
    # array whose dtype is not in native byte order, and anything read straight
    # out of a FITS file is big-endian, so `bkg_map=<array from a FITS>` failed
    # with "Error interpreting argument to FlatSky as an abstract array".
    # Placed before `background_dec` is copied off below, so the deconvolved
    # copy -- which goes into the same jitted FlatSky at `_compute_sky_totals`
    # -- is converted too.
    if convolution_mode == 'GPU' and np.ndim(background) > 0:
        background = jnp.array(background)

    # Load residual_2D if not already loaded
    # This ensures residual_2D is always defined for downstream use
    if residualdata_2D_ is not None:
        residual_2D = residualdata_2D_
    else:
        try:
            residual_2D = pf.getdata(residualname)
        except:
            # Fallback: use background as residual if file doesn't exist
            residual_2D = background
        else:
            residual_2D = _ensure_2D_image(residual_2D,
                                           name=f'residual ({residualname})',
                                           logger=logger)

    # Create deconvolved background copy BEFORE any convolution
    # This preserves the unconvolved version for deconvolved model components
    #
    # The scalar test is `np.ndim(...) == 0`, not isinstance(int, float):
    # mad_std on a float32 image returns np.float32, which is NOT a subclass of
    # Python float (only np.float64 is), so the isinstance test let a 0-d value
    # through to fftconvolve and it died on "in1 and in2 should have the same
    # dimensionality". That path was unreachable before sky_mode existed.
    if np.ndim(background) == 0:
        # FlatSky_level case - background is a scalar, not an array
        background_dec = background
    else:
        # Array case - create copy before potential convolution
        background_dec = background.copy()

    # Convolve background with PSF if needed
    # Only convolve if: (1) PSF exists, (2) map not pre-convolved, (3) background is array
    if is_bkg_map_conv is False and PSF_CONV:
        if logger is not None:
            logger.debug(f" ==> RMS map is not convolved, convolving with PSF now.")

        if np.ndim(background) > 0:
            # Only convolve if background is an array, not a scalar
            if convolution_mode == 'GPU':
                background = _fftconvolve_jax(background, PSF_DATA)
            elif convolution_mode == 'CPU':
                background = scipy.signal.fftconvolve(background, PSF_DATA_raw, 'same')

    return background, background_dec, residual_2D, FlatSky_level


def _setup_coordinate_grid(data_2D, convolution_mode):
    """
    Create coordinate meshgrid for model evaluation.

    Parameters
    ----------
    data_2D : array
        Image data (for shape)
    convolution_mode : str
        'GPU' or 'CPU'

    Returns
    -------
    xy : tuple of arrays
        Coordinate meshgrid
    size : tuple
        Image shape
    """
    size = data_2D.shape
    if convolution_mode == 'GPU':
        x, y = jnp.meshgrid(jnp.arange((size[1])), jnp.arange((size[0])))
        xy = jnp.stack([x, y], axis=0)
    else:
        xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))

    return xy, size


# ============================================================================
# HELPER FUNCTIONS - Model Construction
# ============================================================================

def _create_minimizer_functions(nfunctions, xy, data_2D, data_2D_gpu, background,
                                weights_map,
                                residual_2D, PSF_DATA, PSF_DATA_raw, mask_for_fit,
                                convolution_mode, add_background_to_model=True):
    """
    Create the residual functions for minimization.

    Parameters
    ----------
    nfunctions : int
        Number of Sersic components
    xy : tuple of arrays
        Coordinate meshgrid
    data_2D : array
        Image data (CPU)
    data_2D_gpu : array or None
        Image data (GPU)
    background : array
        Background map
    residual_2D : array
        Residual data
    PSF_DATA : array
        PSF data in computation format
    PSF_DATA_raw : array
        PSF data in raw format
    mask_for_fit : array or None
        Fitting mask
    convolution_mode : str
        'GPU' or 'CPU'

    Returns
    -------
    min_residual_2D : callable or None
        CPU residual function
    min_residual_2D_GPU : callable or None
        GPU residual function
    build_model : callable or None
        Model building function (for GPU)
    func : callable or None
        Parameter splitting function (for GPU)
    """

    # CPU version
    def min_residual_2D(params):
        dict_model = {}
        model = 0
        for i in range(1, nfunctions + 1):
            model = model + sersic2D(xy, params['f' + str(i) + '_x0'],
                                     params['f' + str(i) + '_y0'],
                                     params['f' + str(i) + '_PA'],
                                     params['f' + str(i) + '_ell'],
                                     params['f' + str(i) + '_n'],
                                     params['f' + str(i) + '_In'],
                                     params['f' + str(i) + '_Rn'],
                                     params['f' + str(i) + '_cg'],
                                    #  params['f' + str(i) + '_Rtrunc'],
                                    #  params['f' + str(i) + '_delta_r']
                                     )
        # print(model.shape)
        # model = model + FlatSky_cpu(FlatSky_level, params['s_a'])*background
        # model = model + FlatSky_cpu(background, params['s_a'])
        MODEL_2D_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same')
        if add_background_to_model:
            MODEL_2D_conv = MODEL_2D_conv + FlatSky_cpu(background, params['s_a'])
        # residual = (data_2D - MODEL_2D_conv)*background
        residual = (data_2D - MODEL_2D_conv)*(weights_map)
        return np.ravel(residual)

    # GPU version
    try:
        # @partial(jit, static_argnums=1)
        @jit
        def func(x):
            return jnp.split(x, nfunctions)

        @jit
        def build_model(xy, param_matrix):
            model = 0
            for model_params in param_matrix:
                model = model + sersic2D_GPU(xy, model_params[0],
                                             model_params[1],
                                             model_params[2],
                                             model_params[3],
                                             model_params[4],
                                             model_params[5],
                                             model_params[6],
                                             model_params[7],
                                            #  model_params[8],
                                            #  model_params[9]
                                             )
            return model

    except Exception as e:
        print("JIT compilation failed, running without JIT.")
        def func(x):
            return jnp.split(x, nfunctions)

        def build_model(xy, param_matrix):
            model = 0
            for model_params in param_matrix:
                model = model + sersic2D_GPU(xy, model_params[0],
                                             model_params[1],
                                             model_params[2],
                                             model_params[3],
                                             model_params[4],
                                             model_params[5],
                                             model_params[6],
                                             model_params[7],
                                            #  model_params[8],
                                            #  model_params[9]
                                             )
            return model

    def min_residual_2D_GPU(params):
        model = 0
        for i in range(1, nfunctions + 1):
            model = model + sersic2D_GPU(xy,
                                         params['f' + str(i) + '_x0'].value,
                                         params['f' + str(i) + '_y0'].value,
                                         params['f' + str(i) + '_PA'].value,
                                         params['f' + str(i) + '_ell'].value,
                                         params['f' + str(i) + '_n'].value,
                                         params['f' + str(i) + '_In'].value,
                                         params['f' + str(i) + '_Rn'].value,
                                         params['f' + str(i) + '_cg'].value,
                                        #  params['f' + str(i) + '_Rtrunc'].value,
                                        #  params['f' + str(i) + '_delta_r'].value
                                         )

        # # param_matrix = extract_params(params)
        # param_matrix = func(jnp.array(list(params.valuesdict().values()))[:-1])
        # model = build_model(xy,param_matrix)

        # The one background convention: the data is never modified, the sky
        # rides on the model as s_a*B with s_a fitted. `add_background_to_model`
        # is False only when fit_background=False, i.e. no sky term at all.
        # See do_fit2D.
        MODEL_2D_conv = _fftconvolve_jax(model, PSF_DATA)
        if add_background_to_model:
            MODEL_2D_conv = MODEL_2D_conv + FlatSky(background, params['s_a'].value)

        # MODEL_2D_conv = _fftconvolve_jax(model+
        #                                  FlatSky(background,params['s_a'].value),
        #                                  PSF_DATA)
        # residual = ((data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])/
        #             (1000*(abs(residual_2D[mask_for_fit])+1.0e-6)))
        # residual = ((data_2D_gpu - MODEL_2D_conv)*background)[mask_for_fit]
        # weights = 1/((background[mask_for_fit])/data_2D_gpu[mask_for_fit])
        # weightned_residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit]) * jnp.sqrt(weights)
        # return np.asarray(weightned_residual).copy()

        # residual = ((jnp.log10(data_2D_gpu+0.1) - jnp.log10(MODEL_2D_conv+0.1))*(weights_map))[mask_for_fit]
        # residual = ((jnp.log10(data_2D_gpu - MODEL_2D_conv + 1))*(weights_map))[mask_for_fit]
        residual = ((data_2D_gpu-MODEL_2D_conv)*(weights_map))[mask_for_fit]
        # return np.asarray(jnp.nansum(residual**2)).copy()
        return np.asarray(residual).copy()

    if convolution_mode == 'CPU':
        return min_residual_2D, None, build_model, func
    else:
        return None, min_residual_2D_GPU, build_model, func


def _setup_jax_convolution(convolution_mode):
    """
    Setup JAX convolution function if in GPU mode.

    Parameters
    ----------
    convolution_mode : str
        'GPU' or 'CPU'

    Returns
    -------
    jax_convolve : callable or None
        JIT-compiled convolution function
    """
    if convolution_mode == 'GPU':
        @jit
        def convolve_on_gpu(image, psf):
            """
            This was before jax.scipy implementing fftconvolve.
            It provides the same result, at the same speed.

            This function also accepts PSFs with a different shape of the image.

            """
            # Calculate the new padded shape
            padded_shape = (image.shape[0] + psf.shape[0] - 1,
                            image.shape[1] + psf.shape[1] - 1)

            # Pad both image and psf to the new shape
            pad_shape = [(0, ts - s) for s, ts in zip(image.shape, padded_shape)]
            image_padded = jnp.pad(image, pad_shape, mode='constant')
            pad_shape = [(0, ts - s) for s, ts in zip(psf.shape, padded_shape)]
            psf_padded = jnp.pad(psf, pad_shape, mode='constant')
            # psf_padded = pad_for_convolution(psf, padded_shape)
            image_fft = jnp.fft.fft2(image_padded)
            psf_fft = jnp.fft.fft2(psf_padded)

            conv_fft = image_fft * psf_fft

            # Get the real part of the inverse FFT and crop to the original image size
            result_full = jnp.real(jnp.fft.ifft2(conv_fft))
            return result_full[psf.shape[0] // 2:image.shape[0] + psf.shape[0] // 2,
            psf.shape[1] // 2:image.shape[1] + psf.shape[1] // 2]

        jax_convolve = jit(convolve_on_gpu)
        return jax_convolve
    return None


# ============================================================================
# HELPER FUNCTIONS - Optimization
# ============================================================================

def _perform_optimization(mini, params, method1, method2, convolution_mode,
                          contrain_nelder, workers, max_nfev, tr_solver,
                          regularize, x_scale, f_scale, ftol, xtol, gtol,
                          verbose, loss, maxiter, maxfev, xatol, fatol,
                          return_all, disp, de_options, parameters_mini_init):
    """
    Execute the two-stage optimization process.

    Parameters
    ----------
    mini : lmfit.Minimizer
        Minimizer object
    params : lmfit.Parameters
        Initial parameters
    method1 : str
        First optimization method
    method2 : str
        Second optimization method
    convolution_mode : str
        'GPU' or 'CPU'
    contrain_nelder : bool
        Whether to constrain Nelder-Mead parameters
    workers : int
        Number of parallel workers
    max_nfev : int
        Maximum function evaluations
    tr_solver : str
        Trust region solver
    regularize : bool
        Whether to regularize
    x_scale : str or array
        Parameter scaling
    f_scale : float
        Function scaling
    ftol : float
        Function tolerance
    xtol : float
        Parameter tolerance
    gtol : float
        Gradient tolerance
    verbose : int
        Verbosity level
    loss : str
        Loss function
    maxiter : int
        Maximum iterations for Nelder-Mead
    maxfev : int
        Maximum function evaluations for Nelder-Mead
    xatol : float
        Absolute parameter tolerance for Nelder-Mead
    fatol : float
        Absolute function tolerance for Nelder-Mead
    return_all : bool
        Return all Nelder-Mead results
    disp : bool
        Display convergence messages
    de_options : dict
        Options for differential evolution
    parameters_mini_init : lmfit.Parameters or None
        Initial parameters from previous run

    Returns
    -------
    result : lmfit.MinimizerResult
        Final optimization result
    result_1 : lmfit.MinimizerResult
        First stage optimization result
    result_extra : None
        Placeholder for additional results
    """
    # initial minimization.
    if verbose > 0:
        print(' >> Using', method1, ' solver for first optimisation run... ')
    # take parameters from previous run, and re-optimize them.
    #     method2 = 'ampgo'#'least_squares'
    #     method2 = 'least_squares'
    result_extra = None

    if method1 == 'nelder':
        # very robust, but takes time....
        #         print(' >> Using', method1,' solver for first optimisation run... ')
        result_1 = mini.minimize(method='nelder',
                                 #                                  xatol = 1e-12, fatol = 1e-12, disp=True,
                                 #                                  adaptive = True,max_nfev = 30000,
                                 options={'maxiter': maxiter, 'maxfev': maxfev,
                                          'xatol': xatol, 'fatol': fatol,
                                          'return_all': return_all,
                                          'disp': disp}
                                 )

    elif method1 == 'least_squares':
        # faster, but usually not good for first run.
        # if results_previous_run is not None:
        if verbose > 0:
            print(' >> Using', tr_solver, 'for tr solver, with regularize set to', regularize,
                  ' Loss is', loss, '.')

        if parameters_mini_init is not None:
            if verbose > 0:
                print(f'  ++==>> Using initial mini parameters from a previous run.')
            # try:
            result_1 = mini.minimize(method='least_squares',
                                     params=parameters_mini_init,
                                     max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
                                     tr_solver=tr_solver,
                                     tr_options={'regularize': regularize,
                                                 },
                                     ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
                                     loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
        else:
            result_1 = mini.minimize(method='least_squares',
                                     max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
                                     tr_solver=tr_solver,
                                     tr_options={'regularize': regularize,
                                                 #                                              'min_delta':1e-14, 'eta':0.05,
                                                 #                                              'xtol':1e-14, 'gtol':1e-14,
                                                 #                                              'ftol':1e-14
                                                 },
                                     ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
                                     loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)

    elif method1 == 'differential_evolution':
        # de is giving some issues, I do not know why.
        result_1 = mini.minimize(method='differential_evolution',
                                 options={'disp': True, 'workers': workers,
                                          'max_nfev': max_nfev, 'vectorized': True,
                                          'strategy': 'randtobest1bin',
                                          'mutation': (0.5, 1.5),
                                          'recombination': [0.2, 0.9],
                                          'init': 'random', 'tol': 0.00001,
                                          'updating': 'deferred',
                                          'popsize': 600})
        # result_1 = mini.minimize(method='differential_evolution', popsize=600,
        #                          disp=True,  # init = 'random',
        #                          # mutation=(0.5, 1.5), recombination=[0.2, 0.9],
        #                          max_nfev=20000,
        #                          workers=1, updating='deferred', vectorized=True)
    else:
        raise ValueError(f"Unknown method1: {method1}")

    if verbose > 0:
        print(' >> Using', method2, ' solver for second optimisation run... ')

    second_run_params = result_1.params
    if (contrain_nelder == True) and (method2 == 'nelder'):
        """
        It seems that least_squares is ignoring the best-parameters provided by
        Nelder-mead, which means that it is lookig the parameter space far away
        from the optimised Nelder-Mead ones.

        So, with this condition, we force a much smaller searching region, but
        it assumes that Nelder opt was good (which is not always true).

        YOU MUST CHECK YOUR RESULTS!!!!

        """
        print('Constraining Nelder-Mead Parameters for method', method2)
        params_constrained = constrain_nelder_mead_params(result_1.params,
                                                          max_factor=1.03,
                                                          min_factor=0.97)
        # UPDATE THE SECOND RUN PARAMETERS TO BE THE CONSTRAINED ONES.
        second_run_params = params_constrained

    if method2 == 'nelder':
        result = mini.minimize(method='nelder', params=second_run_params,
                               options={'maxiter': maxiter, 'maxfev': maxfev,
                                        'xatol': xatol, 'fatol': fatol,
                                        'disp': disp})

    elif method2 == 'ampgo':
        # ampgo is not workin well/ takes so long ???
        result = mini.minimize(method='ampgo', params=second_run_params,
                               maxfunevals=10000, totaliter=30, disp=True,
                               maxiter=5, glbtol=1e-8)

    elif method2 == 'least_squares':
        # faster, usually converges and provide errors.
        # Very robust if used in second opt from first opt parameters.
        result = mini.minimize(method='least_squares',
                               params=second_run_params,
                               max_nfev=max_nfev,
                               tr_solver=tr_solver,
                               tr_options={'regularize': regularize,
                                           #                                            'min_delta': 1e-14, 'eta': 0.05,
                                           #                                            'xtol': 1e-14, 'gtol': 1e-14,
                                           #                                            'ftol': 1e-14
                                           },
                               x_scale=x_scale, f_scale=f_scale,
                               ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
                               loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)

    elif method2 == 'differential_evolution':
        # result = mini.minimize(method='differential_evolution',
        #                        params=second_run_params,
        #                        options={'maxiter': 30000, 'workers': -1,
        #                                 'tol': 0.001, 'vectorized': True,
        #                                 'strategy': 'randtobest1bin',
        #                                 'updating': 'deferred', 'disp': True,
        #                                 'seed': 1}
        #                        )
        result = mini.minimize(method='differential_evolution',
                               params=second_run_params,
                               options=de_options
                               )
    else:
        raise ValueError(f"Unknown method2: {method2}")

    return result, result_1, result_extra


# ============================================================================
# HELPER FUNCTIONS - Output Generation
# ============================================================================

def _generate_model_components(params, ncomponents, xy, size, background,
                               background_dec, PSF_CONV, PSF_DATA, PSF_DATA_raw,
                               convolution_mode, is_bkg_map_conv):
    """
    Generate individual and total model components.

    Parameters
    ----------
    params : lmfit.Parameters
        Optimized parameters
    ncomponents : int
        Number of Sersic components
    xy : tuple of arrays
        Coordinate meshgrid
    size : tuple
        Image shape
    background : array
        Background map (possibly convolved)
    background_dec : array
        Deconvolved background map
    PSF_CONV : bool
        Whether PSF convolution is enabled
    PSF_DATA : array
        PSF data in computation format
    PSF_DATA_raw : array
        PSF data in raw format
    convolution_mode : str
        'GPU' or 'CPU'
    is_bkg_map_conv : bool
        Whether RMS map is already convolved

    Returns
    -------
    model_dict : dict
        Dictionary containing all model arrays
    flat_sky_total : array
        Total flat sky background (convolved)
    flat_sky_total_dec : array
        Total flat sky background (deconvolved)
    bkg_comp_i : array
        Background component for individual models
    bkg_comp_i_dec : array
        Deconvolved background component for individual models
    """
    model_temp = Model(sersic2D_GPU)
    xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
    model = 0
    model_dict = {}

    if convolution_mode == 'GPU':
        flat_sky_total = np.asarray(FlatSky(background, params['s_a'].value))
        flat_sky_total_dec = np.asarray(FlatSky(background_dec, params['s_a'].value))
    elif convolution_mode == 'CPU':
        flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
        flat_sky_total_dec = FlatSky_cpu(background_dec, params['s_a'].value)
    else:
        # Fallback
        flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
        flat_sky_total_dec = FlatSky_cpu(background_dec, params['s_a'].value)

    # With no background map, `background` is the scalar mad_std of the image and
    # s_a*B is a flat pedestal -- a 0-d array. Everything downstream (the saved
    # conv_bkg/deconv_bkg FITS, the radial profiles) expects an image, and
    # pf.writeto rejects a 0-d array outright. Broadcast once, here.
    flat_sky_total = np.broadcast_to(np.asarray(flat_sky_total),
                                     (size[0], size[1])).copy()
    flat_sky_total_dec = np.broadcast_to(np.asarray(flat_sky_total_dec),
                                         (size[0], size[1])).copy()

    bkg_comp_i = flat_sky_total.copy()
    bkg_comp_i_dec = flat_sky_total_dec.copy()
    bkg_sign = +1.0

    for i in range(1, ncomponents + 1):
        model_temp = sersic2D_GPU(xy, params['f' + str(i) + '_x0'].value,
                                  params['f' + str(i) + '_y0'].value,
                                  params['f' + str(i) + '_PA'].value,
                                  params['f' + str(i) + '_ell'].value,
                                  params['f' + str(i) + '_n'].value,
                                  params['f' + str(i) + '_In'].value,
                                  params['f' + str(i) + '_Rn'].value,
                                  params['f' + str(i) + '_cg'].value,
                                #   params['f' + str(i) + '_Rtrunc'].value,
                                #   params['f' + str(i) + '_delta_r'].value
                                  )

        model = model + model_temp
        # to each individual component, add the bkg map.
        model_dict['model_c' + str(i)] = np.asarray(model_temp + bkg_sign * bkg_comp_i_dec)

        if PSF_CONV == True:
            if convolution_mode == 'GPU':
                # model_dict['model_c' + str(i) + '_conv'] = np.asarray(jax_convolve(model_temp, PSF_DATA)).copy()
                # model_dict['model_c' + str(i) + '_conv'] = (
                #     np.asarray(_fftconvolve_jax(model_temp,PSF_DATA).copy()+bkg_comp_i))
                # to each individual component, add the bkg map.
                # model_dict['model_c' + str(i) + '_conv'] = (
                #     np.asarray(_fftconvolve_jax(model_temp+bkg_comp_i,PSF_DATA).copy()))

                model_dict['model_c' + str(i) + '_conv'] = (
                                                               np.asarray(_fftconvolve_jax(model_temp,
                                                                                           PSF_DATA).copy())) + bkg_sign * bkg_comp_i

            elif convolution_mode == 'CPU':
                # model_dict['model_c' + str(i) + '_conv'] = (
                #         scipy.signal.fftconvolve(model_temp+bkg_comp_i, PSF_DATA_raw,'same'))
                model_dict['model_c' + str(i) + '_conv'] = (
                        scipy.signal.fftconvolve(model_temp, PSF_DATA_raw, 'same') + bkg_sign * bkg_comp_i
                )
            else:
                model_dict['model_c' + str(i) + '_conv'] = (
                        scipy.signal.fftconvolve(model_temp, PSF_DATA_raw, 'same') + bkg_sign * bkg_comp_i
                )
        else:
            model_dict['model_c' + str(i) + '_conv'] = model_temp + bkg_sign * bkg_comp_i

    #     model = model
    model_dict['model_total_dec'] = np.asarray(model + bkg_sign * flat_sky_total_dec)  # +FlatSky_cpu(background,
    # params['s_a'].value)

    if PSF_CONV == True:
        # model_dict['model_total_conv'] = scipy.signal.fftconvolve(model,
        #                                                           PSF_DATA_raw,
        #                                                           'same')  # + FlatSky(FlatSky_level, params['s_a'])
        if convolution_mode == 'GPU':
            # model_dict['model_total_conv'] = np.asarray(jax_convolve(model,
            #                                                          PSF_DATA)).copy()
            # model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + FlatSky_cpu(background, params['s_a'].value
            # model_conv = _fftconvolve_jax(model+flat_sky_total,PSF_DATA).copy()
            model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + bkg_sign * flat_sky_total
        elif convolution_mode == 'CPU':
            model_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + bkg_sign * flat_sky_total
        else:
            model_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + bkg_sign * flat_sky_total
        model_dict['model_total_conv'] = model_conv
    else:
        model_dict['model_total_conv'] = model + bkg_sign * flat_sky_total

    return model_dict, flat_sky_total, flat_sky_total_dec, bkg_comp_i, bkg_comp_i_dec


def _compute_residuals(data_2D, model_dict, flat_sky_total, is_bkg_map_conv,
                       PSF_CONV, PSF_DATA, convolution_mode):
    """
    Compute residual images and background components.

    Parameters
    ----------
    data_2D : array
        Original image data
    model_dict : dict
        Dictionary of model components
    flat_sky_total : array
        Total flat sky background (convolved)
    is_bkg_map_conv : bool
        Whether RMS map is already convolved
    PSF_CONV : bool
        Whether PSF convolution is enabled
    PSF_DATA : array
        PSF data
    convolution_mode : str
        'GPU' or 'CPU'

    Returns
    -------
    model_dict : dict
        Updated model dictionary with residuals
    """
    # model_dict['best_residual'] = data_2D - model_dict['model_total']
    # bkg_comp_total
    model_dict['model_total_conv'] = np.asarray(model_dict['model_total_conv'])

    if is_bkg_map_conv is False and PSF_CONV is True:
        if convolution_mode == 'GPU':
            model_dict['conv_bkg'] = np.asarray(_fftconvolve_jax(flat_sky_total, PSF_DATA).copy())
        else:
            model_dict['conv_bkg'] = np.asarray(flat_sky_total)
    else:
        model_dict['conv_bkg'] = np.asarray(flat_sky_total)

    # The residual is what the minimiser saw: data minus the FULL model, and the
    # full model includes s_a*B. The is_bkg_map_conv branch used to add
    # flat_sky_total back in, which left a whole background sitting in the saved
    # residual map and in the RESIDUAL curve of every diagnostic plot.
    model_dict['best_residual_conv'] = (np.asarray(data_2D)
                                        - model_dict['model_total_conv'] + model_dict['conv_bkg'])


    return model_dict


def _save_fitting_outputs(imagename, model_dict, ncomponents, special_name,
                          save_name_append, flat_sky_total_dec):
    """
    Save all FITS files and generate output file lists.

    Parameters
    ----------
    imagename : str
        Base image filename
    model_dict : dict
        Dictionary containing all model arrays
    ncomponents : int
        Number of Sersic components
    special_name : str
        Special identifier for output files
    save_name_append : str
        Additional string to append to filenames
    flat_sky_total_dec : array
        Deconvolved flat sky background

    Returns
    -------
    image_results_conv : list
        List of convolved output filenames
    image_results_deconv : list
        List of deconvolved output filenames
    bkg_images : list
        List of background image filenames
    total_image_results_conv : list
        List of total convolved model filenames
    total_image_results_deconv : list
        List of total deconvolved model filenames
    """
    image_results_conv = []
    image_results_deconv = []
    total_image_results_conv = []
    total_image_results_deconv = []
    bkg_images = []

    # Save individual components
    for i in range(1, ncomponents + 1):
        # Convolved component
        conv_filename = (imagename.replace('.fits', '') +
                         "_" + "model_component_" + str(i) +
                         special_name + save_name_append + '.fits')
        pf.writeto(conv_filename, model_dict['model_c' + str(i) + '_conv'],
                   overwrite=True)
        copy_header(imagename, conv_filename, conv_filename)
        image_results_conv.append(conv_filename)

        # Deconvolved component
        dec_filename = (imagename.replace('.fits', '') +
                        "_" + "dec_model_component_" + str(i) +
                        special_name + save_name_append + '.fits')
        pf.writeto(dec_filename, model_dict['model_c' + str(i)],
                   overwrite=True)
        copy_header(imagename, dec_filename, dec_filename)
        image_results_deconv.append(dec_filename)

    # Save total convolved model
    conv_model_filename = (imagename.replace('.fits', '') +
                           "_" + "conv_model" + special_name + save_name_append + '.fits')
    pf.writeto(conv_model_filename, model_dict['model_total_conv'], overwrite=True)
    copy_header(imagename, conv_model_filename, conv_model_filename)
    total_image_results_conv.append(conv_model_filename)
    image_results_conv.append(conv_model_filename)

    # Save total deconvolved model
    dec_model_filename = (imagename.replace('.fits', '') +
                          "_" + "dec_model" + special_name + save_name_append + '.fits')
    pf.writeto(dec_model_filename, model_dict['model_total_dec'], overwrite=True)
    copy_header(imagename, dec_model_filename, dec_model_filename)
    total_image_results_deconv.append(dec_model_filename)
    image_results_deconv.append(dec_model_filename)

    # Save residual
    residual_filename = (imagename.replace('.fits', '') +
                         "_" + "residual" + special_name + save_name_append + ".fits")
    pf.writeto(residual_filename, model_dict['best_residual_conv'], overwrite=True)
    copy_header(imagename, residual_filename, residual_filename)
    image_results_conv.append(residual_filename)

    # Save deconvolved background
    model_dict['deconv_bkg'] = np.asarray(flat_sky_total_dec)
    deconv_bkg_filename = (imagename.replace('.fits', '') +
                           "_" + "deconv_bkg" + special_name + save_name_append + '.fits')
    pf.writeto(deconv_bkg_filename, model_dict['deconv_bkg'], overwrite=True)
    copy_header(imagename, deconv_bkg_filename, deconv_bkg_filename)
    bkg_images.append(deconv_bkg_filename)

    # Save convolved background
    conv_bkg_filename = (imagename.replace('.fits', '') +
                         "_" + "conv_bkg" + special_name + save_name_append + '.fits')
    pf.writeto(conv_bkg_filename, model_dict['conv_bkg'], overwrite=True)
    copy_header(imagename, conv_bkg_filename, conv_bkg_filename)
    bkg_images.append(conv_bkg_filename)

    # # initial minimization.
    # method1 = 'differential_evolution'
    # print(' >> Using', method1, ' solver for first optimisation run... ')
    # # take parameters from previous run, and re-optimize them.
    # #     method2 = 'ampgo'#'least_squares'
    # method2 = 'least_squares'

    # # save mini results (full) to a pickle file.
    # with open(imagename.replace('.fits',
    #                             '_' + 'fit' +
    #                             special_name + save_name_append + '.pickle'),
    #           "wb") as f:
    #     pickle.dump(result, f)

    # with open(imagename.replace('.fits',
    #                             '_' + 'fit' +
    #                             special_name + save_name_append + '_modeldict.pickle'),
    #           "wb") as f:
    #     pickle.dump(model_dict, f)

    return (image_results_conv, image_results_deconv, bkg_images,
            total_image_results_conv, total_image_results_deconv)





def do_fit2D(imagename, params_values_init_IMFIT=None, ncomponents=None,
             init_constraints=None, data_2D_=None, residualdata_2D_=None,
             residualname=None, which_residual='shuffled', observation_type='radio',
             init_params=0.25, final_params=4.0, constrained=True,
             fix_n=None,
             fix_value_n=False,
             fix_max_value_n=False, 
             fix_min_value_n=False,
             fix_max_value_Rn=False, 
             fix_min_value_Rn=False,
             dr_fix=None,
             fix_x0_y0=False, psf_name=None, convolution_mode='GPU',
             convolve_cutout=False, cut_size=512, self_bkg=False,
             bkg_map=None, rms_map=None, 
             use_weights = False,
             sky_mode=None,
             fit_background=None,
             is_background_subtracted=False,
             sky_scale_bounds=None,
             background_mode=None,
             fit_background_scale=None,
             rms_convention='sigma',
             weight_mode='inverse_variance',
             scale_covar=True,
             is_bkg_map_conv=False,
             fix_geometry=None, force_circular=None, trunc=False,
             contrain_nelder=False, workers=6, mask_region=None,
             Npsf=0,  
             special_name='', method1='least_squares', method2='least_squares',
             reduce_fcn='neglogcauchy', loss="cauchy", tr_solver="exact", x_scale='jac',
             ftol=1e-10, xtol=1e-10, gtol=1e-10, verbose=0, max_nfev=200000,
             regularize=True, f_scale=1.0,
             maxiter=30000, maxfev=30000, xatol=1e-12,
             fatol=1e-12, return_all=True, disp=True,
             de_options=None, parameters_mini_init=None,
             save_name_append='', logger=None):
    """
    Perform a Robust and Fast Multi-Sersic Decomposition with GPU acceleration.
    tr_solver:

    Parameters
    ----------
    imagename: str
        Name of the image to be fitted.
    params_values_init_IMFIT: list
        Initial parameters values for the model.
    ncomponents: int
        Number of components to be fitted.
    init_constraints: dict
        Initial constraints for the model.
    data_2D_: 2D array
        Image to be fitted.
    residualdata_2D_: 2D array
        Residual image to be fitted.
    residualname: str
        Name of the residual image to be fitted.
    which_residual: str
        Which residual to be used for the fitting.
        Options: 'shuffled' or 'natural'.
    init_params: float
        Initial parameters for the model.
    final_params: float
        Final parameters for the model.
    constrained: bool
        If True, use initial constraints for the model.
    fix_n: bool
        If True, fix the Sersic index of the model.
    fix_value_n: float
        If True, fix the Sersic index of the model to this value.
    dr_fix: float
        If True, fix the centre position of the model.
    fix_x0_y0: bool
        If True, fix the centre position of the model.
    psf_name: str
        Name of the PSF image to be used for the convolution.
    convolution_mode: str
        If 'GPU', use GPU acceleration for the convolution.
    convolve_cutout: bool
        If True, convolve the image cutout with the PSF.
    cut_size: int
        Size of the cutout image.
    bkg_map: 2D array
        Background map to be used for the fitting.
    is_bkg_map_conv: bool
        If True, the background map is already convolved.
    self_bkg: bool
        If True, use the image background as the residual background.
    rms_map: 2D array
        RMS map to be used for the fitting.
    use_weights : bool
        Whether to use weights for fitting. If True, rms_map can be provided, otherwise it will be computed internally.
    sky_mode : str or None
        How the sky is modelled. The image is never modified in any of them; the
        sky always rides on the MODEL as `+ s_a * B`.

            'none'  no sky term at all. s_a is frozen at 0. Use on data that is
                    already sky-subtracted AND trusted.
            'flat'  fit a single flat pedestal. B is a unit constant, so `s_a`
                    is the sky level in image units and its sign is unambiguous:
                    positive means the image carries a pedestal, negative means
                    the model's outer wings overshoot or the data was
                    over-subtracted upstream.
            'map'   fit the supplied `bkg_map` (or the residual) scaled by s_a.
                    Same mechanics as 'flat', but B has a shape.

        None resolves from the superseded arguments, then from
        `is_background_subtracted`: subtracted -> 'none', raw -> 'flat'. Radio
        with nothing passed keeps its historical behaviour untouched.
    fit_background : bool or None
        The single background switch. The data handed to the minimiser is always
        the image as read; the background rides on the model as `s_a * bkg` with
        `s_a` a fitted parameter, so

            True   ->  residual = ( I - conv(M) - s_a*bkg ) * w
            False  ->  residual = ( I - conv(M) ) * w, with s_a frozen at 0

        None (default) resolves to True, or to whatever the superseded
        `background_mode` / `fit_background_scale` arguments imply.
    is_background_subtracted : bool
        Declares that the sky was already removed from this image by an upstream
        step. It does NOT change what data enters the fit -- it only sets the
        starting point for `s_a` (0 instead of 1), because the correct answer for
        a cleaned image is "no residual sky". Morphen cannot infer this from the
        file, so it has to be stated.
    sky_scale_bounds : tuple or None
        `(init, min, max)` for `s_a`, overriding the defaults. Useful to opt a
        radio fit into the new prior, which is otherwise left untouched.
    background_mode : str or None
        Superseded, still accepted. 'none' -> fit_background=False;
        'add' -> True; 'subtract' -> True with a warning, since it no longer
        pre-subtracts anything.
    fit_background_scale : bool or None
        Superseded and ignored: the amplitude is always fitted when
        `fit_background` is on. Passing True still switches the term on.
    fix_geometry: bool
        If True, fix the geometry of the model.
    contrain_nelder: bool
        If True, constrain the Nelder-Mead optimised parameters.
    workers: int
        Number of workers to be used for the fitting.
    mask_region: 2D array
        Mask to be used for the fitting.
    special_name: str
        Special name to be used for the output files.
    method1: str
        Method to be used for the fitting.
    method2: str
        Method to be used for the fitting.
    reduce_fcn: str

    loss: str

    tr_solver: str

    x_scale: str

    ftol: float

    xtol: float

    gtol: float

    verbose: int

    max_nfev: int

    regularize: bool
        If True, regularize the model.
    f_scale: float

    maxiter: int

    maxfev: int

    xatol: float

    fatol: float

    return_all: bool

    disp: bool

    de_options: dict

    save_name_append: str

    logger: logger

    Npsf : float, optional
        Number of PSF FWHM to exclude from central region during fitting.
        If Npsf=0 (default), entire image is fitted.
        If Npsf=2, central region with radius = 2*PSF_FWHM is excluded.
        Only applied when psf_name is provided.

    returns
    -------
    result: dict
        Dictionary containing the results of the fitting.


    """
    # ========================================================================
    # INITIALIZATION
    # ========================================================================

    # Check JAX availability and adjust convolution mode if needed
    try:
        from jax import jit
    except:
        convolution_mode = 'CPU'

    if de_options is None:
        de_options = {'disp': True, 'workers': 6,
                      'max_nfev': 20000, 'vectorized': True,
                      # 'strategy': 'randtobest1bin',
                      'mutation': (0.5, 1.5),
                      'recombination': [0.2, 0.9],
                      'init': 'random', 'tol': 0.00001,
                      'updating': 'deferred',
                      'popsize': 600}

    startTime = time.time()
    try:
        logger.info(f"Fitting image: {imagename}")
    except:
        pass

    # ========================================================================
    # DATA LOADING AND PREPARATION
    # ========================================================================

    # Every user-supplied array is normalised to 2D up front, so that a 4D
    # (NAXIS=4) input -- what `t_cutout_2D_radec` writes for a 4D parent image
    # -- cannot leak into the model grid. See `_ensure_2D_image`.
    data_2D_ = _ensure_2D_image(data_2D_, name='data_2D_', logger=logger)
    residualdata_2D_ = _ensure_2D_image(residualdata_2D_,
                                        name='residualdata_2D_', logger=logger)
    mask_region = _ensure_2D_image(mask_region, name='mask_region', logger=logger)
    bkg_map = _ensure_2D_image(bkg_map, name='bkg_map', logger=logger)
    rms_map = _ensure_2D_image(rms_map, name='rms_map', logger=logger)

    # Load image data
    _data_2D, _data_2D_gpu = _initialize_image_data(imagename, data_2D_,
                                                  convolution_mode, logger)

    # Load PSF data
    PSF_CONV, PSF_DATA, PSF_DATA_raw = _initialize_psf_data(psf_name,
                                                            convolution_mode, logger)

    # Centres for the optional PSF exclusion holes: one per model component, from
    # the source-extraction positions. Falls back to the brightest pixel when
    # there are no per-component positions (a free fit with no init_constraints),
    # which is what this used to do for *every* case -- leaving every component
    # but the brightest with its PSF-dominated core still in the fit.
    psf_centres = []
    if init_constraints is not None:
        for _i in range(1, int(init_constraints.get('ncomps', 0)) + 1):
            if f'c{_i}_x0c' in init_constraints:
                psf_centres.append((init_constraints[f'c{_i}_x0c'],
                                    init_constraints[f'c{_i}_y0c']))
    if not psf_centres:
        if mask_region is None:
            psf_centres = [nd.maximum_position(_data_2D)[::-1]]
        else:
            psf_centres = [nd.maximum_position(_data_2D * mask_region)[::-1]]
    # Prepare mask
    mask_for_fit = _prepare_mask(mask_region, PSF_DATA_raw, Npsf,
                                 _data_2D.shape,
                                 psf_centres,
                                 convolution_mode, logger)

    # ------------------------------------------------------------------
    # How the sky enters the fit. There is exactly ONE convention:
    #
    #     residual = ( I - conv(M) - s_a * B ) * w
    #
    # The data handed to the minimiser is ALWAYS the image as read. The sky
    # lives on the MODEL side, scaled by the fitted parameter `s_a`, and `I - B`
    # exists only in the diagnostic plots. `sky_mode` says what B is:
    #
    #     'none'   no term at all; s_a frozen at 0
    #     'flat'   B = 1, a unit pedestal, so s_a IS the sky in image units
    #     'map'    B = the supplied bkg_map, or the residual
    #
    # This mirrors IMFIT/GALFIT (subtract-and-no-term, or fit-a-sky-component,
    # never both) and pysersic's typed sky choice. Before this, the non-radio
    # path SUBTRACTED B from the data and ALSO added s_a*B to the model, so the
    # background was counted (1 + s_a) times with s_a in [-1e-3, +1.0] starting
    # at 0.99 -- between 1 and 2 times the sky, never less, never zero.
    #
    # --- resolution ------------------------------------------------------------
    # An explicit sky_mode wins. Otherwise the superseded arguments fold into
    # one, then `is_background_subtracted` supplies the default. They all stay
    # accepted so existing notebooks keep running.
    _bkg_available = (bkg_map is not None
                      or bool(self_bkg)
                      or residualname is not None
                      or residualdata_2D_ is not None)

    _bkg_arg_given = (sky_mode is not None
                      or fit_background is not None
                      or background_mode is not None
                      or fit_background_scale is not None
                      or sky_scale_bounds is not None)

    if sky_mode is not None and sky_mode not in ('none', 'flat', 'map'):
        raise ValueError(
            f"sky_mode must be 'none', 'flat', 'map' or None, got {sky_mode!r}.")

    if background_mode is not None:
        if background_mode not in ('subtract', 'none', 'add'):
            raise ValueError(
                f"background_mode must be 'subtract', 'none', 'add' or None, "
                f"got {background_mode!r}.")
        if background_mode in ('subtract', 'add') and not _bkg_available:
            raise ValueError(
                f"background_mode='{background_mode}' needs a background map; "
                f"none was provided (bkg_map is None).")
        if fit_background is None:
            fit_background = (background_mode != 'none')
        if background_mode == 'subtract' and logger is not None:
            logger.warning(
                " !!++==> background_mode='subtract' no longer subtracts the "
                "background from the data. The single convention is "
                "model + s_a*bkg fitted against the unmodified image; "
                "'subtract' now just means 'fit that amplitude, starting at 1'. "
                "Use sky_mode='map' instead.")

    if fit_background_scale is not None and logger is not None:
        logger.warning(
            " !!++==> fit_background_scale is superseded: the sky amplitude "
            "s_a is always fitted when the sky term is on. The argument is "
            "ignored.")
    if fit_background is None and fit_background_scale:
        fit_background = True

    if sky_mode is None:
        if fit_background is False:
            sky_mode = 'none'
        elif fit_background is True:
            # A map if one is available, a flat pedestal otherwise.
            sky_mode = 'map' if _bkg_available else 'flat'
        elif observation_type == 'radio':
            # Nothing asked for on the radio path: keep the historical
            # behaviour exactly, s_a prior included.
            sky_mode = 'map' if _bkg_available else 'flat'
        elif is_background_subtracted:
            sky_mode = 'none'
        else:
            sky_mode = 'flat' if not _bkg_available else 'map'

    if sky_mode == 'map' and not _bkg_available:
        if logger is not None:
            logger.warning(
                " !!++==> sky_mode='map' but no background map or residual was "
                "given; falling back to sky_mode='flat' (a fitted pedestal).")
        sky_mode = 'flat'

    # `sky_mode='flat'` is reachable without a bkg_map: it short-circuits the
    # which_residual='user' requirement rather than raising.
    _prep_sky_mode = 'flat' if (sky_mode == 'flat' and _bkg_arg_given) else None

    # Prepare background and residual
    background, background_dec, residual_2D, FlatSky_level = \
        _prepare_background_and_residual(
            residualname, residualdata_2D_, which_residual,
            self_bkg, bkg_map, imagename, _data_2D,
            PSF_CONV, PSF_DATA, PSF_DATA_raw, convolution_mode,
            is_bkg_map_conv, logger, sky_mode=_prep_sky_mode
        )
    _bkg_is_array = np.ndim(background) > 0
    _flat_pedestal = (_prep_sky_mode == 'flat')
    if _flat_pedestal:
        # A constant convolved with a normalised PSF is the same constant; the
        # only thing convolution adds is a taper at the frame edge, from the
        # zero padding. The MINIMISER used the true constant (the flat branch
        # returns it unconvolved), so convolving it again downstream would make
        # the saved conv_bkg and the plotted bkg curve disagree with what was
        # actually fitted. Declare it already convolved.
        is_bkg_map_conv = True

    add_background_to_model = (sky_mode != 'none')

    data_2D = _data_2D
    data_2D_gpu = _data_2D_gpu



    if use_weights is True:
        if rms_map is None:
            if logger is not None:
                logger.debug(f" ==> Generating RMS map from data.")
            rms_map = estimate_RMS_map(data_2D)
            rms_convention = 'sigma'
        else:
            if logger is not None:
                logger.debug(f" ==> Using provided RMS map "
                             f"(convention: {rms_convention}).")

        # Whatever convention the map arrives in, reduce it to a sigma first.
        _rms = np.asarray(rms_map, dtype=float)
        if rms_convention == 'sigma':
            sigma_map = _rms
        elif rms_convention == 'variance':
            sigma_map = np.sqrt(_rms)
        elif rms_convention in ('weight', 'invvar'):
            with np.errstate(divide='ignore', invalid='ignore'):
                sigma_map = 1.0 / np.sqrt(_rms)
        else:
            raise ValueError(
                f"rms_convention must be 'sigma', 'variance', 'weight' or "
                f"'invvar', got {rms_convention!r}.")

        # Pixels with no data (zero weight) become infinite sigma; keep them
        # finite but negligible so they cannot produce NaNs in the residual.
        _bad = ~np.isfinite(sigma_map) | (sigma_map <= 0)
        if np.any(_bad):
            _good = sigma_map[~_bad]
            _fill = (np.nanmedian(_good) * 1e6) if _good.size else 1.0
            sigma_map = np.where(_bad, _fill, sigma_map)

        if weight_mode == 'inverse_variance':
            # LMFIT squares this residual, so w = 1/sigma gives the statistically
            # correct chi2 = sum( (d - m)^2 / sigma^2 ).
            _weights_map = 1.0 / sigma_map
        elif weight_mode == 'legacy':
            # Historical behaviour: w = 1/sqrt(sigma), i.e. chi2 = sum(
            # (d - m)^2 / sigma ) -- the square root of the correct weight, which
            # under-weights noisy pixels relative to quiet ones. Kept only so
            # earlier results can be reproduced.
            _weights_map = 1.0 / np.sqrt(sigma_map)
            if logger is not None:
                logger.warning(" !!++==> weight_mode='legacy' applies the SQUARE "
                               "ROOT of the correct inverse-variance weight. "
                               "Use 'inverse_variance' for a proper chi2.")
        else:
            raise ValueError(
                f"weight_mode must be 'inverse_variance' or 'legacy', got "
                f"{weight_mode!r}.")

        # The normalisation only sets the overall scale of the residual, and with
        # scale_covar=True (LMFIT's default) it cannot affect the parameter
        # errors: the covariance is rescaled so that reduced chi2 = 1 regardless.
        #
        # With scale_covar=False it matters completely. There the residual has to
        # be in units of sigma for reduced chi2 to mean anything and for the
        # errors to reflect the real noise -- so the normalisation is skipped and
        # w = 1/sigma is used as-is. Normalising anyway would leave redchi orders
        # of magnitude from 1 and inflate every error bar by the same factor.
        if scale_covar:
            weights_map = _weights_map / np.nanmean(_weights_map)
        else:
            if logger is not None:
                logger.debug(" ==> scale_covar=False: using unnormalised "
                             "1/sigma weights so that redchi is absolute.")
                if loss not in (None, 'linear'):
                    # scipy's robust losses act on f^2/f_scale^2, so they are NOT
                    # scale invariant. Un-normalising the weights changes the
                    # residual scale by orders of magnitude, which switches a
                    # loss that was effectively linear into an aggressively
                    # robust one and moves the best fit. Observed on J0014:
                    # Rn went 7.54 -> 9.69 from this alone.
                    logger.warning(
                        f" !!++==> loss='{loss}' with scale_covar=False: the "
                        f"robust loss is not scale invariant and the "
                        f"unnormalised residuals are far larger than "
                        f"f_scale={f_scale}, so the loss is now doing real work "
                        f"and the best fit will differ. Use loss='linear', or "
                        f"set f_scale to the typical residual size.")
            weights_map = _weights_map
    else:
        if logger is not None:
            logger.debug(f" ==> Not using weights for fitting.")    
        weights_map = np.ones_like(data_2D)

    # Convert residual to GPU format if needed
    if convolution_mode == 'GPU':
        residual_2D = jnp.array(residual_2D)
        weights_map = jnp.array(weights_map)

    # Setup coordinate grid
    xy, size = _setup_coordinate_grid(data_2D, convolution_mode)

    if convolve_cutout is True:
        """
        WARNING: DO NOT USE FOR NOW!

        Instead of convolving the entire image,
        convolve only a box.
        Can be 10x faster.

        Issue:
        It causes the flat sky level to be much higher than the real value.

        Need further investigation and proper implementation.
        """
        x0c, y0c = int(size[0] / 2), int(size[1] / 2)

    #     FlatSky_level = background#mad_std(data_2D)

    # if convolution_mode == 'GPU':

    # FlatSky_level = mad_std(data_2D)

    # ========================================================================
    # MODEL CONSTRUCTION
    # ========================================================================

    nfunctions = ncomponents

    # Create minimizer functions
    min_residual_2D, min_residual_2D_GPU, build_model, func = \
        _create_minimizer_functions(nfunctions, xy, data_2D, data_2D_gpu,
                                    background, weights_map, residual_2D, PSF_DATA, PSF_DATA_raw,
                                    mask_for_fit, convolution_mode,
                                    add_background_to_model=add_background_to_model)

    # Setup JAX convolution if needed
    jax_convolve = _setup_jax_convolution(convolution_mode)

    # Construct model parameters
    smodel2D, params = construct_model_parameters(
        params_values_init_IMFIT=params_values_init_IMFIT, n_components=nfunctions,
        init_constraints=init_constraints, observation_type=observation_type,
        fix_n=fix_n, fix_value_n=fix_value_n,
        fix_max_value_n=fix_max_value_n,
        fix_min_value_n=fix_min_value_n,
        fix_max_value_Rn=fix_max_value_Rn, fix_min_value_Rn=fix_min_value_Rn,
        fix_x0_y0=fix_x0_y0, dr_fix=dr_fix, fix_geometry=fix_geometry,
        force_circular=force_circular,
        init_params=init_params, final_params=final_params,
        trunc=trunc,verbose=verbose,
        constrained=constrained)

    # ------------------------------------------------------------------
    # The prior on the background amplitude `s_a`.
    #
    # construct_model_parameters sets it from the module-level sky_*_bound
    # constants: init 0.99, min -1e-3, max +1.0. That is a RAW-data prior -- it
    # starts one background away from zero and cannot reach zero from above --
    # so on an image whose sky is already gone the correct answer (s_a -> 0) sits
    # on the lower bound and fights the source parameters on the way down.
    #
    #   'none'                     -> frozen at 0
    #   'flat'                      -> s_a IS the sky in image units, so its
    #                                 prior is set in those units: init at the
    #                                 sigma-clipped median of the image (0 if the
    #                                 sky is declared already gone), bounds
    #                                 init +/- 5 * mad_std(image). That window
    #                                 comfortably straddles zero and reaches well
    #                                 into the negative side, which raw data
    #                                 legitimately wants: a Sersic with large n
    #                                 over-produces the outskirts and a negative
    #                                 pedestal absorbs it (the degeneracy ProFit
    #                                 warns about).
    #   'map'                      -> s_a is a dimensionless multiplier on B;
    #                                 init 1.0 raw / 0.0 subtracted, adaptive
    #                                 symmetric bounds (below)
    #   radio, nothing asked for   -> untouched, so every existing radio fit
    #                                 reproduces exactly
    #
    # For 'map' the bounds are ADAPTIVE because a fixed range is wrong in one direction or
    # the other. `s_a` multiplies a map whose own amplitude varies by orders of
    # magnitude: on a raw image B is the whole sky and s_a ~ 1, but on an
    # already-subtracted image B is a small residual and the correction the data
    # wants can be several times it. J0014 galclean is exactly that case -- the
    # upstream step removed 0.0056 where the cutout's own sky is about 0.001, so
    # the clean image sits at -0.0045 while the background estimated FROM it has
    # a median of only -0.00094. Correcting that needs s_a ~ 5, and a fixed
    # max of 2.0 pins the parameter on its bound.
    #
    # The rule instead: s_a may shift the model by at most one pixel sigma, and
    # never by less than the plain factor-of-two range.
    #
    #     amplitude = max(|median(B)|, mad_std(B))
    #     limit     = max(2.0, sigma_image / amplitude)
    #
    # `sky_scale_bounds=(init, min, max)` overrides all of it.
    if not add_background_to_model:
        # With the background term switched off, `s_a` multiplies nothing: it has
        # no effect on the residual, so the Jacobian column for it is identically
        # zero, the covariance matrix is singular and LMFIT returns stderr=None
        # for EVERY parameter. Freeze it so the fit keeps its error bars, and so
        # it stops being counted as a free parameter in aic/bic.
        params['s_a'].set(value=0.0, vary=False)
    elif sky_scale_bounds is not None:
        _sa_init, _sa_min, _sa_max = sky_scale_bounds
        params['s_a'].set(value=float(_sa_init), min=float(_sa_min),
                          max=float(_sa_max), vary=True)
    elif observation_type == 'radio' and not _bkg_arg_given:
        # Legacy prior, deliberately left alone. Here s_a scales a shuffled noise
        # realisation added to the model rather than a sky pedestal, and its
        # optimum commonly sits on a bound; widening it would move published
        # radio results. Pass sky_scale_bounds to opt in.
        pass
    elif _flat_pedestal:
        # B == 1, so s_a is the sky itself and its prior belongs in data units.
        _img = np.asarray(data_2D, dtype=float)
        _img_sigma = float(mad_std(_img, ignore_nan=True))
        if is_background_subtracted:
            _sa_init = 0.0
        else:
            try:
                _sa_init = float(sigma_clipped_stats(_img, sigma=3.0)[1])
            except Exception:
                _sa_init = float(np.nanmedian(_img))
        if not np.isfinite(_img_sigma) or _img_sigma <= 0:
            _img_sigma = abs(_sa_init) if _sa_init else 1.0
        _sa_lo, _sa_hi = _sa_init - 5.0 * _img_sigma, _sa_init + 5.0 * _img_sigma
        params['s_a'].set(value=_sa_init, min=_sa_lo, max=_sa_hi, vary=True)
        if logger is not None:
            logger.debug(f" ==> s_a (flat sky, image units): init "
                         f"{_sa_init:+.5g}, bounds [{_sa_lo:+.5g}, "
                         f"{_sa_hi:+.5g}] (image sigma {_img_sigma:.5g}).")
    else:
        _sa_init = 0.0 if is_background_subtracted else 1.0
        _bkg_arr = np.asarray(background, dtype=float)
        _bkg_amp = max(abs(float(np.nanmedian(_bkg_arr))),
                       float(mad_std(_bkg_arr, ignore_nan=True))
                       if _bkg_arr.size > 1 else 0.0)
        _img_sigma = float(mad_std(np.asarray(data_2D, dtype=float),
                                   ignore_nan=True))
        if _bkg_amp > 0 and np.isfinite(_img_sigma):
            _sa_lim = max(2.0, _img_sigma / _bkg_amp)
        else:
            _sa_lim = 2.0
        params['s_a'].set(value=_sa_init, min=-_sa_lim, max=_sa_lim, vary=True)
        if logger is not None:
            logger.debug(f" ==> s_a: init {_sa_init:+.3g}, bounds "
                         f"[{-_sa_lim:+.4g}, {_sa_lim:+.4g}] "
                         f"(bkg amplitude {_bkg_amp:.4g}, "
                         f"image sigma {_img_sigma:.4g}).")

    # Create minimizer
    # `scale_covar=True` (LMFIT's default, and what morphen has always used
    # implicitly) rescales the covariance so that reduced chi2 = 1, which makes
    # the reported errors independent of the overall residual scale -- convenient,
    # but it also means redchi carries no goodness-of-fit information. Pass
    # scale_covar=False with a correctly calibrated rms_map to get an absolute
    # reduced chi2 and errors that reflect the real noise.
    if convolution_mode == 'CPU':
        mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
                               nan_policy='omit', reduce_fcn=reduce_fcn,
                               scale_covar=scale_covar)
    elif convolution_mode == 'GPU':
        mini = lmfit.Minimizer(min_residual_2D_GPU, params, max_nfev=200000,
                               nan_policy='omit', reduce_fcn=reduce_fcn,
                               scale_covar=scale_covar)
    else:
        # Fallback to CPU
        mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
                               nan_policy='omit', reduce_fcn=reduce_fcn,
                               scale_covar=scale_covar)

    # ========================================================================
    # OPTIMIZATION
    # ========================================================================

    result, result_1, result_extra = _perform_optimization(
        mini, params, method1, method2, convolution_mode,
        contrain_nelder, workers, max_nfev, tr_solver,
        regularize, x_scale, f_scale, ftol, xtol, gtol,
        verbose, loss, maxiter, maxfev, xatol, fatol,
        return_all, disp, de_options, parameters_mini_init
    )

    params = result.params

    # A background amplitude sitting on its bound means the fit wanted a larger
    # correction than it was allowed, so the sky is still wrong and everything
    # degenerate with it (the Sersic wings, mostly) absorbed the difference.
    # Worth saying out loud rather than leaving in the parameter table.
    #
    # Not on the radio legacy prior, though: there s_a scales a noise
    # realisation over the deliberately narrow window [-1e-3, +1] and lands on
    # one end or the other as a matter of course, so the warning would fire on
    # every radio fit and mean nothing.
    _sa = params['s_a']
    _sa_legacy_radio = (observation_type == 'radio' and not _bkg_arg_given
                        and sky_scale_bounds is None)
    if (_sa.vary and not _sa_legacy_radio
            and _sa.min is not None and _sa.max is not None):
        _span = float(_sa.max) - float(_sa.min)
        if _span > 0 and min(abs(_sa.value - _sa.min),
                             abs(_sa.value - _sa.max)) < 0.01 * _span:
            if _flat_pedestal:
                _why = ("the fitted pedestal hit the edge of its window. Either "
                        "the image has a much larger offset than 5 sigma, or the "
                        "sky is degenerate with the outer Sersic wings")
            else:
                _why = ("the fit wanted a larger sky correction than allowed; "
                        "the background map's amplitude is probably wrong")
            _msg = (f" !!++==> the sky amplitude s_a converged onto its bound "
                    f"({_sa.value:+.4g} in [{_sa.min:+.4g}, {_sa.max:+.4g}], "
                    f"sky_mode='{sky_mode}'). {_why}. Widen it with "
                    f"sky_scale_bounds=(init, min, max), or re-estimate the "
                    f"background.")
            if logger is not None:
                logger.warning(_msg)
            else:
                print('--==>>' + _msg)

    # ========================================================================
    # MODEL GENERATION AND OUTPUT
    # ========================================================================

    # Generate model components
    model_dict, flat_sky_total, flat_sky_total_dec, bkg_comp_i, bkg_comp_i_dec = \
        _generate_model_components(params, ncomponents, xy, size, background,
                                   background_dec, PSF_CONV, PSF_DATA, PSF_DATA_raw,
                                   convolution_mode, is_bkg_map_conv)

    # Compute residuals
    model_dict = _compute_residuals(data_2D, model_dict, flat_sky_total,
                                    is_bkg_map_conv, PSF_CONV, PSF_DATA,
                                    convolution_mode)

    # Save outputs
    (image_results_conv, image_results_deconv, bkg_images,
     total_image_results_conv, total_image_results_deconv) = \
        _save_fitting_outputs(imagename, model_dict, ncomponents, special_name,
                              save_name_append, flat_sky_total_dec)

    # ========================================================================
    # FINAL PROCESSING
    # ========================================================================




    
    # <<testing>>
    # Extract model_temp for return (last component generated)
    model_temp = Model(sersic2D_GPU)

    # Create a PrettyTable object
    table = PrettyTable()
    table.field_names = ["Parameter", "Best-Fit", "3 * std"]

    if verbose > 0:
        # Add rows to the table with formatted values
        for param_name in result.params.keys():
            param = result.params[param_name]
            value = param.value
            stderr = param.stderr
            
            # Calculate 3 * stderr if available, otherwise use None
            if stderr is not None:
                three_sigma = 3 * stderr
                table.add_row([
                    param_name,
                    f"{value:.6f}",
                    f"{three_sigma:.6f}"
                ])
            else:
                table.add_row([
                    param_name,
                    f"{value:.6f}",
                    "None"
                ])

        # Print the table
        print(table)

    exec_time = time.time() - startTime
    if verbose > 0:
        print('Exec time fitting=', exec_time, 's')


    return (result, mini, result_1, result_extra, model_dict, image_results_conv,
            image_results_deconv, bkg_images, smodel2D, model_temp)


# """
# Refactored do_fit2D with modular structure - Version 2
#
# This refactored version maintains complete backward compatibility with the original
# function signature and return values, while organizing the internal logic into
# clear, maintainable helper functions.
#
# IMPORTANT: This module assumes all necessary imports are already available in the
# calling module's namespace. It does NOT import anything itself. Required imports:
# - time, numpy as np, scipy.signal, astropy.io.fits as pf
# - jax.numpy as jnp, jax.jit (when GPU mode is used)
# - lmfit
# - All radio_utils functions: sersic2D, sersic2D_GPU, FlatSky, FlatSky_cpu,
#   _fftconvolve_jax, construct_model_parameters, constrain_nelder_mead_params,
#   copy_header, save_results_csv, load_fits_data
# - All signal_stats functions: mad_std, shuffle_2D, sep_background
# - lmfit.Model
#
# This design maintains consistency with the existing morphen codebase structure.
# """
#
#
# # ============================================================================
# # HELPER FUNCTIONS - Data Initialization
# # ============================================================================
#
# def _initialize_image_data(imagename, data_2D_, convolution_mode, pf, jnp, logger=None):
#     """
#     Load and prepare image data for fitting.
#
#     Parameters
#     ----------
#     imagename : str
#         Path to image file
#     data_2D_ : array or None
#         Pre-loaded image data
#     convolution_mode : str
#         'GPU' or 'CPU'
#     pf : module
#         astropy.io.fits module
#     jnp : module
#         jax.numpy module (or numpy if JAX unavailable)
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     data_2D : array
#         Image data in numpy format
#     data_2D_gpu : array or None
#         Image data in JAX format (if GPU mode)
#     """
#     if data_2D_ is None:
#         data_2D = pf.getdata(imagename)
#     else:
#         data_2D = data_2D_
#
#     # Prepare GPU array if needed
#     if convolution_mode == 'GPU':
#         data_2D_gpu = jnp.array(data_2D)
#     else:
#         data_2D_gpu = None
#
#     return data_2D, data_2D_gpu
#
#
# def _initialize_psf_data(psf_name, convolution_mode, pf, jnp, load_fits_data, logger=None):
#     """
#     Load and prepare PSF data for convolution.
#
#     Parameters
#     ----------
#     psf_name : str or None
#         Path to PSF file
#     convolution_mode : str
#         'GPU' or 'CPU'
#     pf : module
#         astropy.io.fits module
#     jnp : module
#         jax.numpy module (or numpy if JAX unavailable)
#     load_fits_data : callable
#         Function to load FITS data as fallback
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array or None
#         PSF data in appropriate format for computation
#     PSF_DATA_raw : array or None
#         PSF data in raw numpy format
#     """
#     if psf_name is not None:
#         PSF_CONV = True
#         try:
#             PSF_DATA_raw = pf.getdata(psf_name)
#             if len(PSF_DATA_raw.shape) == 4:
#                 PSF_DATA_raw = PSF_DATA_raw[0][0]
#         except:
#             PSF_DATA_raw = load_fits_data(psf_name)
#
#         if convolution_mode == 'GPU':
#             if logger is not None:
#                 logger.debug(f"---------------------------------------")
#                 logger.debug(f" <<< PERFORMING CONVOLUTION WITH JAX >>> ")
#                 logger.debug(f"---------------------------------------")
#             PSF_DATA = jnp.array(PSF_DATA_raw)
#         elif convolution_mode == 'CPU':
#             PSF_DATA = PSF_DATA_raw
#         else:
#             PSF_DATA = PSF_DATA_raw
#
#         # PSF_DATA = pf.getdata(
#         #     imagename.replace('-image.cutout.fits', '-beampsf.cutout.fits'))
#     else:
#         PSF_CONV = False
#         PSF_DATA = None
#         PSF_DATA_raw = None
#
#     return PSF_CONV, PSF_DATA, PSF_DATA_raw
#
#
# def _prepare_mask(mask_region, convolution_mode, jnp, logger=None):
#     """
#     Prepare mask for constrained fitting.
#
#     Parameters
#     ----------
#     mask_region : array or None
#         Mask array
#     convolution_mode : str
#         'GPU' or 'CPU'
#     jnp : module
#         jax.numpy module (or numpy if JAX unavailable)
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     mask_for_fit : array or None
#         Mask in appropriate format
#     """
#     if mask_region is not None:
#         """
#
#         """
#         if logger is not None:
#             logger.debug(f" ==> Using provided mask region to constrain fit. ")
#             logger.warning(f" !!++==> Fitting with a mask is faster, but experimental!! \n"
#                            f"         Use with caution.")
#         # data_2D = data_2D * mask_region
#         if convolution_mode == 'GPU':
#             mask_for_fit = jnp.array(mask_region)
#         else:
#             mask_for_fit = mask_region
#     else:
#         mask_for_fit = None
#
#     return mask_for_fit
#
#
# def _prepare_background_and_residual(residualname, residualdata_2D_, which_residual,
#                                      self_bkg, rms_map, imagename, data_2D,
#                                      PSF_CONV, PSF_DATA, PSF_DATA_raw, convolution_mode,
#                                      is_bkg_map_conv,
#                                      pf, jnp, scipy, mad_std, shuffle_2D, sep_background,
#                                      _fftconvolve_jax, logger=None):
#     """
#     Prepare background and residual arrays for fitting.
#
#     This function handles the complex logic of preparing background maps,
#     residual images, and RMS estimates based on various input configurations
#     optimized for different observation types (radio vs optical).
#
#     Background/Residual Strategy by Observation Type
#     -------------------------------------------------
#
#     RADIO IMAGES (residualname provided, which_residual != 'user'):
#         - Uses residual map from interferometric cleaning process
#         - 'shuffled' mode: Shuffles residual to create noise realization
#           This breaks spatial correlations while preserving noise statistics
#         - 'natural' mode: Uses residual as-is without shuffling
#         - Residual added to model via optimized scaling factor s_a
#         - Ensures flux conservation relative to cleaned image noise
#         - FlatSky_level computed from residual for initial scaling
#
#     OPTICAL IMAGES (which_residual = 'user'):
#         - User provides custom RMS/error/background map via rms_map parameter
#         - Could be photometric error map from pipeline
#         - Could be background estimate from source extraction
#         - Requires explicit rms_map parameter (raises ValueError if None)
#         - Allows pre-convolved error maps (set is_bkg_map_conv=True)
#
#     SELF-BACKGROUND MODE (self_bkg = True, no residual provided):
#         - Estimates background directly from input image using SEP
#         - Shuffles background estimate to create noise realization
#         - Useful when no external background/residual available
#         - Appropriate for well-behaved optical images
#
#     FLAT SKY MODE (default fallback, no inputs):
#         - Estimates single RMS value from entire image
#         - Uses MAD (Median Absolute Deviation) estimator for robustness
#         - Assumes uniform noise across image
#         - Simplest case, background is scalar not array
#
#     PSF Convolution of Background/Error Maps
#     -----------------------------------------
#     The is_bkg_map_conv flag controls convolution behavior:
#
#     - If is_bkg_map_conv = False AND PSF provided:
#       Background/error map is convolved with PSF to match model space
#       Critical for radio residuals that need beam matching
#       Optical error maps from unconvolved data need this
#
#     - If is_bkg_map_conv = True:
#       Background/error already convolved, use as-is
#       Common for optical pipelines that provide PSF-matched errors
#       Avoids double-convolution artifacts
#
#     - If no PSF (PSF_CONV = False):
#       No convolution regardless of is_bkg_map_conv setting
#       Background/error used directly
#
#     Implementation Notes
#     --------------------
#     - background_dec stores the unconvolved version for deconvolved model output
#     - background stores the (possibly convolved) version for fitting
#     - Both are needed to produce conv/deconv model components correctly
#     - Scalar backgrounds (FlatSky_level case) are not convolved
#
#     Parameters
#     ----------
#     residualname : str or None
#         Path to residual image (radio: cleaning residual, optical: error map)
#     residualdata_2D_ : array or None
#         Pre-loaded residual data (avoids re-reading file)
#     which_residual : str
#         Residual processing mode: 'shuffled', 'natural', or 'user'
#     self_bkg : bool
#         Whether to use self-background estimation from image
#     rms_map : array or None
#         User-provided RMS/background/error map (for 'user' mode)
#     imagename : str
#         Path to main image (for self-background estimation)
#     data_2D : array
#         Main image data (for flat sky MAD estimation)
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array or None
#         PSF data in computation format (GPU or CPU)
#     PSF_DATA_raw : array or None
#         PSF data in raw numpy format (for CPU convolution)
#     convolution_mode : str
#         'GPU' or 'CPU'
#     is_bkg_map_conv : bool
#         Whether the RMS/error map is already PSF-convolved
#         (True prevents double-convolution for optical error maps)
#     pf : module
#         astropy.io.fits module
#     jnp : module
#         jax.numpy module
#     scipy : module
#         scipy module
#     mad_std : callable
#         Function to compute MAD-based standard deviation
#     shuffle_2D : callable
#         Function to shuffle 2D array for noise realization
#     sep_background : callable
#         Function to estimate background using SEP
#     _fftconvolve_jax : callable
#         JAX-based FFT convolution function
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     background : array or scalar
#         Background array for fitting (possibly PSF-convolved)
#         Array for map-based backgrounds, scalar for flat sky
#     background_dec : array or scalar
#         Deconvolved/unconvolved background (always unconvolved version)
#         Used for generating deconvolved model components
#     residual_2D : array
#         Residual data array (loaded from file or background fallback)
#     FlatSky_level : float or None
#         Flat sky RMS estimate (None if using map-based background)
#         Used for scaling in some model configurations
#     """
#     FlatSky_level = None
#
#     if residualname is not None and which_residual != 'user':
#         """
#         This is important for radio image fitting.
#
#         It uses the shuffled version of the residual cleaned image
#         originated from the interferometric deconvolution.
#
#         This ensures that the best model created here will be on top
#         of that rms noise so that flux conservation is maximized.
#
#         However, this residual is not added as model + shuffled_residual
#         only, but instead by a multiplication factor,
#         e.g. model + const* shuffled_residual, and const will be minimized
#         as well during the fitting (here, called `s_a`).
#         """
#         if residualdata_2D_ is not None:
#             residual_2D = residualdata_2D_
#         else:
#             residual_2D = pf.getdata(residualname)
#
#         if which_residual == 'shuffled':
#             if logger is not None:
#                 logger.debug(f" ==> Using clean shuffled background for optmization... ")
#             residual_2D_to_use = shuffle_2D(residual_2D)
#         elif which_residual == 'natural':
#             if logger is not None:
#                 logger.debug(f" ==> Using clean background for optmization... ")
#             """
#             if psf_name is not None:
#                 if logger is not None:
#                     logger.debug(f" ====> Deconvolving residual map... ")
#                 residual_2D_to_use, _ = deconvolve_fft(residual_2D,
#                                                             PSF_DATA_raw/PSF_DATA_raw.sum())
#             else:
#                 residual_2D_to_use = residual_2D
#             """
#             residual_2D_to_use = residual_2D
#         else:
#             residual_2D_to_use = residual_2D
#
#         FlatSky_level = mad_std(residual_2D_to_use)
#         #         background = residual_2D #residual_2D_to_use
#         if convolution_mode == 'GPU':
#             background = jnp.array(residual_2D_to_use)
#         else:
#             background = residual_2D_to_use
#
#     else:
#         if which_residual == 'user':
#             if rms_map is None:
#                 print('--==>> A rms map/background mode was selected (user)')
#                 print('       but no rms/background map was provided.')
#                 print('       Please, provide a rms/background map.')
#                 print('||==>> Stopping code now.')
#                 raise ValueError("rms_map should not be None when which_residual='user'")
#             else:
#                 if logger is not None:
#                     logger.debug(f" ==> Using provided RMS map. ")
#                 background_map = rms_map
#                 background = background_map.copy()
#         else:
#             if self_bkg == True:
#                 if logger is not None:
#                     logger.warning(f" ==> No residual/background provided. Using image bkg map... ")
#                 background_map = sep_background(imagename)
#                 background = shuffle_2D(background_map.back())
#             else:
#                 if logger is not None:
#                     logger.warning(f" ==> Using only flat sky for rms bkg.")
#                 FlatSky_level = mad_std(data_2D)
#                 background = FlatSky_level
#
#     # Load residual_2D if not already loaded
#     # This ensures residual_2D is always defined for downstream use
#     if residualdata_2D_ is not None:
#         residual_2D = residualdata_2D_
#     else:
#         try:
#             residual_2D = pf.getdata(residualname)
#         except:
#             # Fallback: use background as residual if file doesn't exist
#             residual_2D = background
#
#     # Create deconvolved background copy BEFORE any convolution
#     # This preserves the unconvolved version for deconvolved model components
#     if isinstance(background, (int, float)):
#         # FlatSky_level case - background is a scalar, not an array
#         background_dec = background
#     else:
#         # Array case - create copy before potential convolution
#         background_dec = background.copy()
#
#     # Convolve background with PSF if needed
#     # Only convolve if: (1) PSF exists, (2) map not pre-convolved, (3) background is array
#     if is_bkg_map_conv is False and PSF_CONV:
#         if logger is not None:
#             logger.debug(f" ==> RMS map is not convolved, convolving with PSF now.")
#
#         if not isinstance(background, (int, float)):
#             # Only convolve if background is an array, not a scalar
#             if convolution_mode == 'GPU':
#                 background = _fftconvolve_jax(background, PSF_DATA)
#             elif convolution_mode == 'CPU':
#                 background = scipy.signal.fftconvolve(background, PSF_DATA_raw, 'same')
#
#     return background, background_dec, residual_2D, FlatSky_level
#
#
# def _setup_coordinate_grid(data_2D, convolution_mode, np, jnp):
#     """
#     Create coordinate meshgrid for model evaluation.
#
#     Parameters
#     ----------
#     data_2D : array
#         Image data (for shape)
#     convolution_mode : str
#         'GPU' or 'CPU'
#     np : module
#         numpy module
#     jnp : module
#         jax.numpy module
#
#     Returns
#     -------
#     xy : tuple of arrays
#         Coordinate meshgrid
#     size : tuple
#         Image shape
#     """
#     size = data_2D.shape
#     if convolution_mode == 'GPU':
#         x, y = jnp.meshgrid(jnp.arange((size[1])), jnp.arange((size[0])))
#         xy = jnp.stack([x, y], axis=0)
#     else:
#         xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
#
#     return xy, size
#
#
# # ============================================================================
# # HELPER FUNCTIONS - Model Construction
# # ============================================================================
#
# def _create_minimizer_functions(nfunctions, xy, data_2D, data_2D_gpu, background,
#                                 residual_2D, PSF_DATA, PSF_DATA_raw, mask_for_fit,
#                                 convolution_mode,
#                                 np, jnp, scipy, jit,
#                                 sersic2D, sersic2D_GPU, FlatSky, FlatSky_cpu,
#                                 _fftconvolve_jax):
#     """
#     Create the residual functions for minimization.
#
#     Parameters
#     ----------
#     nfunctions : int
#         Number of Sersic components
#     xy : tuple of arrays
#         Coordinate meshgrid
#     data_2D : array
#         Image data (CPU)
#     data_2D_gpu : array or None
#         Image data (GPU)
#     background : array
#         Background map
#     residual_2D : array
#         Residual data
#     PSF_DATA : array
#         PSF data in computation format
#     PSF_DATA_raw : array
#         PSF data in raw format
#     mask_for_fit : array or None
#         Fitting mask
#     convolution_mode : str
#         'GPU' or 'CPU'
#     np : module
#         numpy module
#     jnp : module
#         jax.numpy module
#     scipy : module
#         scipy module
#     jit : callable
#         JAX JIT compiler
#     sersic2D : callable
#         CPU Sersic 2D function
#     sersic2D_GPU : callable
#         GPU Sersic 2D function
#     FlatSky : callable
#         GPU flat sky function
#     FlatSky_cpu : callable
#         CPU flat sky function
#     _fftconvolve_jax : callable
#         JAX FFT convolution
#
#     Returns
#     -------
#     min_residual_2D : callable or None
#         CPU residual function
#     min_residual_2D_GPU : callable or None
#         GPU residual function
#     build_model : callable or None
#         Model building function (for GPU)
#     func : callable or None
#         Parameter splitting function (for GPU)
#     """
#
#     # CPU version
#     def min_residual_2D(params):
#         dict_model = {}
#         model = 0
#         for i in range(1, nfunctions + 1):
#             model = model + sersic2D(xy, params['f' + str(i) + '_x0'],
#                                      params['f' + str(i) + '_y0'],
#                                      params['f' + str(i) + '_PA'],
#                                      params['f' + str(i) + '_ell'],
#                                      params['f' + str(i) + '_n'],
#                                      params['f' + str(i) + '_In'],
#                                      params['f' + str(i) + '_Rn'],
#                                      params['f' + str(i) + '_cg'], )
#         # print(model.shape)
#         # model = model + FlatSky_cpu(FlatSky_level, params['s_a'])*background
#         # model = model + FlatSky_cpu(background, params['s_a'])
#         MODEL_2D_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + \
#                         FlatSky_cpu(background, params['s_a'])
#         residual = data_2D - MODEL_2D_conv
#         return np.ravel(residual)
#
#     # GPU version
#     try:
#         # @partial(jit, static_argnums=1)
#         @jit
#         def func(x):
#             return jnp.split(x, nfunctions)
#
#         @jit
#         def build_model(xy, param_matrix):
#             model = 0
#             for model_params in param_matrix:
#                 model = model + sersic2D_GPU(xy, model_params[0],
#                                              model_params[1],
#                                              model_params[2],
#                                              model_params[3],
#                                              model_params[4],
#                                              model_params[5],
#                                              model_params[6],
#                                              model_params[7])
#             return model
#
#     except:
#         def func(x):
#             return jnp.split(x, nfunctions)
#
#         def build_model(xy, param_matrix):
#             model = 0
#             for model_params in param_matrix:
#                 model = model + sersic2D_GPU(xy, model_params[0],
#                                              model_params[1],
#                                              model_params[2],
#                                              model_params[3],
#                                              model_params[4],
#                                              model_params[5],
#                                              model_params[6],
#                                              model_params[7])
#             return model
#
#     def min_residual_2D_GPU(params):
#         model = 0
#         for i in range(1, nfunctions + 1):
#             model = model + sersic2D_GPU(xy,
#                                          params['f' + str(i) + '_x0'].value,
#                                          params['f' + str(i) + '_y0'].value,
#                                          params['f' + str(i) + '_PA'].value,
#                                          params['f' + str(i) + '_ell'].value,
#                                          params['f' + str(i) + '_n'].value,
#                                          params['f' + str(i) + '_In'].value,
#                                          params['f' + str(i) + '_Rn'].value,
#                                          params['f' + str(i) + '_cg'].value)
#
#         # # param_matrix = extract_params(params)
#         # param_matrix = func(jnp.array(list(params.valuesdict().values()))[:-1])
#         # model = build_model(xy,param_matrix)
#
#         MODEL_2D_conv = _fftconvolve_jax(model, PSF_DATA) + FlatSky(background, params['s_a'].value)
#
#         # MODEL_2D_conv = _fftconvolve_jax(model+
#         #                                  FlatSky(background,params['s_a'].value),
#         #                                  PSF_DATA)
#         # residual = ((data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])/
#         #             (1000*(abs(residual_2D[mask_for_fit])+1.0e-6)))
#         residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])
#         return np.asarray(residual).copy()
#         # weights = 1/((background[mask_for_fit])/data_2D_gpu[mask_for_fit])
#         # weightned_residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit]) * jnp.sqrt(weights)
#         # return np.asarray(weightned_residual).copy()
#
#     if convolution_mode == 'CPU':
#         return min_residual_2D, None, build_model, func
#     else:
#         return None, min_residual_2D_GPU, build_model, func
#
#
# def _setup_jax_convolution(convolution_mode, jnp, jit):
#     """
#     Setup JAX convolution function if in GPU mode.
#
#     Parameters
#     ----------
#     convolution_mode : str
#         'GPU' or 'CPU'
#     jnp : module
#         jax.numpy module
#     jit : callable
#         JAX JIT compiler
#
#     Returns
#     -------
#     jax_convolve : callable or None
#         JIT-compiled convolution function
#     """
#     if convolution_mode == 'GPU':
#         @jit
#         def convolve_on_gpu(image, psf):
#             """
#             This was before jax.scipy implementing fftconvolve.
#             It provides the same result, at the same speed.
#
#             This function also accepts PSFs with a different shape of the image.
#
#             """
#             # Calculate the new padded shape
#             padded_shape = (image.shape[0] + psf.shape[0] - 1,
#                             image.shape[1] + psf.shape[1] - 1)
#
#             # Pad both image and psf to the new shape
#             pad_shape = [(0, ts - s) for s, ts in zip(image.shape, padded_shape)]
#             image_padded = jnp.pad(image, pad_shape, mode='constant')
#             pad_shape = [(0, ts - s) for s, ts in zip(psf.shape, padded_shape)]
#             psf_padded = jnp.pad(psf, pad_shape, mode='constant')
#             # psf_padded = pad_for_convolution(psf, padded_shape)
#             image_fft = jnp.fft.fft2(image_padded)
#             psf_fft = jnp.fft.fft2(psf_padded)
#
#             conv_fft = image_fft * psf_fft
#
#             # Get the real part of the inverse FFT and crop to the original image size
#             result_full = jnp.real(jnp.fft.ifft2(conv_fft))
#             return result_full[psf.shape[0] // 2:image.shape[0] + psf.shape[0] // 2,
#             psf.shape[1] // 2:image.shape[1] + psf.shape[1] // 2]
#
#         jax_convolve = jit(convolve_on_gpu)
#         return jax_convolve
#     return None
#
#
# # ============================================================================
# # HELPER FUNCTIONS - Optimization
# # ============================================================================
#
# def _perform_optimization(mini, params, method1, method2, convolution_mode,
#                           contrain_nelder, workers, max_nfev, tr_solver,
#                           regularize, x_scale, f_scale, ftol, xtol, gtol,
#                           verbose, loss, maxiter, maxfev, xatol, fatol,
#                           return_all, disp, de_options, parameters_mini_init,
#                           constrain_nelder_mead_params):
#     """
#     Execute the two-stage optimization process.
#
#     Parameters
#     ----------
#     mini : lmfit.Minimizer
#         Minimizer object
#     params : lmfit.Parameters
#         Initial parameters
#     method1 : str
#         First optimization method
#     method2 : str
#         Second optimization method
#     convolution_mode : str
#         'GPU' or 'CPU'
#     contrain_nelder : bool
#         Whether to constrain Nelder-Mead parameters
#     workers : int
#         Number of parallel workers
#     max_nfev : int
#         Maximum function evaluations
#     tr_solver : str
#         Trust region solver
#     regularize : bool
#         Whether to regularize
#     x_scale : str or array
#         Parameter scaling
#     f_scale : float
#         Function scaling
#     ftol : float
#         Function tolerance
#     xtol : float
#         Parameter tolerance
#     gtol : float
#         Gradient tolerance
#     verbose : int
#         Verbosity level
#     loss : str
#         Loss function
#     maxiter : int
#         Maximum iterations for Nelder-Mead
#     maxfev : int
#         Maximum function evaluations for Nelder-Mead
#     xatol : float
#         Absolute parameter tolerance for Nelder-Mead
#     fatol : float
#         Absolute function tolerance for Nelder-Mead
#     return_all : bool
#         Return all Nelder-Mead results
#     disp : bool
#         Display convergence messages
#     de_options : dict
#         Options for differential evolution
#     parameters_mini_init : lmfit.Parameters or None
#         Initial parameters from previous run
#     constrain_nelder_mead_params : callable
#         Function to constrain Nelder-Mead parameters
#
#     Returns
#     -------
#     result : lmfit.MinimizerResult
#         Final optimization result
#     result_1 : lmfit.MinimizerResult
#         First stage optimization result
#     result_extra : None
#         Placeholder for additional results
#     """
#     # initial minimization.
#     print(' >> Using', method1, ' solver for first optimisation run... ')
#     # take parameters from previous run, and re-optimize them.
#     #     method2 = 'ampgo'#'least_squares'
#     #     method2 = 'least_squares'
#     result_extra = None
#
#     if method1 == 'nelder':
#         # very robust, but takes time....
#         #         print(' >> Using', method1,' solver for first optimisation run... ')
#         result_1 = mini.minimize(method='nelder',
#                                  #                                  xatol = 1e-12, fatol = 1e-12, disp=True,
#                                  #                                  adaptive = True,max_nfev = 30000,
#                                  options={'maxiter': maxiter, 'maxfev': maxfev,
#                                           'xatol': xatol, 'fatol': fatol,
#                                           'return_all': return_all,
#                                           'disp': disp}
#                                  )
#
#     elif method1 == 'least_squares':
#         # faster, but usually not good for first run.
#         # if results_previous_run is not None:
#         print(' >> Using', tr_solver, 'for tr solver, with regularize set to', regularize,
#               ' Loss is', loss, '.')
#
#         if parameters_mini_init is not None:
#             print(f'  ++==>> Using initial mini parameters from a previous run.')
#             # try:
#             result_1 = mini.minimize(method='least_squares',
#                                      params=parameters_mini_init,
#                                      max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
#                                      tr_solver=tr_solver,
#                                      tr_options={'regularize': regularize,
#                                                  },
#                                      ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                      loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#         else:
#             result_1 = mini.minimize(method='least_squares',
#                                      max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
#                                      tr_solver=tr_solver,
#                                      tr_options={'regularize': regularize,
#                                                  #                                              'min_delta':1e-14, 'eta':0.05,
#                                                  #                                              'xtol':1e-14, 'gtol':1e-14,
#                                                  #                                              'ftol':1e-14
#                                                  },
#                                      ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                      loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#
#     elif method1 == 'differential_evolution':
#         # de is giving some issues, I do not know why.
#         result_1 = mini.minimize(method='differential_evolution',
#                                  options={'disp': True, 'workers': workers,
#                                           'max_nfev': max_nfev, 'vectorized': True,
#                                           'strategy': 'randtobest1bin',
#                                           'mutation': (0.5, 1.5),
#                                           'recombination': [0.2, 0.9],
#                                           'init': 'random', 'tol': 0.00001,
#                                           'updating': 'deferred',
#                                           'popsize': 600})
#         # result_1 = mini.minimize(method='differential_evolution', popsize=600,
#         #                          disp=True,  # init = 'random',
#         #                          # mutation=(0.5, 1.5), recombination=[0.2, 0.9],
#         #                          max_nfev=20000,
#         #                          workers=1, updating='deferred', vectorized=True)
#     else:
#         raise ValueError(f"Unknown method1: {method1}")
#
#     print(' >> Using', method2, ' solver for second optimisation run... ')
#
#     second_run_params = result_1.params
#     if (contrain_nelder == True) and (method2 == 'nelder'):
#         """
#         It seems that least_squares is ignoring the best-parameters provided by
#         Nelder-mead, which means that it is lookig the parameter space far away
#         from the optimised Nelder-Mead ones.
#
#         So, with this condition, we force a much smaller searching region, but
#         it assumes that Nelder opt was good (which is not always true).
#
#         YOU MUST CHECK YOUR RESULTS!!!!
#
#         """
#         print('Constraining Nelder-Mead Parameters for method', method2)
#         params_constrained = constrain_nelder_mead_params(result_1.params,
#                                                           max_factor=1.03,
#                                                           min_factor=0.97)
#         # UPDATE THE SECOND RUN PARAMETERS TO BE THE CONSTRAINED ONES.
#         second_run_params = params_constrained
#
#     if method2 == 'nelder':
#         result = mini.minimize(method='nelder', params=second_run_params,
#                                options={'maxiter': maxiter, 'maxfev': maxfev,
#                                         'xatol': xatol, 'fatol': fatol,
#                                         'disp': disp})
#
#     elif method2 == 'ampgo':
#         # ampgo is not workin well/ takes so long ???
#         result = mini.minimize(method='ampgo', params=second_run_params,
#                                maxfunevals=10000, totaliter=30, disp=True,
#                                maxiter=5, glbtol=1e-8)
#
#     elif method2 == 'least_squares':
#         # faster, usually converges and provide errors.
#         # Very robust if used in second opt from first opt parameters.
#         result = mini.minimize(method='least_squares',
#                                params=second_run_params,
#                                max_nfev=max_nfev,
#                                tr_solver=tr_solver,
#                                tr_options={'regularize': regularize,
#                                            #                                            'min_delta': 1e-14, 'eta': 0.05,
#                                            #                                            'xtol': 1e-14, 'gtol': 1e-14,
#                                            #                                            'ftol': 1e-14
#                                            },
#                                x_scale=x_scale, f_scale=f_scale,
#                                ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#
#     elif method2 == 'differential_evolution':
#         # result = mini.minimize(method='differential_evolution',
#         #                        params=second_run_params,
#         #                        options={'maxiter': 30000, 'workers': -1,
#         #                                 'tol': 0.001, 'vectorized': True,
#         #                                 'strategy': 'randtobest1bin',
#         #                                 'updating': 'deferred', 'disp': True,
#         #                                 'seed': 1}
#         #                        )
#         result = mini.minimize(method='differential_evolution',
#                                params=second_run_params,
#                                options=de_options
#                                )
#     else:
#         raise ValueError(f"Unknown method2: {method2}")
#
#     return result, result_1, result_extra
#
#
# # ============================================================================
# # HELPER FUNCTIONS - Output Generation
# # ============================================================================
#
# def _generate_model_components(params, ncomponents, xy, size, background,
#                                background_dec, PSF_CONV, PSF_DATA, PSF_DATA_raw,
#                                convolution_mode, is_bkg_map_conv,
#                                np, scipy, Model, sersic2D_GPU, FlatSky, FlatSky_cpu,
#                                _fftconvolve_jax):
#     """
#     Generate individual and total model components.
#
#     Parameters
#     ----------
#     params : lmfit.Parameters
#         Optimized parameters
#     ncomponents : int
#         Number of Sersic components
#     xy : tuple of arrays
#         Coordinate meshgrid
#     size : tuple
#         Image shape
#     background : array
#         Background map (possibly convolved)
#     background_dec : array
#         Deconvolved background map
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array
#         PSF data in computation format
#     PSF_DATA_raw : array
#         PSF data in raw format
#     convolution_mode : str
#         'GPU' or 'CPU'
#     is_bkg_map_conv : bool
#         Whether RMS map is already convolved
#     np : module
#         numpy module
#     scipy : module
#         scipy module
#     Model : class
#         lmfit.Model class
#     sersic2D_GPU : callable
#         GPU Sersic function
#     FlatSky : callable
#         GPU flat sky function
#     FlatSky_cpu : callable
#         CPU flat sky function
#     _fftconvolve_jax : callable
#         JAX FFT convolution
#
#     Returns
#     -------
#     model_dict : dict
#         Dictionary containing all model arrays
#     flat_sky_total : array
#         Total flat sky background (convolved)
#     flat_sky_total_dec : array
#         Total flat sky background (deconvolved)
#     bkg_comp_i : array
#         Background component for individual models
#     bkg_comp_i_dec : array
#         Deconvolved background component for individual models
#     """
#     model_temp = Model(sersic2D_GPU)
#     xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
#     model = 0
#     model_dict = {}
#
#     if convolution_mode == 'GPU':
#         flat_sky_total = np.asarray(FlatSky(background, params['s_a'].value))
#         flat_sky_total_dec = np.asarray(FlatSky(background_dec, params['s_a'].value))
#     elif convolution_mode == 'CPU':
#         flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
#         flat_sky_total_dec = FlatSky_cpu(background_dec, params['s_a'].value)
#     else:
#         # Fallback
#         flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
#         flat_sky_total_dec = FlatSky_cpu(background_dec, params['s_a'].value)
#
#     bkg_comp_i = flat_sky_total.copy()
#     bkg_comp_i_dec = flat_sky_total_dec.copy()
#
#     for i in range(1, ncomponents + 1):
#         model_temp = sersic2D_GPU(xy, params['f' + str(i) + '_x0'].value,
#                                   params['f' + str(i) + '_y0'].value,
#                                   params['f' + str(i) + '_PA'].value,
#                                   params['f' + str(i) + '_ell'].value,
#                                   params['f' + str(i) + '_n'].value,
#                                   params['f' + str(i) + '_In'].value,
#                                   params['f' + str(i) + '_Rn'].value,
#                                   params['f' + str(i) + '_cg'].value)
#
#         model = model + model_temp
#         # to each individual component, add the bkg map.
#         model_dict['model_c' + str(i)] = np.asarray(model_temp + bkg_comp_i_dec)
#
#         if PSF_CONV == True:
#             if convolution_mode == 'GPU':
#                 # model_dict['model_c' + str(i) + '_conv'] = np.asarray(jax_convolve(model_temp, PSF_DATA)).copy()
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #     np.asarray(_fftconvolve_jax(model_temp,PSF_DATA).copy()+bkg_comp_i))
#                 # to each individual component, add the bkg map.
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #     np.asarray(_fftconvolve_jax(model_temp+bkg_comp_i,PSF_DATA).copy()))
#
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                                                                np.asarray(_fftconvolve_jax(model_temp,
#                                                                                            PSF_DATA).copy())) + bkg_comp_i
#
#             elif convolution_mode == 'CPU':
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #         scipy.signal.fftconvolve(model_temp+bkg_comp_i, PSF_DATA_raw,'same'))
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                         scipy.signal.fftconvolve(model_temp, PSF_DATA_raw, 'same') + bkg_comp_i
#                 )
#             else:
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                         scipy.signal.fftconvolve(model_temp, PSF_DATA_raw, 'same') + bkg_comp_i
#                 )
#         else:
#             model_dict['model_c' + str(i) + '_conv'] = model_temp + bkg_comp_i
#
#     #     model = model
#     model_dict['model_total_dec'] = np.asarray(model + flat_sky_total_dec)  # +FlatSky_cpu(background,
#     # params['s_a'].value)
#
#     if PSF_CONV == True:
#         # model_dict['model_total_conv'] = scipy.signal.fftconvolve(model,
#         #                                                           PSF_DATA_raw,
#         #                                                           'same')  # + FlatSky(FlatSky_level, params['s_a'])
#         if convolution_mode == 'GPU':
#             # model_dict['model_total_conv'] = np.asarray(jax_convolve(model,
#             #                                                          PSF_DATA)).copy()
#             # model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + FlatSky_cpu(background, params['s_a'].value
#             # model_conv = _fftconvolve_jax(model+flat_sky_total,PSF_DATA).copy()
#             model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + flat_sky_total
#         elif convolution_mode == 'CPU':
#             model_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + flat_sky_total
#         else:
#             model_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + flat_sky_total
#
#         model_dict['model_total_conv'] = model_conv
#     else:
#         model_dict['model_total_conv'] = model + flat_sky_total
#
#     return model_dict, flat_sky_total, flat_sky_total_dec, bkg_comp_i, bkg_comp_i_dec
#
#
# def _compute_residuals(data_2D, model_dict, flat_sky_total, is_bkg_map_conv,
#                        PSF_CONV, PSF_DATA, convolution_mode,
#                        np, _fftconvolve_jax):
#     """
#     Compute residual images and background components.
#
#     Parameters
#     ----------
#     data_2D : array
#         Original image data
#     model_dict : dict
#         Dictionary of model components
#     flat_sky_total : array
#         Total flat sky background (convolved)
#     is_bkg_map_conv : bool
#         Whether RMS map is already convolved
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array
#         PSF data
#     convolution_mode : str
#         'GPU' or 'CPU'
#     np : module
#         numpy module
#     _fftconvolve_jax : callable
#         JAX FFT convolution
#
#     Returns
#     -------
#     model_dict : dict
#         Updated model dictionary with residuals
#     """
#     # model_dict['best_residual'] = data_2D - model_dict['model_total']
#     # bkg_comp_total
#     model_dict['model_total_conv'] = np.asarray(model_dict['model_total_conv'])
#
#     if is_bkg_map_conv is False and PSF_CONV is True:
#         model_dict['best_residual_conv'] = np.asarray(data_2D) - model_dict['model_total_conv']
#         if convolution_mode == 'GPU':
#             model_dict['conv_bkg'] = np.asarray(_fftconvolve_jax(flat_sky_total, PSF_DATA).copy())
#         else:
#             model_dict['conv_bkg'] = np.asarray(flat_sky_total)
#     else:
#         model_dict['best_residual_conv'] = np.asarray(data_2D) - model_dict['model_total_conv'] + flat_sky_total
#         model_dict['conv_bkg'] = np.asarray(flat_sky_total)
#
#     return model_dict
#
#
# def _save_fitting_outputs(imagename, model_dict, ncomponents, special_name,
#                           save_name_append, flat_sky_total_dec,
#                           pf, copy_header):
#     """
#     Save all FITS files and generate output file lists.
#
#     Parameters
#     ----------
#     imagename : str
#         Base image filename
#     model_dict : dict
#         Dictionary containing all model arrays
#     ncomponents : int
#         Number of Sersic components
#     special_name : str
#         Special identifier for output files
#     save_name_append : str
#         Additional string to append to filenames
#     flat_sky_total_dec : array
#         Deconvolved flat sky background
#     pf : module
#         astropy.io.fits module
#     copy_header : callable
#         Function to copy FITS headers
#
#     Returns
#     -------
#     image_results_conv : list
#         List of convolved output filenames
#     image_results_deconv : list
#         List of deconvolved output filenames
#     bkg_images : list
#         List of background image filenames
#     total_image_results_conv : list
#         List of total convolved model filenames
#     total_image_results_deconv : list
#         List of total deconvolved model filenames
#     """
#     image_results_conv = []
#     image_results_deconv = []
#     total_image_results_conv = []
#     total_image_results_deconv = []
#     bkg_images = []
#
#     # Save individual components
#     for i in range(1, ncomponents + 1):
#         # Convolved component
#         conv_filename = (imagename.replace('.fits', '') +
#                          "_" + "model_component_" + str(i) +
#                          special_name + save_name_append + '.fits')
#         pf.writeto(conv_filename, model_dict['model_c' + str(i) + '_conv'],
#                    overwrite=True)
#         copy_header(imagename, conv_filename, conv_filename)
#         image_results_conv.append(conv_filename)
#
#         # Deconvolved component
#         dec_filename = (imagename.replace('.fits', '') +
#                         "_" + "dec_model_component_" + str(i) +
#                         special_name + save_name_append + '.fits')
#         pf.writeto(dec_filename, model_dict['model_c' + str(i)],
#                    overwrite=True)
#         copy_header(imagename, dec_filename, dec_filename)
#         image_results_deconv.append(dec_filename)
#
#     # Save total convolved model
#     conv_model_filename = (imagename.replace('.fits', '') +
#                            "_" + "conv_model" + special_name + save_name_append + '.fits')
#     pf.writeto(conv_model_filename, model_dict['model_total_conv'], overwrite=True)
#     copy_header(imagename, conv_model_filename, conv_model_filename)
#     total_image_results_conv.append(conv_model_filename)
#     image_results_conv.append(conv_model_filename)
#
#     # Save total deconvolved model
#     dec_model_filename = (imagename.replace('.fits', '') +
#                           "_" + "dec_model" + special_name + save_name_append + '.fits')
#     pf.writeto(dec_model_filename, model_dict['model_total_dec'], overwrite=True)
#     copy_header(imagename, dec_model_filename, dec_model_filename)
#     total_image_results_deconv.append(dec_model_filename)
#     image_results_deconv.append(dec_model_filename)
#
#     # Save residual
#     residual_filename = (imagename.replace('.fits', '') +
#                          "_" + "residual" + special_name + save_name_append + ".fits")
#     pf.writeto(residual_filename, model_dict['best_residual_conv'], overwrite=True)
#     copy_header(imagename, residual_filename, residual_filename)
#     image_results_conv.append(residual_filename)
#
#     # Save deconvolved background
#     model_dict['deconv_bkg'] = np.asarray(flat_sky_total_dec)
#     deconv_bkg_filename = (imagename.replace('.fits', '') +
#                            "_" + "deconv_bkg" + special_name + save_name_append + '.fits')
#     pf.writeto(deconv_bkg_filename, model_dict['deconv_bkg'], overwrite=True)
#     copy_header(imagename, deconv_bkg_filename, deconv_bkg_filename)
#     bkg_images.append(deconv_bkg_filename)
#
#     # Save convolved background
#     conv_bkg_filename = (imagename.replace('.fits', '') +
#                          "_" + "conv_bkg" + special_name + save_name_append + '.fits')
#     pf.writeto(conv_bkg_filename, model_dict['conv_bkg'], overwrite=True)
#     copy_header(imagename, conv_bkg_filename, conv_bkg_filename)
#     bkg_images.append(conv_bkg_filename)
#
#     # # initial minimization.
#     # method1 = 'differential_evolution'
#     # print(' >> Using', method1, ' solver for first optimisation run... ')
#     # # take parameters from previous run, and re-optimize them.
#     # #     method2 = 'ampgo'#'least_squares'
#     # method2 = 'least_squares'
#
#     # # save mini results (full) to a pickle file.
#     # with open(imagename.replace('.fits',
#     #                             '_' + 'fit' +
#     #                             special_name + save_name_append + '.pickle'),
#     #           "wb") as f:
#     #     pickle.dump(result, f)
#
#     # with open(imagename.replace('.fits',
#     #                             '_' + 'fit' +
#     #                             special_name + save_name_append + '_modeldict.pickle'),
#     #           "wb") as f:
#     #     pickle.dump(model_dict, f)
#
#     return (image_results_conv, image_results_deconv, bkg_images,
#             total_image_results_conv, total_image_results_deconv)
#
#
# # ============================================================================
# # MAIN FUNCTION
# # ============================================================================
#
# def do_fit2D(imagename, params_values_init_IMFIT=None, ncomponents=None,
#              init_constraints=None, data_2D_=None, residualdata_2D_=None,
#              residualname=None, which_residual='shuffled', observation_type='radio',
#              init_params=0.25, final_params=4.0, constrained=True,
#              fix_n=True, fix_value_n=False,
#              fix_max_value_Rn=False, fix_min_value_Rn=False,
#              dr_fix=3,
#              fix_x0_y0=False, psf_name=None, convolution_mode='GPU',
#              convolve_cutout=False, cut_size=512, self_bkg=False,
#              rms_map=None, is_bkg_map_conv=False,
#              fix_geometry=True, contrain_nelder=False, workers=6, mask_region=None,
#              special_name='', method1='least_squares', method2='least_squares',
#              reduce_fcn='neglogcauchy', loss="cauchy", tr_solver="exact", x_scale='jac',
#              ftol=1e-10, xtol=1e-10, gtol=1e-10, verbose=2, max_nfev=200000,
#              regularize=True, f_scale=1.0,
#              maxiter=30000, maxfev=30000, xatol=1e-12,
#              fatol=1e-12, return_all=True, disp=True,
#              de_options=None, parameters_mini_init=None,
#              save_name_append='', logger=None):
#     """
#     Perform a Robust and Fast Multi-Sersic Decomposition with GPU acceleration.
#     tr_solver:
#
#     Parameters
#     ----------
#     imagename: str
#         Name of the image to be fitted.
#     params_values_init_IMFIT: list
#         Initial parameters values for the model.
#     ncomponents: int
#         Number of components to be fitted.
#     init_constraints: dict
#         Initial constraints for the model.
#     data_2D_: 2D array
#         Image to be fitted.
#     residualdata_2D_: 2D array
#         Residual image to be fitted.
#     residualname: str
#         Name of the residual image to be fitted.
#     which_residual: str
#         Which residual to be used for the fitting.
#         Options: 'shuffled' or 'natural'.
#     init_params: float
#         Initial parameters for the model.
#     final_params: float
#         Final parameters for the model.
#     constrained: bool
#         If True, use initial constraints for the model.
#     fix_n: bool
#         If True, fix the Sersic index of the model.
#     fix_value_n: float
#         If True, fix the Sersic index of the model to this value.
#     dr_fix: float
#         If True, fix the centre position of the model.
#     fix_x0_y0: bool
#         If True, fix the centre position of the model.
#     psf_name: str
#         Name of the PSF image to be used for the convolution.
#     convolution_mode: str
#         If 'GPU', use GPU acceleration for the convolution.
#     convolve_cutout: bool
#         If True, convolve the image cutout with the PSF.
#     cut_size: int
#         Size of the cutout image.
#     self_bkg: bool
#         If True, use the image background as the residual background.
#     rms_map: 2D array
#         RMS map to be used for the fitting.
#     fix_geometry: bool
#         If True, fix the geometry of the model.
#     contrain_nelder: bool
#         If True, constrain the Nelder-Mead optimised parameters.
#     workers: int
#         Number of workers to be used for the fitting.
#     mask_region: 2D array
#         Mask to be used for the fitting.
#     special_name: str
#         Special name to be used for the output files.
#     method1: str
#         Method to be used for the fitting.
#     method2: str
#         Method to be used for the fitting.
#     reduce_fcn: str
#
#     loss: str
#
#     tr_solver: str
#
#     x_scale: str
#
#     ftol: float
#
#     xtol: float
#
#     gtol: float
#
#     verbose: int
#
#     max_nfev: int
#
#     regularize: bool
#         If True, regularize the model.
#     f_scale: float
#
#     maxiter: int
#
#     maxfev: int
#
#     xatol: float
#
#     fatol: float
#
#     return_all: bool
#
#     disp: bool
#
#     de_options: dict
#
#     save_name_append: str
#
#     logger: logger
#
#
#     returns
#     -------
#     result: dict
#         Dictionary containing the results of the fitting.
#
#
#     """
#     # Check JAX availability and adjust convolution mode if needed
#     try:
#         from jax import jit
#     except:
#         convolution_mode = 'CPU'
#
#     # ========================================================================
#     # INITIALIZATION
#     # ========================================================================
#
#     if de_options is None:
#         de_options = {'disp': True, 'workers': 6,
#                       'max_nfev': 20000, 'vectorized': True,
#                       # 'strategy': 'randtobest1bin',
#                       'mutation': (0.5, 1.5),
#                       'recombination': [0.2, 0.9],
#                       'init': 'random', 'tol': 0.00001,
#                       'updating': 'deferred',
#                       'popsize': 600}
#
#     startTime = time.time()
#     try:
#         logger.info(f"Fitting image: {imagename}")
#     except:
#         pass
#
#     # ========================================================================
#     # DATA LOADING AND PREPARATION
#     # ========================================================================
#
#     # Load image data
#     data_2D, data_2D_gpu = _initialize_image_data(imagename, data_2D_,
#                                                   convolution_mode, pf, jnp, logger)
#
#     # Load PSF data
#     PSF_CONV, PSF_DATA, PSF_DATA_raw = _initialize_psf_data(psf_name,
#                                                             convolution_mode,
#                                                             pf, jnp, load_fits_data, logger)
#
#     # Prepare mask
#     mask_for_fit = _prepare_mask(mask_region, convolution_mode, jnp, logger)
#
#     # Prepare background and residual
#     background, background_dec, residual_2D, FlatSky_level = \
#         _prepare_background_and_residual(
#             residualname, residualdata_2D_, which_residual,
#             self_bkg, rms_map, imagename, data_2D,
#             PSF_CONV, PSF_DATA, PSF_DATA_raw, convolution_mode,
#             is_bkg_map_conv,
#             pf, jnp, scipy, mad_std, shuffle_2D, sep_background,
#             _fftconvolve_jax, logger
#         )
#
#     # Convert residual to GPU format if needed
#     if convolution_mode == 'GPU':
#         residual_2D = jnp.array(residual_2D)
#
#     # Setup coordinate grid
#     xy, size = _setup_coordinate_grid(data_2D, convolution_mode, np, jnp)
#
#     if convolve_cutout is True:
#         """
#         WARNING: DO NOT USE FOR NOW!
#
#         Instead of convolving the entire image,
#         convolve only a box.
#         Can be 10x faster.
#
#         Issue:
#         It causes the flat sky level to be much higher than the real value.
#
#         Need further investigation and proper implementation.
#         """
#         x0c, y0c = int(size[0] / 2), int(size[1] / 2)
#
#     #     FlatSky_level = background#mad_std(data_2D)
#
#     # if convolution_mode == 'GPU':
#
#     # FlatSky_level = mad_std(data_2D)
#
#     # ========================================================================
#     # MODEL CONSTRUCTION
#     # ========================================================================
#
#     nfunctions = ncomponents
#
#     # Create minimizer functions
#     min_residual_2D, min_residual_2D_GPU, build_model, func = \
#         _create_minimizer_functions(nfunctions, xy, data_2D, data_2D_gpu,
#                                     background, residual_2D, PSF_DATA, PSF_DATA_raw,
#                                     mask_for_fit, convolution_mode,
#                                     np, jnp, scipy, jit,
#                                     sersic2D, sersic2D_GPU, FlatSky, FlatSky_cpu,
#                                     _fftconvolve_jax)
#
#     # Setup JAX convolution if needed
#     jax_convolve = _setup_jax_convolution(convolution_mode, jnp, jit)
#
#     # Construct model parameters
#     smodel2D, params = construct_model_parameters(
#         params_values_init_IMFIT=params_values_init_IMFIT, n_components=nfunctions,
#         init_constraints=init_constraints, observation_type=observation_type,
#         fix_n=fix_n, fix_value_n=fix_value_n,
#         fix_max_value_Rn=fix_max_value_Rn, fix_min_value_Rn=fix_min_value_Rn,
#         fix_x0_y0=fix_x0_y0, dr_fix=dr_fix, fix_geometry=fix_geometry,
#         init_params=init_params, final_params=final_params,
#         constrained=constrained)
#
#     # Create minimizer
#     if convolution_mode == 'CPU':
#         mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#     elif convolution_mode == 'GPU':
#         mini = lmfit.Minimizer(min_residual_2D_GPU, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#     else:
#         # Fallback to CPU
#         mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#
#     # ========================================================================
#     # OPTIMIZATION
#     # ========================================================================
#
#     result, result_1, result_extra = _perform_optimization(
#         mini, params, method1, method2, convolution_mode,
#         contrain_nelder, workers, max_nfev, tr_solver,
#         regularize, x_scale, f_scale, ftol, xtol, gtol,
#         verbose, loss, maxiter, maxfev, xatol, fatol,
#         return_all, disp, de_options, parameters_mini_init,
#         constrain_nelder_mead_params
#     )
#
#     params = result.params
#
#     # ========================================================================
#     # MODEL GENERATION AND OUTPUT
#     # ========================================================================
#
#     # Generate model components
#     model_dict, flat_sky_total, flat_sky_total_dec, bkg_comp_i, bkg_comp_i_dec = \
#         _generate_model_components(params, ncomponents, xy, size, background,
#                                    background_dec, PSF_CONV, PSF_DATA, PSF_DATA_raw,
#                                    convolution_mode, is_bkg_map_conv,
#                                    np, scipy, Model, sersic2D_GPU, FlatSky, FlatSky_cpu,
#                                    _fftconvolve_jax)
#
#     # Compute residuals
#     model_dict = _compute_residuals(data_2D, model_dict, flat_sky_total,
#                                     is_bkg_map_conv, PSF_CONV, PSF_DATA,
#                                     convolution_mode,
#                                     np, _fftconvolve_jax)
#
#     # Save outputs
#     (image_results_conv, image_results_deconv, bkg_images,
#      total_image_results_conv, total_image_results_deconv) = \
#         _save_fitting_outputs(imagename, model_dict, ncomponents, special_name,
#                               save_name_append, flat_sky_total_dec,
#                               pf, copy_header)
#
#     # ========================================================================
#     # FINAL PROCESSING
#     # ========================================================================
#
#     exec_time = time.time() - startTime
#     print('Exec time fitting=', exec_time, 's')
#
#     # save results to csv file.
#     try:
#         save_results_csv(result_mini=result,
#                          save_name=image_results_conv[-2].replace('.fits', ''),
#                          ext='.csv',
#                          save_corr=True, save_params=True)
#     except:
#         print('Error Saving Results to a csv file!!!')
#         pass
#
#     # Extract model_temp for return (last component generated)
#     model_temp = Model(sersic2D_GPU)
#
#     return (result, mini, result_1, result_extra, model_dict, image_results_conv,
#             image_results_deconv, bkg_images, smodel2D, model_temp)


# # ============================================================================
# # HELPER FUNCTIONS - Data Initialization
# # ============================================================================
#
# def _initialize_image_data(imagename, data_2D_, convolution_mode, logger=None):
#     """
#     Load and prepare image data for fitting.
#
#     Parameters
#     ----------
#     imagename : str
#         Path to image file
#     data_2D_ : array or None
#         Pre-loaded image data
#     convolution_mode : str
#         'GPU' or 'CPU'
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     data_2D : array
#         Image data in numpy format
#     data_2D_gpu : array or None
#         Image data in JAX format (if GPU mode)
#     """
#     if data_2D_ is None:
#         data_2D = pf.getdata(imagename)
#     else:
#         data_2D = data_2D_
#
#     # Prepare GPU array if needed
#     if convolution_mode == 'GPU':
#         data_2D_gpu = jnp.array(data_2D)
#     else:
#         data_2D_gpu = None
#
#     return data_2D, data_2D_gpu
#
#
# def _initialize_psf_data(psf_name, convolution_mode, logger=None):
#     """
#     Load and prepare PSF data for convolution.
#
#     Parameters
#     ----------
#     psf_name : str or None
#         Path to PSF file
#     convolution_mode : str
#         'GPU' or 'CPU'
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array or None
#         PSF data in appropriate format for computation
#     PSF_DATA_raw : array or None
#         PSF data in raw numpy format
#     """
#     if psf_name is not None:
#         PSF_CONV = True
#         try:
#             PSF_DATA_raw = pf.getdata(psf_name)
#             if len(PSF_DATA_raw.shape) == 4:
#                 PSF_DATA_raw = PSF_DATA_raw[0][0]
#         except:
#             PSF_DATA_raw = load_fits_data(psf_name)
#
#         if convolution_mode == 'GPU':
#             if logger is not None:
#                 logger.debug(f"---------------------------------------")
#                 logger.debug(f" <<< PERFORMING CONVOLUTION WITH JAX >>> ")
#                 logger.debug(f"---------------------------------------")
#             PSF_DATA = jnp.array(PSF_DATA_raw)
#         elif convolution_mode == 'CPU':
#             PSF_DATA = PSF_DATA_raw
#         else:
#             PSF_DATA = PSF_DATA_raw
#
#         # PSF_DATA = pf.getdata(
#         #     imagename.replace('-image.cutout.fits', '-beampsf.cutout.fits'))
#     else:
#         PSF_CONV = False
#         PSF_DATA = None
#         PSF_DATA_raw = None
#
#     return PSF_CONV, PSF_DATA, PSF_DATA_raw
#
#
# def _prepare_mask(mask_region, convolution_mode, logger=None):
#     """
#     Prepare mask for constrained fitting.
#
#     Parameters
#     ----------
#     mask_region : array or None
#         Mask array
#     convolution_mode : str
#         'GPU' or 'CPU'
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     mask_for_fit : array or None
#         Mask in appropriate format
#     """
#     if mask_region is not None:
#         """
#
#         """
#         if logger is not None:
#             logger.debug(f" ==> Using provided mask region to constrain fit. ")
#             logger.warning(f" !!++==> Fitting with a mask is faster, but experimental!! \n"
#                            f"         Use with caution.")
#         # data_2D = data_2D * mask_region
#         if convolution_mode == 'GPU':
#             mask_for_fit = jnp.array(mask_region)
#         else:
#             mask_for_fit = mask_region
#     else:
#         mask_for_fit = None
#
#     return mask_for_fit
#
#
# def _prepare_background_and_residual(residualname, residualdata_2D_, which_residual,
#                                      self_bkg, rms_map, imagename, data_2D,
#                                      PSF_CONV, PSF_DATA, PSF_DATA_raw, convolution_mode,
#                                      is_bkg_map_conv, logger=None):
#     """
#     Prepare background and residual arrays for fitting.
#
#     This function handles the complex logic of preparing background maps,
#     residual images, and RMS estimates based on various input configurations.
#
#     Parameters
#     ----------
#     residualname : str or None
#         Path to residual image
#     residualdata_2D_ : array or None
#         Pre-loaded residual data
#     which_residual : str
#         Type of residual: 'shuffled', 'natural', or 'user'
#     self_bkg : bool
#         Whether to use self-background estimation
#     rms_map : array or None
#         User-provided RMS/background map
#     imagename : str
#         Path to main image
#     data_2D : array
#         Main image data
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array or None
#         PSF data in computation format
#     PSF_DATA_raw : array or None
#         PSF data in raw format
#     convolution_mode : str
#         'GPU' or 'CPU'
#     is_bkg_map_conv : bool
#         Whether the RMS map is already convolved
#     logger : logger, optional
#         Logger instance
#
#     Returns
#     -------
#     background : array
#         Background array for fitting (possibly convolved)
#     background_dec : array
#         Deconvolved background array
#     residual_2D : array
#         Residual data
#     FlatSky_level : float or None
#         Flat sky level estimate
#     """
#     # from signal_stats import mad_std, shuffle_2D, sep_background
#
#     FlatSky_level = None
#
#     if residualname is not None and which_residual != 'user':
#         """
#         This is important for radio image fitting.
#
#         It uses the shuffled version of the residual cleaned image
#         originated from the interferometric deconvolution.
#
#         This ensures that the best model created here will be on top
#         of that rms noise so that flux conservation is maximized.
#
#         However, this residual is not added as model + shuffled_residual
#         only, but instead by a multiplication factor,
#         e.g. model + const* shuffled_residual, and const will be minimized
#         as well during the fitting (here, called `s_a`).
#         """
#         if residualdata_2D_ is not None:
#             residual_2D = residualdata_2D_
#         else:
#             residual_2D = pf.getdata(residualname)
#
#         if which_residual == 'shuffled':
#             if logger is not None:
#                 logger.debug(f" ==> Using clean shuffled background for optmization... ")
#             residual_2D_to_use = shuffle_2D(residual_2D)
#         elif which_residual == 'natural':
#             if logger is not None:
#                 logger.debug(f" ==> Using clean background for optmization... ")
#             """
#             if psf_name is not None:
#                 if logger is not None:
#                     logger.debug(f" ====> Deconvolving residual map... ")
#                 residual_2D_to_use, _ = deconvolve_fft(residual_2D,
#                                                             PSF_DATA_raw/PSF_DATA_raw.sum())
#             else:
#                 residual_2D_to_use = residual_2D
#             """
#             residual_2D_to_use = residual_2D
#         else:
#             residual_2D_to_use = residual_2D
#
#         FlatSky_level = mad_std(residual_2D_to_use)
#         #         background = residual_2D #residual_2D_to_use
#         if convolution_mode == 'GPU':
#             background = jnp.array(residual_2D_to_use)
#         else:
#             background = residual_2D_to_use
#
#     else:
#         if which_residual == 'user':
#             if rms_map is None:
#                 print('--==>> A rms map/background mode was selected (user)')
#                 print('       but no rms/background map was provided.')
#                 print('       Please, provide a rms/background map.')
#                 print('||==>> Stopping code now.')
#                 raise ValueError("rms_map should not be None")
#             else:
#                 if logger is not None:
#                     logger.debug(f" ==> Using provided RMS map. ")
#                 background_map = rms_map
#                 background = background_map.copy()
#         else:
#             if self_bkg == True:
#                 if logger is not None:
#                     logger.warning(f" ==> No residual/background provided. Using image bkg map... ")
#                 background_map = sep_background(imagename)
#                 background = shuffle_2D(background_map.back())
#             else:
#                 if logger is not None:
#                     logger.warning(f" ==> Using only flat sky for rms bkg.")
#                 FlatSky_level = mad_std(data_2D)
#                 background = FlatSky_level
#
#     # Load residual_2D if not already loaded
#     if residualdata_2D_ is not None:
#         residual_2D = residualdata_2D_
#     else:
#         try:
#             residual_2D = pf.getdata(residualname)
#         except:
#             residual_2D = background
#
#     # Create deconvolved background copy before convolution
#     if isinstance(background, (int, float)):
#         # FlatSky_level case - background is a scalar
#         background_dec = background
#     else:
#         background_dec = background.copy()
#
#     # Convolve background if needed and PSF is available
#     if is_bkg_map_conv is False and PSF_CONV:
#         if logger is not None:
#             logger.debug(f" ==> RMS map is not convolved.")
#
#         if not isinstance(background, (int, float)):
#             # Only convolve if background is an array
#             if convolution_mode == 'GPU':
#                 background = _fftconvolve_jax(background, PSF_DATA)
#             elif convolution_mode == 'CPU':
#                 background = scipy.signal.fftconvolve(background, PSF_DATA, 'same')
#
#     return background, background_dec, residual_2D, FlatSky_level
#
#
# def _setup_coordinate_grid(data_2D, convolution_mode):
#     """
#     Create coordinate meshgrid for model evaluation.
#
#     Parameters
#     ----------
#     data_2D : array
#         Image data (for shape)
#     convolution_mode : str
#         'GPU' or 'CPU'
#
#     Returns
#     -------
#     xy : tuple of arrays
#         Coordinate meshgrid
#     size : tuple
#         Image shape
#     """
#     size = data_2D.shape
#     if convolution_mode == 'GPU':
#         x, y = jnp.meshgrid(jnp.arange((size[1])), jnp.arange((size[0])))
#         xy = jnp.stack([x, y], axis=0)
#     else:
#         xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
#
#     return xy, size
#
#
# # ============================================================================
# # HELPER FUNCTIONS - Model Construction
# # ============================================================================
#
# def _create_minimizer_functions(nfunctions, xy, data_2D, data_2D_gpu, background,
#                                 residual_2D, PSF_DATA, PSF_DATA_raw, mask_for_fit,
#                                 convolution_mode):
#     """
#     Create the residual functions for minimization.
#
#     Parameters
#     ----------
#     nfunctions : int
#         Number of Sersic components
#     xy : tuple of arrays
#         Coordinate meshgrid
#     data_2D : array
#         Image data (CPU)
#     data_2D_gpu : array or None
#         Image data (GPU)
#     background : array
#         Background map
#     residual_2D : array
#         Residual data
#     PSF_DATA : array
#         PSF data in computation format
#     PSF_DATA_raw : array
#         PSF data in raw format
#     mask_for_fit : array or None
#         Fitting mask
#     convolution_mode : str
#         'GPU' or 'CPU'
#
#     Returns
#     -------
#     min_residual_2D : callable or None
#         CPU residual function
#     min_residual_2D_GPU : callable or None
#         GPU residual function
#     build_model : callable or None
#         Model building function (for GPU)
#     func : callable or None
#         Parameter splitting function (for GPU)
#     """
#     # CPU version
#     def min_residual_2D(params):
#         dict_model = {}
#         model = 0
#         for i in range(1, nfunctions + 1):
#             model = model + sersic2D(xy, params['f' + str(i) + '_x0'],
#                                      params['f' + str(i) + '_y0'],
#                                      params['f' + str(i) + '_PA'],
#                                      params['f' + str(i) + '_ell'],
#                                      params['f' + str(i) + '_n'],
#                                      params['f' + str(i) + '_In'],
#                                      params['f' + str(i) + '_Rn'],
#                                      params['f' + str(i) + '_cg'], )
#         # print(model.shape)
#         # model = model + FlatSky_cpu(FlatSky_level, params['s_a'])*background
#         # model = model + FlatSky_cpu(background, params['s_a'])
#         MODEL_2D_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + \
#                         FlatSky_cpu(background, params['s_a'])
#         residual = data_2D - MODEL_2D_conv
#         return np.ravel(residual)
#
#     # GPU version
#     try:
#         # @partial(jit, static_argnums=1)
#         @jit
#         def func(x):
#             return jnp.split(x, nfunctions)
#
#         @jit
#         def build_model(xy, param_matrix):
#             model = 0
#             for model_params in param_matrix:
#                 model = model + sersic2D_GPU(xy, model_params[0],
#                                              model_params[1],
#                                              model_params[2],
#                                              model_params[3],
#                                              model_params[4],
#                                              model_params[5],
#                                              model_params[6],
#                                              model_params[7])
#             return model
#
#     except:
#         def func(x):
#             return jnp.split(x, nfunctions)
#
#         def build_model(xy, param_matrix):
#             model = 0
#             for model_params in param_matrix:
#                 model = model + sersic2D_GPU(xy, model_params[0],
#                                              model_params[1],
#                                              model_params[2],
#                                              model_params[3],
#                                              model_params[4],
#                                              model_params[5],
#                                              model_params[6],
#                                              model_params[7])
#             return model
#
#     def min_residual_2D_GPU(params):
#         model = 0
#         for i in range(1, nfunctions + 1):
#             model = model + sersic2D_GPU(xy,
#                                          params['f' + str(i) + '_x0'].value,
#                                          params['f' + str(i) + '_y0'].value,
#                                          params['f' + str(i) + '_PA'].value,
#                                          params['f' + str(i) + '_ell'].value,
#                                          params['f' + str(i) + '_n'].value,
#                                          params['f' + str(i) + '_In'].value,
#                                          params['f' + str(i) + '_Rn'].value,
#                                          params['f' + str(i) + '_cg'].value)
#
#         # # param_matrix = extract_params(params)
#         # param_matrix = func(jnp.array(list(params.valuesdict().values()))[:-1])
#         # model = build_model(xy,param_matrix)
#
#         MODEL_2D_conv = _fftconvolve_jax(model, PSF_DATA) + FlatSky(background, params['s_a'].value)
#
#         # MODEL_2D_conv = _fftconvolve_jax(model+
#         #                                  FlatSky(background,params['s_a'].value),
#         #                                  PSF_DATA)
#         # residual = ((data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])/
#         #             (1000*(abs(residual_2D[mask_for_fit])+1.0e-6)))
#         residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])
#         return np.asarray(residual).copy()
#         # weights = 1/((background[mask_for_fit])/data_2D_gpu[mask_for_fit])
#         # weightned_residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit]) * jnp.sqrt(weights)
#         # return np.asarray(weightned_residual).copy()
#
#     if convolution_mode == 'CPU':
#         return min_residual_2D, None, build_model, func
#     else:
#         return None, min_residual_2D_GPU, build_model, func
#
#
# def _setup_jax_convolution(convolution_mode):
#     """
#     Setup JAX convolution function if in GPU mode.
#
#     Parameters
#     ----------
#     convolution_mode : str
#         'GPU' or 'CPU'
#
#     Returns
#     -------
#     jax_convolve : callable or None
#         JIT-compiled convolution function
#     """
#     if convolution_mode == 'GPU':
#         @jit
#         def convolve_on_gpu(image, psf):
#             """
#             This was before jax.scipy implementing fftconvolve.
#             It provides the same result, at the same speed.
#
#             This function also accepts PSFs with a different shape of the image.
#
#             """
#             # Calculate the new padded shape
#             padded_shape = (image.shape[0] + psf.shape[0] - 1,
#                             image.shape[1] + psf.shape[1] - 1)
#
#             # Pad both image and psf to the new shape
#             pad_shape = [(0, ts - s) for s, ts in zip(image.shape, padded_shape)]
#             image_padded = jnp.pad(image, pad_shape, mode='constant')
#             pad_shape = [(0, ts - s) for s, ts in zip(psf.shape, padded_shape)]
#             psf_padded = jnp.pad(psf, pad_shape, mode='constant')
#             # psf_padded = pad_for_convolution(psf, padded_shape)
#             image_fft = jnp.fft.fft2(image_padded)
#             psf_fft = jnp.fft.fft2(psf_padded)
#
#             conv_fft = image_fft * psf_fft
#
#             # Get the real part of the inverse FFT and crop to the original image size
#             result_full = jnp.real(jnp.fft.ifft2(conv_fft))
#             return result_full[psf.shape[0] // 2:image.shape[0] + psf.shape[0] // 2,
#             psf.shape[1] // 2:image.shape[1] + psf.shape[1] // 2]
#
#         jax_convolve = jit(convolve_on_gpu)
#         return jax_convolve
#     return None
#
#
# # ============================================================================
# # HELPER FUNCTIONS - Optimization
# # ============================================================================
#
# def _perform_optimization(mini, params, method1, method2, convolution_mode,
#                           contrain_nelder, workers, max_nfev, tr_solver,
#                           regularize, x_scale, f_scale, ftol, xtol, gtol,
#                           verbose, loss, maxiter, maxfev, xatol, fatol,
#                           return_all, disp, de_options, parameters_mini_init):
#     """
#     Execute the two-stage optimization process.
#
#     Parameters
#     ----------
#     mini : lmfit.Minimizer
#         Minimizer object
#     params : lmfit.Parameters
#         Initial parameters
#     method1 : str
#         First optimization method
#     method2 : str
#         Second optimization method
#     convolution_mode : str
#         'GPU' or 'CPU'
#     contrain_nelder : bool
#         Whether to constrain Nelder-Mead parameters
#     workers : int
#         Number of parallel workers
#     max_nfev : int
#         Maximum function evaluations
#     tr_solver : str
#         Trust region solver
#     regularize : bool
#         Whether to regularize
#     x_scale : str or array
#         Parameter scaling
#     f_scale : float
#         Function scaling
#     ftol : float
#         Function tolerance
#     xtol : float
#         Parameter tolerance
#     gtol : float
#         Gradient tolerance
#     verbose : int
#         Verbosity level
#     loss : str
#         Loss function
#     maxiter : int
#         Maximum iterations for Nelder-Mead
#     maxfev : int
#         Maximum function evaluations for Nelder-Mead
#     xatol : float
#         Absolute parameter tolerance for Nelder-Mead
#     fatol : float
#         Absolute function tolerance for Nelder-Mead
#     return_all : bool
#         Return all Nelder-Mead results
#     disp : bool
#         Display convergence messages
#     de_options : dict
#         Options for differential evolution
#     parameters_mini_init : lmfit.Parameters or None
#         Initial parameters from previous run
#
#     Returns
#     -------
#     result : lmfit.MinimizerResult
#         Final optimization result
#     result_1 : lmfit.MinimizerResult
#         First stage optimization result
#     result_extra : None
#         Placeholder for additional results
#     """
#
#     # initial minimization.
#     print(' >> Using', method1, ' solver for first optimisation run... ')
#     # take parameters from previous run, and re-optimize them.
#     #     method2 = 'ampgo'#'least_squares'
#     #     method2 = 'least_squares'
#     result_extra = None
#
#     if method1 == 'nelder':
#         # very robust, but takes time....
#         #         print(' >> Using', method1,' solver for first optimisation run... ')
#         result_1 = mini.minimize(method='nelder',
#                                  #                                  xatol = 1e-12, fatol = 1e-12, disp=True,
#                                  #                                  adaptive = True,max_nfev = 30000,
#                                  options={'maxiter': maxiter, 'maxfev': maxfev,
#                                           'xatol': xatol, 'fatol': fatol,
#                                           'return_all': return_all,
#                                           'disp': disp}
#                                  )
#
#     elif method1 == 'least_squares':
#         # faster, but usually not good for first run.
#         # if results_previous_run is not None:
#         print(' >> Using', tr_solver, 'for tr solver, with regularize set to', regularize,
#               ' Loss is', loss, '.')
#
#         if parameters_mini_init is not None:
#             print(f'  ++==>> Using initial mini parameters from a previous run.')
#             # try:
#             result_1 = mini.minimize(method='least_squares',
#                                      params=parameters_mini_init,
#                                      max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
#                                      tr_solver=tr_solver,
#                                      tr_options={'regularize': regularize,
#                                                  },
#                                      ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                      loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#         else:
#             result_1 = mini.minimize(method='least_squares',
#                                      max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
#                                      tr_solver=tr_solver,
#                                      tr_options={'regularize': regularize,
#                                                  #                                              'min_delta':1e-14, 'eta':0.05,
#                                                  #                                              'xtol':1e-14, 'gtol':1e-14,
#                                                  #                                              'ftol':1e-14
#                                                  },
#                                      ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                      loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#
#     elif method1 == 'differential_evolution':
#         # de is giving some issues, I do not know why.
#         result_1 = mini.minimize(method='differential_evolution',
#                                  options={'disp': True, 'workers': workers,
#                                           'max_nfev': max_nfev, 'vectorized': True,
#                                           'strategy': 'randtobest1bin',
#                                           'mutation': (0.5, 1.5),
#                                           'recombination': [0.2, 0.9],
#                                           'init': 'random', 'tol': 0.00001,
#                                           'updating': 'deferred',
#                                           'popsize': 600})
#         # result_1 = mini.minimize(method='differential_evolution', popsize=600,
#         #                          disp=True,  # init = 'random',
#         #                          # mutation=(0.5, 1.5), recombination=[0.2, 0.9],
#         #                          max_nfev=20000,
#         #                          workers=1, updating='deferred', vectorized=True)
#     else:
#         raise ValueError(f"Unknown method1: {method1}")
#
#     print(' >> Using', method2, ' solver for second optimisation run... ')
#
#     second_run_params = result_1.params
#     if (contrain_nelder == True) and (method2 == 'nelder'):
#         """
#         It seems that least_squares is ignoring the best-parameters provided by
#         Nelder-mead, which means that it is lookig the parameter space far away
#         from the optimised Nelder-Mead ones.
#
#         So, with this condition, we force a much smaller searching region, but
#         it assumes that Nelder opt was good (which is not always true).
#
#         YOU MUST CHECK YOUR RESULTS!!!!
#
#         """
#         print('Constraining Nelder-Mead Parameters for method', method2)
#         params_constrained = constrain_nelder_mead_params(result_1.params,
#                                                           max_factor=1.03,
#                                                           min_factor=0.97)
#         # UPDATE THE SECOND RUN PARAMETERS TO BE THE CONSTRAINED ONES.
#         second_run_params = params_constrained
#
#     if method2 == 'nelder':
#         result = mini.minimize(method='nelder', params=second_run_params,
#                                options={'maxiter': maxiter, 'maxfev': maxfev,
#                                         'xatol': xatol, 'fatol': fatol,
#                                         'disp': disp})
#
#     elif method2 == 'ampgo':
#         # ampgo is not workin well/ takes so long ???
#         result = mini.minimize(method='ampgo', params=second_run_params,
#                                maxfunevals=10000, totaliter=30, disp=True,
#                                maxiter=5, glbtol=1e-8)
#
#     elif method2 == 'least_squares':
#         # faster, usually converges and provide errors.
#         # Very robust if used in second opt from first opt parameters.
#         result = mini.minimize(method='least_squares',
#                                params=second_run_params,
#                                max_nfev=max_nfev,
#                                tr_solver=tr_solver,
#                                tr_options={'regularize': regularize,
#                                            #                                            'min_delta': 1e-14, 'eta': 0.05,
#                                            #                                            'xtol': 1e-14, 'gtol': 1e-14,
#                                            #                                            'ftol': 1e-14
#                                            },
#                                x_scale=x_scale, f_scale=f_scale,
#                                ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#
#     elif method2 == 'differential_evolution':
#         # result = mini.minimize(method='differential_evolution',
#         #                        params=second_run_params,
#         #                        options={'maxiter': 30000, 'workers': -1,
#         #                                 'tol': 0.001, 'vectorized': True,
#         #                                 'strategy': 'randtobest1bin',
#         #                                 'updating': 'deferred', 'disp': True,
#         #                                 'seed': 1}
#         #                        )
#         result = mini.minimize(method='differential_evolution',
#                                params=second_run_params,
#                                options=de_options
#                                )
#     else:
#         raise ValueError(f"Unknown method2: {method2}")
#
#     return result, result_1, result_extra
#
#
# # ============================================================================
# # HELPER FUNCTIONS - Output Generation
# # ============================================================================
#
# def _generate_model_components(params, ncomponents, xy, size, background,
#                                background_dec, PSF_CONV, PSF_DATA, PSF_DATA_raw,
#                                convolution_mode, is_bkg_map_conv):
#     """
#     Generate individual and total model components.
#
#     Parameters
#     ----------
#     params : lmfit.Parameters
#         Optimized parameters
#     ncomponents : int
#         Number of Sersic components
#     xy : tuple of arrays
#         Coordinate meshgrid
#     size : tuple
#         Image shape
#     background : array
#         Background map (possibly convolved)
#     background_dec : array
#         Deconvolved background map
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array
#         PSF data in computation format
#     PSF_DATA_raw : array
#         PSF data in raw format
#     convolution_mode : str
#         'GPU' or 'CPU'
#     is_bkg_map_conv : bool
#         Whether RMS map is already convolved
#
#     Returns
#     -------
#     model_dict : dict
#         Dictionary containing all model arrays
#     flat_sky_total : array
#         Total flat sky background (convolved)
#     flat_sky_total_dec : array
#         Total flat sky background (deconvolved)
#     bkg_comp_i : array
#         Background component for individual models
#     bkg_comp_i_dec : array
#         Deconvolved background component for individual models
#     """
#
#     model_temp = Model(sersic2D_GPU)
#     xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
#     model = 0
#     model_dict = {}
#
#     if convolution_mode == 'GPU':
#         flat_sky_total = np.asarray(FlatSky(background, params['s_a'].value))
#         flat_sky_total_dec = np.asarray(FlatSky(background_dec, params['s_a'].value))
#     elif convolution_mode == 'CPU':
#         flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
#         flat_sky_total_dec = FlatSky_cpu(background_dec, params['s_a'].value)
#     else:
#         # Fallback
#         flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
#         flat_sky_total_dec = FlatSky_cpu(background_dec, params['s_a'].value)
#
#     bkg_comp_i = flat_sky_total.copy()
#     bkg_comp_i_dec = flat_sky_total_dec.copy()
#
#     for i in range(1, ncomponents + 1):
#         model_temp = sersic2D_GPU(xy, params['f' + str(i) + '_x0'].value,
#                                   params['f' + str(i) + '_y0'].value,
#                                   params['f' + str(i) + '_PA'].value,
#                                   params['f' + str(i) + '_ell'].value,
#                                   params['f' + str(i) + '_n'].value,
#                                   params['f' + str(i) + '_In'].value,
#                                   params['f' + str(i) + '_Rn'].value,
#                                   params['f' + str(i) + '_cg'].value)
#
#         model = model + model_temp
#         # to each individual component, add the bkg map.
#         model_dict['model_c' + str(i)] = np.asarray(model_temp + bkg_comp_i_dec)
#
#         if PSF_CONV == True:
#             if convolution_mode == 'GPU':
#                 # model_dict['model_c' + str(i) + '_conv'] = np.asarray(jax_convolve(model_temp, PSF_DATA)).copy()
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #     np.asarray(_fftconvolve_jax(model_temp,PSF_DATA).copy()+bkg_comp_i))
#                 # to each individual component, add the bkg map.
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #     np.asarray(_fftconvolve_jax(model_temp+bkg_comp_i,PSF_DATA).copy()))
#
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                                                                np.asarray(_fftconvolve_jax(model_temp,
#                                                                                            PSF_DATA).copy())) + bkg_comp_i
#
#             elif convolution_mode == 'CPU':
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #         scipy.signal.fftconvolve(model_temp+bkg_comp_i, PSF_DATA_raw,'same'))
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                         scipy.signal.fftconvolve(model_temp, PSF_DATA_raw, 'same') + bkg_comp_i
#                 )
#             else:
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                         scipy.signal.fftconvolve(model_temp, PSF_DATA_raw, 'same') + bkg_comp_i
#                 )
#         else:
#             model_dict['model_c' + str(i) + '_conv'] = model_temp + bkg_comp_i
#
#     #     model = model
#     model_dict['model_total_dec'] = np.asarray(model + flat_sky_total_dec)  # +FlatSky_cpu(background,
#     # params['s_a'].value)
#
#     if PSF_CONV == True:
#         # model_dict['model_total_conv'] = scipy.signal.fftconvolve(model,
#         #                                                           PSF_DATA_raw,
#         #                                                           'same')  # + FlatSky(FlatSky_level, params['s_a'])
#         if convolution_mode == 'GPU':
#             # model_dict['model_total_conv'] = np.asarray(jax_convolve(model,
#             #                                                          PSF_DATA)).copy()
#             # model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + FlatSky_cpu(background, params['s_a'].value
#             # model_conv = _fftconvolve_jax(model+flat_sky_total,PSF_DATA).copy()
#             model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + flat_sky_total
#         elif convolution_mode == 'CPU':
#             model_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + flat_sky_total
#         else:
#             model_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw, 'same') + flat_sky_total
#
#         model_dict['model_total_conv'] = model_conv
#     else:
#         model_dict['model_total_conv'] = model + flat_sky_total
#
#     return model_dict, flat_sky_total, flat_sky_total_dec, bkg_comp_i, bkg_comp_i_dec
#
#
# def _compute_residuals(data_2D, model_dict, flat_sky_total, is_bkg_map_conv,
#                        PSF_CONV, PSF_DATA, convolution_mode):
#     """
#     Compute residual images and background components.
#
#     Parameters
#     ----------
#     data_2D : array
#         Original image data
#     model_dict : dict
#         Dictionary of model components
#     flat_sky_total : array
#         Total flat sky background (convolved)
#     is_bkg_map_conv : bool
#         Whether RMS map is already convolved
#     PSF_CONV : bool
#         Whether PSF convolution is enabled
#     PSF_DATA : array
#         PSF data
#     convolution_mode : str
#         'GPU' or 'CPU'
#
#     Returns
#     -------
#     model_dict : dict
#         Updated model dictionary with residuals
#     """
#
#     # model_dict['best_residual'] = data_2D - model_dict['model_total']
#     # bkg_comp_total
#     model_dict['model_total_conv'] = np.asarray(model_dict['model_total_conv'])
#
#     if is_bkg_map_conv is False and PSF_CONV is True:
#         model_dict['best_residual_conv'] = np.asarray(data_2D) - model_dict['model_total_conv']
#         if convolution_mode == 'GPU':
#             model_dict['conv_bkg'] = np.asarray(_fftconvolve_jax(flat_sky_total, PSF_DATA).copy())
#         else:
#             model_dict['conv_bkg'] = np.asarray(flat_sky_total)
#     else:
#         model_dict['best_residual_conv'] = np.asarray(data_2D) - model_dict['model_total_conv'] + flat_sky_total
#         model_dict['conv_bkg'] = np.asarray(flat_sky_total)
#
#     return model_dict
#
#
# def _save_fitting_outputs(imagename, model_dict, ncomponents, special_name,
#                           save_name_append, flat_sky_total_dec):
#     """
#     Save all FITS files and generate output file lists.
#
#     Parameters
#     ----------
#     imagename : str
#         Base image filename
#     model_dict : dict
#         Dictionary containing all model arrays
#     ncomponents : int
#         Number of Sersic components
#     special_name : str
#         Special identifier for output files
#     save_name_append : str
#         Additional string to append to filenames
#     flat_sky_total_dec : array
#         Deconvolved flat sky background
#
#     Returns
#     -------
#     image_results_conv : list
#         List of convolved output filenames
#     image_results_deconv : list
#         List of deconvolved output filenames
#     bkg_images : list
#         List of background image filenames
#     total_image_results_conv : list
#         List of total convolved model filenames
#     total_image_results_deconv : list
#         List of total deconvolved model filenames
#     """
#
#     image_results_conv = []
#     image_results_deconv = []
#     total_image_results_conv = []
#     total_image_results_deconv = []
#     bkg_images = []
#
#     # Save individual components
#     for i in range(1, ncomponents + 1):
#         # Convolved component
#         conv_filename = (imagename.replace('.fits', '') +
#                          "_" + "model_component_" + str(i) +
#                          special_name + save_name_append + '.fits')
#         pf.writeto(conv_filename, model_dict['model_c' + str(i) + '_conv'],
#                    overwrite=True)
#         copy_header(imagename, conv_filename, conv_filename)
#         image_results_conv.append(conv_filename)
#
#         # Deconvolved component
#         dec_filename = (imagename.replace('.fits', '') +
#                         "_" + "dec_model_component_" + str(i) +
#                         special_name + save_name_append + '.fits')
#         pf.writeto(dec_filename, model_dict['model_c' + str(i)],
#                    overwrite=True)
#         copy_header(imagename, dec_filename, dec_filename)
#         image_results_deconv.append(dec_filename)
#
#     # Save total convolved model
#     conv_model_filename = (imagename.replace('.fits', '') +
#                            "_" + "conv_model" + special_name + save_name_append + '.fits')
#     pf.writeto(conv_model_filename, model_dict['model_total_conv'], overwrite=True)
#     copy_header(imagename, conv_model_filename, conv_model_filename)
#     total_image_results_conv.append(conv_model_filename)
#     image_results_conv.append(conv_model_filename)
#
#     # Save total deconvolved model
#     dec_model_filename = (imagename.replace('.fits', '') +
#                           "_" + "dec_model" + special_name + save_name_append + '.fits')
#     pf.writeto(dec_model_filename, model_dict['model_total_dec'], overwrite=True)
#     copy_header(imagename, dec_model_filename, dec_model_filename)
#     total_image_results_deconv.append(dec_model_filename)
#     image_results_deconv.append(dec_model_filename)
#
#     # Save residual
#     residual_filename = (imagename.replace('.fits', '') +
#                          "_" + "residual" + special_name + save_name_append + ".fits")
#     pf.writeto(residual_filename, model_dict['best_residual_conv'], overwrite=True)
#     copy_header(imagename, residual_filename, residual_filename)
#     image_results_conv.append(residual_filename)
#
#     # Save deconvolved background
#     model_dict['deconv_bkg'] = np.asarray(flat_sky_total_dec)
#     deconv_bkg_filename = (imagename.replace('.fits', '') +
#                            "_" + "deconv_bkg" + special_name + save_name_append + '.fits')
#     pf.writeto(deconv_bkg_filename, model_dict['deconv_bkg'], overwrite=True)
#     copy_header(imagename, deconv_bkg_filename, deconv_bkg_filename)
#     bkg_images.append(deconv_bkg_filename)
#
#     # Save convolved background
#     conv_bkg_filename = (imagename.replace('.fits', '') +
#                          "_" + "conv_bkg" + special_name + save_name_append + '.fits')
#     pf.writeto(conv_bkg_filename, model_dict['conv_bkg'], overwrite=True)
#     copy_header(imagename, conv_bkg_filename, conv_bkg_filename)
#     bkg_images.append(conv_bkg_filename)
#
#     # # initial minimization.
#     # method1 = 'differential_evolution'
#     # print(' >> Using', method1, ' solver for first optimisation run... ')
#     # # take parameters from previous run, and re-optimize them.
#     # #     method2 = 'ampgo'#'least_squares'
#     # method2 = 'least_squares'
#
#     # # save mini results (full) to a pickle file.
#     # with open(imagename.replace('.fits',
#     #                             '_' + 'fit' +
#     #                             special_name + save_name_append + '.pickle'),
#     #           "wb") as f:
#     #     pickle.dump(result, f)
#
#     # with open(imagename.replace('.fits',
#     #                             '_' + 'fit' +
#     #                             special_name + save_name_append + '_modeldict.pickle'),
#     #           "wb") as f:
#     #     pickle.dump(model_dict, f)
#
#     return (image_results_conv, image_results_deconv, bkg_images,
#             total_image_results_conv, total_image_results_deconv)
#
#
# # ============================================================================
# # MAIN FUNCTION
# # ============================================================================
#
# def do_fit2D(imagename, params_values_init_IMFIT=None, ncomponents=None,
#              init_constraints=None, data_2D_=None, residualdata_2D_=None,
#              residualname=None, which_residual='shuffled', observation_type='radio',
#              init_params=0.25, final_params=4.0, constrained=True,
#              fix_n=True, fix_value_n=False,
#              fix_max_value_Rn=False, fix_min_value_Rn=False,
#              dr_fix=3,
#              fix_x0_y0=False, psf_name=None, convolution_mode='GPU',
#              convolve_cutout=False, cut_size=512, self_bkg=False,
#              rms_map=None, is_bkg_map_conv=False,
#              fix_geometry=True, contrain_nelder=False, workers=6, mask_region=None,
#              special_name='', method1='least_squares', method2='least_squares',
#              reduce_fcn='neglogcauchy', loss="cauchy", tr_solver="exact", x_scale='jac',
#              ftol=1e-10, xtol=1e-10, gtol=1e-10, verbose=2, max_nfev=200000,
#              regularize=True, f_scale=1.0,
#              maxiter=30000, maxfev=30000, xatol=1e-12,
#              fatol=1e-12, return_all=True, disp=True,
#              de_options=None, parameters_mini_init=None,
#              save_name_append='', logger=None):
#     """
#     Perform a Robust and Fast Multi-Sersic Decomposition with GPU acceleration.
#     tr_solver:
#
#     Parameters
#     ----------
#     imagename: str
#         Name of the image to be fitted.
#     params_values_init_IMFIT: list
#         Initial parameters values for the model.
#     ncomponents: int
#         Number of components to be fitted.
#     init_constraints: dict
#         Initial constraints for the model.
#     data_2D_: 2D array
#         Image to be fitted.
#     residualdata_2D_: 2D array
#         Residual image to be fitted.
#     residualname: str
#         Name of the residual image to be fitted.
#     which_residual: str
#         Which residual to be used for the fitting.
#         Options: 'shuffled' or 'natural'.
#     init_params: float
#         Initial parameters for the model.
#     final_params: float
#         Final parameters for the model.
#     constrained: bool
#         If True, use initial constraints for the model.
#     fix_n: bool
#         If True, fix the Sersic index of the model.
#     fix_value_n: float
#         If True, fix the Sersic index of the model to this value.
#     dr_fix: float
#         If True, fix the centre position of the model.
#     fix_x0_y0: bool
#         If True, fix the centre position of the model.
#     psf_name: str
#         Name of the PSF image to be used for the convolution.
#     convolution_mode: str
#         If 'GPU', use GPU acceleration for the convolution.
#     convolve_cutout: bool
#         If True, convolve the image cutout with the PSF.
#     cut_size: int
#         Size of the cutout image.
#     self_bkg: bool
#         If True, use the image background as the residual background.
#     rms_map: 2D array
#         RMS map to be used for the fitting.
#     fix_geometry: bool
#         If True, fix the geometry of the model.
#     contrain_nelder: bool
#         If True, constrain the Nelder-Mead optimised parameters.
#     workers: int
#         Number of workers to be used for the fitting.
#     mask_region: 2D array
#         Mask to be used for the fitting.
#     special_name: str
#         Special name to be used for the output files.
#     method1: str
#         Method to be used for the fitting.
#     method2: str
#         Method to be used for the fitting.
#     reduce_fcn: str
#
#     loss: str
#
#     tr_solver: str
#
#     x_scale: str
#
#     ftol: float
#
#     xtol: float
#
#     gtol: float
#
#     verbose: int
#
#     max_nfev: int
#
#     regularize: bool
#         If True, regularize the model.
#     f_scale: float
#
#     maxiter: int
#
#     maxfev: int
#
#     xatol: float
#
#     fatol: float
#
#     return_all: bool
#
#     disp: bool
#
#     de_options: dict
#
#     save_name_append: str
#
#     logger: logger
#
#
#     returns
#     -------
#     result: dict
#         Dictionary containing the results of the fitting.
#
#
#     """
#     # ========================================================================
#     # INITIALIZATION
#     # ========================================================================
#
#     # Check JAX availability and adjust convolution mode if needed
#     try:
#         from jax import jit
#     except:
#         convolution_mode = 'CPU'
#
#     if de_options is None:
#         de_options = {'disp': True, 'workers': 6,
#                       'max_nfev': 20000, 'vectorized': True,
#                       # 'strategy': 'randtobest1bin',
#                       'mutation': (0.5, 1.5),
#                       'recombination': [0.2, 0.9],
#                       'init': 'random', 'tol': 0.00001,
#                       'updating': 'deferred',
#                       'popsize': 600}
#
#     startTime = time.time()
#     try:
#         logger.info(f"Fitting image: {imagename}")
#     except:
#         pass
#
#     # ========================================================================
#     # DATA LOADING AND PREPARATION
#     # ========================================================================
#
#     # Load image data
#     data_2D, data_2D_gpu = _initialize_image_data(imagename, data_2D_,
#                                                   convolution_mode, logger)
#
#     # Load PSF data
#     PSF_CONV, PSF_DATA, PSF_DATA_raw = _initialize_psf_data(psf_name,
#                                                             convolution_mode, logger)
#
#     # Prepare mask
#     mask_for_fit = _prepare_mask(mask_region, convolution_mode, logger)
#
#     # Prepare background and residual
#     background, background_dec, residual_2D, FlatSky_level = \
#         _prepare_background_and_residual(
#             residualname, residualdata_2D_, which_residual,
#             self_bkg, rms_map, imagename, data_2D,
#             PSF_CONV, PSF_DATA, PSF_DATA_raw, convolution_mode,
#             is_bkg_map_conv, logger
#         )
#
#     # Convert residual to GPU format if needed
#     if convolution_mode == 'GPU':
#         residual_2D = jnp.array(residual_2D)
#
#     # Setup coordinate grid
#     xy, size = _setup_coordinate_grid(data_2D, convolution_mode)
#
#     if convolve_cutout is True:
#         """
#         WARNING: DO NOT USE FOR NOW!
#
#         Instead of convolving the entire image,
#         convolve only a box.
#         Can be 10x faster.
#
#         Issue:
#         It causes the flat sky level to be much higher than the real value.
#
#         Need further investigation and proper implementation.
#         """
#         x0c, y0c = int(size[0] / 2), int(size[1] / 2)
#
#     #     FlatSky_level = background#mad_std(data_2D)
#
#     # if convolution_mode == 'GPU':
#
#     # FlatSky_level = mad_std(data_2D)
#
#     # ========================================================================
#     # MODEL CONSTRUCTION
#     # ========================================================================
#
#     nfunctions = ncomponents
#
#     # Create minimizer functions
#     min_residual_2D, min_residual_2D_GPU, build_model, func = \
#         _create_minimizer_functions(nfunctions, xy, data_2D, data_2D_gpu,
#                                     background, residual_2D, PSF_DATA, PSF_DATA_raw,
#                                     mask_for_fit, convolution_mode)
#
#     # Setup JAX convolution if needed
#     jax_convolve = _setup_jax_convolution(convolution_mode)
#
#     # Construct model parameters
#
#     smodel2D, params = construct_model_parameters(
#         params_values_init_IMFIT=params_values_init_IMFIT, n_components=nfunctions,
#         init_constraints=init_constraints, observation_type=observation_type,
#         fix_n=fix_n, fix_value_n=fix_value_n,
#         fix_max_value_Rn=fix_max_value_Rn, fix_min_value_Rn=fix_min_value_Rn,
#         fix_x0_y0=fix_x0_y0, dr_fix=dr_fix, fix_geometry=fix_geometry,
#         init_params=init_params, final_params=final_params,
#         constrained=constrained)
#
#     # Create minimizer
#     if convolution_mode == 'CPU':
#         mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#     elif convolution_mode == 'GPU':
#         mini = lmfit.Minimizer(min_residual_2D_GPU, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#     else:
#         # Fallback to CPU
#         mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#
#     # ========================================================================
#     # OPTIMIZATION
#     # ========================================================================
#
#     result, result_1, result_extra = _perform_optimization(
#         mini, params, method1, method2, convolution_mode,
#         contrain_nelder, workers, max_nfev, tr_solver,
#         regularize, x_scale, f_scale, ftol, xtol, gtol,
#         verbose, loss, maxiter, maxfev, xatol, fatol,
#         return_all, disp, de_options, parameters_mini_init
#     )
#
#     params = result.params
#
#     # ========================================================================
#     # MODEL GENERATION AND OUTPUT
#     # ========================================================================
#
#     # Generate model components
#     model_dict, flat_sky_total, flat_sky_total_dec, bkg_comp_i, bkg_comp_i_dec = \
#         _generate_model_components(params, ncomponents, xy, size, background,
#                                    background_dec, PSF_CONV, PSF_DATA, PSF_DATA_raw,
#                                    convolution_mode, is_bkg_map_conv)
#
#     # Compute residuals
#     model_dict = _compute_residuals(data_2D, model_dict, flat_sky_total,
#                                     is_bkg_map_conv, PSF_CONV, PSF_DATA,
#                                     convolution_mode)
#
#     # Save outputs
#     (image_results_conv, image_results_deconv, bkg_images,
#      total_image_results_conv, total_image_results_deconv) = \
#         _save_fitting_outputs(imagename, model_dict, ncomponents, special_name,
#                               save_name_append, flat_sky_total_dec)
#
#     # ========================================================================
#     # FINAL PROCESSING
#     # ========================================================================
#
#     exec_time = time.time() - startTime
#     print('Exec time fitting=', exec_time, 's')
#
#     # save results to csv file.
#     try:
#         save_results_csv(result_mini=result,
#                          save_name=image_results_conv[-2].replace('.fits', ''),
#                          ext='.csv',
#                          save_corr=True, save_params=True)
#     except:
#         print('Error Saving Results to a csv file!!!')
#         pass
#
#     # Extract model_temp for return (last component generated)
#     model_temp = Model(sersic2D_GPU)
#
#     return (result, mini, result_1, result_extra, model_dict, image_results_conv,
#             image_results_deconv, bkg_images, smodel2D, model_temp)



# def do_fit2D(imagename, params_values_init_IMFIT=None, ncomponents=None,
#              init_constraints=None, data_2D_=None, residualdata_2D_=None,
#              residualname=None,which_residual='shuffled',observation_type = 'radio',
#              init_params=0.25, final_params=4.0, constrained=True,
#              fix_n=True, fix_value_n=False,
#              fix_max_value_Rn=False,fix_min_value_Rn = False,
#              dr_fix=3,
#              fix_x0_y0=False, psf_name=None, convolution_mode='GPU',
#              convolve_cutout=False, cut_size=512, self_bkg=False,
#              rms_map=None, is_bkg_map_conv = False,
#              fix_geometry=True, contrain_nelder=False, workers=6,mask_region = None,
#              special_name='', method1='least_squares', method2='least_squares',
#              reduce_fcn='neglogcauchy',loss="cauchy",tr_solver="exact",x_scale = 'jac',
#              ftol=1e-10, xtol=1e-10, gtol=1e-10, verbose=2,max_nfev=200000,
#              regularize  = True, f_scale = 1.0,
#              maxiter = 30000, maxfev = 30000, xatol = 1e-12,
#              fatol = 1e-12, return_all = True, disp = True,
#              de_options=None,parameters_mini_init=None,
#              save_name_append='',logger=None):
#     """
#     Perform a Robust and Fast Multi-Sersic Decomposition with GPU acceleration.
#     tr_solver:
#
#     Parameters
#     ----------
#     imagename: str
#         Name of the image to be fitted.
#     params_values_init_IMFIT: list
#         Initial parameters values for the model.
#     ncomponents: int
#         Number of components to be fitted.
#     init_constraints: dict
#         Initial constraints for the model.
#     data_2D_: 2D array
#         Image to be fitted.
#     residualdata_2D_: 2D array
#         Residual image to be fitted.
#     residualname: str
#         Name of the residual image to be fitted.
#     which_residual: str
#         Which residual to be used for the fitting.
#         Options: 'shuffled' or 'natural'.
#     init_params: float
#         Initial parameters for the model.
#     final_params: float
#         Final parameters for the model.
#     constrained: bool
#         If True, use initial constraints for the model.
#     fix_n: bool
#         If True, fix the Sersic index of the model.
#     fix_value_n: float
#         If True, fix the Sersic index of the model to this value.
#     dr_fix: float
#         If True, fix the centre position of the model.
#     fix_x0_y0: bool
#         If True, fix the centre position of the model.
#     psf_name: str
#         Name of the PSF image to be used for the convolution.
#     convolution_mode: str
#         If 'GPU', use GPU acceleration for the convolution.
#     convolve_cutout: bool
#         If True, convolve the image cutout with the PSF.
#     cut_size: int
#         Size of the cutout image.
#     self_bkg: bool
#         If True, use the image background as the residual background.
#     rms_map: 2D array
#         RMS map to be used for the fitting.
#     fix_geometry: bool
#         If True, fix the geometry of the model.
#     contrain_nelder: bool
#         If True, constrain the Nelder-Mead optimised parameters.
#     workers: int
#         Number of workers to be used for the fitting.
#     mask_region: 2D array
#         Mask to be used for the fitting.
#     special_name: str
#         Special name to be used for the output files.
#     method1: str
#         Method to be used for the fitting.
#     method2: str
#         Method to be used for the fitting.
#     reduce_fcn: str
#
#     loss: str
#
#     tr_solver: str
#
#     x_scale: str
#
#     ftol: float
#
#     xtol: float
#
#     gtol: float
#
#     verbose: int
#
#     max_nfev: int
#
#     regularize: bool
#         If True, regularize the model.
#     f_scale: float
#
#     maxiter: int
#
#     maxfev: int
#
#     xatol: float
#
#     fatol: float
#
#     return_all: bool
#
#     disp: bool
#
#     de_options: dict
#
#     save_name_append: str
#
#     logger: logger
#
#
#     returns
#     -------
#     result: dict
#         Dictionary containing the results of the fitting.
#
#
#     """
#     try:
#         from jax import jit
#     except:
#         convolution_mode = 'CPU'
#
#     if de_options is None:
#         de_options = {'disp': True, 'workers': 6,
#                       'max_nfev': 20000, 'vectorized': True,
#                       # 'strategy': 'randtobest1bin',
#                       'mutation': (0.5, 1.5),
#                       'recombination': [0.2, 0.9],
#                       'init': 'random', 'tol': 0.00001,
#                       'updating': 'deferred',
#                       'popsize': 600}
#
#     startTime = time.time()
#     try:
#         logger.info(f"Fitting image: {imagename}")
#     except:
#         pass
#     FlatSky_level = None
#     if data_2D_ is None:
#         data_2D = pf.getdata(imagename)
#     else:
#         data_2D = data_2D_
#
#     if mask_region is not None:
#         """
#
#         """
#         logger.debug(f" ==> Using provided mask region to constrain fit. ")
#         logger.warning(f" !!++==> Fitting with a mask is faster, but experimental!! \n"
#                        f"         Use with caution.")
#         # data_2D = data_2D * mask_region
#         if convolution_mode == 'GPU':
#             mask_for_fit = jnp.array(mask_region)
#         else:
#             mask_for_fit = mask_region
#     else:
#         mask_for_fit = None
#
#     if convolution_mode == 'GPU':
#         data_2D_gpu = jnp.array(data_2D)
#
#     if psf_name is not None:
#         PSF_CONV = True
#         try:
#             PSF_DATA_raw = pf.getdata(psf_name)
#             if len(PSF_DATA_raw.shape) == 4:
#                 PSF_DATA_raw = PSF_DATA_raw[0][0]
#         except:
#             PSF_DATA_raw = load_fits_data(psf_name)
#
#         if convolution_mode == 'GPU':
#             if logger is not None:
#                 logger.debug(f"---------------------------------------")
#                 logger.debug(f" <<< PERFORMING CONVOLUTION WITH JAX >>> ")
#                 logger.debug(f"---------------------------------------")
#             PSF_DATA = jnp.array(PSF_DATA_raw)
#
#         if convolution_mode == 'CPU':
#             PSF_DATA = PSF_DATA_raw
#         # PSF_DATA = pf.getdata(
#         #     imagename.replace('-image.cutout.fits', '-beampsf.cutout.fits'))
#     else:
#         PSF_CONV = False
#         PSF_DATA = None
#
#     if residualname is not None and which_residual != 'user':
#         """
#         This is important for radio image fitting.
#
#         It uses the shuffled version of the residual cleaned image
#         originated from the interferometric deconvolution.
#
#         This ensures that the best model created here will be on top
#         of that rms noise so that flux conservation is maximized.
#
#         However, this residual is not added as model + shuffled_residual
#         only, but instead by a multiplication factor,
#         e.g. model + const* shuffled_residual, and const will be minimized
#         as well during the fitting (here, called `s_a`).
#         """
#         if residualdata_2D_ is not None:
#             residual_2D = residualdata_2D_
#         else:
#             residual_2D = pf.getdata(residualname)
#
#         if which_residual == 'shuffled':
#             if logger is not None:
#                 logger.debug(f" ==> Using clean shuffled background for optmization... ")
#             residual_2D_to_use = shuffle_2D(residual_2D)
#         if which_residual == 'natural':
#             if logger is not None:
#                 logger.debug(f" ==> Using clean background for optmization... ")
#             """
#             if psf_name is not None:
#                 if logger is not None:
#                     logger.debug(f" ====> Deconvolving residual map... ")
#                 residual_2D_to_use, _ = deconvolve_fft(residual_2D,
#                                                             PSF_DATA_raw/PSF_DATA_raw.sum())
#             else:
#                 residual_2D_to_use = residual_2D
#             """
#             residual_2D_to_use = residual_2D
#
#         FlatSky_level = mad_std(residual_2D_to_use)
#         #         background = residual_2D #residual_2D_to_use
#         if convolution_mode == 'GPU':
#             background = jnp.array(residual_2D_to_use)
#         else:
#             background = residual_2D_to_use
#
#     else:
#         if which_residual == 'user':
#             if rms_map is None:
#                 print('--==>> A rms map/background mode was selected (user)')
#                 print('       but no rms/background map was provided.')
#                 print('       Please, provide a rms/background map.')
#                 print('||==>> Stopping code now.')
#                 raise ValueError("rms_map should not be None")
#             else:
#                 logger.debug(f" ==> Using provided RMS map. ")
#                 background_map=rms_map
#                 background = background_map.copy()
#         else:
#             if self_bkg == True:
#                 if logger is not None:
#                     logger.warning(f" ==> No residual/background provided. Using image bkg map... ")
#                 background_map = sep_background(imagename)
#                 background = shuffle_2D(background_map.back())
#             else:
#                 if logger is not None:
#                     logger.warning(f" ==> Using only flat sky for rms bkg.")
#                 FlatSky_level = mad_std(data_2D)
#                 background = FlatSky_level
#
#     if residualdata_2D_ is not None:
#         residual_2D = residualdata_2D_
#     else:
#         try:
#             residual_2D = pf.getdata(residualname)
#         except:
#             residual_2D = background
#
#     background_dec = background.copy()
#     if is_bkg_map_conv is False:
#         logger.debug(f" ==> RMS map is not convolved.")
#
#         if convolution_mode == 'GPU':
#             background = _fftconvolve_jax(background,PSF_DATA)
#         if convolution_mode == 'CPU':
#             background = scipy.signal.fftconvolve(background, PSF_DATA, 'same')
#
#
#     size = data_2D.shape
#     if convolution_mode == 'GPU':
#         x,y = jnp.meshgrid(jnp.arange((size[1])), jnp.arange((size[0])))
#         xy = jnp.stack([x, y], axis=0)
#     else:
#         xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
#
#     if convolve_cutout is True:
#         """
#         WARNING: DO NOT USE FOR NOW!
#
#         Instead of convolving the entire image,
#         convolve only a box.
#         Can be 10x faster.
#
#         Issue:
#         It causes the flat sky level to be much higher than the real value.
#
#         Need further investigation and proper implementation.
#         """
#         x0c, y0c = int(size[0] / 2), int(size[1] / 2)
#
#     #     FlatSky_level = background#mad_std(data_2D)
#
#
#     # if convolution_mode == 'GPU':
#
#     # FlatSky_level = mad_std(data_2D)
#     if convolution_mode == 'GPU':
#         residual_2D = jnp.array(residual_2D)
#
#     nfunctions = ncomponents
#
#     def min_residual_2D(params):
#         dict_model = {}
#         model = 0
#         for i in range(1, nfunctions + 1):
#             model = model + sersic2D(xy, params['f' + str(i) + '_x0'],
#                                      params['f' + str(i) + '_y0'],
#                                      params['f' + str(i) + '_PA'],
#                                      params['f' + str(i) + '_ell'],
#                                      params['f' + str(i) + '_n'],
#                                      params['f' + str(i) + '_In'],
#                                      params['f' + str(i) + '_Rn'],
#                                      params['f' + str(i) + '_cg'], )
#         # print(model.shape)
#         # model = model + FlatSky_cpu(FlatSky_level, params['s_a'])*background
#         # model = model + FlatSky_cpu(background, params['s_a'])
#         MODEL_2D_conv = scipy.signal.fftconvolve(model, PSF_DATA, 'same') + \
#                         FlatSky_cpu(background, params['s_a'])
#         residual = data_2D - MODEL_2D_conv
#         return np.ravel(residual)
#
#     def convert_params_to_numpy(_params):
#         return list(_params)
#
#     try:
#         # @partial(jit, static_argnums=1)
#         @jit
#         def func(x):
#             return jnp.split(x, nfunctions)
#
#         @jit
#         def build_model(xy,param_matrix):
#             model = 0
#             for model_params in param_matrix:
#                 model = model + sersic2D_GPU(xy, model_params[0],
#                                             model_params[1],
#                                             model_params[2],
#                                             model_params[3],
#                                             model_params[4],
#                                             model_params[5],
#                                             model_params[6],
#                                             model_params[7])
#             return model
#
#     except:
#         def func(x):
#             return jnp.split(x, nfunctions)
#
#
#         def build_model(xy,param_matrix):
#             model = 0
#             for model_params in param_matrix:
#                 model = model + sersic2D_GPU(xy, model_params[0],
#                                             model_params[1],
#                                             model_params[2],
#                                             model_params[3],
#                                             model_params[4],
#                                             model_params[5],
#                                             model_params[6],
#                                             model_params[7])
#             return model
#
#     def min_residual_2D_GPU(params):
#         model = 0
#         for i in range(1, nfunctions + 1):
#             model = model + sersic2D_GPU(xy,
#                                          params['f' + str(i) + '_x0'].value,
#                                          params['f' + str(i) + '_y0'].value,
#                                          params['f' + str(i) + '_PA'].value,
#                                          params['f' + str(i) + '_ell'].value,
#                                          params['f' + str(i) + '_n'].value,
#                                          params['f' + str(i) + '_In'].value,
#                                          params['f' + str(i) + '_Rn'].value,
#                                          params['f' + str(i) + '_cg'].value)
#
#         # # param_matrix = extract_params(params)
#         # param_matrix = func(jnp.array(list(params.valuesdict().values()))[:-1])
#         # model = build_model(xy,param_matrix)
#
#         MODEL_2D_conv = _fftconvolve_jax(model,PSF_DATA) + FlatSky(background,params['s_a'].value)
#
#         # MODEL_2D_conv = _fftconvolve_jax(model+
#         #                                  FlatSky(background,params['s_a'].value),
#         #                                  PSF_DATA)
#         # residual = ((data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])/
#         #             (1000*(abs(residual_2D[mask_for_fit])+1.0e-6)))
#         residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])
#         return np.asarray(residual).copy()
#         # weights = 1/((background[mask_for_fit])/data_2D_gpu[mask_for_fit])
#         # weightned_residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit]) * jnp.sqrt(weights)
#         # return np.asarray(weightned_residual).copy()
#
#     if convolution_mode == 'GPU':
#         @jit
#         def convolve_on_gpu(image, psf):
#             """
#             This was before jax.scipy implementing fftconvolve.
#             It provides the same result, at the same speed.
#
#             This function also accepts PSFs with a different shape of the image.
#
#             """
#             # Calculate the new padded shape
#             padded_shape = (image.shape[0] + psf.shape[0] - 1,
#                             image.shape[1] + psf.shape[1] - 1)
#
#             # Pad both image and psf to the new shape
#             pad_shape = [(0, ts - s) for s, ts in zip(image.shape, padded_shape)]
#             image_padded = jnp.pad(image, pad_shape, mode='constant')
#             pad_shape = [(0, ts - s) for s, ts in zip(psf.shape, padded_shape)]
#             psf_padded = jnp.pad(psf, pad_shape, mode='constant')
#             # psf_padded = pad_for_convolution(psf, padded_shape)
#             image_fft = jnp.fft.fft2(image_padded)
#             psf_fft = jnp.fft.fft2(psf_padded)
#
#             conv_fft = image_fft * psf_fft
#
#             # Get the real part of the inverse FFT and crop to the original image size
#             result_full = jnp.real(jnp.fft.ifft2(conv_fft))
#             return result_full[psf.shape[0] // 2:image.shape[0] + psf.shape[0] // 2,
#                    psf.shape[1] // 2:image.shape[1] + psf.shape[1] // 2]
#
#         jax_convolve = jit(convolve_on_gpu)
#
#
#     smodel2D, params = construct_model_parameters(
#         params_values_init_IMFIT=params_values_init_IMFIT, n_components=nfunctions,
#         init_constraints=init_constraints,observation_type=observation_type,
#         fix_n=fix_n, fix_value_n=fix_value_n,
#         fix_max_value_Rn = fix_max_value_Rn, fix_min_value_Rn = fix_min_value_Rn,
#         fix_x0_y0=fix_x0_y0, dr_fix=dr_fix, fix_geometry=fix_geometry,
#         init_params=init_params, final_params=final_params,
#         constrained=constrained)
#
#     if convolution_mode == 'CPU':
#         mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#     if convolution_mode == 'GPU':
#         mini = lmfit.Minimizer(min_residual_2D_GPU, params, max_nfev=200000,
#                                nan_policy='omit', reduce_fcn=reduce_fcn)
#
#     # initial minimization.
#
#     print(' >> Using', method1, ' solver for first optimisation run... ')
#     # take parameters from previous run, and re-optimize them.
#     #     me``rthod2 = 'ampgo'#'least_squares'
#     #     method2 = 'least_squares'
#     result_extra = None
#     if method1 == 'nelder':
#         # very robust, but takes time....
#         #         print(' >> Using', method1,' solver for first optimisation run... ')
#         result_1 = mini.minimize(method='nelder',
#                                  #                                  xatol = 1e-12, fatol = 1e-12, disp=True,
#                                  #                                  adaptive = True,max_nfev = 30000,
#                                  options={'maxiter': maxiter, 'maxfev': maxfev,
#                                           'xatol': xatol, 'fatol': fatol,
#                                           'return_all': return_all,
#                                           'disp': disp}
#                                  )
#
#
#     if method1 == 'least_squares':
#         # faster, but usually not good for first run.
#         # if results_previous_run is not None:
#         print(' >> Using',tr_solver,'for tr solver, with regularize set to',regularize,
#               ' Loss is',loss,'.')
#
#         if parameters_mini_init is not None:
#             print(f'  ++==>> Using initial mini parameters from a previous run.')
#             # try:
#             result_1 = mini.minimize(method='least_squares',
#                                      params=parameters_mini_init,
#                                     max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
#                                     tr_solver=tr_solver,
#                                     tr_options={'regularize': regularize,
#                                                 },
#                                     ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                     loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#         else:
#             result_1 = mini.minimize(method='least_squares',
#                                     max_nfev=max_nfev, x_scale=x_scale, f_scale=f_scale,
#                                     tr_solver=tr_solver,
#                                     tr_options={'regularize': regularize,
#     #                                              'min_delta':1e-14, 'eta':0.05,
#     #                                              'xtol':1e-14, 'gtol':1e-14,
#     #                                              'ftol':1e-14
#                                                 },
#                                     ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                     loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#
#     if method1 == 'differential_evolution':
#         # de is giving some issues, I do not know why.
#         result_1 = mini.minimize(method='differential_evolution',
#                                  options={'disp': True, 'workers': workers,
#                                           'max_nfev': max_nfev, 'vectorized': True,
#                                           'strategy': 'randtobest1bin',
#                                           'mutation': (0.5, 1.5),
#                                           'recombination': [0.2, 0.9],
#                                           'init': 'random', 'tol': 0.00001,
#                                           'updating': 'deferred',
#                                           'popsize': 600})
#         # result_1 = mini.minimize(method='differential_evolution', popsize=600,
#         #                          disp=True,  # init = 'random',
#         #                          # mutation=(0.5, 1.5), recombination=[0.2, 0.9],
#         #                          max_nfev=20000,
#         #                          workers=1, updating='deferred', vectorized=True)
#
#     print(' >> Using', method2, ' solver for second optimisation run... ')
#
#     second_run_params = result_1.params
#     if (contrain_nelder == True) and (method2 == 'nelder'):
#         """
#         It seems that least_squares is ignoring the best-parameters provided by
#         Nelder-mead, which means that it is lookig the parameter space far away
#         from the optimised Nelder-Mead ones.
#
#         So, with this condition, we force a much smaller searching region, but
#         it assumes that Nelder opt was good (which is not always true).
#
#         YOU MUST CHECK YOUR RESULTS!!!!
#
#         """
#         print('Constraining Nelder-Mead Parameters for method', method2)
#         params_constrained = constrain_nelder_mead_params(result_1.params,
#                                                           max_factor=1.03,
#                                                           min_factor=0.97)
#         # UPDATE THE SECOND RUN PARAMETERS TO BE THE CONSTRAINED ONES.
#         second_run_params = params_constrained
#
#     if method2 == 'nelder':
#         result = mini.minimize(method='nelder', params=second_run_params,
#                                options={'maxiter': maxiter, 'maxfev': maxfev,
#                                         'xatol': xatol, 'fatol': fatol,
#                                         'disp': disp})
#
#     if method2 == 'ampgo':
#         # ampgo is not workin well/ takes so long ???
#         result = mini.minimize(method='ampgo', params=second_run_params,
#                                maxfunevals=10000, totaliter=30, disp=True,
#                                maxiter=5, glbtol=1e-8)
#
#     if method2 == 'least_squares':
#         # faster, usually converges and provide errors.
#         # Very robust if used in second opt from first opt parameters.
#         result = mini.minimize(method='least_squares',
#                                params=second_run_params,
#                                max_nfev=max_nfev,
#                                tr_solver=tr_solver,
#                                tr_options={'regularize': regularize,
# #                                            'min_delta': 1e-14, 'eta': 0.05,
# #                                            'xtol': 1e-14, 'gtol': 1e-14,
# #                                            'ftol': 1e-14
#                                           },
#                                x_scale=x_scale,  f_scale=f_scale,
#                                ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
#                                loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)
#
#     if method2 == 'differential_evolution':
#         # result = mini.minimize(method='differential_evolution',
#         #                        params=second_run_params,
#         #                        options={'maxiter': 30000, 'workers': -1,
#         #                                 'tol': 0.001, 'vectorized': True,
#         #                                 'strategy': 'randtobest1bin',
#         #                                 'updating': 'deferred', 'disp': True,
#         #                                 'seed': 1}
#         #                        )
#         result = mini.minimize(method='differential_evolution',
#                                params=second_run_params,
#                                options=de_options
#                                )
#
#     params = result.params
#
#     model_temp = Model(sersic2D)
#     xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
#     model = 0
#     model_dict = {}
#     image_results_conv = []
#     image_results_deconv = []
#     total_image_results_deconv = []
#     total_image_results_conv = []
#     bkg_images = []
#     if convolution_mode == 'GPU':
#         flat_sky_total = np.asarray(FlatSky(background, params['s_a'].value))
#         flat_sky_total_dec = np.asarray(FlatSky(background_dec, params['s_a'].value))
#     if convolution_mode == 'CPU':
#         flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
#         flat_sky_total_dec = FlatSky_cpu(background_dec, params['s_a'].value)
#     bkg_comp_i = flat_sky_total.copy()
#     bkg_comp_i_dec = flat_sky_total_dec.copy()
#     for i in range(1, ncomponents + 1):
#         model_temp = sersic2D_GPU(xy, params['f' + str(i) + '_x0'].value,
#                               params['f' + str(i) + '_y0'].value,
#                               params['f' + str(i) + '_PA'].value,
#                               params['f' + str(i) + '_ell'].value,
#                               params['f' + str(i) + '_n'].value,
#                               params['f' + str(i) + '_In'].value,
#                               params['f' + str(i) + '_Rn'].value,
#                               params['f' + str(i) + '_cg'].value)
#
#         model = model + model_temp
#         #to each individual component, add the bkg map.
#         model_dict['model_c' + str(i)] = np.asarray(model_temp+bkg_comp_i_dec)
#         if PSF_CONV == True:
#             if convolution_mode == 'GPU':
#                 # model_dict['model_c' + str(i) + '_conv'] = np.asarray(jax_convolve(model_temp, PSF_DATA)).copy()
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #     np.asarray(_fftconvolve_jax(model_temp,PSF_DATA).copy()+bkg_comp_i))
#                 # to each individual component, add the bkg map.
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #     np.asarray(_fftconvolve_jax(model_temp+bkg_comp_i,PSF_DATA).copy()))
#
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                     np.asarray(_fftconvolve_jax(model_temp,PSF_DATA).copy())) + bkg_comp_i
#
#             if convolution_mode == 'CPU':
#                 # model_dict['model_c' + str(i) + '_conv'] = (
#                 #         scipy.signal.fftconvolve(model_temp+bkg_comp_i, PSF_DATA_raw,'same'))
#                 model_dict['model_c' + str(i) + '_conv'] = (
#                         scipy.signal.fftconvolve(model_temp, PSF_DATA_raw,'same') + bkg_comp_i
#                         )
#
#
#         else:
#             model_dict['model_c' + str(i) + '_conv'] = model_temp+bkg_comp_i
#
#         pf.writeto(imagename.replace('.fits', '') +
#                    "_" + "model_component_" + str(i) +
#                    special_name + save_name_append + '.fits',
#                    model_dict['model_c' + str(i) + '_conv'], overwrite=True)
#         copy_header(imagename, imagename.replace('.fits', '') +
#                     "_" + "model_component_" + str(i) +
#                     special_name + save_name_append + '.fits',
#                     imagename.replace('.fits', '') +
#                     "_" + "model_component_" + str(
#                         i) + special_name + save_name_append + '.fits')
#         pf.writeto(imagename.replace('.fits', '') +
#                    "_" + "dec_model_component_" + str(i) +
#                    special_name + save_name_append + '.fits',
#                    model_dict['model_c' + str(i)], overwrite=True)
#         copy_header(imagename, imagename.replace('.fits', '') +
#                     "_" + "dec_model_component_" + str(i) +
#                     special_name + save_name_append + '.fits',
#                     imagename.replace('.fits', '') +
#                     "_" + "dec_model_component_" + str(i) +
#                     special_name + save_name_append + '.fits')
#
#         image_results_conv.append(imagename.replace('.fits', '') +
#                                   "_" + "model_component_" + str(i) +
#                                   special_name + save_name_append + '.fits')
#         image_results_deconv.append(imagename.replace('.fits', '') +
#                                     "_" + "dec_model_component_" + str(i) +
#                                     special_name + save_name_append + '.fits')
#
#     #     model = model
#     model_dict['model_total_dec'] = np.asarray(model+flat_sky_total_dec) # +FlatSky_cpu(background,
#     # params['s_a'].value)
#
#     if PSF_CONV == True:
#         # model_dict['model_total_conv'] = scipy.signal.fftconvolve(model,
#         #                                                           PSF_DATA_raw,
#         #                                                           'same')  # + FlatSky(FlatSky_level, params['s_a'])
#         if convolution_mode == 'GPU':
#             # model_dict['model_total_conv'] = np.asarray(jax_convolve(model,
#             #                                                          PSF_DATA)).copy()
#             # model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + FlatSky_cpu(background, params['s_a'].value
#             # model_conv = _fftconvolve_jax(model+flat_sky_total,PSF_DATA).copy()
#             model_conv = _fftconvolve_jax(model,PSF_DATA).copy() + flat_sky_total
#         if convolution_mode == 'CPU':
#             model_conv = scipy.signal.fftconvolve(model, PSF_DATA_raw,'same') + flat_sky_total
#             model_dict['model_total_conv'] = model_conv
#     else:
#         model_dict['model_total_conv'] = model + flat_sky_total
#
#
#     # model_dict['best_residual'] = data_2D - model_dict['model_total']
#     # bkg_comp_total
#     model_dict['model_total_conv'] = np.asarray(model_conv)
#     if is_bkg_map_conv is False and PSF_CONV is True:
#         model_dict['best_residual_conv'] = np.asarray(data_2D) - model_dict['model_total_conv']
#         model_dict['conv_bkg'] = np.asarray(_fftconvolve_jax(flat_sky_total,PSF_DATA).copy())
#     else:
#         model_dict['best_residual_conv'] = np.asarray(data_2D) - model_dict['model_total_conv'] + flat_sky_total
#         model_dict['conv_bkg'] = np.asarray(flat_sky_total)
#
#     model_dict['deconv_bkg'] = np.asarray(flat_sky_total_dec)
#
#
#
#     pf.writeto(imagename.replace('.fits', '') +
#                "_" + "conv_model" + special_name + save_name_append + '.fits',
#                model_dict['model_total_conv'], overwrite=True)
#
#     total_image_results_conv.append(imagename.replace('.fits', '') +
#                "_" + "conv_model" + special_name + save_name_append + '.fits')
#
#     pf.writeto(imagename.replace('.fits', '') +
#                "_" + "dec_model" + special_name + save_name_append + '.fits',
#                model_dict['model_total_dec'], overwrite=True)
#
#     total_image_results_deconv.append(imagename.replace('.fits', '') +
#                "_" + "dec_model" + special_name + save_name_append + '.fits')
#
#
#     pf.writeto(imagename.replace('.fits', '') +
#                "_" + "residual" + special_name + save_name_append + ".fits",
#                model_dict['best_residual_conv'], overwrite=True)
#     copy_header(imagename, imagename.replace('.fits', '') +
#                 "_" + "conv_model" + special_name + save_name_append + '.fits',
#                 imagename.replace('.fits', '') +
#                 "_" + "conv_model" + special_name + save_name_append + '.fits')
#     copy_header(imagename, imagename.replace('.fits', '') +
#                 "_" + "dec_model" + special_name + save_name_append + '.fits',
#                 imagename.replace('.fits', '') +
#                 "_" + "dec_model" + special_name + save_name_append + '.fits')
#     copy_header(imagename, imagename.replace('.fits', '') +
#                 "_" + "residual" + special_name + save_name_append + '.fits',
#                 imagename.replace('.fits', '') +
#                 "_" + "residual" + special_name + save_name_append + '.fits')
#
#     pf.writeto(imagename.replace('.fits', '') +
#                "_" + "dec_model" + special_name + save_name_append + '.fits',
#                model_dict['model_total_dec'], overwrite=True)
#     copy_header(imagename, imagename.replace('.fits', '') +
#                 "_" + "dec_model" + special_name + save_name_append + '.fits',
#                 imagename.replace('.fits', '') +
#                 "_" + "dec_model" + special_name + save_name_append + '.fits')
#
#     pf.writeto(imagename.replace('.fits', '') +
#                "_" + "deconv_bkg" + special_name + save_name_append + '.fits',
#                model_dict['deconv_bkg'], overwrite=True)
#     copy_header(imagename, imagename.replace('.fits', '') +
#                 "_" + "deconv_bkg" + special_name + save_name_append + '.fits',
#                 imagename.replace('.fits', '') +
#                 "_" + "deconv_bkg" + special_name + save_name_append + '.fits')
#
#     pf.writeto(imagename.replace('.fits', '') +
#                "_" + "conv_bkg" + special_name + save_name_append + '.fits',
#                model_dict['conv_bkg'], overwrite=True)
#     copy_header(imagename, imagename.replace('.fits', '') +
#                 "_" + "conv_bkg" + special_name + save_name_append + '.fits',
#                 imagename.replace('.fits', '') +
#                 "_" + "conv_bkg" + special_name + save_name_append + '.fits')
#     bkg_images.append(imagename.replace('.fits', '') +
#                "_" + "deconv_bkg" + special_name + save_name_append + '.fits')
#     bkg_images.append(imagename.replace('.fits', '') +
#                "_" + "conv_bkg" + special_name + save_name_append + '.fits')
#
#
#     # # initial minimization.
#     # method1 = 'differential_evolution'
#     # print(' >> Using', method1, ' solver for first optimisation run... ')
#     # # take parameters from previous run, and re-optimize them.
#     # #     method2 = 'ampgo'#'least_squares'
#     # method2 = 'least_squares'
#
#
#
#     image_results_conv.append(imagename.replace('.fits', '') +
#                               "_" + "conv_model" +
#                               special_name + save_name_append + '.fits')
#     image_results_deconv.append(imagename.replace('.fits', '') +
#                                 "_" + "dec_model" +
#                                 special_name + save_name_append + '.fits')
#     image_results_conv.append(imagename.replace('.fits', '') +
#                               "_" + "residual" +
#                               special_name + save_name_append + ".fits")
#
#     # # save mini results (full) to a pickle file.
#     # with open(imagename.replace('.fits',
#     #                             '_' + 'fit' +
#     #                             special_name + save_name_append + '.pickle'),
#     #           "wb") as f:
#     #     pickle.dump(result, f)
#
#     # with open(imagename.replace('.fits',
#     #                             '_' + 'fit' +
#     #                             special_name + save_name_append + '_modeldict.pickle'),
#     #           "wb") as f:
#     #     pickle.dump(model_dict, f)
#
#
#     exec_time = time.time() - startTime
#     print('Exec time fitting=', exec_time, 's')
#
#
#
#     # save results to csv file.
#     try:
#         save_results_csv(result_mini=result,
#                          save_name=image_results_conv[-2].replace('.fits', ''),
#                          ext='.csv',
#                          save_corr=True, save_params=True)
#     except:
#         print('Error Saving Results to a csv file!!!')
#         pass
#
#     return (result, mini, result_1, result_extra, model_dict, image_results_conv,
#             image_results_deconv, bkg_images, smodel2D, model_temp)


def return_and_save_model(mini_results, imagename, ncomponents, background=0.0,
                          save_results=False,save_name_append=''):
    params = mini_results.params
    data_2D = load_fits_data(imagename)
    model_temp = Model(sersic2D)
    model = 0
    PSF_CONV = True
    size = load_fits_data(imagename).shape
    FlatSky_level = mad_std(data_2D)
    xy = np.meshgrid(np.arange((size[0])), np.arange((size[1])))
    model_dict = {}
    image_results_conv = []
    image_results_deconv = []
    for i in range(1, ncomponents + 1):
        model_temp = sersic2D(xy, params['f' + str(i) + '_x0'],
                              params['f' + str(i) + '_y0'],
                              params['f' + str(i) + '_PA'],
                              params['f' + str(i) + '_ell'],
                              params['f' + str(i) + '_n'],
                              params['f' + str(i) + '_In'],
                              params['f' + str(i) + '_Rn'],
                              params['f' + str(i) + '_cg'],) + \
                     background / ncomponents + FlatSky(FlatSky_level,params['s_a']) / ncomponents
        # print(model_temp[0])
        model = model + model_temp
        # print(model)
        model_dict['model_c' + str(i)] = model_temp

        if PSF_CONV == True:
            model_dict['model_c' + str(i) + '_conv'] = scipy.signal.fftconvolve(
                model_temp, PSF_DATA,
                'same')  # + FlatSky(FlatSky_level, params['s_a'])/ncomponents
        else:
            model_dict['model_c' + str(i) + '_conv'] = model_temp

        if save_results is True:
            pf.writeto(imagename.replace('.fits', '') + "_" + str(
                ncomponents) + "C_model_component_" + str(
                i) + special_name + save_name_append + '.fits',
                       model_dict['model_c' + str(i) + '_conv'], overwrite=True)
            copy_header(imagename, imagename.replace('.fits', '') + "_" + str(
                ncomponents) + "C_model_component_" + str(
                i) + special_name + save_name_append + '.fits',
                        imagename.replace('.fits', '') + "_" + str(
                            ncomponents) + "C_model_component_" + str(
                            i) + special_name + save_name_append + '.fits')
            pf.writeto(imagename.replace('.fits', '') + "_" + str(
                ncomponents) + "C_dec_model_component_" + str(
                i) + special_name + save_name_append + '.fits',
                       model_dict['model_c' + str(i)], overwrite=True)
            copy_header(imagename, imagename.replace('.fits', '') + "_" + str(
                ncomponents) + "C_dec_model_component_" + str(
                i) + special_name + save_name_append + '.fits',
                        imagename.replace('.fits', '') + "_" + str(
                            ncomponents) + "C_dec_model_component_" + str(
                            i) + special_name + save_name_append + '.fits')

            image_results_conv.append(imagename.replace('.fits', '') + "_" + str(
                ncomponents) + "C_model_component_" + str(
                i) + special_name + save_name_append + '.fits')
            image_results_deconv.append(
                imagename.replace('.fits', '') + "_" + str(
                    ncomponents) + "C_dec_model_component_" + str(
                    i) + special_name + save_name_append + '.fits')

    #     model = model
    model_dict['model_total'] = model  # + FlatSky(FlatSky_level, params['s_a'])

    if PSF_CONV == True:
        model_dict['model_total_conv'] = scipy.signal.fftconvolve(model,
                                                                  PSF_DATA,
                                                                  'same')  # + FlatSky(FlatSky_level, params['s_a'])
    else:
        model_dict['model_total_conv'] = model

    model_dict['best_residual'] = data_2D - model_dict['model_total']
    model_dict['best_residual_conv'] = data_2D - model_dict['model_total_conv']

    if save_results == True:
        pf.writeto(imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_model" + special_name + save_name_append + '.fits',
                   model_dict['model_total_conv'], overwrite=True)
        pf.writeto(imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_residual" + special_name + save_name_append + ".fits",
                   model_dict['best_residual_conv'], overwrite=True)
        copy_header(imagename, imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_model" + special_name + save_name_append + '.fits',
                    imagename.replace('.fits', '') + "_" + str(
                        ncomponents) + "C_model" + special_name + save_name_append + '.fits')
        copy_header(imagename, imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_residual" + special_name + save_name_append + '.fits',
                    imagename.replace('.fits', '') + "_" + str(
                        ncomponents) + "C_residual" + special_name + save_name_append + '.fits')

        pf.writeto(imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_dec_model" + special_name + save_name_append + '.fits',
                   model_dict['model_total'], overwrite=True)
        copy_header(imagename, imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_dec_model" + special_name + save_name_append + '.fits',
                    imagename.replace('.fits', '') + "_" + str(
                        ncomponents) + "C_dec_model" + special_name + save_name_append + '.fits')
        # initial minimization.

        image_results_conv.append(imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_model" + special_name + save_name_append + '.fits')
        image_results_deconv.append(imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_dec_model" + special_name + save_name_append + '.fits')
        image_results_conv.append(imagename.replace('.fits', '') + "_" + str(
            ncomponents) + "C_residual" + special_name + save_name_append + ".fits")
        with open(imagename.replace('.fits', '_' + str(
                ncomponents) + 'C_de_fit' + special_name + save_name_append + '.pickle'),
                  "wb") as f:
            pickle.dump(mini_results, f)
    return (model_dict, image_results_conv, image_results_deconv)



def run_image_fitting(imagelist, residuallist, sources_photometries,
                      n_components, comp_ids=None, mask=None, mask_for_fit=None,
                      use_mask_for_fit=False,
                      indices=None, masks_deblended=None,
                      region_grow_to_ref=True,
                      which_residual='shuffled',
                      save_name_append='', z=None, aspect=None,
                      convolution_mode='GPU', workers=6,
                      method1='least_squares', method2='least_squares',
                      loss="cauchy", tr_solver="exact",
                      init_params=0.25, final_params=4.0,sigma=6,
                      fix_n=None,
                      fix_x0_y0=None,
                      fix_value_n=None,
                      fix_max_value_n=None,
                      fix_min_value_n=None,
                      fix_max_value_Rn=None,
                      fix_min_value_Rn=None,
                      fix_geometry=None,
                      force_circular=None,
                      dr_fix=None,
                      trunc=None,
                      logger=None,
                      parameters_mini_init = None,
                      self_bkg=False, bkg_map=None, rms_map=None, use_weights=False,
                      is_bkg_map_conv=False,
                      verbose=0):
    """
    Support function to run the image fitting to a image or to a list of images.

    Note. This function was implemented to  help with my own research, but it may be useable in
    some contexts. It is not a general function, and it is not well documented.

    What it does:
    For a multi-component source, this function will handle information individually for
    each component and store information in a dictionary.

    Component IDs
    -------------
    `comp_ids` and `ext_ids` are lists of INTEGERS, 1-indexed, and so is the
    `comp_ID` column of every frame returned here. Strings are accepted on input
    and converted. The resolved split is written back onto `sources_photometries`
    as `comp_ids` / `ext_ids`, which is how the caller learns what the automatic
    compact/diffuse classification decided.

    Per-region apertures
    --------------------
    `region_grow_to_ref` (default True) controls how the `region_data` /
    `region_diffuse` apertures are built from `masks_deblended`. Source extraction
    detects at a higher SNR threshold than `mask_region` is drawn at, so the
    deblended cores are much tighter than the reference aperture; growing them by a
    fixed number of dilation steps left most of `mask_region` unclaimed and the
    region fluxes could not sum to the whole-source `total` row. With the default
    the regions instead grow until they stop expanding, partitioning `mask_region`
    between them. `region_area_completeness` / `region_flux_completeness` on each
    region row say how much of the reference aperture the union actually covered.
    Pass False to recover the older fixed one-step growth.

    Returns
    -------
    A tuple whose last element is the long-format `decomp_table`; see
    `decomp_rows` for its schema.

    KNOWN LIMITATION: the two per-component property frames are built from
    `list_individual_*_props[0]` -- only the FIRST image of `imagelist`. Both
    drivers pass a single-image list, so this is currently never hit, but a
    genuine multi-image call silently loses the rest. The long `decomp_table`
    does not have this problem: it accumulates over every image in the loop.
    """

    results_fit = []
    lmfit_results = []
    lmfit_results_1st_pass = []
    errors_fit = []
    data_properties_list = []
    models = []
    list_results_compact_conv_morpho = []   # store the morphological properties
                                            # of the sum of all convolved compact components
    list_results_compact_deconv_morpho = [] # store the morphological properties
                                            # of the sum of all deconvolved compact components
    list_results_ext_conv_morpho = []       # store the morphological properties
                                            # of the sum of all convolved extended components
    list_results_ext_deconv_morpho = []     # store the morphological properties
                                            # of the sum of all deconvolved extended components
    list_individual_deconv_props = []       # store the morphological properties
                                            # of each deconvolved component.
    list_individual_conv_props = []         # store the morphological properties
                                            # of each convolved component.
    decomp_rows_all = []                    # long-format rows, one per measured
                                            # entity (see `decomp_rows`), so the
                                            # six wide frames below can be
                                            # concatenated across images and
                                            # frequencies instead of being taken
                                            # apart by hand.


    for i in range(len(imagelist)):
        #         model_dict_results = {}
        # try:
        crop_image = imagelist[i]
        crop_residual = residuallist[i]
        print('  ++==>>  Fitting', os.path.basename(crop_image))
        #             dict_results['#imagename'] = crop_image
        data_2D = load_fits_data(crop_image)
        res_2D = load_fits_data(crop_residual)
        rms_std_data = mad_std(data_2D)
        rms_std_res = mad_std(res_2D)
        if verbose >= 2:
            print('rms data = ', rms_std_data * 1e6,
                    '; rms res = ', rms_std_res * 1e6,
                    '; ratio = ', rms_std_res / rms_std_data)

        sigma_level = 3
        vmin = 3
        # i = 0 #to be used in indices[0], e.g. first component
        # omaj, omin, _, _, _ = beam_shape(crop_image)
        # dilation_size = int(
        # np.sqrt(omaj * omin) / (2 * get_cell_size(crop_image)))

        if verbose >= 1:
            do_PLOT = True
            PLOT = True
            SAVE = True
            show_figure = True
        else:
            do_PLOT = False
            PLOT = False
            SAVE = False
            show_figure = False


        if use_mask_for_fit == True:
            print('++==>> USING MASK FOR FITTTTTTTTTTTTTTTT ')
            if mask_for_fit is None:
                print('++==>> CALCULATING MASK FOR FIT IS NONE')
                _, mask_for_fit = mask_dilation(crop_image,
                                                rms=rms_std_res,
                                                sigma=sigma, dilation_size=None,
                                                iterations=6, PLOT=True)
            else:
                print('++==>> USING PROVIDED MASK FOR FIT IS NOT NONE')
                mask_for_fit = mask_for_fit
                plt.figure()
                plt.imshow(mask_for_fit,origin='lower')
                plt.show()
                plt.clf()
                plt.close()
        
        # mask_region = mask_for_fit * mask_region_i
        if mask is not None:
            _, mask_region_i = mask_dilation(crop_image,
                                            rms=rms_std_res,
                                            sigma=6.0, dilation_size=None,
                                            iterations=4, PLOT=True)
            mask_region = mask * mask_region_i
            # mask_region = mask.copy()
            # last_level = 2.5
        else:
            _, mask_region_i = mask_dilation(crop_image,
                                            rms=rms_std_res,
                                            sigma=6.0, dilation_size=None,
                                            iterations=2, PLOT=True)
            mask_region = mask_region_i.copy()
            # last_level = 1.0
        last_level = 3.0
        data_properties, _, _ = \
            measures(imagename=crop_image,
                     residualname=crop_residual, z=z,
                     sigma_mask=6.0,
                     last_level=last_level, 
                     vmin_factor=1.0,
                     dilation_size=None,
                     results_final={},
                     rms=rms_std_res,
                     apply_mask=True,
                     iterations = 2,
                     mask=mask_region,
                     do_PLOT=do_PLOT, SAVE=SAVE,
                     do_petro = True,
                     show_figure=show_figure,
                     do_measurements='partial',
                     add_save_name='_data',
                     verbose=verbose)
        data_properties_list.append(data_properties)
        # psf_image_size = dilation_size*6
        # psf_image_size = (2 * psf_image_size) // 2 +1

        psf_zise = int(get_beam_size_px(crop_image)[0])
        psf_image_size = int(psf_zise * 10)
        # psf_image_size = int(data_2D.shape[0])
        if verbose >= 2:
            print('++==>> PSF BEAM SIZE is >=> ', psf_zise)
            print('++==>> PSF IMAGE SIZE is ', psf_image_size)

        # psf_image_size = int(load_fits_data(crop_image).shape[0])

        # creates a psf from the beam shape.
        psf_name = tcreate_beam_psf(crop_image, size=(psf_image_size, psf_image_size),
                                    aspect=aspect,
                                    # aspect=None,
                                    )  # ,app_name='_'+str(psf_image_size)+'x'+str(psf_image_size)+'')

        result_mini, mini, result_1, result_extra, model_dict, \
            image_results_conv, image_results_deconv, bkg_images, \
            smodel2D, model_temp = \
                do_fit2D(imagename=crop_image,
                            residualname=crop_residual,
                            which_residual=which_residual,
                            init_constraints=sources_photometries,
                            psf_name=psf_name,
                            params_values_init_IMFIT=None,# imfit_conf_values[0:-1],
                            #fix_n = False,fix_x0_y0=[False,False,False],
                            ncomponents=n_components, constrained=True,
                            fix_n=fix_n,
                            mask_region=mask_for_fit,
                            fix_value_n=fix_value_n,
                            fix_max_value_n=fix_max_value_n,
                            fix_min_value_n=fix_min_value_n,
                            fix_max_value_Rn=fix_max_value_Rn,
                            fix_min_value_Rn=fix_min_value_Rn,
                            fix_x0_y0=fix_x0_y0,
                            dr_fix=dr_fix,
                            self_bkg=self_bkg, bkg_map=bkg_map, 
                            is_bkg_map_conv=is_bkg_map_conv,
                            rms_map=rms_map, use_weights=use_weights,
                            convolution_mode=convolution_mode,
                            fix_geometry=fix_geometry,
                            force_circular=force_circular,
                            trunc=trunc,
                            workers=workers,
                            method1=method1, method2=method2,
                            loss=loss, tr_solver=tr_solver,
                            init_params=init_params, final_params=final_params,
                            parameters_mini_init = parameters_mini_init,
                            save_name_append=save_name_append,logger=logger,verbose=verbose)

        # This block was disabled because it passed `crop_image`, a filename,
        # where fit2D_norm_metric needs an array -- so the radio path had no
        # zeta_norm at all while the general path (morphen.py) did. Both sides
        # carry the background: model_total_conv includes the fitted s_a*bkg,
        # so the data must be the unmodified image.
        try:
            model_dict['zeta_norm'] = fit2D_norm_metric(
                data=load_fits_data(crop_image),
                model=model_dict['model_total_conv'],
                mask_region=mask_region,
                background_level=rms_std_res)
        except Exception as e:
            print(f'  !!>> Could not compute zeta_norm: {e}')

        # print(result_mini.params)
        models.append(model_dict)
        lmfit_results.append(result_mini.params)
        lmfit_results_1st_pass.append(result_1.params)
        special_name = save_name_append
        bkg_deconv = bkg_images[0]
        bkg_conv = bkg_images[1]

        # _, mask_dilated_new = mask_dilation_from_mask(crop_image,
        #                                               mask_region,
        #                                               rms=rms_std_res,
        #                                               iterations=10)
        # print(image_results_deconv)
        # print(image_results_deconv[:-1])
        rms_bkg_deconv = mad_std(load_fits_data(bkg_deconv))
        deconv_model_properties = compute_model_properties(model_list=image_results_deconv[:-1],
                                                           which_model='deconv',
                                                           residualname=crop_residual,
                                                           rms=rms_std_res,
                                                           z = z,
                                                           mask_region = mask_region,
                                                           verbose=verbose
                                                           )
        rms_bkg_conv = mad_std(load_fits_data(bkg_conv))
        conv_model_properties = compute_model_properties(model_list=image_results_conv[:-2],
                                                         which_model='conv',
                                                         residualname=crop_residual,
                                                         rms=rms_std_res,
                                                         z = z,
                                                         mask_region = mask_region,
                                                         verbose=verbose
                                                         )
        list_individual_deconv_props.append(deconv_model_properties)
        list_individual_conv_props.append(conv_model_properties)

        deconv_props = pd.DataFrame(deconv_model_properties).T
        conv_props = pd.DataFrame(conv_model_properties).T
        try:
            class_results = evaluate_compactness(deconv_props, conv_props)
        except:
            #NEED-A-FIX
            class_results = {}
            pass
        # for l in class_results.keys():
        #     ID = 1
        #     deconv_props.loc[l, 'comp_ID'] = ID


        deconv_props.to_csv(image_results_deconv[-1].replace('.fits','_component_properties.csv'),
                            header=True,index=False)
        conv_props.to_csv(image_results_conv[-2].replace('.fits','_component_properties.csv'),
                          header=True, index=False)

        # comp_ids = []
        # print('*************************************')
        # print(class_results)
        # Component IDs are INTEGERS from here down. Notebooks pass them as
        # strings (`comp_ids = ['1']`), so normalise once, into a NEW list --
        # `comp_ids` used to default to a mutable `[]` that the auto-classifier
        # below appended to in place, so a second call in the same session
        # inherited the previous source's IDs.
        comp_ids = [int(c) for c in (comp_ids or [])]
        try:
            if comp_ids == []:
                ID = 1
                for key in class_results.keys():
                    if class_results[key]['final_class'] == 'C':
                        comp_ids.append(ID)
                    ID = ID + 1
            if comp_ids == []:
                comp_ids = [1]
        except:
            #NEED-A-FIX
            pass

        all_comps_ids = np.arange(1, n_components + 1)
        mask_compact_ids = np.isin(all_comps_ids, np.asarray(comp_ids))
        ext_ids = [int(e) for e in all_comps_ids[~mask_compact_ids]]

        # Hand the resolved split back to the caller. `sources_photometries` is
        # the dict the driver already owns and passed in by reference, and it is
        # already where the component bookkeeping lives (`component_types`,
        # `component_overrides`, `nIDs`). The auto-classified IDs used to reach
        # the driver only because `comp_ids` defaulted to a mutable `[]` that
        # this function appended to in place -- an accident that also leaked one
        # source's IDs into the next call.
        sources_photometries['comp_ids'] = list(comp_ids)
        sources_photometries['ext_ids'] = list(ext_ids)
        print('  ++>> Total component IDs modelled:', all_comps_ids)
        print('  ++>> IDs attributed to compact structures:', comp_ids)
        print('  ++>> IDs attributed to extended structures:', ext_ids)

        compact_model = 0
        extended_model = 0
        compact_model_deconv = 0
        extended_model_deconv = 0
        for lc in comp_ids:
            compact_model = (compact_model +
                                model_dict[f'model_c{lc}_conv']-model_dict['conv_bkg'])
            compact_model_deconv = (compact_model_deconv +
                                    model_dict[f'model_c{lc}']-model_dict['deconv_bkg'])
        # if ext_ids is not None:
        if ext_ids == []:
            extended_model = 0
            extended_model_deconv = 0
            nfunctions = 1
        else:
            for le in ext_ids:
                extended_model = (extended_model +
                                    model_dict[f'model_c{le}_conv']-model_dict['conv_bkg'])
                extended_model_deconv = (extended_model_deconv +
                                            model_dict[f'model_c{le}']-model_dict['deconv_bkg'])
                nfunctions = None
                # extended_model = extended_model - model_dict['conv_bkg'] * (len(ext_ids) - 1)
                # extended_model_deconv = extended_model_deconv - model_dict['deconv_bkg'] * (
                #             len(ext_ids) - 1)
            extended_model = extended_model + model_dict['conv_bkg']
            extended_model_deconv = extended_model_deconv + model_dict['deconv_bkg']

        # compact_model = compact_model - model_dict['conv_bkg']*(len(comp_ids)-1)
        # compact_model_deconv = compact_model_deconv - model_dict['deconv_bkg']*(len(comp_ids)-1)
        compact_model = compact_model + model_dict['conv_bkg']
        compact_model_deconv = compact_model_deconv + model_dict['deconv_bkg']
        extended_data = load_fits_data(crop_image) - compact_model
        
        if ext_ids != []:
            pf.writeto(crop_image.replace('.fits', '') +
                    "_" + "dec_ext_model" + save_name_append + ".fits",
                    extended_model_deconv, overwrite=True)
            copy_header(crop_image, crop_image.replace('.fits', '') +
                        "_" + "dec_ext_model" + save_name_append + ".fits",
                        crop_image.replace('.fits', '') +
                        "_" + "dec_ext_model" + save_name_append + ".fits")

            pf.writeto(crop_image.replace('.fits', '') +
                       "_" + "ext_model" + save_name_append + ".fits",
                       extended_model, overwrite=True)
            copy_header(crop_image, crop_image.replace('.fits', '') +
                        "_" + "ext_model" + save_name_append + ".fits",
                        crop_image.replace('.fits', '') +
                        "_" + "ext_model" + save_name_append + ".fits")

        pf.writeto(crop_image.replace('.fits', '') +
                "_" + "dec_compact" + save_name_append + ".fits",
                compact_model_deconv, overwrite=True)
        copy_header(crop_image, crop_image.replace('.fits', '') +
                    "_" + "dec_compact" + save_name_append + ".fits",
                    crop_image.replace('.fits', '') +
                    "_" + "dec_compact" + save_name_append + ".fits")

        
        exteded_file_name = crop_image.replace('.fits', '') + \
                            special_name + '_extended.fits'
        pf.writeto(exteded_file_name,extended_data + load_fits_data(bkg_conv),overwrite=True)
        copy_header(crop_image,exteded_file_name)
        
        compact_file_name = crop_image.replace('.fits', '') + \
                            special_name + '_conv_compact.fits'
        pf.writeto(compact_file_name,compact_model,overwrite=True)
        copy_header(crop_image,compact_file_name)
        

        """#testing (FOR SKA simulations)"""
        exteded_file_name_nobkg = crop_image.replace('.fits', '') + \
                            special_name + '_extended_nobkg.fits'
        pf.writeto(exteded_file_name_nobkg,extended_data-0*res_2D,overwrite=True)
        copy_header(crop_image,exteded_file_name_nobkg)
        
        # compact_file_name_nobkg = crop_image.replace('.fits', '') + \
        #                     special_name + '_conv_compact_nobkg.fits'
        # pf.writeto(compact_file_name_nobkg,compact_model-res_2D,overwrite=True)
        # copy_header(crop_image,compact_file_name_nobkg)
        
        decomp_results = plot_decomp_results(imagename=crop_image,
                                                compact=compact_model,
                                                extended_model=extended_model,
                                                bkg_image = bkg_conv,
                                                rms=rms_std_res,
                                                nfunctions=nfunctions,
                                                # Measure through the aperture
                                                # the minimiser used, and give
                                                # the printed block the fitted
                                                # parameters and the component
                                                # split it needs.
                                                mask=mask_region,
                                                result_mini=result_mini,
                                                comp_ids=comp_ids,
                                                ext_ids=ext_ids,
                                                model_total=model_dict['model_total_conv'],
                                                zeta_norm=model_dict.get('zeta_norm'),
                                                special_name=special_name)
        
        decomp_results['compact_model_image'] = compact_file_name
        decomp_results['extended_model_image'] = exteded_file_name
        print('bkg_conv')
        print(bkg_conv)
        plot_fit_results(crop_image, model_dict, image_results_conv,
                            sources_photometries,
                            bkg_image = load_fits_data(bkg_conv),
                            crop=False, box_size=200,
                            mask=mask_region,
                            vmax_factor=0.3, vmin_factor=1.0)
        # plt.xlim(0,3)
        plot_slices(load_fits_data(crop_image), load_fits_data(crop_residual), model_dict,
                    image_results_conv[-2], sources_photometries)
        
        parameter_results = result_mini.params.valuesdict().copy()

        try:
            for param in result_mini.params.valuesdict().keys():
                parameter_results[param+'_err'] = result_mini.params[param].stderr
        except:
            pass

        
        parameter_results['#imagename'] = os.path.basename(crop_image)
        parameter_results['residualname'] = os.path.basename(crop_residual)
        parameter_results['beam_size_px'] = psf_zise

        # print('++++++++++++++++++++++++++++++++++++++++')
        # print('++++++++++++++++++++++++++++++++++++++++')
        # print('++++++++++++++++++++++++++++++++++++++++')
        # print(compact_model)
        # print(mask)
        # print(mask_region)
        # print('++++++++++++++++++++++++++++++++++++++++')
        # print('++++++++++++++++++++++++++++++++++++++++')
        # print('++++++++++++++++++++++++++++++++++++++++')

        # _rms_model = mad_std(compact_model)
        # print('**************************')
        # print('**************************')
        # print('RMS MODEL COMPACT CONV:', _rms_model)
        # if _rms_model < 1e-6:
        #     rms_model = mad_std(compact_model[compact_model>1e-6])
        # else:
        #     rms_model = _rms_model


        iterations = 2
        # rms_model = mad_std(compact_model) + rms_std_res
        """testing"""
        rms_model = rms_std_res
        rms_compact_conv = rms_bkg_conv # * len(comp_ids)
        
        _, mask_region_conv_comp = mask_dilation(compact_model,
                                        rms=rms_model,
                                        # rms=rms_compact_conv,
                                        sigma=5.0,
                                        dilation_size=get_dilation_size(crop_image),
                                        # dilation_size=2,
                                        iterations=2,
                                        PLOT=PLOT,
                                        special_name=' compact conv ')#*mask_region_i

        """#testing"""
        if np.nansum(mask_region_conv_comp) == 0 and mask_region is not None:
            mask_region_conv_comp = mask_region
        """---"""

        # print('++++ Computing properties of convolved compact model.')

        # if np.sum(mask * mask_region_conv_comp) < np.sum(mask_region_conv_comp):
        #     _mask = mask_dilated_new * mask_region_conv_comp
        # else:
        #     _mask = mask


        results_compact_conv_morpho, _, _ = \
            measures(imagename=crop_image,
                     residualname=crop_residual, z=z,
                     sigma_mask=5.0,
                     last_level=1.5, vmin_factor=1.0,
                     data_2D=compact_model,
                     dilation_size=None,
                     results_final={},
                     rms=rms_model,
                    #  rms=rms_compact_conv,
                     apply_mask=False, do_PLOT=do_PLOT, SAVE=SAVE,
                     do_petro = True,
                     show_figure=show_figure,
                     mask_component=mask_region_conv_comp,
                     mask=mask_region, do_measurements='partial',
                     add_save_name='_compact_conv',verbose=verbose)

        list_results_compact_conv_morpho.append(results_compact_conv_morpho)



        # _rms_model = mad_std(compact_model_deconv)
        # print('**************************')
        # print('**************************')
        # print('RMS MODEL COMPACT DECONV:', _rms_model)
        # if _rms_model < 1e-6:
        #     rms_model = mad_std(compact_model_deconv[compact_model_deconv>1e-6])
        # else:
        #     rms_model = _rms_model

        rms_model = mad_std(compact_model_deconv) + rms_std_res
        rms_compact_deconv = rms_bkg_deconv # / len(comp_ids)
        _, mask_region_deconv_comp = mask_dilation(compact_model_deconv,
                                        rms=rms_model,
                                        # rms=rms_compact_deconv,
                                        sigma=5.0,
                                        dilation_size=2,
                                        # dilation_size=2,
                                        iterations=2,
                                        PLOT=PLOT,
                                        special_name=' compact deconv ')#*mask_region_i

        # if np.sum(mask * mask_region_deconv_comp) < np.sum(mask_region_deconv_comp):
        #     _mask = mask_dilated_new * mask_region_deconv_comp
        # else:
        #     _mask = mask

        # print('++++ Computing properties of deconvolved compact model.')
        try:
            results_compact_deconv_morpho, _, _ = \
                measures(imagename=crop_image,
                         residualname=crop_residual, z=z,
                         sigma_mask=5.0,
                         last_level=1.5, vmin_factor=1.0,
                         data_2D=compact_model_deconv,
                         dilation_size=None,
                         results_final={},
                         rms=rms_model,
                        #  rms=rms_compact_deconv,
                         apply_mask=False, do_PLOT=do_PLOT, SAVE=SAVE,
                         do_petro = True,
                         show_figure=show_figure,
                         mask_component=mask_region_deconv_comp,
                         mask=mask_region, do_measurements='partial',
                         add_save_name='_compact_deconv')

            list_results_compact_deconv_morpho.append(results_compact_deconv_morpho)
        except:
            empty_results = {key: np.nan for key in results_compact_conv_morpho.keys()}
            list_results_compact_deconv_morpho.append(empty_results)

        if nfunctions == 1:
            """
            Consider that the single component fitted represents a 
            compact component. Hence, extended emission is considered
            to be only the residual after removing that component. 
            """
            try:
                results_ext_conv_morpho, _, _ = \
                    measures(imagename=crop_image,
                             residualname=crop_residual, z=z,
                             sigma_mask=6.0,
                             last_level=1.5, vmin_factor=1.0,
                             data_2D=(load_fits_data(crop_image) - compact_model + load_fits_data(bkg_conv)) * mask_region,
                             dilation_size=None,
                             results_final={},
                             rms=rms_std_res,
                             apply_mask=False, do_PLOT=do_PLOT, SAVE=SAVE,
                             do_petro=False,
                             show_figure=show_figure,
                             # mask_component=mask_region_deconv_comp,
                             mask=mask, do_measurements='partial',
                             add_save_name='_extended_conv')

                list_results_ext_conv_morpho.append(results_ext_conv_morpho)
                results_ext_deconv_morpho = results_ext_conv_morpho
                list_results_ext_deconv_morpho.append(results_ext_deconv_morpho)
            except:
                empty_results = {key: np.nan for key in results_compact_conv_morpho.keys()}
                list_results_ext_conv_morpho.append(empty_results)
                list_results_ext_deconv_morpho.append(empty_results)


        else:
            # rms_model = mad_std(extended_model) + rms_std_res
            """testing"""
            rms_model = rms_std_res
            rms_ext_conv = rms_bkg_conv # * len(ext_ids)
            # if _rms_model < 1e-6:
            #     rms_model = mad_std(extended_model[extended_model>1e-6])
            # else:
            #     rms_model = _rms_model
            # NOTE: `mask_region_conv_ext` is NOT used by the active `measures`
            # call below -- the only reference to it is the commented-out
            # alternative further down, so the conv `diffuse_sum` row is measured
            # over ALL of `mask_region`, while the deconv one (further below) is
            # measured over `mask_region_deconv_ext`. Two apertures for the same
            # physical quantity. Left as-is deliberately: changing either would
            # move existing `diffuse_sum` numbers, which is a separate decision
            # needing its own re-validation. The call is also kept because it
            # draws a diagnostic figure when PLOT=True.
            _, mask_region_conv_ext = mask_dilation(extended_model,
                                rms=rms_model,
                                # rms=rms_ext_conv,
                                sigma=5.0,
                                dilation_size=get_dilation_size(crop_image),
                                # dilation_size=2,
                                iterations = 2,
                                PLOT=PLOT,
                                special_name=' extended conv ')#*mask_region_i
            # print('++++ Computing properties of convolved extended model.')

            # if np.sum(mask * mask_region_conv_ext) < np.sum(mask_region_conv_ext):
            #     _mask = mask_dilated_new * mask_region_conv_ext
            # else:
            #     _mask = mask
            try:
                results_ext_conv_morpho, _, _ = \
                    measures(imagename=crop_image,
                             residualname=crop_residual, z=z,
                             sigma_mask=6.0,
                             last_level=1.5, vmin_factor=1.0,
                             data_2D=(load_fits_data(crop_image) - compact_model + res_2D) * mask_region,
                             dilation_size=None,
                             results_final={},
                             rms=rms_std_res,
                             apply_mask=False, do_PLOT=do_PLOT, SAVE=SAVE,
                             do_petro=True,
                             show_figure=show_figure,
                            #  mask_component=mask_region_deconv_comp,
                             mask=mask_region, do_measurements='partial',
                             add_save_name='_extended_conv')
                # results_ext_conv_morpho, _,_ = \
                #     measures(imagename=crop_image,
                #              residualname=crop_residual, z=z,
                #              sigma_mask=5.0,
                #              last_level=1.5, vmin_factor=1.0,
                #              data_2D=extended_model,
                #              dilation_size=None,
                #              results_final={},
                #              rms=rms_model,
                #             #  rms=rms_ext_conv,
                #              apply_mask=False, do_PLOT=do_PLOT, SAVE=SAVE,
                #              do_petro=True,
                #              show_figure=show_figure,
                #              mask_component=mask_region_conv_ext,
                #              mask=mask_region, do_measurements='partial',
                #              add_save_name='_extended_conv')
                list_results_ext_conv_morpho.append(results_ext_conv_morpho)
            except:
                empty_results = {key: np.nan for key in results_compact_conv_morpho.keys()}
                list_results_ext_conv_morpho.append(empty_results)

            # rms_model = mad_std(extended_model_deconv) + rms_std_res
            """testing"""
            rms_model = rms_std_res
            rms_ext_deconv = rms_bkg_deconv # / len(ext_ids)
            # print('**************************')
            # print('**************************')
            # print('RMS MODEL EXTENDED DECONV:', _rms_model)
            # if _rms_model < 1e-6:
            #     rms_model = mad_std(extended_model_deconv[extended_model_deconv>1e-6])
            # else:
            #     rms_model = _rms_model

            _, mask_region_deconv_ext = mask_dilation(extended_model_deconv,
                                rms=rms_model,
                                # rms = rms_ext_deconv,
                                sigma=5.0,
                                dilation_size=get_dilation_size(crop_image),
                                # dilation_size=2,
                                iterations = 2,
                                PLOT=PLOT,
                                special_name=' extended deconv ')#*mask_region_i

            # if np.sum(mask * mask_region_deconv_ext) < np.sum(mask_region_deconv_ext):
            #     _mask = mask_dilated_new * mask_region_deconv_ext
            # else:
            #     _mask = mask

            try:
                # print('++++ Computing properties of deconvolved extended model.')
                results_ext_deconv_morpho, _,_ = \
                    measures(imagename=crop_image,
                             residualname=crop_residual, z=z,
                             sigma_mask=5.0,
                             last_level=1.0, vmin_factor=1.0,
                             data_2D=extended_model_deconv,
                             dilation_size=None,
                             results_final={},
                             rms=rms_model,
                            #  rms=rms_ext_deconv,
                             apply_mask=False, do_PLOT=do_PLOT, SAVE=SAVE,
                             do_petro=True,
                             show_figure=show_figure,
                             mask_component=mask_region_deconv_ext,
                             mask=mask_region, do_measurements='partial',
                             add_save_name='_extended_deconv')

                list_results_ext_deconv_morpho.append(results_ext_deconv_morpho)
            except:
                empty_results = {key: np.nan for key in results_compact_conv_morpho.keys()}
                list_results_ext_deconv_morpho.append(empty_results)

        all_results = {**parameter_results, **decomp_results}
        results_fit.append(all_results)

        """
        Per-region decomposition.

        The compact/diffuse split above is global: one `data - compact_model`
        over the whole aperture. With several detected regions that says nothing
        about which region the diffuse emission belongs to, and nothing about
        whether a given region even had a compact component fitted to it.

        Everything needed is already recorded. `c{j}_parent` (written by
        `prepare_fit` / `add_extra_component`) maps every component -- added ones
        included -- back to its detected region, and `SE.indices` / `SE.masks`
        give the regions themselves. The apertures come from
        `grow_region_masks`, the same non-overlapping growth
        `structural_morphology` uses, so region fluxes from the two paths are
        directly comparable.
        """
        _freq_i = None
        try:
            _freq_i = getfreqs([crop_image])[0]
        except:
            pass

        _row = lambda props, kind, **kw: decomp_rows(
            props, kind, imagename=crop_image, freq=_freq_i, **kw)

        decomp_rows_all += _row(data_properties, 'total', domain='data')
        for _props in deconv_model_properties.values():
            _cid = int(_props.get('comp_ID', 0))
            decomp_rows_all += _row(_props, 'component', domain='deconv',
                                    comp_ID=_cid,
                                    region_ID=int(sources_photometries.get(
                                        f'c{_cid}_parent', 0) or 0),
                                    is_compact=_cid in comp_ids)
        for _props in conv_model_properties.values():
            _cid = int(_props.get('comp_ID', 0))
            decomp_rows_all += _row(_props, 'component', domain='conv',
                                    comp_ID=_cid,
                                    region_ID=int(sources_photometries.get(
                                        f'c{_cid}_parent', 0) or 0),
                                    is_compact=_cid in comp_ids)
        decomp_rows_all += _row(results_compact_conv_morpho, 'compact_sum',
                                domain='conv')
        decomp_rows_all += _row(list_results_compact_deconv_morpho[-1],
                                'compact_sum', domain='deconv')
        decomp_rows_all += _row(list_results_ext_conv_morpho[-1], 'diffuse_sum',
                                domain='conv')
        decomp_rows_all += _row(list_results_ext_deconv_morpho[-1],
                                'diffuse_sum', domain='deconv')
        decomp_rows_all += _row(decomp_results, 'decomposition', domain='conv')

        if indices is not None and masks_deblended is not None:
            _n_regions = len(indices)
            _region_of = {}
            for _j in range(1, n_components + 1):
                _region_of[_j] = int(sources_photometries.get(f'c{_j}_parent',
                                                              _j) or _j)
            try:
                _dil = get_dilation_size(crop_image)
            except:
                _dil = 2
            # `grow_to_ref` fills `mask_region` instead of taking a fixed number
            # of steps out from each core. The cores come from source extraction,
            # which detects at a higher SNR threshold than `mask_region` is drawn
            # at, so a fixed growth leaves most of the reference aperture
            # unclaimed and the region rows cannot sum to `total`.
            _region_masks = grow_region_masks(data_2D, masks_deblended,
                                              mask_region,
                                              dilation_size=_dil,
                                              iterations=1,
                                              grow_to_ref=region_grow_to_ref)
            # Settle the final apertures BEFORE measuring anything, so the
            # coverage numbers below describe the masks actually used rather
            # than the ones `grow_region_masks` proposed.
            _ref_bool = np.asarray(mask_region).astype(bool)
            _final_masks = []
            for _r in range(1, _n_regions + 1):
                _m_r = np.asarray(_region_masks[_r - 1]).astype(bool)
                if np.nansum(_m_r) == 0:
                    # An empty `mask_component` does NOT measure nothing -- it
                    # falls through to `mask` and quietly measures most of the
                    # source, which is how a region with no aperture reported ~90%
                    # of the total flux. The compact path guards the same way (see
                    # `mask_region_conv_comp` above). Fall back to this region's
                    # own deblended core, not to all of `mask_region`, so one bad
                    # region cannot claim the whole source.
                    _m_r = (np.asarray(masks_deblended[_r - 1]).astype(bool)
                            & _ref_bool)
                _final_masks.append(_m_r)
            # How much of the reference aperture the regions between them
            # actually claimed. A shortfall is uncovered `mask_region` area --
            # pixels disconnected from every deblended core -- and is the direct
            # explanation for sum(region flux) < total_flux_mask. Same diagnostic
            # `structural_morphology` reports as `subregions_*_completeness`.
            # NOTE: with the fixed growth these masks can overlap where an empty
            # aperture fell back to an overlapping extraction ellipse, so the
            # union is a lower bound on what was measured and the region fluxes
            # can then over-count. `grow_to_ref` has no such case: it partitions.
            _covered = np.zeros(_ref_bool.shape, dtype=bool)
            _sum_areas = 0
            for _m in _final_masks:
                _covered |= _m
                _sum_areas += int(np.nansum(_m))
            # Coverage alone cannot tell a clean partition from apertures that
            # overlap: an empty aperture falling back to an extraction ellipse can
            # cover the reference mask completely AND double-count most of it.
            # This is the term that distinguishes them -- 0.0 means disjoint,
            # which `grow_to_ref` guarantees. Anything above 0 says the region
            # fluxes over-count by roughly that fraction.
            try:
                _union_area = float(np.nansum(_covered))
                _overlap_fraction = ((_sum_areas - _union_area) / _union_area
                                     if _union_area > 0 else np.nan)
            except:
                _overlap_fraction = np.nan
            try:
                _ref_area = float(np.nansum(_ref_bool))
                _area_completeness = (float(np.nansum(_covered)) / _ref_area
                                      if _ref_area > 0 else np.nan)
            except:
                _area_completeness = np.nan
            try:
                _beam_area_i = beam_area2(crop_image)
                _ref_flux = data_properties.get('total_flux_mask', np.nan)
                _covered_flux = float(np.nansum(
                    load_fits_data(crop_image) * _covered)) / _beam_area_i
                _flux_completeness = (_covered_flux / _ref_flux
                                      if _ref_flux not in (0, None)
                                      and not np.isnan(_ref_flux) else np.nan)
            except:
                _flux_completeness = np.nan
            _bkg_conv_data = load_fits_data(bkg_conv)
            for _r in range(1, _n_regions + 1):
                _mask_r = _final_masks[_r - 1]
                _comps_r = [j for j in range(1, n_components + 1)
                            if _region_of[j] == _r]
                _compact_r = [j for j in _comps_r if j in comp_ids]
                if _compact_r:
                    _model_r = 0
                    for _j in _compact_r:
                        _model_r = (_model_r + model_dict[f'model_c{_j}_conv']
                                    - model_dict['conv_bkg'])
                    _model_r = _model_r + model_dict['conv_bkg']
                    # Same construction as the global diffuse measurement
                    # above (data - compact + residual) * mask_region: the
                    # residual goes back because subtracting a noiseless model
                    # also removes noise that belongs to the data.
                    _diffuse_r = ((load_fits_data(crop_image) - _model_r
                                   + res_2D) * mask_region)
                else:
                    # Nothing compact was fitted here, so the whole region IS
                    # diffuse -- and with nothing subtracted there is no removed
                    # noise to add back. `region_diffuse` then equals
                    # `region_data` exactly, which is the honest statement.
                    _diffuse_r = load_fits_data(crop_image)
                for _kind, _data_r in (('region_data',
                                        load_fits_data(crop_image)),
                                       ('region_diffuse', _diffuse_r)):
                    try:
                        _props_r, _, _ = measures(
                            imagename=crop_image, residualname=crop_residual,
                            z=z, sigma_mask=6.0, last_level=1.5,
                            vmin_factor=1.0, data_2D=_data_r,
                            dilation_size=None, results_final={},
                            rms=rms_std_res, apply_mask=False,
                            do_PLOT=do_PLOT, SAVE=SAVE, do_petro=True,
                            show_figure=show_figure,
                            mask_component=_mask_r, mask=mask_region,
                            do_measurements='partial',
                            add_save_name=f'_{_kind}_{_r}')
                    except:
                        _props_r = {key: np.nan
                                    for key in results_compact_conv_morpho.keys()}
                    decomp_rows_all += _row(
                        _props_r, _kind, domain='data', region_ID=_r,
                        has_compact=bool(_compact_r),
                        n_comps_region=len(_comps_r),
                        region_area_completeness=_area_completeness,
                        region_flux_completeness=_flux_completeness,
                        region_overlap_fraction=_overlap_fraction)

        # except:
            # print('Error on fitting', os.path.basename(crop_image))
            # errors_fit.append(crop_image)

    return (pd.DataFrame(results_fit),result_mini,mini,
            lmfit_results, lmfit_results_1st_pass, errors_fit, models,
            pd.DataFrame(data_properties_list),
            pd.DataFrame(list_results_compact_conv_morpho),
            pd.DataFrame(list_results_compact_deconv_morpho),
            pd.DataFrame(list_results_ext_conv_morpho),
            pd.DataFrame(list_results_ext_deconv_morpho),
            # `from_dict(orient='index')` rather than `pd.DataFrame(...).T`: the
            # transpose produces the same index labels and values but erases
            # per-column dtypes, so the integer `comp_ID` written by
            # `compute_model_properties` came back as 1.0 while the string one
            # written elsewhere came back as '1'. See the `[0]` note in the
            # docstring for the multi-image caveat.
            pd.DataFrame.from_dict(list_individual_deconv_props[0],
                                   orient='index'),
            pd.DataFrame.from_dict(list_individual_conv_props[0],
                                   orient='index'),
            image_results_conv, image_results_deconv,bkg_images,
            class_results,
            compact_model,
            assemble_decomp_table(decomp_rows_all))
