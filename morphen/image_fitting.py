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
    return 2. * n - 1. / 3. + 0 * ((4. / 405.) * n) + ((46. / 25515.) * n ** 2.0)

def sersic2D(xy, x0, y0, PA, ell, n, In, Rn,cg=0.0):
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
    return (model)

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
        return 2. * n - 1. / 3. + 0 * ((4. / 405.) * n) + ((46. / 25515.) * n ** 2.0)
except:
    def bn(n):
        """
        bn function from Cioti .... (1997);
        Used to define the relation between Rn (half-light radii) and total
        luminosity

        Parameters:
            n: sersic index
        """
        return 2. * n - 1. / 3. + 0 * ((4. / 405.) * n) + ((46. / 25515.) * n ** 2.0)

try:
    @jit
    def sersic2D_GPU(xy, x0=256, y0=256, PA=10, ell=0.9,
                    n=1.0, In=0.1, Rn=10.0, cg=0.0):
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
        return (model)
except:
    def sersic2D_GPU(xy, x0=256, y0=256, PA=10, ell=0.9,
                    n=1.0, In=0.1, Rn=10.0, cg=0.0):
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
        return (model)

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
                               constrained=True, fix_n=False, fix_value_n=False,
                               fix_x0_y0=False,dr_fix = None,fix_geometry=True,
                               fix_max_value_Rn = False,
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
    dr = 10


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
                                value=0.5, min=0.49, max=0.51)
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
                                value=_PA, min=_PA - 60,
                                max=_PA + 60)
                        else:
                            smodel2D.set_param_hint(
                                'f' + str(i + 1) + '_' + param,
                                value=eval(param), min=-50.0,
                                max=190.0)
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

        smodel2D.set_param_hint('s_a', value=1, min=0.99, max=1.01)
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
                for j in range(init_constraints['ncomps']):
                    if fix_n is not False:
                        fix_n_j = fix_n[j]
                        fix_value_n_j = fix_value_n[j]
                    else:
                        fix_n_j = False
                        fix_value_n_j = 1.0 # will be skipped.

                    if fix_max_value_Rn is not False:
                        fix_max_value_Rn_j = fix_max_value_Rn[j]
                    else:
                        fix_max_value_Rn_j = False
                        

                    if fix_x0_y0 is not False:
                        fix_x0_y0_j = fix_x0_y0[j]
                        dr_fix_j = dr_fix[j]
                    else:
                        fix_x0_y0_j = False
                        dr_fix_j = False

                    if fix_geometry is not False:
                        fix_geometry_j = fix_geometry[j]
                    else:
                        fix_geometry_j = False

                    jj = str(j + 1)
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
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=0.5, min=0.3, max=8.0)

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
                            dO = 90
                            _PA = init_constraints['c' + jj + '_PA']
                            PA_max = _PA + dO
                            PA_min = _PA - dO
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=_PA, min=PA_min, max=PA_max)
                        if param == 'ell':
                            ell = 1 - init_constraints['c' + jj + '_q']
                            ell_min = ell * 0.2
                            #                         if ell + dell <= 1.0:
                            # if ell * 1.0 <= 0.5:
                            if ell <= 0.15:
                                ell_max = 0.15
                            else:
                                ell_max = 0.75
                                
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
                            if fix_geometry_j == True:
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=0.0, min=-0.01, max=0.01)
                            else:
                                print('Using general elliptical geometry during '
                                      'fitting... may take longer.')
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=0.0, min=-2.0, max=2.0)

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
                                I50_max = I50 * 100
                                I50_min = I50 * 0.1
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=I50, min=I50_min, max=I50_max)
                            else:
                                I50_max = I50 * 25
                                I50_min = I50 * 0.05
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=I50, min=I50_min, max=I50_max)
                        if param == 'Rn':
                            if fix_max_value_Rn_j is not False:
                                """
                                Fix the value of Rn to a maximum value.
                                This is useful if one wants to fit a delta function.
                                That is, Rn=1.0 ~ delta functio, 1 pixel component.
                                """
                                print(f" ++==>> Limiting {param}_{jj} to {fix_max_value_Rn_j}")
                                smodel2D.set_param_hint(
                                    'f' + str(j + 1) + '_' + param,
                                    value=fix_max_value_Rn_j*0.9, 
                                    min=0.1, 
                                    max=fix_max_value_Rn_j)
                                # smodel2D.set_param_hint(
                                #     'f' + str(j + 1) + '_' + param,
                                #     value=fix_max_value_Rn_j, 
                                #     min=fix_max_value_Rn_j-0.01, 
                                #     max=fix_max_value_Rn_j+0.01)
                                
                            else:
                                R50 = init_constraints['c' + jj + '_R50']
                                # R50_max = R50 * 4.0
                                # R50_max = init_constraints['c' + jj + '_Rp']
                                if observation_type == 'radio':
                                    R50_max = R50 * 1.0
                                    R50_min = R50 * 0.01 #should be small.
                                    smodel2D.set_param_hint(
                                        'f' + str(j + 1) + '_' + param,
                                        value=R50, min=R50_min, max=R50_max)
                                else:
                                    R50_max = R50 * 2.0
                                    R50_min = R50 * 0.5 #should be small.
                                    smodel2D.set_param_hint(
                                        'f' + str(j + 1) + '_' + param,
                                        value=R50, min=R50_min, max=R50_max)

                        if param == 'x0':
                            if fix_x0_y0_j is not False:
                                """
                                Fix centre position by no more than dr_fix.
                                """
                                x0c = init_constraints['c' + jj + '_x0c']
                                x0_max = x0c + dr_fix_j
                                x0_min = x0c - dr_fix_j
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
                            print('Limiting ', param)
                            smodel2D.set_param_hint(
                                'f' + str(j + 1) + '_' + param,
                                value=y0c,
                                min=y0_min,
                                max=y0_max)

            smodel2D.set_param_hint('s_a', value=1, min=0.99, max=1.01)
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
                smodel2D.set_param_hint('s_a', value=1, min=0.99, max=1.01)
                # smodel2D.set_param_hint('s_a', value=1, min=0.3, max=6.0)
            except:
                print('Please, if not providing initial parameters file,')
                print('provide basic information for the source.')
                return (ValueError)

    params = smodel2D.make_params()
    # print(smodel2D.param_hints)
    
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

def add_extra_component(petro_properties, copy_from_id):
    """
    Create another component from a dictionary (petro_properties) having
    photometric properties for N detected components in an image.

    Params
    ------
    petro_properties: dict
        Contain parameters of a number o N components obtained
        by a petrosian analysis of all detected sources.
        Example (these are actually the keys() from the dictionary):
        ['c1_PA', 'c1_q', 'c1_area', 'c1_Re',
        'c1_x0c', 'c1_y0c', 'c1_label', 'c1_R50',
        'c1_Snu', 'c1_Rp', 'c1_Rpidx', 'c1_rlast',
        'c1_I50']
    copy_from: int
        From which component copy parameters from.
        This is useful, for example, the source has two components detected,
        1 compact and the other a structure that can not be modelled by a single
        sersic function. Then, we need one function to model the compact structure,
        but 2 sersic functions to model the other structure.

        Assume that we have a blob surrounded by a disky emission ( detected
        as one source). Both are placed on the same region, on top of each other
        (e.g. from optical, we can call for example] a bulge and
        a disk). We need two functions to model this region.

        So, if component i=1 is the blob (or the bulge) we copy the parameters from it and
        create a second component. We just have to ajust some of the parameters.
        E.g. the effective radius of this new component, is in principle, larger than the original component.
        As well, the effective intensity will be smaller because we are adding a component
        further away from the centre. Other quantities, however, are uncertain, such as the Sersic index, position angle
        etc, but may be (or not!) close to those of component i.

    """

    from collections import OrderedDict
    dict_keys = list(petro_properties.keys())
    unique_list = list(OrderedDict.fromkeys(
        [elem.split('_')[1] for elem in dict_keys if '_' in elem]))
    #     print(unique_list)
    print(unique_list)
    petro_properties_copy = petro_properties.copy()
    new_comp_id = petro_properties['ncomps'] + 1
    for k in range(len(unique_list)):
        #         print(unique_list[k])
        # do not change anything for other parameters.
        petro_properties_copy['c' + str(new_comp_id) + '_' + unique_list[k]] = \
        petro_properties_copy['c' + str(copy_from_id) + '_' + unique_list[k]]
        if unique_list[k] == 'R50':
            # multiply the R50 value by a factor, e.g., 5.0
            # factor = 5.0
            factor = np.min([petro_properties['cg_Rp'],2*petro_properties['cg_Rp'] / abs(petro_properties['cg_Rp']-petro_properties_copy['c' + str(copy_from_id) + '_' + unique_list[k]])])
            print(factor )
            petro_properties_copy[
                'c' + str(new_comp_id) + '_' + unique_list[k]] = \
            petro_properties_copy[
                'c' + str(copy_from_id) + '_' + unique_list[k]] * factor
        if unique_list[k] == 'I50':
            # multiply the I50 value by a factor, e.g., 0.2
            factor = 0.01
            petro_properties_copy[
                'c' + str(new_comp_id) + '_' + unique_list[k]] = \
            petro_properties_copy[
                'c' + str(copy_from_id) + '_' + unique_list[k]] * factor
        if unique_list[k] == 'Rp':
            petro_properties_copy[
                'c' + str(new_comp_id) + '_' + unique_list[k]] = petro_properties['cg_Rp']
                    
    # update number of components
    petro_properties_copy['ncomps'] = petro_properties_copy['ncomps'] + 1
    return (petro_properties_copy)



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
    

 

def phot_source_ext(imagename, residual=None, sigma=1.0, iterations=2, dilation_size=None,
                    deblend_nthresh=5, deblend_cont=1e-6, maskthresh=0.0,
                    gain=1, filter_kernel=None, mask=None,
                    segmentation_map=False, clean_param=1.0, clean=True,
                    minarea=100, minarea_factor=1, npixels = None,
                    filter_type='matched', 
                    sort_by='distance', threshold_mode='sigma',
                    bw=64, bh=64, fw=3, fh=3, ell_size_factor=None,
                    apply_mask=False, sigma_mask=6,
                    show_bkg_map=False, show_detection=False,
                    SE_ref=None):
    """
    Simple source extraction algorithm (using SEP https://sep.readthedocs.io/en/v1.1.x/).


    """

    data_2D = load_fits_data(imagename)
    if len(data_2D.shape) == 4:
        data_2D = data_2D[0][0]
    # m, s = np.mean(data_2D), np.std(data_2D)
    if residual is not None:
        m, s = np.mean(data_2D), mad_std(residual)
    else:
        m, s = np.mean(data_2D), mad_std(data_2D)
    # bkg = 0.0
    if apply_mask:
        _, mask = mask_dilation(data_2D, sigma=sigma_mask, iterations=iterations,
                                rms=s,
                                PLOT=True,show_figure=True,
                                dilation_size=dilation_size)


    bkg = sep.Background(data_2D, mask=mask, bw=bw, bh=bh, fw=fw, fh=fh)
    # print(bkg.globalback)
    # print(bkg.globalrms)

    bkg_image = bkg.back()
    bkg_rms = bkg.rms()
    
    data_sub = data_2D - bkg_image

    if show_bkg_map == True:
        # plt.figure()
        # # display bkg map.
        # plt.imshow(data_2D, interpolation='nearest', cmap='gray', vmin=3*s,
        #            vmax=0.2*np.max(data_2D), origin='lower')
        # plt.colorbar()
        # plt.close()
        # plt.show()
        # plt.figure()
        # plt.imshow(bkg_image)
        # plt.colorbar()
        # plt.show()
        # plt.close()
        
        fig = plt.figure(figsize=(18, 6))
        ax0 = fig.add_subplot(1, 3, 2)
        ax0 = eimshow(data_sub,ax=ax0,fig=fig,
                            plot_colorbar=True,
                            # cbar_orientation = 'horizontal',
                            show_axis='off',
                            vmin_factor=0.5,
                            cbar_n_points = 3,
                            add_contours=True)
        ax0.set_title(f'data - bkg')
        
        ax1 = fig.add_subplot(1, 3, 1)
        ax1 = eimshow(data_2D,ax=ax1,fig=fig,
                            plot_colorbar=True,
                            # cbar_orientation = 'horizontal',
                            show_axis='off',
                            vmin_factor=0.5,
                            cbar_n_points=3,
                            plot_title = f'data',
                            add_contours=True)
        ax2 = fig.add_subplot(1, 3, 3)
        ax2 = eimshow(bkg_image,ax=ax2,fig=fig,
                            plot_colorbar=True,
                            # cbar_orientation = 'horizontal',
                            show_axis='off',
                            # vmin_factor=0.5,
                            # vmax_factor=1.0,
                            cbar_n_points=3,
                            plot_title = f'bkg',
                            add_contours=False)
        plt.show()
        plt.clf()
        plt.close()

    
    if mask is not None:
        data_sub = data_sub * mask
    else:
        data_sub = data_sub
        
    if (threshold_mode == 'residual') and residual is not None:
        threshold = sigma * residual
    elif threshold_mode == 'sigma':
        threshold = sigma * s
            
    # else:
    #     mask = None
    # print(data_sub)
    if npixels is None:
        npixels = int(minarea * minarea_factor)
    # print(' INFO: Uinsg min number of pixels of :', npixels)
    cat, segm, seg_maps = make_catalog(image=data_sub,
                                       threshold=threshold,
                                       deblend=True, 
                                       contrast=deblend_cont,
                                       nlevels=deblend_nthresh,
                                       npixels=npixels,
                                       figsize=(20, 20),
                                       plot=show_detection, vmin=1.0 * s)
    # cat, segm, seg_maps = make_catalog(image=data_sub,
    #                                    threshold=sigma * s,
    #                                    deblend=True, 
    #                                    contrast=deblend_cont,
    #                                    nlevels=deblend_nthresh,
    #                                    npixels=npixels,
    #                                    figsize=(20, 20),
    #                                    plot=show_detection, vmin=1.0 * s)
    if sort_by == 'flux':
        indices = list(order_cat(cat, key='segment_flux', reverse=True))
    if sort_by == 'area':
        indices = list(order_cat(cat, key='area', reverse=True))
    if sort_by == 'distance':
        ref_centre = data_2D.shape[0] / 2, data_2D.shape[1] / 2
        if SE_ref is not None:
            coords, distances, indices = \
                sorted_detected_coordinates(SE_ref.objects['xc'],
                                            SE_ref.objects['yc'],
                                            cat.xcentroid, 
                                            cat.ycentroid, 
                                            ref_centre)
            print(indices)
        else:
            
            distances = distances_from_reference(cat.xcentroid, 
                                            cat.ycentroid,
                                            ref_centre)
            indices = np.argsort(distances)
        
    masks_deblended = []
    for k in range(len(indices)):
        # print(k)
        masks_deblended.append(seg_maps == seg_maps.labels[indices[k]])


    # len(objects)
    from matplotlib.patches import Ellipse
    from skimage.draw import ellipse

    # m, s = np.mean(data_sub), np.std(data_sub)
    if show_detection == True:
        fig, ax = plt.subplots(figsize=(10, 10))
        norm = simple_norm(data_sub, stretch='sqrt', asinh_a=0.02, min_cut=s,
                    max_cut=0.2*np.nanmax(data_sub))
        im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
                       norm=norm, origin='lower')
        # im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
        #                vmin=s, vmax=0.2*np.nanmax(data_sub), origin='lower')

    masks_regions = []

    if ell_size_factor is None:
        if mask is not None:
            # ell_size_factor = np.sqrt(np.sum(mask) / (np.pi))/cat[0].equivalent_radius.value
            ell_size_factor = 0.05*np.sqrt(np.sum(mask) / (np.pi))
        else:
            ell_size_factor = 0.5
    
    y, x = np.indices(data_2D.shape[:2])
    
    # objects = cat
    for i in range(len(cat)):
        source = cat[i]
        seg_mask = (seg_maps.data == i + 1)
        e = Ellipse(xy=(source.centroid[0], source.centroid[1]),
                    width=1 * ell_size_factor * source.equivalent_radius.value,
                    height=1 * ell_size_factor * (
                                1 - source.ellipticity.value) * source.equivalent_radius.value,
                    # angle=source.orientation.value * 180. / np.pi
                    angle=source.orientation.value
                    )

        xc = source.centroid[0]
        yc = source.centroid[1]
        a = ell_size_factor * source.equivalent_radius.value
        b = ell_size_factor * (
                    1 - source.ellipticity.value) * source.equivalent_radius.value
        theta = source.orientation.value
        rx = (x - xc) * np.cos(theta) + (y - yc) * np.sin(theta)
        ry = (y - yc) * np.cos(theta) - (x - xc) * np.sin(theta)

        inside = ((rx / a) ** 2 + (ry / b) ** 2) <= 1
        mask_ell = np.zeros_like(data_2D)
        mask_ell[inside] = True
        if show_detection == True:
            e.set_facecolor('none')
            e.set_edgecolor('red')
            ax.add_artist(e)
        masks_regions.append(seg_mask)

    #         plt.savefig('components_SEP.pdf',dpi=300, bbox_inches='tight')
    # flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
    #                                      3.0, err=bkg.globalrms, gain=1.0)
    # for i in range(len(objects)):
    #     print("object {:d}: flux = {:f} +/- {:f}".format(i, flux[i], fluxerr[i]))
    # objects['b'] / objects['a'], np.rad2deg(objects['theta'])

    # sort regions from largest size to smallest size.
    mask_areas = []
    mask_fluxes = []
    for mask_comp in masks_regions:
        area_mask = np.sum(mask_comp)
        sum_mask = np.sum(mask_comp * data_2D)
        mask_areas.append(area_mask)
        mask_fluxes.append(sum_mask)
    mask_areas = np.asarray(mask_areas)
    mask_fluxes = np.asarray(mask_fluxes)
    if sort_by == 'area':
        sorted_indices_desc = np.argsort(mask_areas)[::-1]
        sorted_arr_desc = mask_areas[sorted_indices_desc]
    if sort_by == 'flux':
        sorted_indices_desc = np.argsort(mask_fluxes)[::-1]
        sorted_arr_desc = mask_fluxes[sorted_indices_desc]
    if sort_by == 'distance':
        ref_centre = data_2D.shape[0] / 2, data_2D.shape[1] / 2
        if SE_ref is not None:
            coords, distances, sorted_indices_desc = \
                sorted_detected_coordinates(SE_ref.objects['xc'],
                                            SE_ref.objects['yc'],
                                            cat.xcentroid, 
                                            cat.ycentroid, 
                                            ref_centre)
            # sorted_indices_desc = np.argsort(distances)
            # sorted_arr_desc = distances[sorted_indices_desc]
                
        else:
            distances = distances_from_reference(cat.xcentroid, 
                                                cat.ycentroid,
                                                ref_centre)
            sorted_indices_desc = np.argsort(distances)
            # sorted_arr_desc = distances[sorted_indices_desc]

    objects_sorted = {}
    objects_sorted['xc'] = np.asarray([1] * len(cat))
    objects_sorted['yc'] = np.asarray([1] * len(cat))
    for i in range(len(cat)):
        source = cat[sorted_indices_desc[i]]
        objects_sorted['xc'][i] = source.centroid[0]
        objects_sorted['yc'][i] = source.centroid[1]

    if show_detection == True:
        for i in range(len(cat)):
            source = cat[sorted_indices_desc[i]]
            xc = source.centroid[0]
            yc = source.centroid[1]
            label = str('ID' + str(i + 1))
            
            label_x = xc + 3 * ell_size_factor
            label_y = yc + 20 * ell_size_factor
            
            line_end_y = label_y - 5 
            
            text = Text(label_x, label_y, label, ha='center', va='center', color='red')
            ax.add_artist(text)
            
            ax.plot([xc, label_x], [yc, line_end_y], color='red', alpha=0.4)

        plt.axis('off')
        plt.savefig(imagename + '_SEP_phot.jpg', dpi=300, bbox_inches='tight')
        plt.show()



    if segmentation_map == True:
        return (masks_deblended, sorted_indices_desc, bkg,
                seg_maps, objects_sorted)
    else:
        return (masks_deblended, sorted_indices_desc, bkg,
                objects_sorted)

def prepare_fit(ref_image, ref_res, z, ids_to_add=[1],
                bw=51, bh=51, fw=15, fh=15, sigma=15, ell_size_factor=2.0,
                deblend_cont=1e-7, deblend_nthresh=15,
                minarea=None,minarea_factor=1.0,npixels=None,
                sigma_mask=6,mask=None,dilation_size=None,
                show_detection=True,use_extraction_positions=False,
                clean_param=0.9,clean=True,sort_by='distance',
                apply_mask=False, mask_grow_iterations=3,
                obs_type = 'radio',algorithm='SEP',force_circular=True,
                show_petro_plots=False):
    """
    Prepare the imaging data to be modeled.

    This function runs a source extraction, computes basic petrosian properties
    from the data for each detected source as well as shape morphology
    (e.g. position angle, axis ration, effective intensity and radii).
    """
    crop_image = ref_image
    crop_residual = ref_res
    data_2D = load_fits_data(crop_image)
    if crop_residual is not None:
        residual_2D = load_fits_data(crop_residual)
    else:
        residual_2D = None
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

    # if apply_mask == True:
    #     mask_detection = mask
    # else:
    #     mask_detection = None #np.ones(data_2D.shape)
    # # plt.figure()

    # _, mask = mask_dilation(crop_image, sigma=6, dilation_size=None,
    #                         iterations=2)
    if algorithm == 'SEP':
        masks, indices, bkg, seg_maps, objects = \
            sep_source_ext(crop_image, 
                           residualname=crop_residual,
                           bw=bw,bh=bh,fw=fw, fh=fh,
                           minarea=minarea,
                           minarea_factor=minarea_factor,
                           segmentation_map=True,
                           filter_type='matched',
                           mask = mask_detection,
                           deblend_nthresh=deblend_nthresh,
                           deblend_cont=deblend_cont,
                           clean_param=clean_param,
                           clean=clean,
                           sort_by=sort_by,
                           npixels = npixels,
                           dilation_size=dilation_size,
                           iterations = mask_grow_iterations,
                           sigma=sigma,sigma_mask=sigma_mask,
                           ell_size_factor=ell_size_factor,
                           apply_mask=apply_mask,
                           show_detection=show_detection)
    if algorithm == 'PF':
        masks, indices, bkg, seg_maps, objects = \
            phot_source_ext(crop_image, 
                           residual=residual_2D,
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
                           sigma=sigma,sigma_mask=sigma_mask,
                           iterations = mask_grow_iterations,
                           ell_size_factor=ell_size_factor,
                           apply_mask=apply_mask,
                           show_detection=show_detection)



    sigma_level = 3
    vmin = 3
    # i = 0 #to be used in indices[0], e.g. first component
    sources_photometries = {}  # init dict to store values.
    # if use_extraction_positions == True:
    #     for i in range(len(indices)):
    #         # ii = str(i+1)
    #         positions = np.array([objects['xc'][i], objects['yc'][i]])
    #         mask_component = masks[indices[i]]
    #         data_component = data_2D * mask_component
    #         sources_photometries = compute_petro_source(data_component,
    #                                                     mask_component=mask_component,
    #                                                     sigma_level=1,positions=positions,
    #                                                     i=i, plot=show_petro_plots,
    #                                                     source_props=sources_photometries)

    # else:
    
    
    sources_photometries = compute_petro_source(data_2D = data_2D,
                                                imagename = crop_image,
                                                sigma_level=6,
                                                nlevels=1, contrast=1,
                                                deblend=False, npixels=None,
                                                i='g', mask_component = mask_detection,
                                                plot=show_petro_plots,
                                                source_props=sources_photometries)
    
    for i in tqdm(range(len(indices))):
        # ii = str(i+1)
        mask_component = masks[i]
        data_component = data_2D * mask_component
        sources_photometries = compute_petro_source(data_2D = data_component,
                                                    imagename = crop_image,
                                                    mask_component=mask_component,
                                                    obs_type = obs_type,
                                                    sigma_level=1,nlevels=1, contrast=1,
                                                    deblend=False, npixels=None,
                                                    i=i, plot=show_petro_plots,
                                                    source_props=sources_photometries)


    if obs_type == 'radio':
        """
        PSF image is contained within the header of the original image. 
        A new psf file will be created (as `psf_name`). 
        """
        # omaj, omin, _, _, _ = beam_shape(crop_image)
        # dilation_size = int(
        #     np.sqrt(omaj * omin) / (2 * get_cell_size(crop_image)))
        # psf_image_size = dilation_size*6
        # psf_image_size = (2 * psf_image_size) // 2 +1
        psf_image_size = int(data_2D.shape[0])
        # print('++==>> PSF IMAGE SIZE is', psf_image_size)
        # creates a psf from the beam shape.
        psf_name = tcreate_beam_psf(crop_image, size=(
            psf_image_size, psf_image_size))  # ,app_name='_'+str(psf_image_size)+'x'+str(psf_image_size)+'')
    if obs_type == 'other':
        """
        Provide a psf file.
        """
        psf_name = None

    n_components = len(indices)
    n_IDs = len(indices)
    sources_photometries['ncomps'] = n_components
    sources_photometries['nIDs'] = n_IDs

    
    print("# of structures (IDs) to be fitted =", n_IDs)
    # sources_photometies_new = sources_photometies
    # n_components_new = n_components
    if ids_to_add is not None:
        for id_to_add in ids_to_add:
            sources_photometries = add_extra_component(sources_photometries,
                                                       copy_from_id=id_to_add)
    
    if force_circular:
        print(f"<> Forcing components to be circular.")
        # for i in range(len(indices)):
        for i in range(sources_photometries['ncomps']):
            if i+1 <= n_IDs:
                sources_photometries[f'c{i+1}_q'] = 0.99
            print(f"<> q_{i+1} {sources_photometries[f'c{i+1}_q']}")
    else:
        print(f"<> Leaving components to be elliptical.")
        for i in range(sources_photometries['ncomps']):
            if i+1 > n_IDs:
                sources_photometries[f'c{i+1}_q'] = 0.6
                # sources_photometries[f'c{i+1}_R50'] = sources_photometries[f'cg_Rp']
                # sources_photometries[f'c{i+1}_Rp'] = sources_photometries[f'cg_Rp']
            print('<> q = ',sources_photometries[f'c{i+1}_q'])

    # update variable `n_components`.
    n_components = sources_photometries['ncomps']
    print("# of model components (COMPS) to be fitted =", n_components)
    return (sources_photometries, n_components, n_IDs, masks, indices, objects,
            psf_name, mask, bkg)

def do_fit2D(imagename, params_values_init_IMFIT=None, ncomponents=None,
             init_constraints=None, data_2D_=None, residualdata_2D_=None,
             residualname=None,which_residual='shuffled',observation_type = 'radio',
             init_params=0.25, final_params=4.0, constrained=True,
             fix_n=True, fix_value_n=False, 
             fix_max_value_Rn=False,
             dr_fix=3,
             fix_x0_y0=False, psf_name=None, convolution_mode='GPU',
             convolve_cutout=False, cut_size=512, self_bkg=False, rms_map=None,
             fix_geometry=True, contrain_nelder=False, workers=6,mask_region = None,
             special_name='', method1='least_squares', method2='least_squares',
             reduce_fcn='neglogcauchy',loss="cauchy",tr_solver="exact",x_scale = 'jac',
             ftol=1e-10, xtol=1e-10, gtol=1e-10, verbose=2,max_nfev=200000,
             regularize  = True, f_scale = 1.0,
             maxiter = 30000, maxfev = 30000, xatol = 1e-12,
             fatol = 1e-12, return_all = True, disp = True,
             de_options=None,parameters_mini_init=None,
             save_name_append='',logger=None):
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
    self_bkg: bool
        If True, use the image background as the residual background.
    rms_map: 2D array
        RMS map to be used for the fitting.
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


    returns
    -------
    result: dict
        Dictionary containing the results of the fitting.


    """
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
    FlatSky_level = None
    if data_2D_ is None:
        data_2D = pf.getdata(imagename)
    else:
        data_2D = data_2D_

    if mask_region is not None:
        """
        
        """
        logger.debug(f" ==> Using provided mask region to constrain fit. ")
        logger.warning(f" !!==> Fitting with a mask is experimental! ")
        # data_2D = data_2D * mask_region
        if convolution_mode == 'GPU':
            mask_for_fit = jnp.array(mask_region)
        else:
            mask_for_fit = mask_region
    else:
        mask_for_fit = None

    if convolution_mode == 'GPU':
        data_2D_gpu = jnp.array(data_2D)

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

        if convolution_mode == 'CPU':
            PSF_DATA = PSF_DATA_raw
        # PSF_DATA = pf.getdata(
        #     imagename.replace('-image.cutout.fits', '-beampsf.cutout.fits'))
    else:
        PSF_CONV = False
        PSF_DATA = None

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
            residual_2D = pf.getdata(residualname)

        if which_residual == 'shuffled':
            if logger is not None:
                logger.debug(f" ==> Using clean shuffled background for optmization... ")
            residual_2D_to_use = shuffle_2D(residual_2D)
        if which_residual == 'natural':
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

        FlatSky_level = mad_std(residual_2D_to_use)
        #         background = residual_2D #residual_2D_to_use
        if convolution_mode == 'GPU':
            background = jnp.array(residual_2D_to_use)
        else:
            background = residual_2D_to_use

    else:
        if which_residual == 'user':
            if rms_map is None:
                print('--==>> A rms map/background mode was selected (user)')
                print('       but no rms/background map was provided.')
                print('       Please, provide a rms/background map.')
                print('||==>> Stopping code now.')
                raise ValueError("rms_map should not be None")
            else:
                logger.debug(f" ==> Using provided RMS map. ")
                background_map=rms_map
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

    if residualdata_2D_ is not None:
        residual_2D = residualdata_2D_
    else:
        try:
            residual_2D = pf.getdata(residualname)
        except:
            residual_2D = background

    size = data_2D.shape
    if convolution_mode == 'GPU':
        x,y = jnp.meshgrid(jnp.arange((size[1])), jnp.arange((size[0])))
        xy = jnp.stack([x, y], axis=0)
    else:
        xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))

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
    if convolution_mode == 'GPU':
        residual_2D = jnp.array(residual_2D)

    nfunctions = ncomponents

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
                                     params['f' + str(i) + '_cg'], )
        # print(model.shape)
        # model = model + FlatSky_cpu(FlatSky_level, params['s_a'])*background
        # model = model + FlatSky_cpu(background, params['s_a'])
        MODEL_2D_conv = scipy.signal.fftconvolve(model, PSF_DATA, 'same') + \
                        FlatSky_cpu(background, params['s_a'])
        residual = data_2D - MODEL_2D_conv
        return np.ravel(residual)

    # @partial(jit, static_argnums=2)
    # def build_model(xy,_params,nfunctions):
    #     params = func(_params[:-1])
    #     model = 0
    #     for i in range(0, nfunctions):
    #         params_i = params[i]
    #         print(params_i[5])
    #         model = model + sersic2D_GPU(xy, params_i[0],
    #                                      params_i[1],
    #                                      params_i[2],
    #                                      params_i[3],
    #                                      params_i[4],
    #                                      params_i[5],
    #                                      params_i[6],
    #                                      params_i[7])
    #     return model
    """
    def extract_params(params):
        param_array = jnp.array(
            [params[f'f{i}_{name}'].value for i in range(1, nfunctions + 1) for name in
             ['x0', 'y0', 'PA', 'ell', 'n', 'In', 'Rn', 'cg']])
        param_matrix = np.split(param_array, nfunctions)
        return param_matrix
    # @partial(jit, static_argnums=2)
    # @jit
    # def build_model(xy,param_matrix):
    #     model = 0
    #     for model_params in param_matrix:
    #         model = model + sersic2D_GPU(xy, model_params[0],
    #                                      model_params[1],
    #                                      model_params[2],
    #                                      model_params[3],
    #                                      model_params[4],
    #                                      model_params[5],
    #                                      model_params[6],
    #                                      model_params[7])
    #     return model
    # batched_sersic2D_GPU = vmap(sersic2D_GPU, in_axes=(None, 0))
    # @jit
    # def build_model(xy,param_matrix):
    #     # Use the batched version of sersic2D_GPU to compute all models in parallel
    #     models = batched_sersic2D_GPU(xy, param_matrix)
    #     # Sum the models to get the final model
    #     total_model = models.sum(axis=0)
    #     return total_model


    # @jit
    # def sersic2D_GPU_vectorized(xy, params):
    #     # Assuming sersic2D_GPU is compatible with JAX and can work with batched inputs
    #     # Unpack params directly within the function if necessary
    #     return sersic2D_GPU_new(xy, params)
    #
    # batched_sersic2D_GPU = vmap(sersic2D_GPU_vectorized, in_axes=(None, 0))
    #
    # @jit
    # def build_model(xy, param_matrix):
    #     # Compute all models in parallel
    #     models = batched_sersic2D_GPU(xy, param_matrix)
    #     # Sum the models to get the final model
    #     total_model = models.sum(axis=0)
    #     return total_model
    def min_min_residual_2D_GPU(params):
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
                                         params['f' + str(i) + '_cg'].value)
        # # param_matrix = extract_params(params)
        # param_matrix = func(jnp.array(list(params.valuesdict().values()))[:-1])
        # model = build_model(xy,param_matrix)
        # print(params.values)
        # print(params.valuesdict().values())
        # _params = convert_params_to_numpy(np.asarray(list(params.valuesdict().values())))
        # print(jnp.array(list(params.valuesdict().values())))
        # model = build_model(xy,
        #                     jnp.array(list(params.valuesdict().values())),
        #                     nfunctions)
        # print(model.shape)
        # MODEL_2D_conv = convolve_on_gpu(model, PSF_DATA)
        # MODEL_2D_conv = jax_convolve(model, PSF_DATA)
        # model = model + FlatSky(background, params['s_a'].value)
        # MODEL_2D_conv = _fftconvolve_jax(model, PSF_DATA) + \
        #                 FlatSky(background,params['s_a'].value)
        MODEL_2D_conv = _fftconvolve_jax(model+
                                         FlatSky(background,params['s_a'].value),
                                         PSF_DATA)
        residual = data_2D_gpu - MODEL_2D_conv

        # return np.asarray(residual).copy()
        # return np.asarray(residual).copy().flatten()
        # return np.asarray(residual+0.01*abs(jnp.nanmin(residual))).copy()
        return np.asarray(residual).copy()
    """
    def convert_params_to_numpy(_params):
        return list(_params)

    try:
        # @partial(jit, static_argnums=1)
        @jit
        def func(x):
            return jnp.split(x, nfunctions)

        @jit
        def build_model(xy,param_matrix):
            model = 0
            for model_params in param_matrix:
                model = model + sersic2D_GPU(xy, model_params[0],
                                            model_params[1],
                                            model_params[2],
                                            model_params[3],
                                            model_params[4],
                                            model_params[5],
                                            model_params[6],
                                            model_params[7])
            return model
    
    except:
        # @partial(jit, static_argnums=1)
        
        def func(x):
            return jnp.split(x, nfunctions)

        
        def build_model(xy,param_matrix):
            model = 0
            for model_params in param_matrix:
                model = model + sersic2D_GPU(xy, model_params[0],
                                            model_params[1],
                                            model_params[2],
                                            model_params[3],
                                            model_params[4],
                                            model_params[5],
                                            model_params[6],
                                            model_params[7])
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
                                         params['f' + str(i) + '_cg'].value)

        # # param_matrix = extract_params(params)
        # param_matrix = func(jnp.array(list(params.valuesdict().values()))[:-1])
        # model = build_model(xy,param_matrix)

        MODEL_2D_conv = _fftconvolve_jax(model+
                                         FlatSky(background,params['s_a'].value),
                                         PSF_DATA)
        # residual = ((data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])/
        #             (1000*(abs(residual_2D[mask_for_fit])+1.0e-6)))
        residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit])
        return np.asarray(residual).copy()
        # weights = 1/((background[mask_for_fit])/data_2D_gpu[mask_for_fit])
        # weightned_residual = (data_2D_gpu[mask_for_fit] - MODEL_2D_conv[mask_for_fit]) * jnp.sqrt(weights)
        # return np.asarray(weightned_residual).copy()
    
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


    smodel2D, params = construct_model_parameters(
        params_values_init_IMFIT=params_values_init_IMFIT, n_components=nfunctions,
        init_constraints=init_constraints,observation_type=observation_type,
        fix_n=fix_n, fix_value_n=fix_value_n,
        fix_max_value_Rn = fix_max_value_Rn,
        fix_x0_y0=fix_x0_y0, dr_fix=dr_fix, fix_geometry=fix_geometry,
        init_params=init_params, final_params=final_params,
        constrained=constrained)

    if convolution_mode == 'CPU':
        mini = lmfit.Minimizer(min_residual_2D, params, max_nfev=200000,
                               nan_policy='omit', reduce_fcn=reduce_fcn)
    if convolution_mode == 'GPU':
        mini = lmfit.Minimizer(min_residual_2D_GPU, params, max_nfev=200000,
                               nan_policy='omit', reduce_fcn=reduce_fcn)

    # initial minimization.

    print(' >> Using', method1, ' solver for first optimisation run... ')
    # take parameters from previous run, and re-optimize them.
    #     me``rthod2 = 'ampgo'#'least_squares'
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


    if method1 == 'least_squares':
        # faster, but usually not good for first run.
        # if results_previous_run is not None:
        print(' >> Using',tr_solver,'for tr solver, with regularize set to',regularize,
              ' Loss is',loss,'.')
        
        if parameters_mini_init is not None:
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

    if method1 == 'differential_evolution':
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

    if method2 == 'ampgo':
        # ampgo is not workin well/ takes so long ???
        result = mini.minimize(method='ampgo', params=second_run_params,
                               maxfunevals=10000, totaliter=30, disp=True,
                               maxiter=5, glbtol=1e-8)

    if method2 == 'least_squares':
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
                               x_scale=x_scale,  f_scale=f_scale,
                               ftol=ftol, xtol=xtol, gtol=gtol, verbose=verbose,
                               loss=loss)  # ,f_scale=0.5, max_nfev=5000, verbose=2)

    if method2 == 'differential_evolution':
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

    params = result.params

    model_temp = Model(sersic2D)
    xy = np.meshgrid(np.arange((size[1])), np.arange((size[0])))
    model = 0
    model_dict = {}
    image_results_conv = []
    image_results_deconv = []
    total_image_results_deconv = []
    total_image_results_conv = []
    bkg_images = []
    flat_sky_total = FlatSky_cpu(background, params['s_a'].value)
    bkg_comp_i = flat_sky_total.copy()
    for i in range(1, ncomponents + 1):
        model_temp = sersic2D_GPU(xy, params['f' + str(i) + '_x0'].value,
                              params['f' + str(i) + '_y0'].value,
                              params['f' + str(i) + '_PA'].value,
                              params['f' + str(i) + '_ell'].value,
                              params['f' + str(i) + '_n'].value,
                              params['f' + str(i) + '_In'].value,
                              params['f' + str(i) + '_Rn'].value,
                              params['f' + str(i) + '_cg'].value)

        model = model + model_temp
        #to each individual component, add the bkg map.
        model_dict['model_c' + str(i)] = np.asarray(model_temp+bkg_comp_i).copy()
        if PSF_CONV == True:
            if convolution_mode == 'GPU':
                # model_dict['model_c' + str(i) + '_conv'] = np.asarray(jax_convolve(model_temp, PSF_DATA)).copy()
                # model_dict['model_c' + str(i) + '_conv'] = (
                #     np.asarray(_fftconvolve_jax(model_temp,PSF_DATA).copy()+bkg_comp_i))
                # to each individual component, add the bkg map.
                model_dict['model_c' + str(i) + '_conv'] = (
                    np.asarray(_fftconvolve_jax(model_temp+bkg_comp_i,PSF_DATA).copy()))
            if convolution_mode == 'CPU':
                model_dict['model_c' + str(i) + '_conv'] = (
                        scipy.signal.fftconvolve(model_temp+bkg_comp_i, PSF_DATA_raw,'same'))

        else:
            model_dict['model_c' + str(i) + '_conv'] = model_temp+bkg_comp_i

        pf.writeto(imagename.replace('.fits', '') +
                   "_" + "model_component_" + str(i) +
                   special_name + save_name_append + '.fits',
                   model_dict['model_c' + str(i) + '_conv'], overwrite=True)
        copy_header(imagename, imagename.replace('.fits', '') +
                    "_" + "model_component_" + str(i) +
                    special_name + save_name_append + '.fits',
                    imagename.replace('.fits', '') +
                    "_" + "model_component_" + str(
                        i) + special_name + save_name_append + '.fits')
        pf.writeto(imagename.replace('.fits', '') +
                   "_" + "dec_model_component_" + str(i) +
                   special_name + save_name_append + '.fits',
                   model_dict['model_c' + str(i)], overwrite=True)
        copy_header(imagename, imagename.replace('.fits', '') +
                    "_" + "dec_model_component_" + str(i) +
                    special_name + save_name_append + '.fits',
                    imagename.replace('.fits', '') +
                    "_" + "dec_model_component_" + str(i) +
                    special_name + save_name_append + '.fits')

        image_results_conv.append(imagename.replace('.fits', '') +
                                  "_" + "model_component_" + str(i) +
                                  special_name + save_name_append + '.fits')
        image_results_deconv.append(imagename.replace('.fits', '') +
                                    "_" + "dec_model_component_" + str(i) +
                                    special_name + save_name_append + '.fits')

    #     model = model
    model_dict['model_total_dec'] = np.asarray(model+flat_sky_total) # +FlatSky_cpu(background,
    # params['s_a'].value)

    if PSF_CONV == True:
        # model_dict['model_total_conv'] = scipy.signal.fftconvolve(model,
        #                                                           PSF_DATA_raw,
        #                                                           'same')  # + FlatSky(FlatSky_level, params['s_a'])
        if convolution_mode == 'GPU':
            # model_dict['model_total_conv'] = np.asarray(jax_convolve(model,
            #                                                          PSF_DATA)).copy()
            # model_conv = _fftconvolve_jax(model, PSF_DATA).copy() + FlatSky_cpu(background, params['s_a'].value
            model_conv = _fftconvolve_jax(model+flat_sky_total,PSF_DATA).copy()
        if convolution_mode == 'CPU':
            model_conv = scipy.signal.fftconvolve(model+flat_sky_total, PSF_DATA_raw,'same')
            model_dict['model_total_conv'] = model_conv
    else:
        model_dict['model_total_conv'] = model + flat_sky_total


    # model_dict['best_residual'] = data_2D - model_dict['model_total']
    # bkg_comp_total
    model_dict['model_total_conv'] = np.asarray(model_conv)
    model_dict['best_residual_conv'] = np.asarray(data_2D) - model_dict['model_total_conv']


    model_dict['deconv_bkg'] = np.asarray(flat_sky_total)
    # bkg_comp_total
    model_dict['conv_bkg'] = np.asarray(_fftconvolve_jax(flat_sky_total,PSF_DATA).copy())

    pf.writeto(imagename.replace('.fits', '') +
               "_" + "conv_model" + special_name + save_name_append + '.fits',
               model_dict['model_total_conv'], overwrite=True)

    total_image_results_conv.append(imagename.replace('.fits', '') +
               "_" + "conv_model" + special_name + save_name_append + '.fits')

    pf.writeto(imagename.replace('.fits', '') +
               "_" + "dec_model" + special_name + save_name_append + '.fits',
               model_dict['model_total_dec'], overwrite=True)

    total_image_results_deconv.append(imagename.replace('.fits', '') +
               "_" + "dec_model" + special_name + save_name_append + '.fits')


    pf.writeto(imagename.replace('.fits', '') +
               "_" + "residual" + special_name + save_name_append + ".fits",
               model_dict['best_residual_conv'], overwrite=True)
    copy_header(imagename, imagename.replace('.fits', '') +
                "_" + "conv_model" + special_name + save_name_append + '.fits',
                imagename.replace('.fits', '') +
                "_" + "conv_model" + special_name + save_name_append + '.fits')
    copy_header(imagename, imagename.replace('.fits', '') +
                "_" + "dec_model" + special_name + save_name_append + '.fits',
                imagename.replace('.fits', '') +
                "_" + "dec_model" + special_name + save_name_append + '.fits')
    copy_header(imagename, imagename.replace('.fits', '') +
                "_" + "residual" + special_name + save_name_append + '.fits',
                imagename.replace('.fits', '') +
                "_" + "residual" + special_name + save_name_append + '.fits')

    pf.writeto(imagename.replace('.fits', '') +
               "_" + "dec_model" + special_name + save_name_append + '.fits',
               model_dict['model_total_dec'], overwrite=True)
    copy_header(imagename, imagename.replace('.fits', '') +
                "_" + "dec_model" + special_name + save_name_append + '.fits',
                imagename.replace('.fits', '') +
                "_" + "dec_model" + special_name + save_name_append + '.fits')

    pf.writeto(imagename.replace('.fits', '') +
               "_" + "deconv_bkg" + special_name + save_name_append + '.fits',
               model_dict['deconv_bkg'], overwrite=True)
    copy_header(imagename, imagename.replace('.fits', '') +
                "_" + "deconv_bkg" + special_name + save_name_append + '.fits',
                imagename.replace('.fits', '') +
                "_" + "deconv_bkg" + special_name + save_name_append + '.fits')

    pf.writeto(imagename.replace('.fits', '') +
               "_" + "conv_bkg" + special_name + save_name_append + '.fits',
               model_dict['conv_bkg'], overwrite=True)
    copy_header(imagename, imagename.replace('.fits', '') +
                "_" + "conv_bkg" + special_name + save_name_append + '.fits',
                imagename.replace('.fits', '') +
                "_" + "conv_bkg" + special_name + save_name_append + '.fits')
    bkg_images.append(imagename.replace('.fits', '') +
               "_" + "deconv_bkg" + special_name + save_name_append + '.fits')
    bkg_images.append(imagename.replace('.fits', '') +
               "_" + "conv_bkg" + special_name + save_name_append + '.fits')


    # # initial minimization.
    # method1 = 'differential_evolution'
    # print(' >> Using', method1, ' solver for first optimisation run... ')
    # # take parameters from previous run, and re-optimize them.
    # #     method2 = 'ampgo'#'least_squares'
    # method2 = 'least_squares'



    image_results_conv.append(imagename.replace('.fits', '') +
                              "_" + "conv_model" +
                              special_name + save_name_append + '.fits')
    image_results_deconv.append(imagename.replace('.fits', '') +
                                "_" + "dec_model" +
                                special_name + save_name_append + '.fits')
    image_results_conv.append(imagename.replace('.fits', '') +
                              "_" + "residual" +
                              special_name + save_name_append + ".fits")

    # save mini results (full) to a pickle file.
    with open(imagename.replace('.fits',
                                '_' + 'fit' +
                                special_name + save_name_append + '.pickle'),
              "wb") as f:
        pickle.dump(result, f)

    with open(imagename.replace('.fits',
                                '_' + 'fit' +
                                special_name + save_name_append + '_modeldict.pickle'),
              "wb") as f:
        pickle.dump(model_dict, f)


    exec_time = time.time() - startTime
    print('Exec time fitting=', exec_time, 's')



    # save results to csv file.
    try:
        save_results_csv(result_mini=result,
                         save_name=image_results_conv[-2].replace('.fits', ''),
                         ext='.csv',
                         save_corr=True, save_params=True)
    except:
        print('Error Saving Results to a csv file!!!')
        pass

    return (result, mini, result_1, result_extra, model_dict, image_results_conv,
            image_results_deconv, bkg_images, smodel2D, model_temp)


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
                      n_components, comp_ids=[], mask=None, mask_for_fit=None,
                      use_mask_for_fit=False,
                      which_residual='shuffled',
                      save_name_append='', z=None, aspect=None,
                      convolution_mode='GPU', workers=6,
                      method1='least_squares', method2='least_squares',
                      loss="cauchy", tr_solver="exact",
                      init_params=0.25, final_params=4.0,sigma=6,
                      fix_n=[True, True, True, True, True, True, False],
                      fix_x0_y0 = [True, True, True, True, True, True, True, True],
                      fix_value_n=[0.5, 0.5, 0.5, 1.0], 
                      fix_max_value_Rn = [False, False, False, False],
                      fix_geometry=True,
                      dr_fix=[10, 10, 10, 10, 10, 10, 10, 10],logger=None,
                      parameters_mini_init = None,
                      self_bkg=False, bkg_rms_map=None,verbose=0):
    """
    Support function to run the image fitting to a image or to a list of images.

    Note. This function was implemented to  help with my own research, but it may be useable in
    some contexts. It is not a general function, and it is not well documented.

    What it does:
    For a multi-component source, this function will handle information individually for
    each component and store information in a dictionary.



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

        _, mask_region_i = mask_dilation(crop_image,
                                        rms=rms_std_res,
                                        sigma=5.0, dilation_size=None,
                                        iterations=4, PLOT=True)
        
        mask_region = mask_for_fit.copy() * mask_region_i
        
        data_properties, _, _ = \
            measures(imagename=crop_image,
                     residualname=crop_residual, z=z,
                     sigma_mask=6.0,
                     last_level=1.0, vmin_factor=1.0,
                     dilation_size=None,
                     results_final={},
                     rms=rms_std_res,
                     apply_mask=True,
                     iterations = 3,
                     mask=mask_region_i*mask_region,
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
                            fix_max_value_Rn=fix_max_value_Rn,
                            fix_x0_y0=fix_x0_y0,
                            dr_fix=dr_fix,
                            self_bkg=self_bkg, rms_map=bkg_rms_map,
                            convolution_mode=convolution_mode,
                            fix_geometry=fix_geometry,
                            workers=workers,
                            method1=method1, method2=method2,
                            loss=loss, tr_solver=tr_solver,
                            init_params=init_params, final_params=final_params,
                            parameters_mini_init = parameters_mini_init,
                            save_name_append=save_name_append,logger=logger)

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

        rms_bkg_deconv = mad_std(load_fits_data(bkg_deconv))
        deconv_model_properties = compute_model_properties(model_list=image_results_deconv[:-1],
                                                           which_model='deconv',
                                                           residualname=crop_residual,
                                                           rms=rms_std_res,
                                                           mask_region = mask_region,
                                                           verbose=verbose
                                                           )
        rms_bkg_conv = mad_std(load_fits_data(bkg_conv))
        conv_model_properties = compute_model_properties(model_list=image_results_conv[:-2],
                                                         which_model='conv',
                                                         residualname=crop_residual,
                                                         rms=rms_std_res,
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
        try:
            if comp_ids == []:
                ID = 1
                for key in class_results.keys():
                    if class_results[key]['final_class'] == 'C':
                        comp_ids.append(str(ID))
                    ID = ID + 1
            if comp_ids == []:
                comp_ids = ['1']
        except:
            #NEED-A-FIX
            pass

        all_comps_ids = np.arange(1, n_components + 1).astype('str')
        mask_compact_ids = np.isin(all_comps_ids, np.asarray(comp_ids))
        ext_ids = list(all_comps_ids[~mask_compact_ids])
        print('  ++>> Total component IDs modelled:', all_comps_ids)
        print('  ++>> IDs attributed to compact structures:', comp_ids)
        print('  ++>> IDs attributed to extended structures:', ext_ids)

        compact_model = 0
        extended_model = 0
        compact_model_deconv = 0
        extended_model_deconv = 0
        for lc in comp_ids:
            compact_model = (compact_model +
                                model_dict['model_c' + lc + '_conv']-model_dict['conv_bkg'])
            compact_model_deconv = (compact_model_deconv +
                                    model_dict['model_c' + lc]-model_dict['deconv_bkg'])
        # if ext_ids is not None:
        if ext_ids == []:
            extended_model = 0
            extended_model_deconv = 0
            nfunctions = 1
        else:
            for le in ext_ids:
                extended_model = (extended_model +
                                    model_dict['model_c' + le + '_conv']-model_dict['conv_bkg'])
                extended_model_deconv = (extended_model_deconv +
                                            model_dict['model_c' + le]-model_dict['deconv_bkg'])
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
        pf.writeto(exteded_file_name,extended_data,overwrite=True)
        copy_header(crop_image,exteded_file_name)
        
        compact_file_name = crop_image.replace('.fits', '') + \
                            special_name + '_conv_compact.fits'
        pf.writeto(compact_file_name,compact_model,overwrite=True)
        copy_header(crop_image,compact_file_name)
        
        
        
        decomp_results = plot_decomp_results(imagename=crop_image,
                                                compact=compact_model,
                                                extended_model=extended_model,
                                                rms=rms_std_res,
                                                nfunctions=nfunctions,
                                                special_name=special_name)
        
        decomp_results['compact_model_image'] = compact_file_name
        decomp_results['extended_model_image'] = exteded_file_name

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
        rms_model = mad_std(compact_model) + rms_std_res
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



        # print('++++ Computing properties of convolved compact model.')

        # if np.sum(mask * mask_region_conv_comp) < np.sum(mask_region_conv_comp):
        #     _mask = mask_dilated_new * mask_region_conv_comp
        # else:
        #     _mask = mask


        results_compact_conv_morpho, _, _ = \
            measures(imagename=crop_image,
                     residualname=crop_residual, z=z,
                     sigma_mask=5.0,
                     last_level=1.0, vmin_factor=1.0,
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
                         last_level=1.0, vmin_factor=1.0,
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
                             data_2D=(load_fits_data(crop_image) - compact_model) * mask_region,
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
            rms_model = mad_std(extended_model) + rms_std_res
            rms_ext_conv = rms_bkg_conv # * len(ext_ids)
            # if _rms_model < 1e-6:
            #     rms_model = mad_std(extended_model[extended_model>1e-6])
            # else:
            #     rms_model = _rms_model
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
                             data_2D=(load_fits_data(crop_image) - compact_model) * mask_region,
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

            rms_model = mad_std(extended_model_deconv) + rms_std_res
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
            pd.DataFrame(list_individual_deconv_props[0]).T,
            pd.DataFrame(list_individual_conv_props[0]).T,
            image_results_conv, image_results_deconv,bkg_images,
            class_results,
            compact_model)
