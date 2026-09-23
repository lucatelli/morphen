"""
                                                          ..___|**_
                                                  .|||||||||*+@+*__*++.
                                              _||||.           .*+;].,#_
                                         _|||*_                _    .@@@#@.
                                   _|||||_               .@##@#| _||_
       Morphen                |****_                   .@.,/\..@_.
                             #///#+++*|    .       .@@@;#.,.\@.
                              .||__|**|||||*||*+@#];_.  ;,;_
 Geferson Lucatelli                            +\*_.__|**#
                                              |..      .]]
                                               ;@       @.*.
                                                #|       _;]];|.
                                                 ]_          _+;]@.
                                                 _/_             |]\|    .  _
                                              ...._@* __ .....     ]]+ ..   _
                                                  .. .       . .. .|.|_ ..

"""
__versions__ = ('0.3.1alpha-1', '0.4.0alpha-1', '0.5.0alpha-1', '0.7.0alpha-1', '0.8.0alpha-1','1.0.0alpha-1')
__codenames__ = ('Pelicoto', 'Pelicoto', 'Pelicoto', 'Saurinho', 'Goba', 'Lito')
__package_name__ = 'morphen'
__dates__ =  ('2024 03 25','2024 11 13', '2024 12 18', '2025 11 04', '2026 02', '2026 09')
__version__ = '1.0.0alpha-1'
__codename__ = 'Lito'
__author__ = 'Geferson Lucatelli'
# __coauthors__ = ('Javier Moldon, Rob Beswick, '
                #   'Fabricio Ferrari, Leonardo Ferreira')
__email__ = 'gefersonlucatelli@gmail.com; gefersonlucatelli@furg.br'
__date__ = '2026 09'
print(__doc__)
print('Version',__version__, '('+__codename__+')')
print('By',__author__)
print('Date',__date__)
import os
import sys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from matplotlib.text import Text
from matplotlib.patches import Ellipse
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.text import Text
from matplotlib import rcParams
from matplotlib import colors
from matplotlib.ticker import ScalarFormatter
from matplotlib.offsetbox import AnchoredText
from matplotlib.patches import Ellipse
import matplotlib.figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
import io
import gc
from IPython.display import Image, display

import numpy as np
# np.set_printoptions(precision=4, suppress=True, linewidth=100)
np.set_printoptions(legacy='1.21')
# from sympy import *
import casatasks
from casatasks import *
import casatools
import casaviewer.imview as imview
# from casatools import *
from scipy.ndimage import rotate
import matplotlib.pyplot as plt
import matplotlib as mpl
import astropy.io.fits as pf
from astropy.coordinates import SkyCoord
import astropy.units as u
from casatools import image as IA


import lmfit
from lmfit import Model
from lmfit import Parameters, fit_report, minimize
import emcee
from prettytable import PrettyTable
from joblib import Parallel, delayed


import string
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl_
import glob
from astropy.nddata import Cutout2D
from astropy.wcs import WCS
import sep
import fitsio

from astropy.stats import mad_std
from scipy.ndimage import gaussian_filter
from astropy import visualization
from astropy.visualization import simple_norm
from photutils.segmentation import detect_sources, detect_threshold
from astropy.convolution import convolve, Gaussian2DKernel
from astropy.stats import SigmaClip, sigma_clipped_stats
from photutils.background import Background2D, MedianBackground, MADStdBackgroundRMS, SExtractorBackground, MMMBackground, StdBackgroundRMS, LocalBackground
from skimage.measure import perimeter_crofton
from scipy import ndimage
from scipy.ndimage import morphology
from scipy.ndimage import shift
from scipy.special import huber
from skimage.morphology import disk, square
from skimage.morphology import dilation
from skimage.segmentation import watershed
from skimage.filters import gaussian
from scipy.spatial import ConvexHull
from itertools import combinations


from scipy.optimize import leastsq, fmin, curve_fit
import scipy.ndimage as nd
import scipy
from scipy.stats import circmean, circstd
from scipy.signal import savgol_filter



from astropy.cosmology import FlatLambdaCDM
from astropy import units as u
from astropy import coordinates
import pandas as pd
import sys
import pickle
import time
import corner
import re
from tqdm import tqdm
import itertools

import dynesty
from dynesty import plotting as dyplot
import corner



from scipy import ndimage
from sklearn.neighbors import KNeighborsClassifier


# import pymc3 as pm
# try:
from petrofit import make_radius_list
from petrofit import Petrosian
from petrofit import source_photometry
from petrofit import make_catalog, plot_segments
from petrofit import plot_segment_residual
from petrofit import order_cat
# except:
#     pass

# --- compat patch: petrofit 0.6.0 `radial_photometry` calls float() on the
# result of photutils Aperture.do_photometry(), which numpy >= 2.0 rejects
# for any array with ndim > 0 (even single-element ones), and which
# photutils >= 3.0 always returns as an ndarray rather than a scalar.
# Delete this once petrofit ships a fix for numpy 2.x / photutils 3.x.
import petrofit.photometry as _pf_photometry
import petrofit.segmentation as _pf_segmentation


def _scalar(x):
    return float(np.asarray(x).reshape(-1)[0])


def _patched_radial_photometry(
    image,
    position,
    r_list,
    error=None,
    mask=None,
    elong=1.0,
    theta=0.0,
    plot=False,
    vmin=0,
    vmax=None,
    method="exact",
):
    flux_arr = []
    error_arr = []
    area_arr = []

    if plot:
        ax = plt.gca()
        plt.imshow(image, vmin=vmin, vmax=image.mean() * 10 if vmax is None else vmax)
        ax.set_title("Image and Aperture Radii")
        ax.set_xlabel("Pixels")
        ax.set_ylabel("Pixels")

    mask = ~mask if mask is not None else None
    for i, r in enumerate(r_list):
        aperture = _pf_photometry.radial_elliptical_aperture(
            position, r, elong=elong, theta=theta
        )

        photometric_value, photometric_err = aperture.do_photometry(
            data=image, error=error, mask=mask, method=method
        )
        aperture_area, aperture_area_err = aperture.do_photometry(
            data=np.ones_like(image), error=None, mask=mask, method=method
        )

        aperture_area = _scalar(np.round(aperture_area, 6))
        photometric_value = _scalar(np.round(photometric_value, 6))
        photometric_err = (
            _scalar(np.round(photometric_err, 6)) if photometric_err.size > 0 else np.nan
        )

        if np.isnan(photometric_value):
            raise Exception("Nan photometric_value")

        if plot:
            aperture.plot(plt.gca(), color="w", alpha=0.5)

        flux_arr.append(photometric_value)
        area_arr.append(aperture_area)
        error_arr.append(photometric_err)

    return np.array(flux_arr), np.array(area_arr), np.array(error_arr)


_pf_photometry.radial_photometry = _patched_radial_photometry
_pf_segmentation.radial_photometry = _patched_radial_photometry
import copy
# from copy import copy
import astropy.io.fits as fits
import matplotlib.ticker as mticker
import coloredlogs
import logging
import warnings
warnings.filterwarnings('ignore', module='photutils')
warnings.filterwarnings('ignore', module='astropy')

from functools import partial
try:
    import jax
    from jax import jit, vmap
    from jax.numpy.fft import fft2, ifft2, fftshift
    import jax.numpy as jnp
    import jax.scipy as jscipy
except:
    print('Jax was not imported/installed correctly, Sersic Fitting Will FAIL! ')
    print('Jax/GPU Libraries not imported.')
    pass



# sys.path.append('../../scripts/analysis_scripts/')
sys.path.append('./analysis_scripts/')

morphen_path = os.path.dirname(os.path.abspath(__file__))
print(f' > {__package_name__} path: {morphen_path}')
# libs_path = os.path.join(current_dir, "config.py")


def exec_module_with_tracking(module_path, module_name, namespace):
    """
    Execute a Python file with proper file tracking for better error messages.
    
    Parameters
    ----------
    module_path : str
        Full path to the .py file
    module_name : str
        Name of the module (for display)
    namespace : dict
        Namespace to execute in (usually globals())
    """
    with open(module_path, 'r') as f:
        code_string = f.read()
    
    # Compile with the actual filename - this preserves file info in tracebacks
    compiled_code = compile(code_string, module_path, 'exec')
    
    # Execute in the provided namespace
    exec(compiled_code, namespace)


modules_to_load = [
    ('config.py', 'config'),
    ('utils.py', 'utils'),
    ('data_io.py', 'data_io'),
    ('image_alignment.py', 'image_alignment'),
    ('fit_ellipse.py', 'fit_ellipse'),
    ('cosmo.py', 'cosmo'),
    ('image_fitting.py', 'image_fitting'),
    ('image_morphometry.py', 'image_morphometry'),
    ('image_photometry.py', 'image_photometry'),
    ('plotting.py', 'plotting'),
    ('alignment_viz.py', 'alignment_viz'),
    ('radio_sed.py', 'radio_sed'),
    ('radio_utils.py', 'radio_utils'),
    ('signal_stats.py', 'signal_stats'),
    ('source_extraction.py', 'source_extraction'),
    ('field_extraction.py', 'field_extraction'),
    ('utils_analysis.py', 'utils_analysis'),
    ('testing_deploy.py', 'testing_deploy'),
]

for module_file, module_name in modules_to_load:
    module_path = os.path.join(morphen_path, module_file)
    if os.path.exists(module_path):
        exec_module_with_tracking(module_path, module_name, globals())
    else:
        print(f"    Warning: {module_file} not found")


reset_rc_params()
import fit_ellipse
from fit_ellipse import fit_ellipse_to_galaxy