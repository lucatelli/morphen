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
__versions__ = ('0.3.1alpha-1', '0.4.0alpha-1', '0.5.0alpha-1')
__codenames__ = ('Pelicoto', 'Saurinho', '')
__package_name__ = 'morphen'
__version__ = '0.5.0alpha-2'
__codename__ = 'Saurinho'
__author__ = 'Geferson Lucatelli'
__email__ = 'geferson.lucatelli@postgrad.manchester.ac.uk; gefersonlucatelli@gmail.com'
__date__ = '2025 04 01'
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


import numpy as np
from sympy import *
import casatasks
from casatasks import *
import casatools
import casaviewer.imview as imview
# from casatools import *
from scipy.ndimage import rotate
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
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
from astropy.convolution import Gaussian2DKernel
from skimage.measure import perimeter_crofton
from scipy import ndimage
from scipy.ndimage import morphology
from scipy.ndimage import shift
from scipy.special import huber
from skimage.morphology import disk, square
from skimage.morphology import dilation
from scipy.spatial import ConvexHull
from itertools import combinations


from scipy.optimize import leastsq, fmin, curve_fit
import scipy.ndimage as nd
import scipy
from scipy.stats import circmean, circstd
from scipy.signal import savgol_filter


from astropy.cosmology import FlatLambdaCDM
import numpy as np
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



import numpy as np
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



"""
This is a temporary fix to import all modules into the main namespace (morphen).
This will change in the future. Something like this:

    >>> from config import *
    >>> from fit_ellipse import *
    >>> from utils import *
    >>> from data_io import *

"""
exec(open(f"{morphen_path}/config.py").read())
exec(open(f"{morphen_path}/utils.py").read())
exec(open(f"{morphen_path}/data_io.py").read())
exec(open(f"{morphen_path}/fit_ellipse.py").read())
exec(open(f"{morphen_path}/cosmo.py").read())
exec(open(f"{morphen_path}/image_fitting.py").read())
exec(open(f"{morphen_path}/image_morphometry.py").read())
exec(open(f"{morphen_path}/image_photometry.py").read())
exec(open(f"{morphen_path}/plotting.py").read())
exec(open(f"{morphen_path}/radio_sed.py").read())
exec(open(f"{morphen_path}/radio_utils.py").read())
exec(open(f"{morphen_path}/signal_stats.py").read())
exec(open(f"{morphen_path}/source_extraction.py").read())
exec(open(f"{morphen_path}/utils_analysis.py").read())
exec(open(f"{morphen_path}/testing_deploy.py").read())


reset_rc_params()

import fit_ellipse
from fit_ellipse import fit_ellipse_to_galaxy
