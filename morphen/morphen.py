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
# print(__doc__)


import argparse
import atexit
import math
import os
import sys
import matplotlib as mpl
import logging
from matplotlib import use as mpluse
# sys.path.append("/mirror/scratch/lucatelli/app/miniconda3/envs/casa6/lib/python3.8/site-packages/")
sys.path.append('./')
sys.path.append('./analysis_scripts/')
# import sys
import mlibs as mlibs
import analysisUtils as au
from analysisUtils import *
import coloredlogs

# F_nu[Jy] = F_lambda[erg/s/cm2/Angstrom] * lambda[Angstrom]**2 * 1e23 / c[Angstrom/s].
# With c in Angstrom/s the constant is 1e23 / 2.99792458e18 = 33356.4095, and it
# MULTIPLIES -- an earlier version of the HST branch divided by a rounded 3.34e4
# instead, which is wrong by 3.34e4**2 = 1.1156e9 (22.62 mag). Cross-check: WFC3/IR
# F160W carries both keywords, and PHOTFLAM * PHOTPLAM**2 * this constant reproduces
# its PHOTFNU (1.52665582e-07) to seven significant figures. The exact value is used
# rather than 3.34e4 because the rounding alone costs 0.0014 mag.
_FLAM_TO_FNU_JY = 1e23 / 2.99792458e18

class config():
    """
    Configuration Class to specify basic parameters.
    """

    def reset_rc_params():
        """
        Global configuration for matplotlib.pyplot
        """
        # global_font_size = 14
        global_font_size = 14
        mpl.rcParams.update({'font.size': global_font_size,
                             'text.usetex': False, 
                             'font.family': 'sans-serif',
                             'mathtext.fontset': 'stix',
                             'font.family': 'sans',
                             'font.weight': 'medium',  
                             'font.family': 'STIXGeneral',
                            #  'text.usetex' : True,
                            #  'font.family' : 'serif',
                            #  'font.serif' : ['Garamond Libre', 'EB Garamond', 'Cormorant Garamond', 'serif'],
                            #  'text.latex.preamble': r'''
                            #     \usepackage{ebgaramond-maths}
                            #     \usepackage{garamondlibre}
                            #     %\usepackage{accanthis}
                            #     \usepackage{amsmath}
                            #     \usepackage{amssymb}
                            #     \usepackage{mathrsfs}
                            #     \DeclareFontFamily{U}{BOONDOX-calo}{\skewchar\font=45}
                            #     \DeclareFontShape{U}{BOONDOX-calo}{m}{n}{<-> s*[1.05] BOONDOX-r-calo}{}
                            #     \DeclareFontShape{U}{BOONDOX-calo}{b}{n}{<-> s*[1.05] BOONDOX-b-calo}{}
                            #     \DeclareMathAlphabet{\mcb}{U}{BOONDOX-calo}{m}{n}
                            #     \SetMathAlphabet{\mcb}{bold}{U}{BOONDOX-calo}{b}{n}
                            #     \DeclareMathAlphabet{\mbcb}{U}{BOONDOX-calo}{b}{n}
                            #     ''',
                             'xtick.labelsize': global_font_size,
                             'figure.figsize': (6, 4),
                             'ytick.labelsize': global_font_size,
                             'axes.labelsize': global_font_size,
                             'xtick.major.width': 1,
                             'ytick.major.width': 1,
                             'axes.linewidth': 1.5,
                             'axes.edgecolor':'orange',
                             'lines.linewidth': 2,
                             'legend.fontsize': global_font_size,
                             'grid.linestyle': '--',
                             # 'grid.color':'black',
                             #  'figure.dpi': 96,
                             'axes.grid.which': 'major',  
                             'axes.grid.axis': 'both', 
                             'axes.spines.right': True,
                             'axes.grid': True,
                             'axes.titlesize' : global_font_size,
                             'legend.framealpha': 1.0
                             })
        # mpl.rcParams['text.usetex'] = True
        # mpl.rcParams['font.family'] = 'serif'
        # mpl.rcParams['font.serif'] = ['Garamond Libre', 'EB Garamond', 'Cormorant Garamond', 'serif']
        # # LaTeX preamble with your specific packages
        # mpl.rcParams['text.latex.preamble'] = r'''
        # \usepackage{ebgaramond-maths}
        # \usepackage{garamondlibre}
        # \usepackage{amsmath}
        # \usepackage{amssymb}
        # '''
        pass

    reset_rc_params()
    sigma=3
    mask_iterations = 1
    show_plots = True
    ext = '.jpg'
    log_file_name = 'logfile.log'

    if "--noshow" in sys.argv:
        mpluse('Agg')

    def __init__(self):
        print("Initializing Morphen")

# class _logging_():
#     def __init__(self,log_file_name):
#         self.log_file_name = log_file_name.replace('.fits','.log')
#         self.start_log()
#
#
#     def start_log(self):
#         self.logger = logging.getLogger(__name__)
#         # Set the log level
#         self.logger.setLevel(logging.DEBUG)
#         # Fancy format
#         log_format = "%(asctime)s - %(levelname)s - %(message)s"
#         # Use colored logs to add color to the log messages
#         coloredlogs.install(level='DEBUG', logger=self.logger, fmt=log_format)
#         # coloredlogs.install(level='CALC', logger=logger, fmt=log_format)
#         # config.log_file_name = config.file_name('.fits', '.log')
#         file_handler = logging.FileHandler(self.log_file_name)
#         file_handler.setLevel(logging.DEBUG)
#         file_handler.setFormatter(logging.Formatter(log_format))
#         self.logger.addHandler(file_handler)
#         self.logger.info("Initializing Logging!")
#         config.loger_file = True

class _logging_():
    logger = logging.getLogger(__name__)
    # Set the log level
    logger.setLevel(logging.DEBUG)
    # Fancy format
    log_format = "%(asctime)s - %(levelname)s - %(message)s"
    # Use colored logs to add color to the log messages
    coloredlogs.install(level='DEBUG', logger=logger, fmt=log_format)
    # coloredlogs.install(level='CALC', logger=logger, fmt=log_format)

    try:
        # print("# Removing previous log file.")
        os.system(f"rm -r {config.log_file_name}")
    except:
        pass
    file_handler = logging.FileHandler(config.log_file_name)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(logging.Formatter(log_format))
    logger.addHandler(file_handler)


    def __init__(self):
        logger.info("Initializing Logging!")

class read_data():
    """
    Read Input Data
    """
    def __init__(self, filename=None,residualname=None,psfname=None,
                 imagelist = [],residuallist=[],
                 is_background_subtracted=False,
                 sky_offset=None,
                 invarfilename=None, wtmapfilename=None,
                 varfilename=None, rmsfilename=None):
        """
        Parameters
        ----------
        invarfilename, wtmapfilename, varfilename, rmsfilename : str, optional
            Auxiliary noise products that ship alongside the science image.
            Surveys describe the same information in four conventions, so each
            one gets its own argument rather than a single file plus a mode
            flag -- the file name then records the provenance:

                invarfilename  inverse variance   sigma = 1 / sqrt(invvar)
                wtmapfilename  weight map         sigma = 1 / sqrt(wt)
                varfilename    variance           sigma = sqrt(var)
                rmsfilename    sigma directly     used as-is

            Which survey gives which:

                Legacy Survey  --invvar (IMAGETYP=INVVAR)   -> invarfilename
                Legacy Survey  --wtmap (Data Lab)           -> wtmapfilename
                HSC            --download-variance          -> varfilename
                Euclid         --rms (ESA/Q2 and ERO)       -> rmsfilename
                JWST           --rms (MAST ERR, DJA 1/sqrt(WHT))
                                                            -> rmsfilename

            A weight map and an inverse-variance map convert identically -- the
            Legacy Survey calls the same quantity by both names -- so the two
            arguments differ only in what they document.

            Whichever is given, `self.rms_data_2D` ends up holding a SIGMA map
            in image units, ready to pass on as `rms_map=` to `do_fit2D` /
            `sersic_multifit_general` (with the default
            `rms_convention='sigma'`) or to `compute_image_properties`. The
            arrays as loaded are kept too, under `invvar_data_2D`,
            `weight_data_2D` and `variance_data_2D`, since
            `compute_image_properties` and `compute_bkg_rms_maps` also take
            those conventions natively.

            If more than one is supplied the most direct wins --
            rms > invvar > variance > weight -- and nothing is combined: a
            survey's weight and variance maps are usually two views of the same
            numbers, so averaging them would understate the noise. The choice is
            reported on `self.rms_source`.

            Nothing is auto-wired: these are loaded and converted, never
            substituted into a fit behind your back, so no existing call
            changes. Pass `rms_map=input_data.rms_data_2D` explicitly.

            The scalar `rms_img` / `rms_res` (mad_std of the image / residual)
            are left alone and keep their current meaning.
        sky_offset : float, str or None, optional
            Add a constant back to the image as it is loaded, turning an
            already-sky-subtracted product back into an un-subtracted one.

            This exists because a cleaned image and a sky-subtracted image are
            two different things that usually arrive glued together. In
            morfometryka the sky is removed BEFORE galclean runs
            (`galclean(self.gal0 - self.skybg, ...)`); galclean itself only
            replaces contaminating sources with draws from the background and
            never subtracts anything. So the pedestal removal and the crowded
            -field cleaning are separable, and the removal is a single scalar:
            on J001431 `raw - galclean` is 0.00315254 on 95% of pixels with
            exactly zero scatter, identical in every annulus, and
            `galclean + skybg` reproduces the raw image to 6e-19.

            Adding it back gives you the cleaning without the sky subtraction --
            which matters when that subtraction over-estimated the sky, as it
            does whenever the galaxy fills its cutout.

                float   add this constant
                'mfmtk' read `skybg` from the sibling <rootname>.mfmtk

            Either way `is_background_subtracted` is forced to False, because
            the image now carries its pedestal again, and the value used is
            recorded on `self.sky_offset`.
        is_background_subtracted : bool, optional
            Declare that this image has ALREADY had its sky removed by an
            upstream tool (morfometryka `galclean`, a survey pipeline, anything).
            This is a property of the data, which is why it is declared here
            rather than at the fit, and it cannot be inferred from the file --
            provenance is not recoverable in general, so it has to be stated.

            Its effect: the fitter will not estimate and subtract a background
            from it. Doing so would be a double subtraction -- removing a
            pedestal that is already gone -- which pushes source-free regions
            negative and biases the outer profile.

            Defaults to False (raw image), so existing calls keep their current
            meaning.
        """

        self.filename = filename
        self.residualname = residualname
        self.psfname = psfname
        self.is_background_subtracted = is_background_subtracted
        self.sky_offset = sky_offset
        self.invarfilename = invarfilename
        self.wtmapfilename = wtmapfilename
        self.varfilename = varfilename
        self.rmsfilename = rmsfilename
        self.print_names()
        self.get_data()
        self.get_info()
        # try:
        #     self.get_info()
        # except:
        #     pass

    def print_names(self):
        if self.filename != None:
            print('++>> Image File:', os.path.basename(self.filename))
        if self.residualname != None:
            print('++>> Residual File:', os.path.basename(self.residualname))
        elif self.residualname == None:
            print('-->> No Residual File was provided.')
        if self.psfname != None:
            print('++>> PSF File:', os.path.basename(self.psfname))
        elif self.psfname == None:
            print('-->> No PSF File was provided.')
        for _label, _name in (('Inverse Variance', self.invarfilename),
                              ('Weight Map', self.wtmapfilename),
                              ('Variance', self.varfilename),
                              ('RMS', self.rmsfilename)):
            if _name is not None:
                print(f'++>> {_label} File:', os.path.basename(_name))

    def get_data(self):
        self.image_data_2D = None
        self.residual_data_2D = None
        self.psf_data_2D = None
        self.rms_img = None
        self.rms_res = None
        self.invvar_data_2D = None
        self.weight_data_2D = None
        self.variance_data_2D = None
        self.rms_data_2D = None
        self.rms_source = 'none'
        
        if self.filename != None:
            self.image_data_2D = mlibs.load_fits_data(self.filename)
            if self.sky_offset is not None:
                if isinstance(self.sky_offset, str):
                    # <rootname>.mfmtk: one '#'-prefixed header line of comma
                    # separated names, then one line of values.
                    _mf = self.filename.replace('_galclean.fits', '.fits')
                    _mf = _mf.replace('.fits', '.mfmtk')
                    with open(_mf) as _fh:
                        _keys, _vals = [l.strip().lstrip('#').split(',')
                                        for l in _fh][:2]
                    self.sky_offset = float(dict(zip(_keys, _vals))
                                            [self.sky_offset
                                             if self.sky_offset != 'mfmtk'
                                             else 'skybg'])
                    print(f'++>> sky_offset read from {os.path.basename(_mf)}: '
                          f'{self.sky_offset:.6g}')
                self.sky_offset = float(self.sky_offset)
                self.image_data_2D = self.image_data_2D + self.sky_offset
                # do_fit2D, the output writer and the plots all re-read the FITS
                # from `filename`, so an in-memory offset alone would never
                # reach the fit. Write the shifted image once and point the
                # object at it, which keeps every downstream consumer consistent
                # without threading an array through them.
                _out = self.filename.replace('.fits', '_skyadd.fits')
                mlibs.pf.writeto(_out, self.image_data_2D, overwrite=True)
                mlibs.copy_header(self.filename, _out, _out)
                self.original_filename = self.filename
                self.filename = _out
                # The pedestal is back, so the image is no longer subtracted.
                self.is_background_subtracted = False
                print(f'++>> Added a sky offset of {self.sky_offset:.6g} back to '
                      f'the image -> {os.path.basename(_out)}; '
                      f'is_background_subtracted set to False.')
            self.rms_img = mlibs.mad_std(self.image_data_2D,ignore_nan=True)
        if self.residualname != None:
            self.residual_data_2D = mlibs.load_fits_data(self.residualname)
            self.rms_res = mlibs.mad_std(self.residual_data_2D,ignore_nan=True)
        else:
            self.rms_res = self.rms_img # Assume same RMS if no residual provided
        if self.psfname != None:
            self.psf_data_2D = mlibs.load_fits_data(self.psfname)
            self.psf_fwhm = mlibs.psf_params(self.psf_data_2D)

        # Auxiliary noise products, each kept in the convention it arrived in.
        for _fname, _attr in ((self.invarfilename, 'invvar_data_2D'),
                              (self.wtmapfilename, 'weight_data_2D'),
                              (self.varfilename, 'variance_data_2D'),
                              (self.rmsfilename, 'rms_data_2D')):
            if _fname is None:
                continue
            _arr = mlibs.load_fits_data(_fname)
            # load_fits_data RETURNS ValueError (the class) instead of raising
            # when it cannot read a file, which would otherwise travel silently
            # all the way to the weights.
            if not isinstance(_arr, mlibs.np.ndarray):
                raise ValueError(f'Could not read {_attr} from '
                                 f'{os.path.basename(_fname)}.')
            if (self.image_data_2D is not None
                    and _arr.shape != mlibs.np.shape(self.image_data_2D)):
                raise ValueError(
                    f'{_attr} shape {_arr.shape} from '
                    f'{os.path.basename(_fname)} does not match the image '
                    f'shape {mlibs.np.shape(self.image_data_2D)}.')
            setattr(self, _attr, _arr)

        if (self.invvar_data_2D is not None or self.weight_data_2D is not None
                or self.variance_data_2D is not None
                or self.rms_data_2D is not None):
            # One shared implementation of the four conversions, precedence and
            # zero-weight handling; see `_sigma_from_noise_inputs`. Zero invvar
            # means 'no data here' and would give an infinite sigma, so those
            # pixels come back finite but huge, i.e. effectively unweighted.
            self.rms_data_2D, self.rms_source = mlibs._sigma_from_noise_inputs(
                rms_map=self.rms_data_2D,
                invvar_map=self.invvar_data_2D,
                variance_map=self.variance_data_2D,
                weight_map=self.weight_data_2D)
            _med = float(mlibs.np.nanmedian(self.rms_data_2D))
            print(f'++>> RMS map available as `.rms_data_2D` (from '
                  f'{self.rms_source}): median sigma = {_med:.6g}')
            if self.rms_img is not None and _med > 0:
                # Not a consistency check -- mad_std over the whole frame is
                # raised by the source itself, so it SHOULD sit above the
                # survey sigma. A ratio far from ~1 in the other direction is
                # the interesting case.
                print(f'     image mad_std = {self.rms_img:.6g} '
                      f'(ratio mad_std / sigma = {self.rms_img / _med:.3g})')
            print('     Pass it on explicitly, e.g. '
                  'rms_map=input_data.rms_data_2D, use_weights=True.')
    
    def get_info(self):
        self.cell_size = mlibs.get_cell_size(self.filename)
        try:
            self.beam_area_px = mlibs.beam_area2(self.filename)
        except:
            self.beam_area_px = 1
        self.get_pixel_scale()
        self.get_jy_conversion_factor()
        self.get_filter_info()
        self.report_units()
    
    def get_pixel_scale(self):
        cell_size_arcsec = mlibs.get_cell_size(self.filename)
        self.pixel_scale = cell_size_arcsec
        self.pixel_scale_x = cell_size_arcsec
        self.pixel_scale_y = cell_size_arcsec
        self.scale_source = "cell size"
        
        
    # def get_pixel_scale(self):
    #     """
    #     Determine the pixel scale from the FITS header in arcseconds.
    #     Order of preference:
    #     1. D00*SCAL (for drizzled HST images)
    #     2. CD matrix calculation (accounts for rotation)
    #     3. IDCSCALE (HST specific)
    #     4. CDELT values (for non-rotated images)
    #     """
    #     try:
    #         with mlibs.pf.open(self.filename) as hdul:
    #             # Try each method in order of preference
                
    #             # Method 1: Use drizzle scale if available (HST/WFC3 drizzled images)
    #             drizzle_scale_keys = [k for k in hdul[0].header.keys() if k.startswith('D') and k.endswith('SCAL')]
    #             if drizzle_scale_keys:
    #                 self.pixel_scale = hdul[0].header[drizzle_scale_keys[0]]
    #                 self.pixel_scale_x = self.pixel_scale_y = self.pixel_scale
    #                 self.scale_source = f"drizzle parameter ({drizzle_scale_keys[0]})"
    #                 return
                
    #             # Method 2: Check for CD matrix (accounts for rotation)
    #             for i, hdu in enumerate(hdul):
    #                 if 'CD1_1' in hdu.header and 'CD1_2' in hdu.header and 'CD2_1' in hdu.header and 'CD2_2' in hdu.header:
    #                     cd1_1 = hdu.header.get('CD1_1', 0)
    #                     cd1_2 = hdu.header.get('CD1_2', 0)
    #                     cd2_1 = hdu.header.get('CD2_1', 0)
    #                     cd2_2 = hdu.header.get('CD2_2', 0)
                        
    #                     # Calculate pixel scales using the proper formula for rotated images
    #                     self.pixel_scale_x = np.sqrt(cd1_1**2 + cd2_1**2) * 3600.0  # Convert from degrees to arcsec
    #                     self.pixel_scale_y = np.sqrt(cd1_2**2 + cd2_2**2) * 3600.0  # Convert from degrees to arcsec
                        
    #                     # Use average of x and y scales for overall scale
    #                     self.pixel_scale = (self.pixel_scale_x + self.pixel_scale_y) / 2.0
    #                     self.scale_source = "CD matrix calculation"
    #                     return
                
    #             # Method 3: Check for IDCSCALE (HST specific)
    #             for hdu in hdul:
    #                 if 'IDCSCALE' in hdu.header:
    #                     self.pixel_scale = hdu.header['IDCSCALE']
    #                     self.pixel_scale_x = self.pixel_scale_y = self.pixel_scale
    #                     self.scale_source = "IDCSCALE"
    #                     return
                
    #             # Method 4: Check for CDELT values (for non-rotated images)
    #             for hdu in hdul:
    #                 if 'CDELT1' in hdu.header and 'CDELT2' in hdu.header:
    #                     self.pixel_scale_x = abs(hdu.header.get('CDELT1', 0)) * 3600.0  # Convert from degrees to arcsec
    #                     self.pixel_scale_y = abs(hdu.header.get('CDELT2', 0)) * 3600.0  # Convert from degrees to arcsec
    #                     self.pixel_scale = (self.pixel_scale_x + self.pixel_scale_y) / 2.0
    #                     self.scale_source = "CDELT values"
    #                     return
                
    #             # If we get here, no scale information was found
    #             self.pixel_scale = None
    #             self.pixel_scale_x = self.pixel_scale_y = None
    #             self.scale_source = None
    #             print("-->> Warning: Could not determine pixel scale from header")
                    
    #     except Exception as e:
    #         print(f"-->> Error determining pixel scale: {e}")
    #         self.pixel_scale = None
    #         self.pixel_scale_x = self.pixel_scale_y = None
    #         self.scale_source = None
        
    # def get_pixel_scale(self):
    #     """
    #     Determine the pixel scale from the FITS header in arcseconds
    #     """
    #     try:
    #         with mlibs.pf.open(self.filename) as hdul:
    #             # Try different methods to get pixel scale
    #             # Method 1: Check for CD matrix (most common)
    #             if 'CD1_1' in hdul[0].header or ('CD1_1' in hdul[1].header if len(hdul) > 1 else False):
    #                 header = hdul[0].header if 'CD1_1' in hdul[0].header else hdul[1].header
    #                 cd1_1 = header.get('CD1_1', 0)
    #                 cd1_2 = header.get('CD1_2', 0)
    #                 cd2_1 = header.get('CD2_1', 0)
    #                 cd2_2 = header.get('CD2_2', 0)
    #                 self.pixel_scale_x = np.sqrt(cd1_1**2 + cd2_1**2) * 3600.0  # Convert from degrees to arcsec
    #                 self.pixel_scale_y = np.sqrt(cd1_2**2 + cd2_2**2) * 3600.0  # Convert from degrees to arcsec
    #                 self.pixel_scale = (self.pixel_scale_x + self.pixel_scale_y) / 2.0
    #                 print(f"Pixel scale determined from CD matrix: {self.pixel_scale:.4f} arcsec/pixel")
                
    #             # Method 2: Check for CDELT
    #             elif 'CDELT1' in hdul[0].header or ('CDELT1' in hdul[1].header if len(hdul) > 1 else False):
    #                 header = hdul[0].header if 'CDELT1' in hdul[0].header else hdul[1].header
    #                 self.pixel_scale_x = abs(header.get('CDELT1', 0)) * 3600.0  # Convert from degrees to arcsec
    #                 self.pixel_scale_y = abs(header.get('CDELT2', 0)) * 3600.0  # Convert from degrees to arcsec
    #                 self.pixel_scale = (self.pixel_scale_x + self.pixel_scale_y) / 2.0
    #                 print(f"Pixel scale determined from CDELT: {self.pixel_scale:.4f} arcsec/pixel")
                
    #             # Method 3: Check for other possible keywords (HST specific)
    #             elif 'IDCSCALE' in hdul[0].header or ('IDCSCALE' in hdul[1].header if len(hdul) > 1 else False):
    #                 header = hdul[0].header if 'IDCSCALE' in hdul[0].header else hdul[1].header
    #                 self.pixel_scale = header.get('IDCSCALE')
    #                 self.pixel_scale_x = self.pixel_scale_y = self.pixel_scale
    #                 print(f"Pixel scale determined from IDCSCALE: {self.pixel_scale:.4f} arcsec/pixel")
                
    #             # Method 4: Use drizzle scale if available (HST/WFC3 specific)
    #             elif 'D001SCAL' in hdul[0].header:
    #                 self.pixel_scale = hdul[0].header.get('D001SCAL') 
    #                 self.pixel_scale_x = self.pixel_scale_y = self.pixel_scale
    #                 print(f"Pixel scale determined from D001SCAL: {self.pixel_scale:.4f} arcsec/pixel")
                
    #             else:
    #                 # Default value if no scale information found
    #                 self.pixel_scale = None
    #                 self.pixel_scale_x = self.pixel_scale_y = None
    #                 print("-->> Warning: Could not determine pixel scale from header")
                
    #     except Exception as e:
    #         print(f"-->> Error determining pixel scale: {e}")
    #         self.pixel_scale = None
    #         self.pixel_scale_x = self.pixel_scale_y = None
            
    # def get_jy_conversion_factor(self):
    #     """
    #     Determine the conversion factor from raw pixel values to Janskys
    #     """
    #     # try:
    #     with mlibs.pf.open(self.filename) as hdul:
    #         # Check for instrument-specific keywords
    #         self.instrument = None
    #         self.bunit = None
            
    #         # First, determine the instrument and the units
    #         for i, hdu in enumerate(hdul):
    #             if 'INSTRUME' in hdu.header:
    #                 self.instrument = hdu.header['INSTRUME'].strip()
                
    #             # Check for BUNIT in this HDU
    #             if 'BUNIT' in hdu.header:
    #                 self.bunit = hdu.header['BUNIT'].strip()
    #                 # For HST data, science data is often in extension, so prioritize that
    #                 if i > 0 and 'SCI' in hdu.name:
    #                     break  # Found units in SCI extension, prioritize this
            
    #         # Now find appropriate conversion factor
    #         self.jy_conversion = None
            
    #         # For HST WFC3 data
    #         if self.instrument == 'WFC3' or self.instrument == 'WFC3  ':
    #             # Look for PHOTFNU (direct conversion to Jy)
    #             for hdu in hdul:
    #                 if 'PHOTFNU' in hdu.header:
    #                     self.jy_conversion = hdu.header['PHOTFNU']
    #                     break
                
    #             # If PHOTFNU not found, try PHOTFLAM with PHOTPLAM
    #             if self.jy_conversion is None:
    #                 for hdu in hdul:
    #                     if 'PHOTFLAM' in hdu.header and 'PHOTPLAM' in hdu.header:
    #                         photflam = hdu.header['PHOTFLAM']
    #                         photplam = hdu.header['PHOTPLAM']
    #                         # Convert from FLAM to Jy: Jy = FLAM * PHOTPLAM² / 3.34e4
    #                         self.jy_conversion = photflam * (photplam**2) / 3.34e4
    #                         break
            
    #         # For Legacy Survey data (nanomaggies)
    #         elif 'Legacy' in str(hdul[0].header.get('SURVEY', '')):
    #             # Standard conversion for nanomaggies
    #             self.jy_conversion = 3.631e-6
            
    #         # For other instruments, provide warning
    #         if self.jy_conversion is None:
    #             print("-->> Warning: Could not determine Jy conversion factor from header")
    #         else:
    #             # Create conversion methods for convenience
    #             self.to_jy = lambda data: data * self.jy_conversion
    #             self.to_mjy = lambda data: data * self.jy_conversion * 1000  # Convert to mJy
                    
    #     # except Exception as e:
    #     #     print(f"-->> Error determining Jy conversion: {e}")
    #     #     self.jy_conversion = None
    
    # def get_filter_info(self):
    #     """
    #     Get information about the filter used for the observation
    #     """
    #     # try:
    #     with mlibs.pf.open(self.filename) as hdul:
    #         self.filter = None
    #         self.filter_wavelength = None
    #         self.filter_frequency = None
            
    #         # Try to get filter name
    #         for hdu in hdul:
    #             if 'FILTER' in hdu.header:
    #                 self.filter = hdu.header['FILTER'].strip()
    #                 break
            
    #         # Try to get filter wavelength info
    #         for hdu in hdul:
    #             if 'PHOTPLAM' in hdu.header:
    #                 # PHOTPLAM is in Angstroms
    #                 self.filter_wavelength = hdu.header['PHOTPLAM'] * 1e-10  # Convert to meters
    #                 # Calculate frequency in Hz
    #                 self.filter_frequency = 2.99792458e8 / self.filter_wavelength
    #                 break
            
    #         # Additional info for HST WFC3
    #         if self.instrument == 'WFC3' or self.instrument == 'WFC3  ':
    #             for hdu in hdul:
    #                 if 'PHOTBW' in hdu.header:
    #                     self.filter_bandwidth = hdu.header['PHOTBW'] * 1e-10  # Convert to meters
    #                     break
                
    #     # except Exception as e:
    #     #     print(f"-->> Error getting filter info: {e}")
    #     #     self.filter = None
    #     #     self.filter_wavelength = None
    #     #     self.filter_frequency = None
    
    # def report_units(self):
    #     """
    #     Print a report of the units and conversion factors
    #     """
    #     print("\n===== Image Units and Conversion Information =====")
    #     print(f"Instrument: {self.instrument}")
    #     print(f"Filter: {self.filter}")
        
    #     if self.filter_wavelength is not None:
    #         print(f"Filter wavelength: {self.filter_wavelength*1e9:.2f} nm")
    #         print(f"Filter frequency: {self.filter_frequency/1e9:.2f} GHz")
        
    #     print(f"Original units: {self.bunit}")
        
    #     if self.jy_conversion is not None:
    #         print(f"Conversion factor to Jy: {self.jy_conversion:.6e}")        
    #     if self.pixel_scale is not None:
    #         print(f"Pixel scale: {self.pixel_scale:.4f} arcsec/pixel")
        
    #     print("=================================================\n")

    def get_jy_conversion_factor(self):
        """
        Determine the conversion factor from raw pixel values to Janskys.
        Supports radio (JY/BEAM), nanomaggies (Legacy Survey / DECam),
        JWST, HST, surface-brightness (MJy/sr) and AB zero-point calibrated
        (HSC/Subaru, Rubin/LSST) data.
        """
        with mlibs.pf.open(self.filename) as hdul:
            self.instrument = None
            self.bunit = None
            self.jy_conversion = None
            # Only set by Strategy 7 below; always defined so callers can test
            # it without hasattr guards.
            self.zero_point = None

            # Collect instrument and best BUNIT across all HDUs;
            # science extensions take priority for HST/JWST multi-extension files.
            for i, hdu in enumerate(hdul):
                if 'INSTRUME' in hdu.header:
                    self.instrument = hdu.header['INSTRUME'].strip()
                if 'BUNIT' in hdu.header:
                    candidate = hdu.header['BUNIT'].strip()
                    if i > 0 and 'SCI' in hdu.name:
                        self.bunit = candidate
                        break
                    if self.bunit is None:
                        self.bunit = candidate

            # Survey provenance, where the file states it. Legacy Survey
            # cutouts carry SURVEY/VERSION but no INSTRUME; HSC carries
            # HSCRERUN ('pdr3_wide'), which is the same kind of information.
            self.survey = None
            self.survey_version = None
            for hdu in hdul:
                if 'SURVEY' in hdu.header and self.survey is None:
                    self.survey = str(hdu.header['SURVEY']).strip()
                if 'VERSION' in hdu.header and self.survey_version is None:
                    self.survey_version = str(hdu.header['VERSION']).strip()
                if 'HSCRERUN' in hdu.header and self.survey_version is None:
                    self.survey_version = str(hdu.header['HSCRERUN']).strip()

            # Neither HSC cutouts nor Legacy Survey cutouts carry
            # INSTRUME/TELESCOP, so fall back on what does identify them.
            if self.instrument is None:
                for hdu in hdul:
                    if any(k in hdu.header
                           for k in ('HSCFILT', 'HSCRERUN', 'HSCPROD')):
                        self.instrument = 'HSC'
                        if self.survey is None:
                            self.survey = 'HSC-SSP'
                        break
            if self.instrument is None and self.survey is not None:
                # 'LegacySurvey' names a programme rather than one camera --
                # its images come from DECam, 90Prime (BASS) and MOSAIC-3
                # (MzLS) depending on band and hemisphere, and the cutout
                # header does not say which. Reporting the survey is the
                # honest answer; `self.survey`/`self.survey_version` carry it
                # verbatim for anyone who needs to be precise.
                self.instrument = self.survey
            if self.instrument is None:
                # Radio images (WSClean/CASA) name the array in TELESCOP and
                # leave INSTRUME out entirely.
                for hdu in hdul:
                    if 'TELESCOP' in hdu.header:
                        _tel = str(hdu.header['TELESCOP']).strip()
                        if _tel:
                            self.instrument = _tel
                            break

            # --- Strategy 1: BUNIT-driven conversion ---
            # Handles the most common cases purely from the header unit keyword.
            if self.bunit is not None:
                bl = self.bunit.lower()

                # Radio: pixel values already in Jy/beam or Jy
                if bl in ('jy/beam', 'jy/bm', 'jy beam-1', 'jy'):
                    self.jy_conversion = 1.0

                # Optical survey: nanomaggies (AB zeropoint 22.5 mag → 3.631 µJy)
                elif bl in ('nanomaggy', 'nanomaggies', 'nmgy', 'nmgys',
                            'nmgy/pix', 'nanomaggy/pix'):
                    self.jy_conversion = 3.631e-6

                # microJansky
                elif bl in ('ujy', 'ujy/beam', 'microjy', 'microjansky',
                            'µjy', 'ujy beam-1'):
                    self.jy_conversion = 1e-6

                # milliJansky
                elif bl in ('mjy/beam', 'millijy', 'millijansky', 'mjy beam-1'):
                    self.jy_conversion = 1e-3

                # nanoJansky (JWST native; may carry a numeric prefix, e.g. "10.0*nanoJansky")
                elif 'nanojansky' in bl or bl in ('njy', 'njy/beam', 'nanojy'):
                    if '*' in bl:
                        try:
                            multiplier = float(bl.split('*')[0].strip())
                            self.jy_conversion = multiplier * 1e-9
                        except ValueError:
                            self.jy_conversion = 1e-9
                    else:
                        self.jy_conversion = 1e-9

                # MJy/sr (Spitzer, WISE, Herschel, some JWST pipeline products)
                # Surface brightness, so it needs the solid angle of a pixel to
                # become a flux per pixel.
                elif 'mjy/sr' in bl or bl in ('mjy sr-1', 'mjansky/sr'):
                    pixel_area_sr = None
                    # The JWST pipeline states the pixel area outright. Prefer
                    # it: deriving the area from `pixel_scale` squares an
                    # average of the two axis scales, which is not an area
                    # unless the pixels are square, and it inherits whichever
                    # extension `get_cell_size` happened to pick.
                    for hdu in hdul:
                        if 'PIXAR_SR' in hdu.header:
                            try:
                                _pa = float(hdu.header['PIXAR_SR'])
                            except (ValueError, TypeError):
                                _pa = None
                            if _pa and _pa > 0:
                                pixel_area_sr = _pa
                            break
                    if (pixel_area_sr is None and getattr(self, 'pixel_scale', None)
                            is not None):
                        # pixel_scale is in arcsec; convert to radians
                        pixel_scale_rad = self.pixel_scale * (3.14159265358979 / (180.0 * 3600.0))
                        pixel_area_sr = pixel_scale_rad ** 2
                    if pixel_area_sr is not None:
                        self.jy_conversion = 1e6 * pixel_area_sr  # MJy/pix → Jy/pix

            # --- Strategy 2: JWST instruments (PHOTFNU is authoritative) ---
            if self.jy_conversion is None and self.instrument in (
                    'NIRCAM', 'NIRISS', 'NIRSPEC', 'MIRI'):
                for hdu in hdul:
                    if 'PHOTFNU' in hdu.header:
                        self.jy_conversion = hdu.header['PHOTFNU']
                        break
                if self.jy_conversion is None:
                    print("Warning: JWST data detected but PHOTFNU keyword not found")

            # --- Strategy 3: HST instruments (PHOTFNU > PHOTFLAM+PHOTPLAM) ---
            if self.jy_conversion is None and self.instrument in (
                    'WFC3', 'WFC3  ', 'ACS', 'NICMOS', 'WFPC2', 'WFC3/IR', 'WFC3/UVIS'):
                for hdu in hdul:
                    if 'PHOTFNU' in hdu.header:
                        self.jy_conversion = hdu.header['PHOTFNU']
                        break
                if self.jy_conversion is None:
                    for hdu in hdul:
                        if 'PHOTFLAM' in hdu.header and 'PHOTPLAM' in hdu.header:
                            photflam = hdu.header['PHOTFLAM']
                            photplam = hdu.header['PHOTPLAM']
                            # F_lambda -> Jy; see _FLAM_TO_FNU_JY. ACS and WFPC2
                            # never write PHOTFNU, so this is their only route.
                            self.jy_conversion = (photflam * (photplam ** 2)
                                                  * _FLAM_TO_FNU_JY)
                            break

            # --- Strategy 4: Legacy Survey cutout service (SURVEY keyword) ---
            if self.jy_conversion is None:
                survey = str(hdul[0].header.get('SURVEY', '')).lower()
                if 'legacy' in survey:
                    self.jy_conversion = 3.631e-6
                    if self.bunit is None:
                        self.bunit = 'nanomaggy'

            # --- Strategy 5: nanomaggies inferred from MAGZERO = 22.5 ---
            if self.jy_conversion is None:
                for hdu in hdul:
                    if 'MAGZERO' in hdu.header:
                        try:
                            if abs(float(hdu.header['MAGZERO']) - 22.5) < 0.1:
                                self.jy_conversion = 3.631e-6
                                if self.bunit is None:
                                    self.bunit = 'nanomaggy'
                        except (ValueError, TypeError):
                            pass
                        break

            # --- Strategy 6: Generic PHOTFNU / PHOTFLAM fallback ---
            # Strategy 3 is gated on INSTRUME, so a header that carries the
            # photometric keywords but names no instrument never reaches it --
            # e.g. HLA colour composites, which are single-HDU products with
            # PHOTFLAM/PHOTPLAM/BUNIT and nothing else. The keywords mean the
            # same thing whoever wrote them, so honour them here rather than
            # giving up.
            if self.jy_conversion is None:
                for hdu in hdul:
                    if 'PHOTFNU' in hdu.header:
                        self.jy_conversion = hdu.header['PHOTFNU']
                        break
            if self.jy_conversion is None:
                for hdu in hdul:
                    if 'PHOTFLAM' in hdu.header and 'PHOTPLAM' in hdu.header:
                        self.jy_conversion = (hdu.header['PHOTFLAM']
                                              * (hdu.header['PHOTPLAM'] ** 2)
                                              * _FLAM_TO_FNU_JY)
                        break

            # --- Strategy 7: AB zero-point keywords (HSC/Subaru, Rubin/LSST) ---
            # Last resort by design: every strategy above is a more direct
            # statement of the units and must not be pre-empted by this one.
            #
            # A zero-point ZP means m_AB = -2.5*log10(counts) + ZP, and
            # m_AB = 0 is 3631 Jy, so one count is 3631 * 10**(-0.4*ZP) Jy.
            #
            # The LSST Science Pipelines (and therefore HSC) do not write ZP
            # directly; they write FLUXMAG0, the count level of a zero-magnitude
            # source, so ZP = 2.5*log10(FLUXMAG0). HSC PDR coadds are calibrated
            # to FLUXMAG0 = 6.3095734448e10, i.e. ZP = 27.0 exactly and
            # 5.7547e-08 Jy (57.5 nJy) per count.
            if self.jy_conversion is None:
                zero_point = None
                for hdu in hdul:
                    if 'FLUXMAG0' in hdu.header:
                        try:
                            fluxmag0 = float(hdu.header['FLUXMAG0'])
                            if fluxmag0 > 0:
                                zero_point = 2.5 * math.log10(fluxmag0)
                                break
                        except (ValueError, TypeError):
                            pass
                if zero_point is None:
                    for hdu in hdul:
                        for key in ('MAGZP', 'PHOTZP', 'ZP', 'ZPTMAG', 'MAGZPT'):
                            if key in hdu.header:
                                try:
                                    zero_point = float(hdu.header[key])
                                except (ValueError, TypeError):
                                    zero_point = None
                                break
                        if zero_point is not None:
                            break
                if zero_point is not None:
                    self.zero_point = zero_point
                    self.jy_conversion = 3631.0 * 10 ** (-0.4 * zero_point)
                    if self.bunit is None:
                        self.bunit = (f'counts (AB zero-point '
                                      f'{zero_point:.4g} mag)')

            # --- Final: build convenience conversion lambdas ---
            if self.jy_conversion is None:
                print("Warning: Could not determine Jy conversion factor from header")
                print(f"Instrument: {self.instrument}, BUNIT: {self.bunit}")
            else:
                # Any Jy-per-pixel-unit scale implies an AB zero-point, since
                # m_AB = 0 is 3631 Jy: ZP = -2.5*log10(jy_conversion / 3631).
                # Strategy 7 reads it straight from the header, and this
                # inverts exactly for those files; for everything else it fills
                # the value in -- 22.5 for nanomaggies, 8.9 for an array
                # already in Jy -- so `to_abmag` works for any calibrated image
                # rather than only the zero-point-calibrated ones.
                if self.zero_point is None:
                    self.zero_point = -2.5 * math.log10(self.jy_conversion
                                                        / 3631.0)
                self.to_jy = lambda data: data * self.jy_conversion
                self.to_mjy = lambda data: data * self.jy_conversion * 1000
                self.to_ujy = lambda data: data * self.jy_conversion * 1e6
                self.to_njy = lambda data: data * self.jy_conversion * 1e9
                self.to_abmag = lambda data: (-2.5 * mlibs.np.log10(data)
                                              + self.zero_point)

        # self.flux_units/self.flux_conversion_factor: the (flux_units,
        # flux_conversion_factor) vocabulary that `eimshow`'s colorbar/plot_rms
        # formatting and `compute_image_properties`' flux reporting understand,
        # derived from self.jy_conversion above. Stored as plain attributes (not
        # just returned from a method) so notebook code can do
        # `eimshow(..., flux_units=input_data.flux_units,
        # flux_conversion_factor=input_data.flux_conversion_factor)` directly,
        # same as any other read_data attribute (cell_size, rms_img, ...).
        # Neither consumer rescales the array itself from these -- they only
        # drive labels and the reported flux numbers -- so flux_units must
        # describe the units the raw pixel array is already in, not the units
        # you want displayed.
        #
        # flux_conversion_factor has ONE meaning everywhere: mJy per native
        # pixel unit, i.e. jy_conversion * 1000. That is exactly the quantity
        # `compute_image_properties` multiplies its raw sums by, so the pair can
        # be handed to it unmodified. `eimshow` reads the factor only in its
        # 'mJy/px' branch -- 'mJy' and 'nanomaggies' carry their own hardcoded
        # factors -- so keeping it populated in every branch costs those two
        # nothing and stops the value from being a trap when it is passed on.
        if self.jy_conversion is None:
            # Unknown units: don't claim a photometric scale we can't back up.
            # 1.0 here is "leave the numbers alone", not a photometric claim.
            self.flux_units = 'any'
            self.flux_conversion_factor = 1.0
        else:
            self.flux_conversion_factor = self.jy_conversion * 1000
            if math.isclose(self.jy_conversion, 1.0, rel_tol=1e-6):
                # Raw array already in Jy(/beam) -- e.g. typical radio images.
                self.flux_units = 'mJy'
            elif math.isclose(self.jy_conversion, 3.631e-6, rel_tol=1e-6):
                # Legacy Survey / DECam nanomaggies -- eimshow has a dedicated,
                # hardcoded factor for this case.
                self.flux_units = 'nanomaggies'
            else:
                # Everything else (HST PHOTFNU/PHOTFLAM, JWST PHOTFNU, MJy/sr,
                # HSC/LSST zero-point counts, ...): the raw array is in some
                # native per-pixel unit with no name of its own, and
                # jy_conversion carries the Jy-per-native-unit scale.
                self.flux_units = 'mJy/px'

    def get_flux_kwargs(self):
        """Convenience dict form of self.flux_units/self.flux_conversion_factor,
        for `**`-unpacking into anything that takes the pair -- `eimshow`,
        `compute_image_properties`, `measures`, `structural_morphology`:

            mlibs.compute_image_properties(..., **input_data.get_flux_kwargs())

        The two must travel together: `flux_units` alone leaves an
        instrument whose native unit has no name ('mJy/px') unscaled, and the
        factor alone is ignored by every consumer outside that branch.
        """
        return dict(flux_units=self.flux_units,
                    flux_conversion_factor=self.flux_conversion_factor)

    # Kept for existing notebooks/CLI; the pair is no longer eimshow-specific.
    get_eimshow_flux_kwargs = get_flux_kwargs

    # Nominal effective wavelengths of the broad bands, in Angstroms (SVO
    # Filter Profile Service). Neither HSC nor Legacy Survey cutouts carry
    # PHOTPLAM, so these tables are the only way to give
    # `filter_wavelength`/`filter_frequency` a value. They are band-nominal,
    # good to roughly a percent, and no substitute for a real filter curve if
    # you are doing precision photometry.
    #
    # HSC: 'i2'/'r2' are the re-coated i and r filters; their effective
    # wavelengths differ from the originals by well under the width of the
    # band, so they map onto the same nominal value here.
    HSC_FILTER_WAVELENGTHS = {'g': 4816.0, 'r': 6234.0, 'r2': 6234.0,
                              'i': 7741.0, 'i2': 7741.0, 'z': 8912.0,
                              'y': 9762.0}

    # Legacy Survey: DECam values (DECaLS/DELVE/DR10 south). The northern
    # 90Prime g,r and MOSAIC-3 z differ by a few percent, and the header does
    # not say which camera took the image, so one table covers both.
    LS_FILTER_WAVELENGTHS = {'g': 4798.0, 'r': 6412.0, 'i': 7822.0,
                             'z': 9158.0, 'y': 9887.0}

    def get_filter_info(self):
        """
        Get information about the filter used for the observation.
        Enhanced to handle JWST filter naming conventions.
        """
        with mlibs.pf.open(self.filename) as hdul:
            self.filter = None
            self.filter_wavelength = None
            self.filter_frequency = None
            self.filter_bandwidth = None
            self.pupil = None

            # Try to get filter name
            for hdu in hdul:
                if 'FILTER' in hdu.header:
                    self.filter = hdu.header['FILTER'].strip()
                    break

            # ACS (and WFPC2) put the filter in a pair of wheel keywords rather
            # than in FILTER, with the unused wheel parked on CLEAR1L/CLEAR2L.
            # Raw pipeline products (*_drc/_drz straight from the archive) and
            # HLA products have no FILTER keyword at all, so without this they
            # report no filter.
            if self.filter is None:
                _cands = []
                for _key in ('FILTER1', 'FILTER2', 'FILTNAM1', 'FILTNAM2'):
                    for hdu in hdul:
                        if _key in hdu.header:
                            _v = str(hdu.header[_key]).strip()
                            if _v and not _v.upper().startswith('CLEAR'):
                                _cands.append(_v)
                            break
                if any(c.lower() == 'detection' for c in _cands):
                    # HLA detection stacks set FILTER1='detection' and leave a
                    # real filter in FILTER2; the stack is not that filter.
                    self.filter = 'detection'
                elif len(_cands) == 1:
                    self.filter = _cands[0]
                elif len(_cands) > 1:
                    # Crossed/ramp configurations really are two filters.
                    self.filter = '+'.join(dict.fromkeys(_cands))

            # HSC: prefer HSCFILT ('HSC-Z') over FILTER and normalise to the
            # plain band letter. FILTER holds 'i2'/'r2' for the re-coated
            # filters, which is a fact about the glass rather than the band, and
            # it makes band-keyed bookkeeping across a multi-band set awkward.
            if self.instrument == 'HSC':
                for hdu in hdul:
                    if 'HSCFILT' in hdu.header:
                        _hf = str(hdu.header['HSCFILT']).strip()
                        self.filter = _hf.split('-')[-1].lower() or _hf
                        break
                else:
                    if self.filter is not None:
                        self.filter = self.filter.lower()

            # Legacy Survey cutouts have no FILTER; the band is in BANDS
            # ('g'), with BAND0 naming the first plane of a multi-band stack.
            if self.filter is None:
                for hdu in hdul:
                    for _key in ('BANDS', 'BAND0'):
                        if _key in hdu.header:
                            _b = str(hdu.header[_key]).strip()
                            # A multi-band stack lists every plane ('grz');
                            # that is not one filter, so leave it verbatim
                            # rather than pretending it is a single band.
                            self.filter = _b.lower() if _b else None
                            break
                    if self.filter is not None:
                        break

            # For JWST, also check for PUPIL element
            if self.instrument in ['NIRCAM', 'NIRISS', 'NIRSPEC', 'MIRI']:
                for hdu in hdul:
                    if 'PUPIL' in hdu.header:
                        self.pupil = hdu.header['PUPIL'].strip()
                        break
            
            # Try to get filter wavelength information
            for hdu in hdul:
                if 'PHOTPLAM' in hdu.header:
                    # PHOTPLAM is in Angstroms
                    self.filter_wavelength = hdu.header['PHOTPLAM'] * 1e-10  # Convert to meters
                    # Calculate frequency in Hz (c = λν)
                    self.filter_frequency = 2.99792458e8 / self.filter_wavelength
                    break
            
            # HSC, Legacy Survey and JWST have no PHOTPLAM, so fall back on the
            # nominal band tables. Kept after the PHOTPLAM loop and guarded on
            # it, so a header that does state its pivot wavelength always wins.
            if self.filter_wavelength is None and self.filter is not None:
                _table = None
                if self.instrument == 'HSC':
                    _table = self.HSC_FILTER_WAVELENGTHS
                elif self.survey is not None and 'legacy' in self.survey.lower():
                    _table = self.LS_FILTER_WAVELENGTHS
                else:
                    # Filter-name keyed (F150W, F160W, ...) -- shared with
                    # radio_utils.get_frequencies so both agree by construction.
                    _table = getattr(mlibs, 'FILTER_PIVOT_WAVELENGTHS', None)
                if _table is not None and self.filter in _table:
                    self.filter_wavelength = _table[self.filter] * 1e-10
                    self.filter_frequency = 2.99792458e8 / self.filter_wavelength

            # Radio images state their frequency outright, on the third WCS
            # axis. Same keywords and unit handling as `get_frequencies` in
            # radio_utils.py, so a radio image reports observing frequency
            # where an optical one reports a band.
            if self.filter_frequency is None:
                for hdu in hdul:
                    if ('CRVAL3' in hdu.header
                            and 'FREQ' in str(hdu.header.get('CTYPE3', '')).upper()):
                        try:
                            _freq = float(hdu.header['CRVAL3'])
                        except (ValueError, TypeError):
                            break
                        _unit = str(hdu.header.get('CUNIT3', 'Hz')).upper()
                        _scale = {'KHZ': 1e3, 'MHZ': 1e6, 'GHZ': 1e9}
                        for _u, _s in _scale.items():
                            if _u in _unit:
                                _freq *= _s
                                break
                        if _freq > 0:
                            self.filter_frequency = _freq
                            self.filter_wavelength = 2.99792458e8 / _freq
                        break

            # Try to get filter bandwidth
            for hdu in hdul:
                if 'PHOTBW' in hdu.header:
                    # PHOTBW is in Angstroms
                    self.filter_bandwidth = hdu.header['PHOTBW'] * 1e-10  # Convert to meters
                    break


    def report_units(self):
        """
        Print a comprehensive report of the units and conversion factors.
        Enhanced to display JWST-specific information.
        """
        print("\n" + "="*60)
        print("Image Units and Conversion Information")
        print("="*60)
        
        # Instrument information
        print(f"Instrument: {self.instrument if self.instrument else 'Unknown'}")
        if getattr(self, 'survey', None):
            _ver = getattr(self, 'survey_version', None)
            print(f"Survey: {self.survey}" + (f" ({_ver})" if _ver else ""))
        print(f"Filter: {self.filter if self.filter else 'Unknown'}")
        
        # JWST-specific: show pupil element if present
        if hasattr(self, 'pupil') and self.pupil:
            print(f"Pupil: {self.pupil}")
        
        # Wavelength and frequency information
        if self.filter_wavelength is not None:
            # Radio wavelengths are centimetres, optical ones are Angstroms --
            # one fixed unit makes the other unreadable.
            if self.filter_wavelength >= 1e-3:
                print(f"Wavelength: {self.filter_wavelength*100:.4g} cm")
                print(f"Frequency: {self.filter_frequency/1e9:.4f} GHz")
            else:
                print(f"Filter wavelength: {self.filter_wavelength*1e9:.2f} nm ({self.filter_wavelength*1e10:.1f} Angstroms)")
                print(f"Filter frequency: {self.filter_frequency/1e9:.2f} GHz")
        
        if hasattr(self, 'filter_bandwidth') and self.filter_bandwidth is not None:
            print(f"Filter bandwidth: {self.filter_bandwidth*1e9:.2f} nm ({self.filter_bandwidth*1e10:.1f} Angstroms)")
        
        # Units and conversion factors
        print(f"\nOriginal units (BUNIT): {self.bunit if self.bunit else 'Not specified'}")

        if getattr(self, 'zero_point', None) is not None:
            print(f"AB zero-point: {self.zero_point:.4f} mag")

        if self.jy_conversion is not None:
            print(f"Conversion factor to Jy: {self.jy_conversion:.6e}")
            # print(f"Conversion factor to mJy: {self.jy_conversion*1000:.6e}")
            # print(f"Conversion factor to microJy: {self.jy_conversion*1e6:.6e}")
            
            # Show example conversion
            # print(f"\nExample: pixel value = 1.0")
            # print(f"  -> {self.jy_conversion:.6e} Jy")
            # print(f"  -> {self.jy_conversion*1000:.6e} mJy")
            # print(f"  -> {self.jy_conversion*1e6:.6e} microJy")
            # print(f"  -> {self.jy_conversion*1e9:.6e} nJy")
        else:
            print("Conversion factor to Jy: Not available")
        
        # Pixel scale information
        if hasattr(self, 'pixel_scale') and self.pixel_scale is not None:
            print(f"\nPixel scale: {self.pixel_scale:.4f} arcsec/pixel")
            # Calculate pixel area for surface brightness conversions
            pixel_area_arcsec2 = self.pixel_scale**2
            print(f"Pixel area: {pixel_area_arcsec2:.6f} arcsec^2")
        
        print("="*60 + "\n")
        
        


class radio_image_analysis():
    def __init__(self, input_data,z = None,do_petro=False,
                 # logger=None,
                 crop=False,box_size=256,
                 apply_mask=True,mask=None,dilation_size = None,
                 sigma_level=3, sigma_mask=6,vmin_factor=3,last_level=3,
                 results=None,mask_component=None,
                 npixels=128,kernel_size=21,fwhm=81,
                 SAVE=True, show_figure=True):
        self.input_data = input_data
        # self.logger = logger
        self.crop = crop
        self.box_size = box_size
        self.do_petro = do_petro
        self.apply_mask = apply_mask
        self.dilation_size = dilation_size
        self.sigma_level = sigma_level
        self.sigma_mask = sigma_mask
        self.mask = mask
        self.mask_component = mask_component
        self.vmin_factor = vmin_factor
        self.last_level = last_level
        self.npixels = npixels
        self.kernel_size = kernel_size
        self.fwhm = fwhm
        self.z = z
        self.results = results
        self.SAVE = SAVE
        self.show_figure = show_figure

        self.image_properties()


    def image_properties(self):
        try:
            self.cell_size = mlibs.get_cell_size(self.input_data.filename)
        except:
            # print('!! WARNING !! Setting cellsize/pixelsize to unity.')
            _logging_.logger.warning("Setting cellsize/pixelsize to unity.")
            self.cell_size = 1

        _logging_.logger.info("Computing image level statistics.")
        if self.z is None:
            _logging_.logger.warning("The redshift of the source was not specified."
                                     "Conversions to physical units will not be "
                                     "performed.")
        # _logging_.file_handler.info("Computing image level statistics.")

        # self.image_level_statistics = \
        #     mlibs.level_statistics(img=self.input_data.filename,
        #                            cell_size=self.cell_size, crop=self.crop,
        #                            sigma = self.sigma_level,
        #                            apply_mask=self.apply_mask,
        #                            results=self.results, SAVE=self.SAVE,
        #                            ext=config.ext,
        #                            show_figure=config.show_plots)
        #
        _logging_.logger.info("Computing image properties.")
        self.levels, self.fluxes, self.Lgrow, self.Lgrow_err, self.Lgrow_norm, self.Lgrow_err_norm, self.radii, \
            self.agrow, self.omask, self.mask, self.results_im_props = \
            mlibs.compute_image_properties(img=self.input_data.filename,
                                           cell_size=self.cell_size,
                                           residual=self.input_data.residualname,
                                           sigma_mask=self.sigma_mask,
                                           dilation_size=self.dilation_size,
                                           crop=self.crop,
                                           iterations=config.mask_iterations,
                                           box_size=self.box_size,
                                           last_level=self.last_level,
                                           mask=self.mask,
                                           apply_mask=self.apply_mask,
                                           vmin_factor=self.vmin_factor,
                                           results=self.results,
                                           show_figure=self.show_figure,
                                           logger=_logging_.logger)
        
        # self.img_stats = \
        #     mlibs.get_image_statistics(imagename=self.input_data.filename,
        #                                residual_name=self.input_data.residualname,
        #                                cell_size=self.cell_size,
        #                                mask_component=None,
        #                                mask=self.mask,
        #                                region='', dic_data=None,
        #                                sigma_mask=self.sigma_mask,
        #                                apply_mask=self.apply_mask,
        #                                fracX=0.15, fracY=0.15)

        # self.image_measures, self.mask, self.omask = \
        #     mlibs.measures(imagename=self.input_data.filename,
        #                    residualname=self.input_data.residualname,
        #                    z=self.z,
        #                    mask_component=self.mask_component,
        #                    sigma_mask=self.sigma_mask,
        #                    last_level=self.last_level,
        #                    vmin_factor=self.vmin_factor,
        #                    plot_catalog=True, 
        #                    data_2D=self.input_data.image_data_2D,
        #                    npixels=self.npixels,fwhm=self.fwhm,
        #                    kernel_size=self.kernel_size,
        #                    dilation_size=self.dilation_size,
        #                    main_feature_index=0,
        #                    results_final={},
        #                    crop=self.crop, box_size=self.box_size,
        #                    iterations=config.mask_iterations,
        #                    fracX=0.15, fracY=0.15,
        #                    deblend=False, bkg_sub=False,
        #                    bkg_to_sub=None, rms=None,
        #                    do_petro=self.do_petro,
        #                    apply_mask=self.apply_mask,
        #                    do_PLOT=True, SAVE=self.SAVE,
        #                    show_figure=self.show_figure,
        #                    mask=self.mask,
        #                    do_measurements='all',
        #                    compute_A=True,
        #                    add_save_name='',logger=_logging_.logger)


class source_extraction():
    """
    Source extraction class, responsible to find relevant regions of emission.

    For now, only SEP and PF+Photutils are implemented. Soon, other algorithms will be added
    in this class as alternatives for source extraction. These are:
        - PyBDSF
        - AstroDendro
        - Photutils (alone)
    """
    def __init__(self, input_data, z=0.05, ids_to_add=None,
                ids_types=None,
                default_component_type=None,
                crop=False, box_size=256,
                apply_mask=False, mask=None, dilation_size=None,
                sigma_level=6, sigma_mask=6, vmin_factor=3, mask_component=None,
                bwf=1, bhf=1, fwf=1, fhf=1,
                segmentation_map=True, filter_type='conv',
                deblend_nthresh=25, deblend_cont=1e-3,
                clean_param=0.5, clean=True,
                minarea_factor=1.0, npixels=None,
                sort_by='distance', first_ID_only=False,
                sigma=6,
                mask_grow_iterations=2,
                ell_size_factor=None,
                obs_type='radio', algorithm='SEP', threshold_mode='sigma',
                # DEPRECATED, see `force_circular` in the fit drivers: declare
                # circularity on the component instead, with `ids_types` /
                # `ids_to_add`. These two only speak about whole groups.
                force_circular=False, force_circular_all=False,
                show_detection=False, show_petro_plots=False,
                show_bkg_map=False,
                SAVE=True, show_figure=True, dry_run=False, SE_ref=None,
                # New robust-specific parameters
                multiscale_levels=3,
                extended_threshold_factor=0.5,
                watershed_connectivity=2,
                merge_threshold=0.5,
                adaptive_background=True,
                preserve_extended=True,
                min_separation=None,
                # New field/mosaic-specific parameters (algorithm='field')
                bkg_box_size=None,
                bkg_filter_size=(3, 3),
                bkg_estimator='median',
                rms_map=None,
                save_products=False,
                products_path=None,
                products=('data_sub', 'bkg', 'rms'),
                nproc=1,
                edge_margin=0,
                reject_edge_sources=False,
                snr_min=None,
                cutout_size_factor=4.0,
                name_prefix='J',
                return_masks=False,
                max_masks=200,
                max_labels_to_plot=500,
                verbose=1):
        """
        Parameters
        ----------
        input_data : str
            Path to the image.
        z : float
            Redshift of the source.
        ids_to_add : dict
            Extra model components to fit on top of what was detected, keyed by
            parent ID (the 1-indexed detected region)::

                ids_to_add = {'1': ['point-like', 'disk'], '2': ['sersic']}

            Purely additive: the detected component for each region stays as
            detected, and each named type appends one further component cloned
            from that region. A region may be given several. Component types come
            from `mlibs.COMPONENT_TYPE_PRESETS` -- 'point-like', 'gaussian',
            'compact', 'disk', 'sersic'.

            The old flat-list form (`['1','2','3']`) is no longer accepted; it only
            ever said *which parent to clone*, never what to make, and every clone
            silently came out a disk. Passing one now raises a `ValueError` showing
            the equivalent dict.
        ids_types : dict, optional
            Re-type a *detected* component without adding anything, e.g.
            `{'1': 'point-like'}` when detected region 1 is unresolved. Untyped
            detections use `default_component_type`.
        default_component_type : str, optional
            Preset for any component not explicitly typed. `None` (the default)
            defers to the observation kind: 'gaussian' (n=0.5) for radio,
            'sersic' (n free) for anything else.
        crop : bool
            If True, crop the image.
        box_size : int
            Size of the box to be cropped.
        apply_mask : bool
            If True, apply the mask to the image.
        mask : array
            Mask to be applied to the image.
        dilation_size : int
            Size of the binary dilation kernel.
        sigma_level : float
            Sigma level for detection.
        sigma_mask : float
            Sigma level for the mask.
        vmin_factor : float
            Factor to be multiplied by the standard deviation to set the minimum
            value for the imshow plot.
        mask_component : array
            Mask to be applied to the image.
        bwf : int
            Box width fraction in terms of the beam size
            for the background estimation.
        bhf : int
            Box height fraction in terms of the beam size
            for the background estimation.
        fwf : int
            Filter width fraction in terms of the beam size
            for the background estimation.
        fhf : int
            Filter height fraction in terms of the beam size
            for the background estimation.
        segmentation_map : bool
            If True, returns the segmentation map.
        filter_type : str
            Type of filter to be used.
        deblend_nthresh : int
            Number of thresholds for deblending.
        deblend_cont : float
            Minimum contrast ratio for deblending.
        clean_param : float
            Cleaning parameter.
        clean : bool
            If True, clean the image.
        minarea_factor : float
            Factor to be multiplied by the minimum area for detection.
            Default is 1.0, i.e. one restoring beam size. Any structure smaller
            than one beam size will not be detected. This is critical if you have
            oversampled data.
        sort_by : str
            Sort the output by flux or area.
        sigma : float
            Sigma level for detection.
        ell_size_factor : int
            Size factor of the ellipse to be drawn in the detected structures.
        show_detection : bool
            If True, show the detection plot.
        show_petro_plots : bool
            If True, show the petrosian plots.
        SAVE : bool
            If True, save plots.
        show_figure : bool
            If True, show the figure.
        dry_run : bool
            If True, do not compute source properties. In a first run, use True
            to inspect how well the source detection was.

        Notes on ``algorithm='field'``
        ------------------------------
        Field/mosaic mode is a different beast from the postage-stamp
        algorithms ('SEP', 'PF', 'astphot', 'robust'). It detects *many* sources
        across a wide image against a 2D RMS map and produces a catalogue
        (``self.catalogue``) rather than per-component photometry; it never
        calls ``prepare_fit``, which is postage-stamp machinery. The intended
        flow is:

            SE = mp.source_extraction(input_data, algorithm='field', sigma=5.0)
            SE.to_csv(); SE.to_regions()          # inspect / QA
            SE.filter_catalogue(snr_min=6.0)
            out = SE.make_cutouts(cutout_size='auto')
            # then analyse out['imagelist'] with the usual per-source tools.

        bkg_box_size : int
            Background2D mesh size in pixels. Default ~10 beams.
        bkg_filter_size : tuple
            Median filter size for the background mesh.
        bkg_estimator : str
            'median', 'sextractor' or 'mmm'.
        rms_map : str or array
            Externally supplied noise map, bypassing RMS estimation.
        save_products : bool
            Write the full-mosaic intermediates (data-bkg, bkg, rms, and
            optionally the segmentation images) to FITS. Also available
            after the fact as ``SE.save_products()``.
        products_path : str
            Where to write them; defaults to the image's own directory.
        products : tuple of str
            Which to write: 'data_sub', 'bkg', 'rms', 'segm', 'segm_parent'.
        nproc : int
            Processes used for deblending a crowded field.
        edge_margin : int
            Border width, in pixels, used to flag (or reject) edge sources.
        reject_edge_sources : bool
            Remove labels touching `edge_margin` instead of only flagging them.
        snr_min : float
            Drop catalogue rows below this peak SNR.
        cutout_size_factor : float
            Multiplier on the parent island size used for `cutout_size='auto'`.
        name_prefix : str
            Prefix for the IAU-style source designations.
        return_masks : bool
            Materialise per-component boolean masks. Off by default: for a field
            with thousands of detections these are thousands of full-size arrays.
        max_masks : int
            Cap on the number of masks built when `return_masks=True`.
        max_labels_to_plot : int
            Cap on labelled sources in the field detection plot.
        """

        self.input_data = input_data
        self.z = z
        self.ids_to_add = ids_to_add
        self.ids_types = ids_types
        # None means "let the fit driver decide" -- radio wants 'gaussian',
        # optical wants 'sersic'. Resolved in `contruct_source_properties`.
        self.default_component_type = default_component_type
        self.crop = crop
        self.box_size = box_size
        self.apply_mask = apply_mask
        self.mask = mask
        self.dilation_size = dilation_size
        self.sigma_level = sigma_level
        self.sigma_mask = sigma_mask
        self.sigma = sigma
        self.minarea_factor = minarea_factor
        self.npixels = npixels
        self.ell_size_factor = ell_size_factor
        self.vmin_factor = vmin_factor
        self.mask_component = mask_component
        self.algorithm = algorithm
        self.SE_ref = SE_ref
        # self.bw = bw
        # self.bh = bh
        # self.fw = fw
        # self.fh = fh
        # if (bw == None) & (bw == None) & (bw == None) & (bw == None):
        try:
            self.bspx, self.aO, self.bO = \
                mlibs.get_beam_size_px(self.input_data.filename)
            print(self.bspx)
            self.bw = self.aO / bwf
            self.bh = self.bO / bhf
            self.fw = self.aO / fwf
            self.fh = self.bO / fhf
        except:
            self.bw = int((input_data.image_data_2D.shape[0]*0.2)/bwf)
            self.bh = int((input_data.image_data_2D.shape[1]*0.2)/bhf)
            self.fw = int((input_data.image_data_2D.shape[0]*0.1)/fwf)
            self.fh = int((input_data.image_data_2D.shape[1]*0.1)/fhf)
        try:
            self.minarea = mlibs.beam_area2(self.input_data.filename)
        except:
            self.minarea = self.input_data.image_data_2D.shape[0]/30
        # self.bw, self.bh, self.fw, self.fh = bw, bh, fw, fh
        self.segmentation_map = segmentation_map
        self.filter_type = filter_type
        self.deblend_nthresh = deblend_nthresh
        self.deblend_cont = deblend_cont
        self.clean_param = clean_param
        self.clean = clean
        self.sort_by = sort_by
        self.first_ID_only = first_ID_only
        self.SAVE = SAVE
        self.show_figure = show_figure
        self.show_detection = show_detection
        self.show_bkg_map = show_bkg_map
        self.show_petro_plots = show_petro_plots
        self.obs_type = obs_type
        self.force_circular = force_circular
        self.force_circular_all = force_circular_all
        self.mask_grow_iterations = mask_grow_iterations
        self.threshold_mode = threshold_mode
        self.multiscale_levels = multiscale_levels
        self.extended_threshold_factor = extended_threshold_factor
        self.watershed_connectivity = watershed_connectivity
        self.merge_threshold = merge_threshold
        self.adaptive_background = adaptive_background
        self.preserve_extended = preserve_extended
        self.min_separation = min_separation
        # field/mosaic mode
        self.bkg_box_size = bkg_box_size
        self.bkg_filter_size = bkg_filter_size
        self.bkg_estimator = bkg_estimator
        self.rms_map = rms_map
        self.save_products_flag = save_products
        self.products_path = products_path
        self.products = products
        self.nproc = nproc
        self.edge_margin = edge_margin
        self.reject_edge_sources = reject_edge_sources
        self.snr_min = snr_min
        self.cutout_size_factor = cutout_size_factor
        self.name_prefix = name_prefix
        self.return_masks = return_masks
        self.max_masks = max_masks
        self.max_labels_to_plot = max_labels_to_plot
        self.verbose = verbose
        self.catalogue = None
        self.cutouts = None

        if self.algorithm == 'field':
            """
            Field/mosaic mode always stops at the catalogue: `prepare_fit` is
            postage-stamp machinery and must not run on a wide-field image.
            `dry_run` keeps its usual meaning of "show me the detection", which
            here is simply the default behaviour.
            """
            if dry_run is True:
                self.show_detection = True
            self.get_sources()
        else:
            if dry_run is True:
                self.show_detection = True
                self.get_sources()

            if dry_run is not True:
                self.contruct_source_properties()

    def get_sources(self):
        if self.algorithm == 'SEP':
            """
            It uses the SEP library to perform source extraction.
            """
            self.masks, self.indices, self.bkg, self.seg_maps, self.objects = \
                mlibs.sep_source_ext(self.input_data.filename,
                                     residualname = self.input_data.residualname,
                               bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                               # filtering options for source detection
                               minarea=self.minarea,
                               minarea_factor=self.minarea_factor,
                               segmentation_map=self.segmentation_map,
                               filter_type=self.filter_type, 
                               mask=self.mask,
                               deblend_nthresh=self.deblend_nthresh,
                               deblend_cont=self.deblend_cont,
                               clean_param=self.clean_param,
                               clean=self.clean,
                               sort_by=self.sort_by,
                               npixels = self.npixels,
                               dilation_size=self.dilation_size,
                               iterations = self.mask_grow_iterations,
                               sigma=self.sigma,sigma_mask=self.sigma_mask,
                               ell_size_factor=self.ell_size_factor,
                               apply_mask=self.apply_mask,
                               show_detection=self.show_detection,
                               show_bkg_map=self.show_bkg_map)
        if self.algorithm == 'PF':
            """
            It uses PetroFit routines which call Photutils functions
            to perform source extraction.
            """
            (self.masks, self.indices, self.bkg,
             self.seg_maps, self.objects, self.cat) = \
                mlibs.phot_source_ext(self.input_data.filename,
                                      residual=self.input_data.residual_data_2D,
                               bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                               # filtering options for source detection
                               segmentation_map=self.segmentation_map,
                               psf_data = self.input_data.psf_data_2D,
                               filter_type=self.filter_type, mask=self.mask,
                               deblend_nthresh=self.deblend_nthresh,
                               deblend_cont=self.deblend_cont,
                               clean_param=self.clean_param,
                               clean=self.clean,
                               sort_by=self.sort_by, first_ID_only=self.first_ID_only,
                               threshold_mode=self.threshold_mode,
                               sigma=self.sigma,sigma_mask=self.sigma_mask,
                               dilation_size=self.dilation_size,
                               iterations = self.mask_grow_iterations,
                               minarea=self.minarea,
                               npixels = self.npixels,
                               minarea_factor = self.minarea_factor,
                               ell_size_factor=self.ell_size_factor,
                               apply_mask=self.apply_mask,
                               show_detection=self.show_detection,
                               show_bkg_map=self.show_bkg_map,
                               SE_ref=self.SE_ref)
                
        if self.algorithm == 'astphot':
            """
            It uses PetroFit routines which call Photutils functions
            to perform source extraction.
            """
            (self.masks, self.indices, self.bkg,
             self.seg_maps, self.objects, self.cat) = \
                mlibs.astphot_source_ext(self.input_data.filename,
                                      residual=self.input_data.residual_data_2D,
                               bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                               # filtering options for source detection
                               segmentation_map=self.segmentation_map,
                               psf_data = self.input_data.psf_data_2D,
                               filter_type=self.filter_type, mask=self.mask,
                               deblend_nthresh=self.deblend_nthresh,
                               deblend_cont=self.deblend_cont,
                               clean_param=self.clean_param,
                               clean=self.clean,
                               sort_by=self.sort_by, first_ID_only=self.first_ID_only,
                               threshold_mode=self.threshold_mode,
                               sigma=self.sigma,sigma_mask=self.sigma_mask,
                               dilation_size=self.dilation_size,
                               iterations = self.mask_grow_iterations,
                               minarea=self.minarea,
                               npixels = self.npixels,
                               minarea_factor = self.minarea_factor,
                               ell_size_factor=self.ell_size_factor,
                               apply_mask=self.apply_mask,
                               show_detection=self.show_detection,
                               show_bkg_map=self.show_bkg_map,
                               SE_ref=self.SE_ref)
        if self.algorithm == 'robust':
            (self.masks, self.indices, self.bkg,
            self.seg_maps, self.objects) = \
                mlibs.robust_source_ext(self.input_data.filename,
                                    residualname=self.input_data.residualname,
                                    bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                                    minarea=self.minarea,
                                    minarea_factor=self.minarea_factor,
                                    segmentation_map=self.segmentation_map,
                                    filter_type=self.filter_type,
                                    mask=self.mask,
                                    deblend_nthresh=self.deblend_nthresh,
                                    deblend_cont=self.deblend_cont,
                                    clean_param=self.clean_param,
                                    clean=self.clean,
                                    sort_by=self.sort_by,
                                    npixels=self.npixels,
                                    dilation_size=self.dilation_size,
                                    iterations=self.mask_grow_iterations,
                                    sigma=self.sigma,
                                    sigma_mask=self.sigma_mask,
                                    ell_size_factor=self.ell_size_factor,
                                    apply_mask=self.apply_mask,
                                    show_detection=self.show_detection,
                                    show_bkg_map=self.show_bkg_map,
                                    # robust-specific parameters
                                    multiscale_levels=self.multiscale_levels,
                                    extended_threshold_factor=self.extended_threshold_factor,
                                    watershed_connectivity=self.watershed_connectivity,
                                    merge_threshold=self.merge_threshold,
                                    adaptive_background=self.adaptive_background,
                                    preserve_extended=self.preserve_extended,
                                    min_separation=self.min_separation)
        if self.algorithm == 'robust_opt':
            (self.masks, self.indices, self.bkg,
            self.seg_maps, self.objects) = \
                mlibs.robust_optical_source_ext(self.input_data.filename,
                                    residualname=self.input_data.residualname,
                                    bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                                    minarea=self.minarea,
                                    minarea_factor=self.minarea_factor,
                                    segmentation_map=self.segmentation_map,
                                    filter_type=self.filter_type,
                                    mask=self.mask,
                                    deblend_nthresh=self.deblend_nthresh,
                                    deblend_cont=self.deblend_cont,
                                    clean_param=self.clean_param,
                                    clean=self.clean,
                                    sort_by=self.sort_by,
                                    npixels=self.npixels,
                                    dilation_size=self.dilation_size,
                                    iterations=self.mask_grow_iterations,
                                    sigma=self.sigma,
                                    sigma_mask=self.sigma_mask,
                                    ell_size_factor=self.ell_size_factor,
                                    apply_mask=self.apply_mask,
                                    show_detection=self.show_detection,
                                    show_bkg_map=self.show_bkg_map,
                                    # robust-specific parameters
                                    multiscale_levels=self.multiscale_levels,
                                    extended_threshold_factor=self.extended_threshold_factor,
                                    watershed_connectivity=self.watershed_connectivity,
                                    merge_threshold=self.merge_threshold,
                                    adaptive_background=self.adaptive_background,
                                    preserve_extended=self.preserve_extended,
                                    # min_separation=self.min_separation
                                    )

        if self.algorithm == 'field':
            """
            Field/mosaic source extraction: detection against a 2D RMS map over
            a wide image, producing a source catalogue rather than per-component
            photometry. See morphen/field_extraction.py.
            """
            (self.masks, self.indices, self.bkg, self.bkg_rms,
             self.seg_maps, self.seg_maps_parent, self.objects,
             self.cat, self.cat_parent, self.catalogue, self.wcs) = \
                mlibs.field_source_ext(self.input_data.filename,
                                       residualname=self.input_data.residualname,
                                       residual=self.input_data.residual_data_2D,
                                       sigma=self.sigma,
                                       minarea=self.minarea,
                                       minarea_factor=self.minarea_factor,
                                       npixels=self.npixels,
                                       deblend_nthresh=self.deblend_nthresh,
                                       deblend_cont=self.deblend_cont,
                                       nproc=self.nproc,
                                       filter_type=self.filter_type,
                                       psf_data=self.input_data.psf_data_2D,
                                       bkg_box_size=self.bkg_box_size,
                                       bkg_filter_size=self.bkg_filter_size,
                                       bkg_estimator=self.bkg_estimator,
                                       rms_map=self.rms_map,
                                       mask=self.mask,
                                       apply_mask=self.apply_mask,
                                       sigma_mask=self.sigma_mask,
                                       dilation_size=self.dilation_size,
                                       iterations=self.mask_grow_iterations,
                                       edge_margin=self.edge_margin,
                                       reject_edge_sources=self.reject_edge_sources,
                                       snr_min=self.snr_min,
                                       sort_by=self.sort_by,
                                       cutout_size_factor=self.cutout_size_factor,
                                       name_prefix=self.name_prefix,
                                       return_masks=self.return_masks,
                                       max_masks=self.max_masks,
                                       show_detection=self.show_detection,
                                       show_bkg_map=self.show_bkg_map,
                                       max_labels_to_plot=self.max_labels_to_plot,
                                       save_plot=self.SAVE,
                                       save_products=self.save_products_flag,
                                       products_path=self.products_path,
                                       products=self.products,
                                       verbose=self.verbose)

    @property
    def data_sub(self):
        """
        Background-subtracted field image, computed on demand.

        Deliberately not cached: on a large mosaic this is another full-size
        float array, and the object already holds the data, the background, the
        RMS map and two segmentation images.
        """
        if self.algorithm != 'field' or getattr(self, 'bkg', None) is None:
            return None
        return mlibs.np.squeeze(self.input_data.image_data_2D) - self.bkg

    def _check_field_mode(self, what):
        """Guard for the methods that only make sense in field/mosaic mode."""
        if self.algorithm != 'field':
            raise RuntimeError(
                f"`{what}` is only available with algorithm='field'. "
                f"This object was built with algorithm='{self.algorithm}', "
                f"which does not produce a field catalogue.")
        if self.catalogue is None or len(self.catalogue) == 0:
            raise RuntimeError('The field catalogue is empty; nothing to do. '
                               'Try lowering `sigma` or `minarea_factor`.')

    def _default_output_name(self, extension):
        base = self.input_data.filename.replace('.fits', '')
        return f'{base}_field_catalogue.{extension}'

    def to_csv(self, filename=None):
        """Save the field catalogue as CSV."""
        self._check_field_mode('to_csv')
        filename = filename or self._default_output_name('csv')
        self.catalogue.to_csv(filename, index=False)
        print(f' ++==>> Catalogue ({len(self.catalogue)} rows) saved to {filename}')
        return filename

    def to_fits(self, filename=None):
        """Save the field catalogue as a FITS binary table."""
        self._check_field_mode('to_fits')
        from astropy.table import Table
        filename = filename or self._default_output_name('fits')
        Table.from_pandas(self.catalogue).write(filename, overwrite=True)
        print(f' ++==>> Catalogue ({len(self.catalogue)} rows) saved to {filename}')
        return filename

    def to_regions(self, filename=None, **kwargs):
        """Write a ds9/CARTA region file for quick visual QA of the detections."""
        self._check_field_mode('to_regions')
        filename = filename or self._default_output_name('reg')
        return mlibs.catalogue_to_regions(self.catalogue, filename, **kwargs)

    def filter_catalogue(self, inplace=True, **kwargs):
        """
        Filter the field catalogue (SNR, flux, area, edges, primary-only).

        See `mlibs.filter_field_catalogue`. By default the SE object's own
        catalogue is replaced, so a subsequent `make_cutouts` acts on the
        filtered set.
        """
        self._check_field_mode('filter_catalogue')
        kwargs.setdefault('image_shape', self.input_data.image_data_2D.shape[-2:])
        filtered = mlibs.filter_field_catalogue(self.catalogue, **kwargs)
        if inplace:
            self.catalogue = filtered
        return filtered

    def show_field(self, **kwargs):
        """
        Field overview in the standard photutils segmentation style: the
        background-subtracted data with Kron apertures, alongside the deblended
        segmentation image. Both panels are subset to whatever is currently in
        `self.catalogue`, so it reflects any filtering already applied.
        """
        self._check_field_mode('show_field')
        kwargs.setdefault('max_labels_to_plot', self.max_labels_to_plot)
        return mlibs.plot_field_detections(self.input_data.filename,
                                           self.catalogue,
                                           segm=self.seg_maps,
                                           cat=self.cat,
                                           data_sub=self.data_sub,
                                           bkg_rms=self.bkg_rms, **kwargs)

    def save_products(self, products=('data_sub', 'bkg', 'rms'),
                      save_path=None, prefix=None, **kwargs):
        """
        Write the full-mosaic intermediates to FITS: background-subtracted
        data, background, RMS map, and optionally the segmentation images
        ('segm', 'segm_parent').

        The RMS map is worth keeping - feeding it back as `rms_map=` on a later
        run skips re-estimating the background over the whole mosaic.
        """
        self._check_field_mode('save_products')
        return mlibs.save_field_products(self.input_data.filename,
                                         bkg_image=self.bkg,
                                         bkg_rms=self.bkg_rms,
                                         data_sub=self.data_sub,
                                         seg_maps=self.seg_maps,
                                         seg_maps_parent=self.seg_maps_parent,
                                         products=products,
                                         save_path=save_path, prefix=prefix,
                                         **kwargs)

    def make_cutouts(self, imagelist=None, residuallist=None, **kwargs):
        """
        Build one cutout per catalogued source, for one or many images.

        Defaults to cutting the image this SE object was built from. The
        returned `imagelist`/`residuallist` are ready to be handed straight to
        the existing per-source tools (`mp.source_extraction(algorithm='PF')`,
        `measures`, `structural_morphology`, ...) - a cutout is exactly the
        postage stamp they expect.

        See `mlibs.cutouts_from_catalogue` for the full parameter list.
        """
        self._check_field_mode('make_cutouts')
        if imagelist is None:
            # Only the object's own image/residual pair may be defaulted. If the
            # caller supplies their own imagelist (several bands, say), the
            # residuals must come with it -- silently pairing every band with
            # this object's single residual would cut the wrong noise image.
            imagelist = [self.input_data.filename]
            if residuallist is None and self.input_data.residualname is not None:
                residuallist = [self.input_data.residualname]
        kwargs.setdefault('cutout_size_factor', self.cutout_size_factor)
        kwargs.setdefault('verbose', self.verbose)

        self.cutouts = mlibs.cutouts_from_catalogue(self.catalogue, imagelist,
                                                    residuallist=residuallist,
                                                    **kwargs)
        self.imagelist = self.cutouts['imagelist']
        self.residuallist = self.cutouts['residuallist']
        return self.cutouts

    def make_multiband_cutouts(self, bands=None, MFS_images=None,
                               MFS_residuals=None, imagelist=None,
                               residuallist=None, freqlist_MFS=None,
                               freqlist=None, band_names=None, **kwargs):
        """
        Cut every catalogued source out of every band, MFS and sub-band alike.

        Either pass a `bands` spec from `mlibs.build_band_spec`, or hand over
        the six lists `mlibs.prepare_data` returns and they will be grouped
        here. Output is one directory per source:

            <workdir>/<iau_name>/<iau_name>-<band>-<token>-{image,residual}.fits

        Astrometric alignment is **off unless you ask for it**: pass
        ``correct_shift=True`` (default False, as on `cutout_2D_radec`), and
        optionally ``align='global'`` to shift every band and sub-band onto the
        reference band's MFS cutout - only valid on a shared pixel grid, and
        the cell sizes are checked - or ``align='band'`` to align each band's
        sub-bands to their own MFS. ``align`` alone does nothing while
        ``correct_shift`` is False.

        See `mlibs.multiband_cutouts_from_catalogue` for the full parameter
        list; everything there is accepted here as a keyword. Per-source
        results land in ``out['sources'][name]`` (``imagelist`` /
        ``residuallist`` for the sub-bands, ``MFS_imagelist`` /
        ``MFS_residuallist`` for the MFS).
        """
        self._check_field_mode('make_multiband_cutouts')
        if bands is None:
            if MFS_images is None:
                raise ValueError('Pass either `bands` (from build_band_spec) '
                                 'or at least `MFS_images`.')
            bands = mlibs.build_band_spec(MFS_images,
                                          MFS_residuals=MFS_residuals,
                                          imagelist=imagelist,
                                          residuallist=residuallist,
                                          freqlist_MFS=freqlist_MFS,
                                          freqlist=freqlist,
                                          band_names=band_names,
                                          verbose=self.verbose)
        kwargs.setdefault('cutout_size_factor', self.cutout_size_factor)
        kwargs.setdefault('verbose', self.verbose)

        self.multiband_cutouts = mlibs.multiband_cutouts_from_catalogue(
            self.catalogue, bands, **kwargs)
        self.manifest = self.multiband_cutouts['manifest']
        return self.multiband_cutouts

    def contruct_source_properties(self):
        # Radio data is dominated by beam-scale blobs, optical by real profiles,
        # so the sensible baseline component differs. Kept as a parameter with a
        # per-observation default rather than hardcoded in `prepare_fit`.
        default_component_type = self.default_component_type
        if default_component_type is None:
            default_component_type = ('gaussian' if self.obs_type == 'radio'
                                      else 'sersic')

        (self.sources_photometries, self.n_components, self.n_IDs,
         self.masks, self.indices, self.objects,
         self.psf_name, self.mask, self.bkg) = \
            mlibs.prepare_fit(self.input_data.filename,
                              self.input_data.residualname,
                              z=self.z,ids_to_add = self.ids_to_add,
                              ids_types = self.ids_types,
                              default_component_type=default_component_type,
                              bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                              sigma=self.sigma, sigma_mask=self.sigma_mask,
                              psf_data = self.input_data.psf_data_2D,
                              apply_mask=self.apply_mask,mask=self.mask,
                              mask_grow_iterations = self.mask_grow_iterations,
                              deblend_nthresh=self.deblend_nthresh,
                              ell_size_factor=self.ell_size_factor,
                              minarea=self.minarea,npixels = self.npixels,
                              minarea_factor=self.minarea_factor,
                              deblend_cont=self.deblend_cont,
                              clean_param=self.clean_param,
                              obs_type=self.obs_type,algorithm=self.algorithm,
                              SE = self.SE_ref,
                              sort_by=self.sort_by,
                              multiscale_levels=self.multiscale_levels, 
                              extended_threshold_factor=self.extended_threshold_factor, 
                              watershed_connectivity=self.watershed_connectivity, 
                              merge_threshold=self.merge_threshold, 
                              adaptive_background=self.adaptive_background, 
                              preserve_extended=self.preserve_extended, 
                              min_separation=self.min_separation,
                              force_circular = self.force_circular,
                              force_circular_all = self.force_circular_all,
                              show_petro_plots=self.show_petro_plots)

    def construct_source_properties(self):
        """
        Correctly-spelled alias of `contruct_source_properties`.

        The original name is missing an 's' and is used throughout the existing
        notebooks, so it stays as the implementation.
        """
        return self.contruct_source_properties()

    def component_types(self):
        """
        The resolved component type of every model component, in fit order.

        Returns
        -------
        list of str
            Length `n_components`; entry `i` is the type of component `i+1`.
            Detected components come first, then those added by `ids_to_add`.
        """
        return list(self.sources_photometries['component_types'])


class evaluate_source_structure():
    """
    This will be designed to evaluate the souce structure
    in order to check its complexity and compute how many model
    components will be required to perform the multi-sersic fitting.
    
    Also, this will compute basic source morphology, in order to 
    quantify which component represents a compact or a extended 
    structure. 
    """
    pass
    

class sersic_multifit_radio():
    """
    Multi-Sersic Fitting Decomposition.

    Perform a semi-automated and robust multi-sersic image decomposition.
    It supports GPU-acceleration using Jax. If no GPU is present, Jax still 
    will benefit from CPU parallel processing. Do not worry, you do not have 
    to change anything, Jax will automatically detect wheter you are runnin on
    CPU or GPU. 

    Basic principles:
        - run a source extraction, to identify relevant emission.
        - compute basic properties for each identified region, such as
          size, intensity, shape and orientation
        - uses that information to construct an object and prepare the
          settings to start the fit
        - compute statistics of the fit
        - if asked, calculates the relative fluxes of each component.
    To improve:
        - run an MCMC on model parameters (optional, takes time)
    To-do:
        - automated evaluation of which components is compact (unresolved)
        and which components is extended
        - automated evaluation of which structure cannot be modelled by a
        single function. For example, a spiral galaxy is reconized as a single
        source, but it can not be modelled by a single function: we require to
        model the bulge/bar/disk, for example.
    """
    def __init__(self, input_data, SE, aspect=None,
                 which_residual='shuffled',
                 fix_geometry = None,
                 force_circular = None,
                 comp_ids = [],
                 region_grow_to_ref = True,
                 default_component_type = 'gaussian',
                 fix_n = None,
                 fix_value_n = None,
                 fix_max_value_Rn = None,
                 fix_min_value_Rn = None,
                 fix_max_value_n = None,
                 fix_min_value_n = None,
                 dr_fix = None,fix_x0_y0=None,
                 trunc =None,
                 sigma=6.0, use_mask_for_fit=False,mask_for_fit=None,
                 mask=None,
                 tr_solver = "exact",loss = 'cauchy',
                 convolution_mode='GPU',method1='least_squares',
                 self_bkg = False, bkg_map = None,
                 rms_map = None, use_weights = False,
                 is_bkg_map_conv = False,
                 method2='least_squares',
                 parameters_mini_init = None,
                 z = 0.01,
                 verbose=0):
        """
        Parameters
        ----------
        input_data : object
            Input data object. See read_data class.
        SE : object
            Source extraction object.
        aspect : float
            Aspect ratio of the image.
        which_residual : str
            Which residual to use for the fitting.
        fix_geometry : bool or dict or list
            If True, fix the geometry (the boxiness parameter `cg`) of the
            components. See the note on per-component arguments below.
        comp_ids : list of str
            Which *final* component indices are the compact ones, as strings,
            e.g. `['1','3']`. Note this indexes the full component list --
            detected components first, then those added by `ids_to_add` -- so it
            has to be revisited whenever the component set changes.
        region_grow_to_ref : bool, optional
            How the per-region `region_data` / `region_diffuse` apertures in
            `decomp_table` are grown out of `SE.masks`. Source extraction detects
            at a higher SNR threshold than the reference mask is drawn at, so the
            deblended cores are much tighter than that mask; the default True
            grows each region until it stops expanding, so the regions partition
            the reference aperture and their fluxes sum to the whole-source
            `total` row. Pass False for the older fixed one-step growth. Check
            `region_area_completeness` on the region rows to see how much of the
            reference aperture was actually covered.
        default_component_type : str
            Preset used for components the source-extraction step did not type.
            Radio data defaults to 'gaussian' (Sersic with n locked to 0.5), which
            suits beam-scale blobs. Use 'sersic' to leave n free.
        force_circular : bool or float or dict or list, optional
            Per-component circularity. `True` pins `ell` to 0 and freezes `PA`
            and `cg` with it; a float caps `ell` at that value instead. Normally
            declared on the component itself, which is clearer because it needs
            no final index::

                ids_types  = {'1': ('point-like', {'force_circular': True})}
                ids_to_add = {'1': [('sersic', {'force_circular': 0.05})]}

            This argument is the by-final-index form of the same thing, and wins
            over what the component declared. Both win over
            `source_extraction`'s deprecated force_circular/force_circular_all.
        fix_n, fix_value_n, fix_max_value_Rn, fix_min_value_Rn, dr_fix,
        fix_x0_y0, fix_geometry, trunc : dict or list or scalar, optional
            Per-component fit controls. **Leave these unset unless you need to
            override a preset**: each component's type already supplies them, and
            the full dict is built automatically at the right length. When you do
            pass something, a dict keyed by component number is the clear form::

                fix_max_value_Rn = {1: 1.0}   # component 1 unresolved

            A positional list (indexed from component 1) and a bare scalar
            (applied to every component) are both still accepted. Overriding a
            value the component's preset *locks* -- n for 'gaussian'/'disk', the
            1 px radius for 'point-like' -- warns and is ignored; switch that
            component to the 'sersic' type to take control of it.
        sigma : float
            Sigma level for detection.
        tr_solver : str
            Solver for the trust region problem.
        convolution_mode : str
            Convolution mode.
        method1 : str
            Method for the first pass of the fit.
        method2 : str
            Method for the second pass of the fit.
        z : float
            Redshift of the source.
        """
        self.input_data = input_data
        self.SE = SE
        self.use_mask_for_fit = use_mask_for_fit
        if self.use_mask_for_fit == True:
            if mask_for_fit is None:
                _logging_.logger.info(f" ++>> Using a mask for fitting was requested, "
                                      f"but no mask was provided. Using the mask from the source "
                                      f"extraction object (SE.mask).")
                self.mask_for_fit = self.SE.mask
            else:
                _logging_.logger.info(f" ++>> Using the provided  mask for fitting.")
                self.mask_for_fit = mask_for_fit
        else:
            _logging_.logger.info(f" ++>> Fitting without a mask.")
            if mask_for_fit is None:
                self.mask_for_fit = None
            else:
                self.mask_for_fit = mask_for_fit

        if mask is None:
            self.mask = self.SE.mask
        else:
            self.mask = mask
            
        self.aspect = aspect
        # Component IDs are integers everywhere downstream (`comp_ID` columns,
        # `ext_ids`, the `model_cN` keys). Notebooks pass them as strings, so
        # normalise once here rather than at each consumer.
        self.comp_ids = [int(c) for c in (comp_ids or [])]
        # How the per-region apertures are grown out of `SE.masks`; see
        # `run_image_fitting`'s "Per-region apertures" note.
        self.region_grow_to_ref = region_grow_to_ref
        self.convolution_mode = convolution_mode
        self.method1 = method1
        self.method2 = method2
        self.tr_solver = tr_solver
        self.loss = loss
        self.parameters_mini_init = parameters_mini_init
        self.z = z
        self.which_residual = which_residual
        self.sigma = sigma
        self.self_bkg = self_bkg
        self.bkg_map = bkg_map
        self.rms_map = rms_map
        self.use_weights = use_weights
        self.is_bkg_map_conv = is_bkg_map_conv
        self.verbose = verbose

        # Every per-component fit control is materialised from each component's
        # preset, at exactly `n_components` entries, with anything the caller
        # supplied layered on top. This replaces the hand-written parallel lists
        # that had to be sized against a component count only known after source
        # extraction -- and that were routinely the wrong length, which went
        # unnoticed because the fit only ever indexed `[j]`.
        self.default_component_type = default_component_type
        _controls = mlibs.build_fit_control_maps(
            self.SE.sources_photometries,
            self.SE.n_components,
            overrides={'fix_n': fix_n,
                       'fix_value_n': fix_value_n,
                       'fix_max_value_n': fix_max_value_n,
                       'fix_min_value_n': fix_min_value_n,
                       'fix_max_value_Rn': fix_max_value_Rn,
                       'fix_min_value_Rn': fix_min_value_Rn,
                       'dr_fix': dr_fix,
                       'fix_x0_y0': fix_x0_y0,
                       'fix_geometry': fix_geometry,
                       'trunc': trunc,
                       'force_circular': force_circular},
            default_component_type=self.default_component_type)

        self.fix_n = _controls['fix_n']
        self.fix_value_n = _controls['fix_value_n']
        self.fix_max_value_n = _controls['fix_max_value_n']
        self.fix_min_value_n = _controls['fix_min_value_n']
        self.fix_max_value_Rn = _controls['fix_max_value_Rn']
        self.fix_min_value_Rn = _controls['fix_min_value_Rn']
        self.dr_fix = _controls['dr_fix']
        self.fix_x0_y0 = _controls['fix_x0_y0']
        self.fix_geometry = _controls['fix_geometry']
        self.trunc = _controls['trunc']
        # Absent from `_controls` unless a component declared it or the caller
        # passed one, which is exactly the "nobody said anything" case that lets
        # `source_extraction`'s deprecated force_circular flags still apply.
        self.force_circular = _controls.get('force_circular')

        self.__sersic_radio()

    def __sersic_radio(self):
        (self.results_fit, self.result_mini, self.mini, self.lmfit_results,
         self.lmfit_results_1st_pass,
         self.errors_fit, self.models, self.data_properties,
         self.results_compact_conv_morpho,
         self.results_compact_deconv_morpho, self.results_ext_conv_morpho,
         self.results_ext_deconv_morpho,
         self.components_deconv_props, self.components_conv_props,
         self.image_results_conv, self.image_results_deconv,self.bkg_images,
         self.class_resuts, self.compact_model,
         self.decomp_table) = \
            mlibs.run_image_fitting(imagelist=[self.input_data.filename],
                                    residuallist=[self.input_data.residualname],
                                    aspect=self.aspect,
                                    which_residual=self.which_residual,
                                    comp_ids=self.comp_ids,# which IDs refers to compact components?
                                    # The detected regions, so the diffuse
                                    # emission can be measured per region and
                                    # not only globally. Same objects
                                    # `structural_morphology` takes.
                                    indices=list(self.SE.indices),
                                    masks_deblended=self.SE.masks,
                                    region_grow_to_ref=self.region_grow_to_ref,
                                    sources_photometries=self.SE.sources_photometries,
                                    n_components=self.SE.n_components,
                                    z=self.z,
                                    convolution_mode=self.convolution_mode,
                                    method1=self.method1,
                                    method2=self.method2,
                                    mask=self.mask,
                                    use_mask_for_fit=self.use_mask_for_fit,
                                    mask_for_fit=self.mask_for_fit,
                                    bkg_map=self.bkg_map,
                                    self_bkg=self.self_bkg,
                                    rms_map=self.rms_map,
                                    use_weights = self.use_weights,
                                    is_bkg_map_conv = self.is_bkg_map_conv,
                                    save_name_append='',
                                    fix_n=self.fix_n,
                                    loss=self.loss,
                                    tr_solver = self.tr_solver,
                                    fix_value_n=self.fix_value_n,
                                    fix_max_value_n=self.fix_max_value_n,
                                    fix_min_value_n=self.fix_min_value_n,
                                    trunc=self.trunc,
                                    fix_max_value_Rn = self.fix_max_value_Rn,
                                    fix_min_value_Rn = self.fix_min_value_Rn,
                                    fix_geometry=self.fix_geometry,  # unstable if  False
                                    force_circular=self.force_circular,
                                    dr_fix=self.dr_fix,
                                    fix_x0_y0 = self.fix_x0_y0,
                                    parameters_mini_init = self.parameters_mini_init,
                                    sigma=self.sigma,
                                    logger=_logging_.logger,verbose=self.verbose)

        # `comp_ids` may have been left empty for `run_image_fitting` to work out
        # from the compact/diffuse classification. Read the resolved split back,
        # so `decompose_emission` and the summary tables see the IDs that were
        # actually used. This used to happen implicitly, by that function
        # appending to the very list object held here.
        self.comp_ids = list(
            self.SE.sources_photometries.get('comp_ids', self.comp_ids))
        self.ext_ids = list(
            self.SE.sources_photometries.get('ext_ids', []))

        # compute sizes
        try:
            self.cell_size = mlibs.get_cell_size(self.input_data.filename)
            self.pix_to_pc = \
                mlibs.pixsize_to_pc(z=self.z, cell_size=self.cell_size)
            self.size_unit = ' pc'
        except:
            self.cell_size = 1.0
            self.pix_to_pc = 1
            self.size_unit = ' px'

        self.beam_size_px = self.results_fit['beam_size_px']
        self.beam_size_pc = self.beam_size_px * self.pix_to_pc

        # 50% core-compact/unresolved deconvolved radii
        self.C50comp_radii_deconv = \
            self.results_compact_deconv_morpho['C50radii'] * self.pix_to_pc
        # 50% core-compact/unresolved convolved radii
        self.C50comp_radii_conv = \
            self.results_compact_conv_morpho['C50radii'] * self.pix_to_pc
        # 95% core-compact/unresolved deconvolved radii
        self.C95comp_radii_deconv = \
            self.results_compact_deconv_morpho['C95radii'] * self.pix_to_pc
        # 95% core-compact/unresolved convolved radii
        self.C95comp_radii_conv = \
            self.results_compact_conv_morpho['C95radii'] * self.pix_to_pc

        # 50% core-compact/unresolved deconvolved radii
        self.C50ext_radii_deconv = \
            self.results_ext_deconv_morpho['C50radii'] * self.pix_to_pc
        # 50% core-compact/unresolved convolved radii
        self.C50ext_radii_conv = \
            self.results_ext_conv_morpho['C50radii'] * self.pix_to_pc
        # 95% core-compact/unresolved deconvolved radii
        self.C95ext_radii_deconv = \
            self.results_ext_deconv_morpho['C95radii'] * self.pix_to_pc
        # 95% core-compact/unresolved convolved radii
        self.C95ext_radii_conv = \
            self.results_ext_conv_morpho['C95radii'] * self.pix_to_pc

        # Rn main core-compact/unresolved component (ID1)
        self.Rn_comp = self.results_fit['f1_Rn'] * self.pix_to_pc
        self.theta2_Rnfit = 2 * self.results_fit['f1_Rn'] * self.pix_to_pc
        self.theta1_Rnfit = (2 * (1-self.results_fit['f1_ell']) * self.results_fit['f1_Rn'] * self.pix_to_pc)
        # self.Rn_comp_err = self.results_fit['f1_Rn'][0].stderr * self.pix_to_pc

        mlibs.print_logger_header(title="Core-Compact Component Sizes",
                            logger=_logging_.logger)
        _logging_.logger.info(f" >=> 1 px = {self.pix_to_pc:.2f} pc")
        _logging_.logger.info(f" >=> Beam Size = "
                              f"{self.beam_size_px[0]:.2f} px")
        _logging_.logger.info(f" >=> Beam Size = "
                              f"{self.beam_size_pc[0]:.2f} {self.size_unit}")        

        _logging_.logger.info(f" >=> Rn Main Compact = "
                              f"{self.Rn_comp[0]:.2f} {self.size_unit}")
        _logging_.logger.info(f" >=> major axis FWHM = "
                              f"{self.theta2_Rnfit[0]:.2f} {self.size_unit}")

        _logging_.logger.info(f" >=> C50 Compact Deconv Radii = "
                              f"{self.C50comp_radii_deconv[0]:.2f} {self.size_unit}")
        _logging_.logger.info(f" >=> C50 Compact Conv Radii = "
                              f"{self.C50comp_radii_conv[0]:.2f} {self.size_unit}")
        _logging_.logger.info(f" >=> C95 Compact Deconv Radii = "
                              f"{self.C95comp_radii_deconv[0]:.2f} {self.size_unit}")
        _logging_.logger.info(f" >=> C95 Compact Conv Radii = "
                              f"{self.C95comp_radii_conv[0]:.2f} {self.size_unit}")

        mlibs.print_logger_header(title="Extended Component Sizes",
                            logger=_logging_.logger)
        _logging_.logger.info(f" >=> C50 Extended Deconv Radii = "
                              f"{self.C50ext_radii_deconv[0]:.2f} {self.size_unit} "
                              f"[flagged={self.results_ext_deconv_morpho['flag50'][0]}]")
        _logging_.logger.info(f" >=> C50 Extended Conv Radii = "
                              f"{self.C50ext_radii_conv[0]:.2f} {self.size_unit} "
                              f"[flagged={self.results_ext_conv_morpho['flag50'][0]}]")
        _logging_.logger.info(f" >=> C95 Extended Deconv Radii = "
                              f"{self.C95ext_radii_deconv[0]:.2f} {self.size_unit}"
                              f"[flagged={self.results_ext_deconv_morpho['flag9095'][0]}]")
        _logging_.logger.info(f" >=> C95 Extended Conv Radii = "
                              f"{self.C95ext_radii_conv[0]:.2f} {self.size_unit}"
                              f"[flagged={self.results_ext_conv_morpho['flag9095'][0]}]")

class sersic_multifit_general():
    """
    Multi-Sersic Fitting Decomposition.

    Perform a semi-automated and robust multi-sersic image decomposition.
    It supports GPU-acceleration using Jax. If no GPU is present, Jax still
    will benefit from CPU parallel processing. Do not worry, you do not have
    to change anything, Jax will automatically detect wheter you are running on
    CPU or GPU.
    
    This class it to help in modelling optical data, but is pure experimental.
    
    Major milestones: 
        - improve source detection
        - improve background estimation
        - improve PSF modelling, especially for JWST data. 
    """

    def __init__(self, input_data, SE,
                 fix_geometry = None,
                 force_circular = None,
                 comp_ids = ['1'],
                 default_component_type = 'sersic',
                 fix_n = None,
                 fix_value_n = None,
                 fix_max_value_Rn = None,
                 fix_min_value_Rn = None,
                 fix_max_value_n = None,
                 fix_min_value_n = None,
                 dr_fix = None, fix_x0_y0=None,
                 trunc = None,
                 constrained=True, self_bkg=False,
                 sigma=6.0, 
                 use_mask_for_fit=False,mask_for_fit=None,mask=None,
                 Npsf=1,
                 bkg_map = None,
                 rms_map = None,
                 use_weights = False,
                 sky_mode = None,
                 fit_background = None,
                 sky_scale_bounds = None,
                 background_mode = None,
                 fit_residual_background = None,
                 fit_background_scale = None,
                 rms_convention = 'sigma',
                 weight_mode = 'inverse_variance',
                 scale_covar = True,
                 loss='cauchy', tr_solver='exact',
                 regularize=True, f_scale=1.0, ftol=1e-10,
                 xtol=1e-10, gtol=1e-10,
                 init_params=0.2, final_params=5.0,
                 which_residual = 'user',
                 is_bkg_map_conv = False,
                 convolution_mode='GPU',method1='least_squares',
                 method2='least_squares',z = 0.01,
                 save_name_append = ''):

        self.input_data = input_data
        self.SE = SE
        if use_mask_for_fit == True:
            if mask_for_fit is None:
                self.mask_for_fit = self.SE.mask
            else:
                self.mask_for_fit = mask_for_fit
        else:
            self.mask_for_fit = None
            
        if mask is None:
            self.mask = self.SE.mask
        else:
            self.mask = mask
            
        # Component IDs are integers everywhere downstream (`comp_ID` columns,
        # `ext_ids`, the `model_cN` keys). Notebooks pass them as strings, so
        # normalise once here rather than at each consumer.
        self.comp_ids = [int(c) for c in (comp_ids or [])]
        self.convolution_mode = convolution_mode
        self.constrained = constrained
        self.Npsf = Npsf
        self.method1 = method1
        self.method2 = method2
        self.init_params = init_params
        self.final_params = final_params
        self.tr_solver = tr_solver
        self.regularize = regularize
        self.f_scale = f_scale
        self.ftol = ftol
        self.xtol = xtol
        self.gtol = gtol
        self.loss = loss
        self.z = z
        self.sigma = sigma
        self.which_residual = which_residual
        self.self_bkg = self_bkg
        self.bkg_map = bkg_map
        self.rms_map = rms_map
        self.use_weights = use_weights
        self.is_bkg_map_conv = is_bkg_map_conv

        if self.self_bkg == True and self.bkg_map is None:
            self.bkg_map = self.SE.bkg
        else:
            self.self_bkg = False

        # See the equivalent block in `sersic_multifit_radio`: the per-component
        # fit controls are built from each component's preset at exactly
        # `n_components` entries, with caller-supplied values layered on top.
        # Optical data defaults to the 'sersic' preset, which locks nothing, so
        # the Sersic index stays free unless a component says otherwise.
        self.default_component_type = default_component_type
        _controls = mlibs.build_fit_control_maps(
            self.SE.sources_photometries,
            self.SE.n_components,
            overrides={'fix_n': fix_n,
                       'fix_value_n': fix_value_n,
                       'fix_max_value_n': fix_max_value_n,
                       'fix_min_value_n': fix_min_value_n,
                       'fix_max_value_Rn': fix_max_value_Rn,
                       'fix_min_value_Rn': fix_min_value_Rn,
                       'dr_fix': dr_fix,
                       'fix_x0_y0': fix_x0_y0,
                       'fix_geometry': fix_geometry,
                       'trunc': trunc,
                       'force_circular': force_circular},
            default_component_type=self.default_component_type)

        self.fix_n = _controls['fix_n']
        self.fix_value_n = _controls['fix_value_n']
        self.fix_max_value_n = _controls['fix_max_value_n']
        self.fix_min_value_n = _controls['fix_min_value_n']
        self.fix_max_value_Rn = _controls['fix_max_value_Rn']
        self.fix_min_value_Rn = _controls['fix_min_value_Rn']
        self.dr_fix = _controls['dr_fix']
        self.fix_x0_y0 = _controls['fix_x0_y0']
        self.fix_geometry = _controls['fix_geometry']
        self.trunc = _controls['trunc']
        # Absent from `_controls` unless a component declared it or the caller
        # passed one, which is exactly the "nobody said anything" case that lets
        # `source_extraction`'s deprecated force_circular flags still apply.
        self.force_circular = _controls.get('force_circular')

        # ------------------------------------------------------------------
        # The sky. One convention, and a typed choice of what the sky IS:
        #
        #     residual = ( I - conv(M) - s_a * B ) * w
        #
        #     sky_mode='none'   no term; s_a frozen at 0
        #     sky_mode='flat'   B = 1, so s_a is the sky in image units
        #     sky_mode='map'    B = bkg_map, s_a a dimensionless multiplier
        #
        # The data is never modified in any of them. This is the IMFIT/GALFIT
        # split -- subtract the sky and fit no term, or fit a sky component on
        # un-subtracted data, never both.
        #
        # `is_background_subtracted` (declared on read_data) does NOT select a
        # data path. It only (a) supplies the default mode and where s_a starts,
        # and (b) guards against incoherent setups.
        already_subtracted = bool(getattr(self.input_data,
                                          'is_background_subtracted', False))

        if sky_mode is not None and sky_mode not in ('none', 'flat', 'map'):
            raise ValueError(
                f"sky_mode must be 'none', 'flat', 'map' or None, "
                f"got {sky_mode!r}.")

        if background_mode == 'subtract' and already_subtracted:
            # Guard A. Kept as a hard error: it means the caller is still
            # thinking in terms of removing a sky that is not there.
            raise ValueError(
                "background_mode='subtract' on an image declared "
                "is_background_subtracted=True is incoherent. The background is "
                "no longer subtracted from the data at all -- it is fitted on "
                "the model as s_a*bkg. Use sky_mode='map' and let s_a go to ~0, "
                "or sky_mode='none' for no sky term.")
        if fit_residual_background is not None:
            print("--==>> `fit_residual_background` is superseded by "
                  "`sky_mode`; the sky amplitude s_a is always fitted when the "
                  "term is on.")
            if fit_residual_background and fit_background is None:
                fit_background = True

        # Resolve the mode here too, so the guards below and the printout can
        # talk about the mode the fit will actually use. do_fit2D repeats the
        # same resolution, and agrees, because it is handed the resolved value.
        if sky_mode is None:
            if fit_background is False:
                sky_mode = 'none'
            elif fit_background is True:
                sky_mode = 'map' if self.bkg_map is not None else 'flat'
            elif background_mode is not None:
                sky_mode = 'none' if background_mode == 'none' else 'map'
            elif already_subtracted:
                sky_mode = 'none'
            else:
                sky_mode = 'map' if self.bkg_map is not None else 'flat'
        print(f"--==>> Sky: sky_mode='{sky_mode}'"
              f"{' (image declared already sky-subtracted)' if already_subtracted else ''}")

        # Guard B. A full sky pedestal offered for an image whose sky is gone.
        if already_subtracted and sky_mode == 'map' and self.bkg_map is not None:
            try:
                _b = mlibs.np.asarray(self.bkg_map, dtype=float)
                _b_med = float(mlibs.np.nanmedian(_b))
                _b_rms = float(mlibs.mad_std(_b[mlibs.np.isfinite(_b)]))
                if _b_rms > 0 and abs(_b_med) > _b_rms:
                    print(f"--==>> WARNING: the image is declared "
                          f"is_background_subtracted=True, but the supplied "
                          f"bkg_map has a pedestal (median {_b_med:.4g}) larger "
                          f"than its own spatial rms ({_b_rms:.4g}). Prefer "
                          f"sky_mode='flat', where s_a is the residual offset in "
                          f"image units and reads directly, or sky_mode='none'.")
            except Exception:
                pass

        # Guard C, the converse. Fitting no sky at all on an image that still
        # has one leaves the pedestal for the Sersic wings to absorb.
        if sky_mode == 'none' and not already_subtracted:
            print("--==>> WARNING: sky_mode='none' on an image NOT declared "
                  "is_background_subtracted=True. If the sky is still in the "
                  "data, the outer Sersic wings will absorb it and n will be "
                  "biased high. Use sky_mode='flat', or declare the image "
                  "subtracted on read_data.")

        self.sky_mode = sky_mode
        self.fit_background = fit_background
        self.sky_scale_bounds = sky_scale_bounds
        self.is_background_subtracted = already_subtracted
        self.background_mode = background_mode
        self.fit_residual_background = fit_residual_background
        self.fit_background_scale = fit_background_scale
        # How a supplied rms_map should be read ('sigma', 'variance', 'weight',
        # 'invvar') and how it becomes a weight. 'inverse_variance' is the
        # statistically correct w = 1/sigma; 'legacy' reproduces the old
        # w = 1/sqrt(sigma). Only affects fits with use_weights=True.
        self.rms_convention = rms_convention
        self.weight_mode = weight_mode
        self.scale_covar = scale_covar

        self.save_name_append = save_name_append

        self.__sersic_general()

    def __sersic_general(self):
        (self.result_mini, self.mini, self.result_1, self.result_extra,
         self.model_dict, self.image_results_conv, self.image_results_deconv, self.bkg_images,
         self.smodel2D,  self.model_temp) = \
            mlibs.do_fit2D(imagename=self.input_data.filename,
                           residualname=self.input_data.residualname,
                           init_constraints=self.SE.sources_photometries,
                           psf_name=self.input_data.psfname,
                           params_values_init_IMFIT=None,
                           ncomponents=self.SE.n_components,
                           constrained=self.constrained,
                           self_bkg=self.self_bkg,
                           bkg_map = self.bkg_map,
                           rms_map = self.rms_map,
                           use_weights = self.use_weights,
                           sky_mode = self.sky_mode,
                           fit_background = self.fit_background,
                           is_background_subtracted = self.is_background_subtracted,
                           sky_scale_bounds = self.sky_scale_bounds,
                           background_mode = self.background_mode,
                           fit_background_scale = self.fit_background_scale,
                           rms_convention = self.rms_convention,
                           weight_mode = self.weight_mode,
                           scale_covar = self.scale_covar,
                           is_bkg_map_conv = self.is_bkg_map_conv,
                           which_residual=self.which_residual,
                           observation_type = self.SE.obs_type,
                           mask_region = self.mask_for_fit,
                           Npsf=self.Npsf,
                           # rms_map=self.SE.bkg,
                           # rms_map=None,
                           fix_n=self.fix_n,
                           fix_value_n=self.fix_value_n,
                           fix_max_value_Rn=self.fix_max_value_Rn,
                           fix_min_value_Rn=self.fix_min_value_Rn,
                           fix_max_value_n = self.fix_max_value_n,
                           fix_min_value_n = self.fix_min_value_n,
                           dr_fix=self.dr_fix, fix_x0_y0 = self.fix_x0_y0,
                           convolution_mode=self.convolution_mode,
                           fix_geometry=self.fix_geometry, 
                           force_circular=self.force_circular,
                           trunc=self.trunc,
                           workers=-1, #not working 
                           method1=self.method1,
                           method2=self.method2,
                           init_params=self.init_params,
                           final_params=self.final_params,
                           loss=self.loss, tr_solver=self.tr_solver,
                           regularize=self.regularize, f_scale=self.f_scale,
                           ftol=self.ftol,
                           xtol=self.xtol, gtol=self.gtol,
                           save_name_append=self.save_name_append,
                           verbose=0,
                           logger=_logging_.logger)

        all_comps_ids = np.arange(1, self.SE.n_components + 1)
        mask_compact_ids = np.isin(all_comps_ids, np.asarray(self.comp_ids))
        ext_ids = [int(e) for e in all_comps_ids[~mask_compact_ids]]

        special_name = ''
        compact_model = 0
        extended_model = 0
        compact_model_deconv = 0
        extended_model_deconv = 0
        psf_fwhm = mlibs.psf_params(mlibs.load_fits_data(self.input_data.psfname))
        print(f" PSF FWHM = {psf_fwhm} px")
        # Same helper the fit mask uses, so the radius drawn on the diagnostic
        # plots is the radius that was actually excluded. These used to be
        # computed independently -- here, in _create_psf_exclusion_mask, and again
        # in plot_fit_results -- and disagreed by up to ~1 px.
        psf_px_size = mlibs.psf_exclusion_radius_px(psf_fwhm, self.Npsf)

        if self.input_data.rms_res == None:
            rms_std_res = self.input_data.rms_img
        else:
            rms_std_res = self.input_data.rms_res

        # Every `model_cN_conv` carries its own copy of the fitted sky (see the
        # `bkg_sign * bkg_comp_i` term where model_dict is built), so summing N
        # of them accumulates N backgrounds. Take the sky out per component and
        # add exactly one back, which is what run_image_fitting already does on
        # the radio side.
        for lc in self.comp_ids:
            compact_model = (compact_model +
                                self.model_dict[f'model_c{lc}_conv']
                                - self.model_dict['conv_bkg'])
            compact_model_deconv = (compact_model_deconv +
                                    self.model_dict[f'model_c{lc}']
                                    - self.model_dict['deconv_bkg'])
        compact_model = compact_model + self.model_dict['conv_bkg']
        compact_model_deconv = compact_model_deconv + self.model_dict['deconv_bkg']
        # if ext_ids is not None:
        if ext_ids == []:
            extended_model = 0
            extended_model_deconv = 0
            nfunctions = 1
        else:
            for le in ext_ids:
                extended_model = (extended_model +
                                    self.model_dict[f'model_c{le}_conv']
                                    - self.model_dict['conv_bkg'])
                extended_model_deconv = (extended_model_deconv +
                                            self.model_dict[f'model_c{le}']
                                            - self.model_dict['deconv_bkg'])
                nfunctions = None
            extended_model = extended_model + self.model_dict['conv_bkg']
            extended_model_deconv = (extended_model_deconv
                                     + self.model_dict['deconv_bkg'])

        # Both sides carry the background: `model_total_conv` includes the fitted
        # s_a*bkg, so the data must be the unmodified image. It used to subtract
        # the raw bkg_map from the data only, which removed a background from one
        # side of the comparison and left it on the other.
        zeta_metric = mlibs.fit2D_norm_metric(
            data = self.input_data.image_data_2D,
            model = self.model_dict['model_total_conv'],
            mask_region = self.mask,
            background_level = self.input_data.rms_res)
        self.model_dict['zeta_norm'] = zeta_metric
        self.decomposition_results = mlibs.plot_decomp_results(
            imagename=self.input_data.filename,
            compact=compact_model,
            extended_model=extended_model,
            data_2D_=self.input_data.image_data_2D,
            # The FITTED sky, s_a*B -- the same array plot_fit_results and
            # plot_fitting_summary are given, so all three agree on what the
            # background was.
            bkg_image=self.model_dict['conv_bkg'],
            rms=rms_std_res,
            nfunctions=nfunctions,
            obs_type=self.SE.obs_type,
            # Measure through the aperture the minimiser used, not a mask
            # re-derived from the data at a different threshold.
            mask=self.mask,
            result_mini=self.result_mini,
            comp_ids=self.comp_ids,
            ext_ids=ext_ids,
            model_total=self.model_dict['model_total_conv'],
            zeta_norm=self.model_dict['zeta_norm'],
            # The PSF image itself, not `psf_px_size`: that one is the
            # Npsf-scaled exclusion-hole radius (Npsf*FWHM/2), which is the
            # right thing to mask but the wrong thing to call a PSF size. The
            # resolved/unresolved test wants FWHM/2. Optical headers carry no
            # beam, so without this `beam_shape` returns its "one pixel is one
            # resolution element" sentinel and the report showed 0.5 px.
            psfname=self.input_data.psfname,
            special_name=special_name)
        # print(self.bkg_images[-1])
        mlibs.plot_fit_results(imagename=self.input_data.filename, 
                               model_dict=self.model_dict,
                               image_results_conv=self.image_results_conv,
                               sources_photometries=self.SE.sources_photometries,
                               result_mini=self.result_mini,
                               # The FITTED background, s_a*bkg, which is what the
                               # minimiser actually used. Passing the raw bkg_map
                               # here made the plot subtract a full, unscaled sky
                               # from a curve the fit had scaled by s_a -- and on
                               # an already-subtracted image, a sky that was never
                               # in the data at all.
                               bkg_image=self.bkg_images[-1],
                               crop=False, box_size=200,
                               mask=self.mask,
                            #    plotlim = int(4.0*mlibs.area_to_radii(np.nansum(self.mask))),
                            #    plotlim = int(4.0*mlibs.area_to_radii(np.nansum(self.mask))),
                               plotlim = 2.5*self.SE.sources_photometries['cg_Rp'],
                               psf_px_size=psf_px_size,
                               obs_type=self.SE.obs_type,
                               vmax_factor=0.3, vmin_factor=1.0)

        # mlibs.plot_slices(data_2D=self.input_data.image_data_2D,
        #                   model_dict=self.model_dict,
        #                   image_results_conv=self.image_results_conv[-2],
        #                   Rp_props=self.SE.sources_photometries,
        #                   residual_2D=None)
        
        mlibs.plot_fitting_summary(imagename = self.input_data.filename, 
                                   result_mini = self.result_mini, 
                                   model_dict = self.model_dict, 
                                   image_results_conv = self.image_results_conv,
                                   ncomponents=self.SE.sources_photometries['ncomps'],
                                   # Same fitted background as plot_fit_results.
                                   bkg_map = self.model_dict['conv_bkg'],
                                   observation_type=self.SE.obs_type,
                                   mask=self.mask,
                                #    plotlim = int(1.5*mlibs.area_to_radii(np.nansum(self.mask))),
                                   plotlim = 2.5*self.SE.sources_photometries['cg_Rp'],
                                   psf_px_size=psf_px_size,
                                   # so each component's radial profile is taken
                                   # about its own centre, not the global peak
                                   sources_photometries=self.SE.sources_photometries,
                                   rms=rms_std_res, 
                                #    cell_size=self.input_data.cell_size, 
                                   cell_size=1.0,
                                   show_figure=True, figsize=(12*0.8, 8*0.8),
                                #    box_size=None, crop=True,
                                #    save_name=None, 
                                   vmax_factor=0.5, vmin_factor=3, add_contours=True)

        # save results to csv file.
        try:
            mlibs.save_results_csv(result_mini=self.result_mini,
                                   other_param_dict={'zeta_norm': self.model_dict['zeta_norm']},
                            #  save_name=image_results_conv[-2].replace('.fits', ''),
                            save_name = self.input_data.filename.replace('.fits', ''),
                            ext='.csv',
                            save_corr=True, save_params=True)
        except:
            print('Error Saving Results to a csv file!!!')
            pass

        self.parameter_results = self.result_mini.params.valuesdict().copy()
        try:
            for param in self.result_mini.params.valuesdict().keys():
                if self.result_mini.params[param].stderr is not None:
                    self.parameter_results[param+'_err'] = self.result_mini.params[param].stderr
                else:
                    self.parameter_results[param+'_err'] = 0.1 * self.parameter_results[param]
        except:
            pass
        self.parameter_results['#imagename'] = os.path.basename(self.input_data.filename)
        # self.parameter_results['zeta_norm'] = self.model_dict['zeta_norm']
        # self.results_fit = {**self.parameter_results, **self.decomposition_results}
        self.results_fit = {**self.parameter_results}
        # self.results_fit.append(all_results)
        pass

class morphometry():
    """
    Core functionalities from Morfometryka. 
    Morfometryka is not publically available yet, 
    so these functions will be added in a later stage. 

    << 2025 Update >> Will not be added, as Morfometryka contains everything.
    Thus will be removed in future versions.
    """
    def __init__(self, input_data, aspect=None):
        self.input_data = input_data
        

    def _concentration(self):
        pass
    def _asymetry(self):
        pass
    def _momentum(self):
        pass
    def _sigma_psi(self):
        pass
    def _entropy(self):
        """
        This will be provided here soon, since it was my Master Thesis research.
        """
        pass
    def _kurvature(self):
        """
        This will be provided here soon, since it was my Master Thesis research.
        """
        pass
        
    pass

class decompose_emission():
    def __init__(self,input_data,SMFR,z=0.01,frequency=None):
        self.input_data = input_data
        self.SMFR = SMFR
        self.z = z
        _logging_.logger.info(f"Imagename: {os.path.basename(self.input_data.filename)}")
        _logging_.logger.info(f"Imagename full path: {self.input_data.filename}")
        if frequency is None:
            try:
                imhd = mlibs.imhead(self.input_data.filename)
                frequency = mlibs.imhd['refval'][2] / 1e9 # GHz
                _logging_.logger.info(f"Using frequency of {frequency:.2f} GHz for "
                                      f"star formation estimate.")
            except:
                try:
                    # GHz
                    frequency = mlibs.get_frequency(self.input_data.filename)
                    # _logging_.logger.info(f"Using frequency of {frequency:.2f} GHz for "
                    #                       f"star formation estimate.")
                except:
                    _logging_.logger.warning('Frequency may be wrong. Please, '
                                             'provide the frequency of the observation.')
                    frequency = 5.0  # GHz

        else:
            frequency = frequency
            
        self.frequency = frequency
        
        self.cell_size = mlibs.get_cell_size(self.input_data.filename)
        self.pix_to_pc = mlibs.pixsize_to_pc(z=self.z,
                                  cell_size=self.cell_size)
        self.beam_size_px = mlibs.get_beam_size_px(self.input_data.filename)[0]
        self.beam_size_pc = self.beam_size_px * self.pix_to_pc
        
        self.compute_compact_properties()
        self.compute_extended_properties()
        self.summary_decomp = {**self.compact_dec_summary,
                  **self.compact_conv_summary,
                  **self.ext_dec_summary,
                  **self.ext_conv_summary
                 }
        self.summary_decomp_df = mlibs.pd.DataFrame([self.summary_decomp])
        self.summary_decomp_df['imagename'] = os.path.basename(self.input_data.filename).replace('.fits','')

        self.decomp_table = self.build_decomp_table()

    def build_decomp_table(self):
        """
        The same derived quantities as `summary_decomp_df`, in long format.

        `summary_decomp_df` is one WIDE row per image, with the component number
        baked into every column name (`f2_Rn`, `comp_2_Snu`). That shape is what
        forced per-ID column juggling whenever several images or frequencies had
        to be combined. Here each component gets a ROW instead, with the ID in a
        `comp_ID` column, so `pd.concat` across images just works and a component
        is selected with a filter:

            t[(t.comp_ID == 2) & (t.domain == 'conv')]

        Nothing is recomputed -- the prefixes are stripped off the very same
        dicts, so the two views cannot disagree. `summary_decomp_df` is unchanged
        and still built first.

        Two rows per component: `domain='deconv'` carries the fitted and
        deconvolved quantities, `domain='conv'` the convolved ones.
        """
        globals_ = ('cell_size', 'pix_to_pc', 'beam_size_px', 'beam_size_pc',
                    'frequency')
        rows = []
        for is_compact, ids, dec, conv in (
                (True, self.SMFR.comp_ids,
                 self.compact_dec_summary, self.compact_conv_summary),
                (False, getattr(self.SMFR, 'ext_ids', []),
                 self.ext_dec_summary, self.ext_conv_summary)):
            for cid in ids:
                for domain, summary in (('deconv', dec), ('conv', conv)):
                    row = {'imagename': os.path.basename(
                               self.input_data.filename).replace('.fits', ''),
                           'freq': self.frequency,
                           'kind': 'component_physical',
                           'comp_ID': int(cid),
                           'region_ID': int(self.SMFR.SE.sources_photometries.get(
                               f'c{cid}_parent', 0) or 0),
                           'domain': domain,
                           'is_compact': is_compact}
                    for key in globals_:
                        if key in summary:
                            row[key] = summary[key]
                    # `f{id}_...` are the fitted parameters, `comp_{id}_...` the
                    # measured ones; both lose the ID, which now lives in its own
                    # column.
                    for prefix in (f'f{cid}_', f'comp_{cid}_'):
                        for key, value in summary.items():
                            if key.startswith(prefix):
                                row[key[len(prefix):]] = value
                    rows.append(row)
        return mlibs.assemble_decomp_table(rows)

    def _summarise_component(self, component_id, dec_summary, conv_summary,
                             label='Compact'):
        """
        Derive the physical quantities of ONE fitted component and write them
        into the two summary dicts.

        This is the body that `compute_compact_properties` and
        `compute_extended_properties` each carried a copy of -- 160 identical
        statements, differing only in which dicts they wrote to and the word in
        the log line. Keeping one copy is what makes the compact and diffuse
        sides provably consistent; they had already drifted apart by two keys
        (see the `s_Snu_A50` note below).

        `label` is used only for the log line. Values and key names are
        unchanged from the compact version.
        """

        mlibs.print_logger_header(title="Decomposition Results",
                logger=_logging_.logger)
        

        
        _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        _logging_.logger.info(f'{label} component ID {component_id}.')
        _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        
        # results from sersic fitting
        dec_summary[f'f{component_id}_Rn'] = \
            self.SMFR.result_mini.params[f'f{component_id}_Rn'].value * self.pix_to_pc
        dec_summary[f'f{component_id}_Rn_err'] = \
            self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr * self.pix_to_pc
        A_Rn,A_Rn_err = mlibs.radii_to_area(self.SMFR.result_mini.params[f'f{component_id}_Rn'].value,
                                            self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr)
        dec_summary[f'f{component_id}_A_Rn'] = \
            mlibs.pix_area_to_kpc_area(A_Rn,self.pix_to_pc)
        dec_summary[f'f{component_id}_A_Rn_err'] = \
            mlibs.pix_area_to_kpc_area(A_Rn_err,self.pix_to_pc)
        
        
        theta_Rn_maj_pc,theta_Rn_maj_pc_err,\
            theta_Rn_min_pc,theta_Rn_min_pc_err,\
            theta_Rn_pc,theta_Rn_pc_err = \
            mlibs.R50_to_fwhm(self.SMFR.result_mini.params[f'f{component_id}_Rn'].value,
                              self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr,
                              q_ratio = (1-self.SMFR.result_mini.params[f'f{component_id}_ell'].value),
                              scale = self.pix_to_pc)
            
        theta_Rn_maj_asec,theta_Rn_maj_asec_err, \
            theta_Rn_min_asec,theta_Rn_min_asec_err, \
            theta_Rn_asec,theta_Rn_asec_err = \
            mlibs.R50_to_fwhm(self.SMFR.result_mini.params[f'f{component_id}_Rn'].value,
                              self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr,
                              q_ratio = (1-self.SMFR.result_mini.params[f'f{component_id}_ell'].value),
                              scale = self.cell_size)
        
        # print('/-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+/')
        # print(np.sqrt(theta_Rn_maj_asec**2 + theta_Rn_min_asec**2))
        # print(np.sqrt(theta_Rn_maj_asec_err**2 + theta_Rn_min_asec_err**2))
        # print('/-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+/')
        theta_R50_maj_pc,theta_R50_maj_pc_err,\
            theta_R50_min_pc,theta_R50_min_pc_err,\
            theta_R50_pc,theta_R50_pc_err = \
            mlibs.R50_to_fwhm(self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii'],
                              self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii_err'],
                              scale=self.pix_to_pc)
        
        theta_R50_maj_asec,theta_R50_maj_asec_err, \
            theta_R50_min_asec,theta_R50_min_asec_err, \
            theta_R50_asec,theta_R50_asec_err = \
            mlibs.R50_to_fwhm(self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii'],
                              self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii_err'],
                              scale=self.cell_size)
        
        
        # Use the `peak_snr` column computed by `compute_image_properties`, which
        # is the peak over the *local* annulus rms. Dividing by `rms_residual`
        # here was a third, map-wide definition of the same quantity, and the
        # Condon (1997)/Fomalont (1999) size error below wants the local one.
        _cprops = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']
        peak_snr = _cprops['peak_snr'] if 'peak_snr' in self.SMFR.components_conv_props.columns \
            else _cprops['peak_of_flux'] / _cprops['rms_residual']


        # theta_Rn_maj_asec_err = self.SMFR.components_conv_props['bmajor'].iloc[0]/(2*peak_snr*theta_Rn_maj_asec)
        # theta_Rn_min_asec_err = self.SMFR.components_conv_props['bminor'].iloc[0]/(2*peak_snr*theta_Rn_min_asec)
        theta_Rn_maj_asec_err = self.SMFR.components_conv_props['bmajor'].iloc[0]/(2*peak_snr)
        theta_Rn_min_asec_err = self.SMFR.components_conv_props['bminor'].iloc[0]/(2*peak_snr)
        theta_Rn_asec_err = np.sqrt(theta_Rn_maj_asec_err**2 + theta_Rn_min_asec_err**2)

        
        Snu = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['total_flux_mask']
        Snu_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['total_flux_error']
        dec_summary[f'comp_{component_id}_Snu'] = Snu
        dec_summary[f'comp_{component_id}_Snu_err'] = Snu_err
        # `Snu` is measured on the CONVOLVED component (it comes straight from
        # `components_conv_props`), so it belongs on the conv row too -- without
        # this the `domain='conv'` rows of `decomp_table` carry sizes but a NaN
        # flux. The value is the same object used for `s_Snu_A50`/`s_Snu_A95`
        # below, so the two views cannot disagree, and `summary_decomp` merges
        # the dicts key-by-key with identical values, leaving it unchanged.
        conv_summary[f'comp_{component_id}_Snu'] = Snu
        conv_summary[f'comp_{component_id}_Snu_err'] = Snu_err
        # The same flux measured on the DECONVOLVED component, as a convergence
        # check: convolution conserves flux, so `Snu_d` and `Snu` should agree to
        # within the aperture and noise. It is deliberately a separate key -- SFR,
        # TB and the surface brightnesses all use `Snu` (the convolved one) and
        # must keep doing so. `_d` follows `R50_d` / `A_R50_d`.
        _dcp = self.SMFR.components_deconv_props.loc[
            f'model_c_deconv_{component_id}_props']
        dec_summary[f'comp_{component_id}_Snu_d'] = _dcp['total_flux_mask']
        dec_summary[f'comp_{component_id}_Snu_d_err'] = _dcp['total_flux_error']
        
        TB_Rn,TB_Rn_err  = mlibs.Tb_source(Snu=Snu,
                                    freq=self.frequency, #in GHz
                                    theta1=theta_Rn_maj_asec, 
                                    theta2=theta_Rn_min_asec,
                                    z=self.z,
                                    theta1_err=theta_Rn_maj_asec_err, 
                                    theta2_err=theta_Rn_min_asec_err,
                                    Snu_err = Snu_err)
        
        TB_R50,TB_R50_err  = mlibs.Tb_source(Snu=Snu,
                                    freq=self.frequency, #in GHz
                                    theta1=theta_R50_maj_asec, 
                                    theta2=theta_R50_min_asec,
                                    z=self.z,
                                    theta1_err=theta_R50_maj_asec_err, 
                                    theta2_err=theta_R50_min_asec_err,
                                    Snu_err = Snu_err)

        
        dec_summary[f'f{component_id}_theta_Rn_maj_pc'] = theta_Rn_maj_pc
        dec_summary[f'f{component_id}_theta_Rn_min_pc'] = theta_Rn_min_pc
        dec_summary[f'f{component_id}_theta_Rn_pc'] = theta_Rn_pc
        dec_summary[f'f{component_id}_theta_Rn_maj_pc_err'] = theta_Rn_maj_pc_err 
        dec_summary[f'f{component_id}_theta_Rn_min_pc_err'] = theta_Rn_min_pc_err 
        dec_summary[f'f{component_id}_theta_Rn_pc_err'] = theta_Rn_pc_err

        dec_summary[f'f{component_id}_theta_Rn_maj_asec'] = theta_Rn_maj_asec
        dec_summary[f'f{component_id}_theta_Rn_min_asec'] = theta_Rn_min_asec
        dec_summary[f'f{component_id}_theta_Rn_asec'] = theta_Rn_asec
        dec_summary[f'f{component_id}_theta_Rn_asec_err'] = theta_Rn_asec_err
        dec_summary[f'f{component_id}_theta_Rn_maj_asec_err'] = theta_Rn_maj_asec_err
        dec_summary[f'f{component_id}_theta_Rn_min_asec_err'] = theta_Rn_min_asec_err
        
        dec_summary[f'f{component_id}_theta_R50_pc'] = theta_R50_pc
        dec_summary[f'f{component_id}_theta_R50_maj_pc'] = theta_R50_maj_pc
        dec_summary[f'f{component_id}_theta_R50_min_pc'] = theta_R50_min_pc
        dec_summary[f'f{component_id}_theta_R50_pc_err'] = theta_R50_pc_err
        dec_summary[f'f{component_id}_theta_R50_maj_pc_err'] = theta_R50_maj_pc_err
        dec_summary[f'f{component_id}_theta_R50_min_pc_err'] = theta_R50_min_pc_err
        
        dec_summary[f'f{component_id}_theta_R50_asec'] = theta_R50_asec
        dec_summary[f'f{component_id}_theta_R50_maj_asec'] = theta_R50_maj_asec
        dec_summary[f'f{component_id}_theta_R50_min_asec'] = theta_R50_min_asec
        dec_summary[f'f{component_id}_theta_R50_asec_err'] = theta_R50_asec_err
        dec_summary[f'f{component_id}_theta_R50_maj_asec_err'] = theta_R50_maj_asec_err
        dec_summary[f'f{component_id}_theta_R50_min_asec_err'] = theta_R50_min_asec_err
        
        
        
        # dec_summary[f'f{component_id}_theta_Rn_maj_pc_err'] = theta_Rn_maj_asec_err * self.pix_to_pc/self.cell_size
        # dec_summary[f'f{component_id}_theta_Rn_min_pc_err'] = theta_Rn_min_asec_err * self.pix_to_pc/self.cell_size

        dec_summary[f'f{component_id}_Tb_Rn'] = TB_Rn
        dec_summary[f'f{component_id}_Tb_Rn_err'] = TB_Rn_err
        dec_summary[f'f{component_id}_Tb_R50'] = TB_R50
        dec_summary[f'f{component_id}_Tb_R50_err'] = TB_R50_err
        
    
        #properties from model images
        dec_summary[f'comp_{component_id}_R50_d'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii'] * self.pix_to_pc
        dec_summary[f'comp_{component_id}_R50_d_err'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii_err'] * self.pix_to_pc
        dec_summary[f'comp_{component_id}_R95_d'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C95radii'] * self.pix_to_pc
        dec_summary[f'comp_{component_id}_R95_d_err'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C95radii_err'] * self.pix_to_pc
        


        A_R50_d = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix50']
        A_R50_d_err = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix50_err']
        
        dec_summary[f'comp_{component_id}_A_R50_d'] = mlibs.pix_area_to_kpc_area(A_R50_d,self.pix_to_pc)
        dec_summary[f'comp_{component_id}_A_R50_d_err'] = mlibs.pix_area_to_kpc_area(A_R50_d_err,self.pix_to_pc)
            
        A_R95_d = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix95']
        A_R95_d_err = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix95_err']

        dec_summary[f'comp_{component_id}_A_R95_d'] = mlibs.pix_area_to_kpc_area(A_R95_d,self.pix_to_pc)
        dec_summary[f'comp_{component_id}_A_R95_d_err'] = mlibs.pix_area_to_kpc_area(A_R95_d_err,self.pix_to_pc)


        dec_summary[f'comp_{component_id}_s_Snu_A50_d'] = 0.5*Snu / dec_summary[f'comp_{component_id}_A_R50_d']
        dec_summary[f'comp_{component_id}_s_Snu_A50_d_err'] = 0.5*Snu_err / dec_summary[f'comp_{component_id}_A_R50_d']
        dec_summary[f'comp_{component_id}_s_Snu_A95_d'] = Snu / dec_summary[f'comp_{component_id}_A_R95_d']
        dec_summary[f'comp_{component_id}_s_Snu_A95_d_err'] = Snu_err / dec_summary[f'comp_{component_id}_A_R95_d']
        dec_summary[f'comp_{component_id}_s_Snu_A_Rn'] = 0.5*Snu / dec_summary[f'f{component_id}_A_Rn']
        dec_summary[f'comp_{component_id}_s_Snu_A_Rn_err'] = 0.5*Snu_err / dec_summary[f'f{component_id}_A_Rn']

        #CONVOLVED QUANTITIES
        conv_summary[f'comp_{component_id}_R50'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C50radii'] * self.pix_to_pc
        conv_summary[f'comp_{component_id}_R50_err'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C50radii_err'] * self.pix_to_pc
        conv_summary[f'comp_{component_id}_R95'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C95radii'] * self.pix_to_pc
        conv_summary[f'comp_{component_id}_R95_err'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C95radii_err'] * self.pix_to_pc



        A_R50 = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix50']
        A_R50_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix50_err']
        
        conv_summary[f'comp_{component_id}_A_R50'] = mlibs.pix_area_to_kpc_area(A_R50,self.pix_to_pc)
        conv_summary[f'comp_{component_id}_A_R50_err'] = mlibs.pix_area_to_kpc_area(A_R50_err,self.pix_to_pc)
            
        A_R95 = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix95']
        A_R95_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix95_err']

        conv_summary[f'comp_{component_id}_A_R95'] = mlibs.pix_area_to_kpc_area(A_R95,self.pix_to_pc)
        conv_summary[f'comp_{component_id}_A_R95_err'] = mlibs.pix_area_to_kpc_area(A_R95_err,self.pix_to_pc)
        # Was written only for the extended components. It is the same
        # quantity in both cases, so it is emitted for every component now;
        # this ADDS two columns to the compact summary and removes none.
        conv_summary[f'comp_{component_id}_s_Snu_A50'] = 0.5*Snu / conv_summary[f'comp_{component_id}_A_R50']
        conv_summary[f'comp_{component_id}_s_Snu_A50_err'] = 0.5*Snu_err / conv_summary[f'comp_{component_id}_A_R50']
        conv_summary[f'comp_{component_id}_s_Snu_A95'] = Snu / conv_summary[f'comp_{component_id}_A_R95']
        conv_summary[f'comp_{component_id}_s_Snu_A95_err'] = Snu_err / conv_summary[f'comp_{component_id}_A_R95']
        
        _logging_.logger.info(f"Snu = {1000*dec_summary[f'comp_{component_id}_Snu']:.2f}" 
                              f" +/- {1000*dec_summary[f'comp_{component_id}_Snu_err']:.2f} mJy")
        _logging_.logger.info(f"Rn = {dec_summary[f'f{component_id}_Rn']:.2f}" 
                              f" +/- {dec_summary[f'f{component_id}_Rn_err']:.2f} pc")
        _logging_.logger.info(f"R50_d = {dec_summary[f'comp_{component_id}_R50_d']:.2f}"
                              f" +/- {dec_summary[f'comp_{component_id}_R50_d_err']:.2f} pc")
        # _logging_.logger.info(f"R95_d = {dec_summary[f'comp_{component_id}_R95_d']:.2f} pc")
        _logging_.logger.info(f"R95_d = {dec_summary[f'comp_{component_id}_R95_d']:.2f}"
                              f" +/- {dec_summary[f'comp_{component_id}_R95_d_err']:.2f} pc")
        
        _logging_.logger.info(f"R50 = {conv_summary[f'comp_{component_id}_R50']:.2f}"
                              f" +/- {conv_summary[f'comp_{component_id}_R50_err']:.2f} pc")
        _logging_.logger.info(f"R95 = {conv_summary[f'comp_{component_id}_R95']:.2f}"
                              f" +/- {conv_summary[f'comp_{component_id}_R95_err']:.2f} pc")            

        _logging_.logger.info(f"Theta_Rn = {dec_summary[f'f{component_id}_theta_Rn_pc']:.2f}" 
                              f" +/- {dec_summary[f'f{component_id}_theta_Rn_pc_err']:.2f} pc")
        _logging_.logger.info(f"Theta_R50 = {dec_summary[f'f{component_id}_theta_R50_pc']:.2f}"
                                f" +/- {dec_summary[f'f{component_id}_theta_R50_pc_err']:.2f} pc")
        _logging_.logger.info(f"Tb_Rn = ({dec_summary[f'f{component_id}_Tb_Rn']:.2f}" 
                              f" +/- {dec_summary[f'f{component_id}_Tb_Rn_err']:.2f}) x 10^5 K")
        _logging_.logger.info(f"Tb_R50 = ({dec_summary[f'f{component_id}_Tb_R50']:.2f}"
                                f" +/- {dec_summary[f'f{component_id}_Tb_R50_err']:.2f}) x 10^5 K")


        _logging_.logger.info(f"A_Rn = ({1e3*dec_summary[f'f{component_id}_A_Rn']:.3f}" 
                              f" +/- {1e3*dec_summary[f'f{component_id}_A_Rn_err']:.3f}) x 1e-3 kpc^2")
        _logging_.logger.info(f"A50_d = ({1e3*dec_summary[f'comp_{component_id}_A_R50_d']:.3f}"
                              f" +/- {1e3*dec_summary[f'comp_{component_id}_A_R50_d_err']:.3f}) x 1e-3 kpc^2")
        _logging_.logger.info(f"A95_d = ({1e3*dec_summary[f'comp_{component_id}_A_R95_d']:.3f}"
                              f" +/- {1e3*dec_summary[f'comp_{component_id}_A_R95_d_err']:.3f}) x 1e-3 kpc^2")
        _logging_.logger.info(f"A50 = ({1e3*conv_summary[f'comp_{component_id}_A_R50']:.3f}"
                              f" +/- {1e3*conv_summary[f'comp_{component_id}_A_R50_err']:.3f}) x 1e-3 kpc^2")
        _logging_.logger.info(f"A95 = ({1e3*conv_summary[f'comp_{component_id}_A_R95']:.3f}"
                              f" +/- {1e3*conv_summary[f'comp_{component_id}_A_R95_err']:.3f}) x 1e-3 kpc^2")
        _logging_.logger.info(f"sSnu_A95_d = {dec_summary[f'comp_{component_id}_s_Snu_A95_d']:.3f}"
                              f" +/- {dec_summary[f'comp_{component_id}_s_Snu_A95_d_err']:.3f} Jy/kpc^2")
        _logging_.logger.info(f"sSnu_A95 = {conv_summary[f'comp_{component_id}_s_Snu_A95']:.3f}"
                              f" +/- {conv_summary[f'comp_{component_id}_s_Snu_A95_err']:.3f} Jy/kpc^2")

        # _logging_.logger.info(f"Rn =  {dec_summary[f'f{component_id}_Rn']}")

    def compute_compact_properties(self):
        """
        Computes the total flux Density of radio compact components.
        """
        compact_dec_summary = {}
        compact_conv_summary = {}
        
        compact_dec_summary['cell_size'] = self.cell_size
        compact_dec_summary['pix_to_pc'] = self.pix_to_pc
        compact_dec_summary['beam_size_px'] = self.beam_size_px
        compact_dec_summary['beam_size_pc'] = self.beam_size_pc
        compact_dec_summary['frequency'] = self.frequency
        
        compact_conv_summary['cell_size'] = self.cell_size
        compact_conv_summary['pix_to_pc'] = self.pix_to_pc
        compact_conv_summary['beam_size_px'] = self.beam_size_px
        compact_conv_summary['beam_size_pc'] = self.beam_size_pc
        compact_conv_summary['frequency'] = self.frequency
        
        
        _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        _logging_.logger.info(f'Total emission.')
        _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        # print(self.SMFR.data_properties['total_flux_mask'])
        A_R50_data = mlibs.pix_area_to_kpc_area(self.SMFR.data_properties['npix50'].values[0],self.pix_to_pc)
        A_R50_data_err = mlibs.pix_area_to_kpc_area(self.SMFR.data_properties['npix50_err'].values[0],self.pix_to_pc)
        A_R95_data = mlibs.pix_area_to_kpc_area(self.SMFR.data_properties['npix95'].values[0],self.pix_to_pc)
        A_R95_data_err = mlibs.pix_area_to_kpc_area(self.SMFR.data_properties['npix95_err'].values[0],self.pix_to_pc)
        

        _logging_.logger.info(f"Snu = {1000*self.SMFR.data_properties['total_flux_mask'].values[0]:.2f}" 
                                f" +/- {1000*self.SMFR.data_properties['total_flux_error'].values[0]:.2f} mJy")
        
        _logging_.logger.info(f"R50 = {self.pix_to_pc * self.SMFR.data_properties['C50radii'].values[0]:.2f}"
                                f" +/- {self.pix_to_pc * self.SMFR.data_properties['C50radii_err'].values[0]:.2f} pc")
        _logging_.logger.info(f"R95 = {self.pix_to_pc * self.SMFR.data_properties['C95radii'].values[0]:.2f}"
                                f" +/- {self.pix_to_pc * self.SMFR.data_properties['C95radii_err'].values[0]:.2f} pc")
        _logging_.logger.info(f"A50 = {A_R50_data:.3f}"
                                f" +/- {A_R50_data_err:.3f} [kpc^2]")
        _logging_.logger.info(f"A95 = {A_R95_data:.3f}"
                                f" +/- {A_R95_data_err:.3f} [kpc^2]")
        _logging_.logger.info(f"A50 = ({1e3*A_R50_data:.3f}"
                                f" +/- {1e3*A_R50_data_err:.3f}) x 1e3 [kpc^2]")
        _logging_.logger.info(f"A95 = ({1e3*A_R95_data:.3f}"
                                f" +/- {1e3*A_R95_data_err:.3f}) x 1e3 [kpc^2]")


        for component_id in self.SMFR.comp_ids:
            self._summarise_component(component_id, compact_dec_summary,
                                      compact_conv_summary,
                                      label='Compact')

        compact_conv_summary[f'Snu_comp_total'] = self.SMFR.results_compact_conv_morpho['total_flux_mask'].values[0]
        compact_conv_summary[f'Snu_comp_total_err'] = self.SMFR.results_compact_conv_morpho['total_flux_error'].values[0]
        _logging_.logger.info(f"Snu Compact Total = {1000*compact_conv_summary[f'Snu_comp_total']:.2f}" 
                                f" +/- {1000*compact_conv_summary[f'Snu_comp_total_err']:.2f} mJy")
         
        self.compact_dec_summary = compact_dec_summary
        self.compact_conv_summary = compact_conv_summary
        
        pass
    
    def compute_extended_properties(self):
        all_comps_ids = np.arange(1, self.SMFR.SE.n_components + 1)
        mask_compact_ids = np.isin(all_comps_ids, np.asarray(self.SMFR.comp_ids))
        self.SMFR.ext_ids = [int(e) for e in all_comps_ids[~mask_compact_ids]]
        
        ext_dec_summary = {}
        ext_conv_summary = {}
        
        ext_dec_summary['cell_size'] = self.cell_size
        ext_dec_summary['pix_to_pc'] = self.pix_to_pc
        ext_dec_summary['beam_size_px'] = self.beam_size_px
        ext_dec_summary['beam_size_pc'] = self.beam_size_pc
        ext_dec_summary['frequency'] = self.frequency
        
        ext_conv_summary['cell_size'] = self.cell_size
        ext_conv_summary['pix_to_pc'] = self.pix_to_pc
        ext_conv_summary['beam_size_px'] = self.beam_size_px
        ext_conv_summary['beam_size_pc'] = self.beam_size_pc
        ext_conv_summary['frequency'] = self.frequency
        
        
        for component_id in self.SMFR.ext_ids:
            self._summarise_component(component_id, ext_dec_summary,
                                      ext_conv_summary, label='Extended')

        # ext_conv_summary[f'Snu_ext_total'] = self.SMFR.results_ext_conv_morpho['total_flux_mask'].values[0]
        # ext_conv_summary[f'Snu_ext_total_err'] = self.SMFR.results_ext_conv_morpho['total_flux_error'].values[0]
        ext_conv_summary[f'Snu_ext_total'] = self.SMFR.data_properties['total_flux_mask'].values[0]-self.SMFR.results_compact_conv_morpho['total_flux_mask'].values[0]
        ext_conv_summary[f'Snu_ext_total_2'] = self.SMFR.results_fit['flux_density_ext']
        ext_conv_summary[f'Snu_ext_total_err'] = self.SMFR.data_properties['total_flux_error'].values[0]
        _logging_.logger.info(f"Snu Extended Total = {1000*ext_conv_summary[f'Snu_ext_total']:.2f}" 
                                f" +/- {1000*ext_conv_summary[f'Snu_ext_total_err']:.2f} mJy")
        
        _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        
        A50_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix50'].values[0],self.pix_to_pc)
        A50_err_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix50_err'].values[0],self.pix_to_pc)
        A95_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix95'].values[0],self.pix_to_pc)
        A95_err_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix95_err'].values[0],self.pix_to_pc)
        
        _logging_.logger.info(f"Snu Extended Data = {1000*self.SMFR.results_ext_conv_morpho[f'total_flux_mask'].values[0]:.2f}" 
                                f" +/- {1000*self.SMFR.results_ext_conv_morpho[f'total_flux_error'].values[0]:.2f} mJy")

        _logging_.logger.info(f"R50 ext data = {self.pix_to_pc*self.SMFR.results_ext_conv_morpho[f'C50radii'].values[0]:.2f}"
                                f" +/- {self.pix_to_pc*self.SMFR.results_ext_conv_morpho[f'C50radii_err'].values[0]:.2f} pc")
        _logging_.logger.info(f"R95 ext data = {self.pix_to_pc*self.SMFR.results_ext_conv_morpho[f'C95radii'].values[0]:.2f}"
                                f" +/- {self.pix_to_pc*self.SMFR.results_ext_conv_morpho[f'C95radii_err'].values[0]:.2f} pc")        
        
        _logging_.logger.info(f"A50 ext data = ({1e3*A50_kpc_ext_conv:.3f}"
                                f" +/- {1e3*A50_err_kpc_ext_conv:.3f}) x 1e-3 kpc^2")
        _logging_.logger.info(f"A95 ext data = ({1e3*A95_kpc_ext_conv:.3f}"
                                f" +/- {1e3*A95_err_kpc_ext_conv:.3f}) x 1e-3 kpc^2")
        
        ext_conv_summary['A50_kpc_ext_data'] = A50_kpc_ext_conv
        ext_conv_summary['A50_kpc_ext_data_err'] = A50_err_kpc_ext_conv
        ext_conv_summary['A95_kpc_ext_data'] = A95_kpc_ext_conv
        ext_conv_summary['A95_kpc_ext_data_err'] = A95_err_kpc_ext_conv
        
        self.ext_dec_summary = ext_dec_summary
        self.ext_conv_summary = ext_conv_summary
        
        pass

class radio_star_formation():
    """
    Compute star-formation estimates from radio emission, given the converstion
    law.
    """

    def __init__(self, input_data, SMFR,SMFR_decomp,z=0.01,
                 calibration_kind='Murphy12',
                 alpha = -0.85, alpha_NT = -0.85, frequency = None):
        self.input_data = input_data
        self.SMFR = SMFR
        self.SMFR_decomp = SMFR_decomp
        self.z = z
        self.frequency = frequency #in GHz
        self.alpha = alpha
        self.alpha_NT = alpha_NT
        self.cell_size = mlibs.get_cell_size(self.input_data.filename)
        self.pix_to_pc = mlibs.pixsize_to_pc(z=self.z,
                                  cell_size=self.cell_size)


        if frequency is None:
            try:
                imhd = mlibs.imhead(self.input_data.filename)
                frequency = mlibs.imhd['refval'][2] / 1e9 # GHz
                _logging_.logger.info(f"Using frequency of {frequency:.2f} GHz for "
                                      f"star formation estimate.")
            except:
                try:
                    # GHz
                    frequency = mlibs.get_frequency(self.input_data.filename)
                    _logging_.logger.info(f"Using frequency of {frequency:.2f} GHz for "
                                          f"star formation estimate.")
                except:
                    _logging_.logger.warning('Frequency may be wrong. Please, '
                                             'provide the frequency of the observation.')
                    frequency = 5.0  # GHz

        else:
            frequency = frequency
            
        self.frequency = frequency
        self.calibration_kind = calibration_kind
        self.compute_SFR()
        
        
        # self.compute_surface_areas_SFR()
        # self.brightness_temperature()


    # def compute_flux_compact(self):
    #     """
    #     Computes the total flux Density of radio compact components.
    #     """
    #     pass

    def compute_SFR(self):
        """
        Computes star-formation for radio extended components.
        """
        ## Compact components
        
        SFR_estimates_extended = {}
        SFR_estimates_compact = {}
        SFR_estimates_extended['imagename'] = os.path.basename(self.input_data.filename).replace('.fits','')
        SFR_estimates_compact['imagename'] = os.path.basename(self.input_data.filename).replace('.fits','')
        
        mlibs.print_logger_header(title="SFR Estimates (extended regions)",
                            logger=_logging_.logger)
        
        self.SFR_compact, self.SFR_compact_err = \
            mlibs.compute_SFR_general(flux=self.SMFR_decomp.compact_conv_summary[f'Snu_comp_total'],
                                        flux_error=self.SMFR_decomp.compact_conv_summary[f'Snu_comp_total_err'],
                                        frequency=self.frequency, z=self.z,
                                        alpha=self.alpha, alpha_NT=self.alpha_NT,
                                        calibration_kind=self.calibration_kind)
            
            
        self.SFR_extended, self.SFR_extended_err = \
            mlibs.compute_SFR_general(flux=self.SMFR_decomp.ext_conv_summary[f'Snu_ext_total'],
                                        flux_error=self.SMFR_decomp.ext_conv_summary[f'Snu_ext_total_err'],
                                        frequency=self.frequency, z=self.z,
                                        alpha=self.alpha, alpha_NT=self.alpha_NT,
                                        calibration_kind=self.calibration_kind)
        SFR_estimates_extended['SFR_ext_total'] = self.SFR_extended
        SFR_estimates_extended['SFR_ext_err'] = self.SFR_extended_err
        SFR_estimates_compact['SFR_compact_total'] = self.SFR_compact
        SFR_estimates_compact['SFR_compact_total_err'] = self.SFR_compact_err
    

        
        for component_id in self.SMFR.ext_ids:
            self.SFR, self.SFR_err = \
                mlibs.compute_SFR_general(flux=self.SMFR_decomp.ext_dec_summary[f'comp_{component_id}_Snu'],
                                          flux_error=self.SMFR_decomp.ext_dec_summary[f'comp_{component_id}_Snu_err'],
                                          frequency=self.frequency, z=self.z,
                                          alpha=self.alpha, alpha_NT=self.alpha_NT,
                                          calibration_kind=self.calibration_kind)
            self.sSFR_95_d = self.SFR/self.SMFR_decomp.ext_dec_summary[f'comp_{component_id}_A_R95_d']
            self.sSFR_95_d_err = self.SFR_err/self.SMFR_decomp.ext_dec_summary[f'comp_{component_id}_A_R95_d']
            self.sSFR_95_d_err_2 = self.sSFR_95_d * np.sqrt((self.SFR_err/self.SFR)**2.0 + (self.SMFR_decomp.ext_dec_summary[f'comp_{component_id}_A_R95_d_err']/self.SMFR_decomp.ext_dec_summary[f'comp_{component_id}_A_R95_d'])**2.0)
            self.sSFR_95 = self.SFR/self.SMFR_decomp.ext_conv_summary[f'comp_{component_id}_A_R95']
            self.sSFR_95_err = self.SFR_err/self.SMFR_decomp.ext_conv_summary[f'comp_{component_id}_A_R95']
            self.sSFR_95_err_2 = self.sSFR_95 * np.sqrt((self.SFR_err/self.SFR)**2.0 + (self.SMFR_decomp.ext_conv_summary[f'comp_{component_id}_A_R95_err']/self.SMFR_decomp.ext_conv_summary[f'comp_{component_id}_A_R95'])**2.0)
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            _logging_.logger.info(f'Extended component ID {component_id}.')
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            _logging_.logger.info(f" ==> SFR = {self.SFR:.2f} +/- {self.SFR_err:.2f} "
                                f"Mo/yr")
            _logging_.logger.info(f" ==> sSFR (deconvolved areas) = {self.sSFR_95_d:.2f} +/- {self.sSFR_95_d_err:.2f} "
                                f"Mo/yr/kpc^2")
            _logging_.logger.info(f" ==> sSFR (convolved areas) = {self.sSFR_95:.2f} +/- {self.sSFR_95_err:.2f} "
                                f"Mo/yr/kpc^2")
            _logging_.logger.info(f" ==> Testing new error computation")
            _logging_.logger.info(f" ==> sSFR (deconvolved areas) = {self.sSFR_95_d:.2f} +/- {self.sSFR_95_d_err_2:.2f} "
                                f"Mo/yr/kpc^2")
            _logging_.logger.info(f" ==> sSFR (convolved areas) = {self.sSFR_95:.2f} +/- {self.sSFR_95_err_2:.2f} "
                                f"Mo/yr/kpc^2")
            

            SFR_estimates_extended[f'comp_{component_id}_SFR'] = self.SFR
            SFR_estimates_extended[f'comp_{component_id}_SFR_err'] = self.SFR_err
            SFR_estimates_extended[f'comp_{component_id}_sSFR_95_d'] = self.sSFR_95_d
            SFR_estimates_extended[f'comp_{component_id}_sSFR_95_d_err'] = self.sSFR_95_d_err
            SFR_estimates_extended[f'comp_{component_id}_sSFR_95'] = self.sSFR_95
            SFR_estimates_extended[f'comp_{component_id}_sSFR_95_err'] = self.sSFR_95

        _logging_.logger.info(f" ==> SFR (extended total) = {self.SFR_extended:.2f} +/- {self.SFR_extended_err:.2f} "
                            f"Mo/yr")

        mlibs.print_logger_header(title="SFR Estimates (compact regions)",
                            logger=_logging_.logger)
        
        _logging_.logger.warning(f'The SFRs below are hypothetical. These would be the SFR on compact/unresolved regions.\n'
                                f'                              If the region is a starburst, then these SFR have a meaning (nuclear compact star formation).\n'
                                f'                              If the region is a AGN, then the SFR is unrealistic and should not\n'
                                f'                              be interpreted as real star formation.')
        
        for component_id in self.SMFR.comp_ids:
            self.SFR, self.SFR_err = \
                mlibs.compute_SFR_general(flux=self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_Snu'],
                                          flux_error=self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_Snu_err'],
                                          frequency=self.frequency, z=self.z,
                                          alpha=self.alpha, alpha_NT=self.alpha_NT,
                                          calibration_kind=self.calibration_kind)
            self.sSFR_95_d = self.SFR/self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_A_R95_d']
            self.sSFR_95_d_err = self.SFR_err/self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_A_R95_d']
            self.sSFR_95_d_err_2 = self.sSFR_95_d * np.sqrt((self.SFR_err/self.SFR)**2.0 + (self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_A_R95_d_err']/self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_A_R95_d'])**2.0)
            self.sSFR_95 = self.SFR/self.SMFR_decomp.compact_conv_summary[f'comp_{component_id}_A_R95']
            self.sSFR_95_err = self.SFR_err/self.SMFR_decomp.compact_conv_summary[f'comp_{component_id}_A_R95']
            self.sSFR_95_err_2 = self.sSFR_95 * np.sqrt((self.SFR_err/self.SFR)**2.0 + (self.SMFR_decomp.compact_conv_summary[f'comp_{component_id}_A_R95_err']/self.SMFR_decomp.compact_conv_summary[f'comp_{component_id}_A_R95'])**2.0)

            
            
            beam_area_kpc = mlibs.pix_area_to_kpc_area(self.input_data.beam_area_px,self.pix_to_pc)

            # self.sSFR_95_corr = self.SFR/np.sqrt(self.SMFR_decomp.compact_conv_summary[f'comp_{component_id}_A_R95']*beam_area_kpc)
            # self.sSFR_95_corr_err = self.SFR_err/np.sqrt(self.SMFR_decomp.compact_conv_summary[f'comp_{component_id}_A_R95']*beam_area_kpc)
            
            # self.sSFR_95_d_corr = self.SFR/np.sqrt(self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_A_R95_d']*beam_area_kpc)
            # self.sSFR_95_d_corr_err = self.SFR_err/np.sqrt(self.SMFR_decomp.compact_dec_summary[f'comp_{component_id}_A_R95_d']*beam_area_kpc)
            
            self.sSFR_95_corr = self.SFR/np.sqrt(beam_area_kpc*beam_area_kpc)
            self.sSFR_95_corr_err = self.SFR_err/np.sqrt(beam_area_kpc*beam_area_kpc)
            
            self.sSFR_95_d_corr = self.SFR/np.sqrt(beam_area_kpc*beam_area_kpc)
            self.sSFR_95_d_corr_err = self.SFR_err/(beam_area_kpc*beam_area_kpc)
            
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            _logging_.logger.info(f'Compact component ID {component_id}.')
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            _logging_.logger.info(f" ==> SFR = {self.SFR:.2f} +/- {self.SFR_err:.2f} "
                                f"Mo/yr")
            _logging_.logger.info(f" ==> sSFR (deconvolved areas) = {self.sSFR_95_d:.2f} +/- {self.sSFR_95_d_err:.2f} "
                                f"Mo/yr/kpc^2")
            _logging_.logger.info(f" ==> sSFR (convolved areas) = {self.sSFR_95:.2f} +/- {self.sSFR_95_err:.2f} "
                                f"Mo/yr/kpc^2")
            _logging_.logger.info(f" ==> Testing new error computation")
            _logging_.logger.info(f" ==> sSFR (deconvolved areas) = {self.sSFR_95_d:.2f} +/- {self.sSFR_95_d_err_2:.2f} "
                                f"Mo/yr/kpc^2")
            _logging_.logger.info(f" ==> sSFR (convolved areas) = {self.sSFR_95:.2f} +/- {self.sSFR_95_err_2:.2f} "
                                f"Mo/yr/kpc^2")

            _logging_.logger.info(f" ==> Beam Corrected sSFR")
            _logging_.logger.info(f'Beam Area px {self.input_data.beam_area_px}.')
            _logging_.logger.info(f'Beam Area kpc {beam_area_kpc}.')
            _logging_.logger.info(f" ==> sSFR-corr (deconvolved areas) = {self.sSFR_95_d_corr:.2f} +/- {self.sSFR_95_d_corr_err:.2f} "
                                f"Mo/yr/kpc^2")
            _logging_.logger.info(f" ==> sSFR-corr (convolved areas) = {self.sSFR_95_corr:.2f} +/- {self.sSFR_95_corr_err:.2f} "
                                f"Mo/yr/kpc^2")

            SFR_estimates_compact[f'comp_{component_id}_SFR'] = self.SFR
            SFR_estimates_compact[f'comp_{component_id}_SFR_err'] = self.SFR_err
            SFR_estimates_compact[f'comp_{component_id}_sSFR_95_d'] = self.sSFR_95_d
            SFR_estimates_compact[f'comp_{component_id}_sSFR_95_d_err'] = self.sSFR_95_d_err
            SFR_estimates_compact[f'comp_{component_id}_sSFR_95'] = self.sSFR_95
            SFR_estimates_compact[f'comp_{component_id}_sSFR_95_err'] = self.sSFR_95

            
            
        _logging_.logger.info(f" ==> SFR (compact total) = {self.SFR_compact:.2f} +/- {self.SFR_compact_err:.2f} "
                            f"Mo/yr")
        
        self.SFR_estimates_extended = SFR_estimates_extended
        self.SFR_estimates_compact = SFR_estimates_compact


    def compute_surface_areas_SFR(self):
        """
        From the Sersic fitting results, compute the convolved and deconvolved areas. 
        These will be used to determine the surface density star formation rates.

        """
        def get_area_in_kpc(df,region,pix_to_pc):
            area_region = (df[region] * df['beam_area'] * (pix_to_pc**2.0)/(1000**2.0))
            return(area_region)

        # The 50% deconvolved area for core-compact/unresolved components
        self.A50_kpc_comp_deconv = get_area_in_kpc(self.SMFR.results_compact_deconv_morpho,
                                             region='A50',
                                             pix_to_pc = self.pix_to_pc)
        # The 50% convolved area for core-compact/unresolved components
        self.A50_kpc_comp_conv = get_area_in_kpc(self.SMFR.results_compact_conv_morpho,
                                           region='A50',
                                           pix_to_pc = self.pix_to_pc)
        # The 50% deconvolved area for diffuse components
        self.A50_kpc_ext_deconv = get_area_in_kpc(self.SMFR.results_ext_deconv_morpho,
                                            region='A50',
                                            pix_to_pc = self.pix_to_pc)
        # The 50% convolved area for diffuse components
        self.A50_kpc_ext_conv = get_area_in_kpc(self.SMFR.results_ext_conv_morpho,
                                          region='A50',
                                          pix_to_pc = self.pix_to_pc)
        # The 95% deconvolved area for core-compact/unresolved components
        self.A95_kpc_comp_deconv = get_area_in_kpc(self.SMFR.results_compact_deconv_morpho,
                                             region='A95',
                                             pix_to_pc = self.pix_to_pc)
        # The 95% convolved area for core-compact/unresolved components
        self.A95_kpc_comp_conv = get_area_in_kpc(self.SMFR.results_compact_conv_morpho,
                                           region='A95',
                                           pix_to_pc = self.pix_to_pc)
        # The 95% deconvolved area for diffuse components
        self.A95_kpc_ext_deconv = get_area_in_kpc(self.SMFR.results_ext_deconv_morpho,
                                            region='A95',
                                            pix_to_pc = self.pix_to_pc)
        # The 95% convolved area for diffuse components
        self.A95_kpc_ext_conv = get_area_in_kpc(self.SMFR.results_ext_conv_morpho,
                                          region='A95',
                                          pix_to_pc = self.pix_to_pc)

        self.sSFR50_deconv_ext = self.SFR_ext / self.A50_kpc_ext_deconv
        self.sSFR50_conv_ext = self.SFR_ext / self.A50_kpc_ext_conv
        self.sSFR95_deconv_ext = self.SFR_ext / self.A95_kpc_ext_deconv
        self.sSFR95_conv_ext = self.SFR_ext / self.A95_kpc_ext_conv


        mlibs.print_logger_header(title="sSFR Estimates",
                            logger=_logging_.logger)
        _logging_.logger.info(f" >=> Deconvolved A50 sSFR = "
                              f"{self.sSFR50_deconv_ext[0]:.2f} Mo/(yr kpc^2)")
        _logging_.logger.info(f" >=> Convolved A50 sSFR = "
                              f"{self.sSFR50_conv_ext[0]:.2f} Mo/(yr kpc^2)")
        _logging_.logger.info(f" >=> Deconvolved A95 sSFR = "
                              f"{self.sSFR95_deconv_ext[0]:.2f} Mo/(yr kpc^2)")
        _logging_.logger.info(f" >=> Convolved A95 sSFR = "
                              f"{self.sSFR95_conv_ext[0]:.2f} Mo/(yr kpc^2)")


        

    def brightness_temperature(self):
        """
        To Do: Loop over all model components, and compute TB individually.

        """
        self.theta1_Rnfit = 2 * self.cell_size * self.SMFR.results_fit['f1_Rn']
        self.theta2_Rnfit = (2 * (1 - self.SMFR.results_fit['f1_ell']) * self.cell_size *
                             self.SMFR.results_fit['f1_Rn'])

        self.TB_Rnfit  = mlibs.Tb_source(Snu=self.SMFR.results_compact_deconv_morpho[
            'total_flux_mask'],
                                   freq=self.frequency,
                                   theta1=self.theta1_Rnfit, theta2=self.theta2_Rnfit,
                                   z=self.z)

        self.theta1_R50fit = 2 * self.cell_size * self.SMFR.results_compact_deconv_morpho[
            'C50radii']
        self.theta2_R50fit = (2 * (self.SMFR.results_compact_deconv_morpho['qm']) *
                              self.cell_size *
                             self.SMFR.results_compact_deconv_morpho['C50radii'])
        self.TB_R50fit  = mlibs.Tb_source(Snu=self.SMFR.results_compact_deconv_morpho['total_flux_mask'],
                                   freq=self.frequency,
                                   theta1=self.theta1_R50fit, theta2=self.theta2_R50fit,
                                   z=self.z)
        mlibs.print_logger_header(title="Brightness Temperature",
                            logger=_logging_.logger)
        _logging_.logger.info(f" ==> TB Rn Fit = {self.TB_Rnfit[0]:.2f} e5 K")
        _logging_.logger.info(f" ==> TB R50 Fit = {self.TB_R50fit[0]:.2f} e5 K")



class make_plots():
    """
    """
    pass

class save_results():
    """
    """
    pass


class wsclean_imaging():
    """
    """
    pass


class casa_imaging():
    """
    """
    pass

class selfcalibration():
    """
    """
    pass


if __name__ == '__main__':

    # --- (DEV) New-style CLI mode scaffolding ---------------------------
    # Registry of new run modes: each mode is a boolean CLI flag mapped to a
    # handler function taking the parsed `args`. To add a mode later, write
    # one handler + one `@register_mode(...)` line -- nothing else here, and
    # nothing in the legacy `-mode` chain below, needs to change.
    MODES = {}

    def register_mode(flag_name):
        def _decorator(func):
            MODES[flag_name] = func
            return func
        return _decorator

    def _parse_kwargs(kwargs_list):
        """Turn ['vmax_factor=0.3', 'CM=viridis', "scalebar_length=2*mlibs.u.kpc"]
        into a dict. Each item is tried, in order, as:
          1. a Python literal (ast.literal_eval) -- numbers/bools/tuples/strings.
          2. a Python expression (eval), with `mlibs` (and therefore mlibs.np,
             mlibs.u, mlibs.find_z_NED, etc.) in scope -- covers unit-bearing
             values and one-off lookups.
          3. otherwise, kept as a plain string.

        Note: `eval` here runs whatever expression you pass on the CLI, same
        as pasting it into the notebook -- fine for this single-user research
        tool, but don't wrap this in anything that accepts untrusted input.

        Shell quoting: each key=value pair must arrive as ONE argv token, so
        anything containing spaces or shell metacharacters (`*`, `(`, `)`,
        quotes, ...) needs to be quoted as a whole, e.g.:
            --kwargs "scalebar_length=2 * mlibs.u.kpc" \\
                     "source_distance=mlibs.find_z_NED('UGC02369', return_luminosity_distance=True)[1]"
        Without the quotes, the shell splits on the spaces (so 'scalebar_length=2',
        '*', 'mlibs.u.kpc' arrive as three separate items -- and expands the bare
        '*' as a filename glob) before this function ever sees it.
        """
        import ast
        parsed = {}
        for item in kwargs_list or []:
            key, sep, raw_value = item.partition('=')
            if not sep:
                raise ValueError(
                    f"--kwargs item {item!r} has no '=' -- if this was meant to be "
                    f"part of a value with spaces, quote the whole key=value pair.")
            try:
                parsed[key] = ast.literal_eval(raw_value)
            except (ValueError, SyntaxError):
                try:
                    parsed[key] = eval(raw_value)
                except Exception as e:
                    raise ValueError(
                        f"Could not parse --kwargs value for {key!r}: {raw_value!r} ({e})")
        return parsed

    @register_mode('plot_image')
    def mode_plot_image(args):
        input_data = read_data(filename=args.filename,
                                residualname=args.residualname,
                                psfname=args.psfname)
        # Header-derived flux_units/flux_conversion_factor are defaults only:
        # eimshow uses flux_conversion_factor only inside its 'mJy/px' branch,
        # so flux_units and flux_conversion_factor must be set together and
        # consistently (see read_data.get_eimshow_flux_kwargs) -- passing the
        # factor alone, as before, was silently ignored for any BUNIT other
        # than the default flux_units='mJy'. Any matching key in --kwargs
        # still overrides these (dict.update keeps the last value, and
        # extra_kwargs is applied last).
        call_kwargs = dict(cell_size=input_data.cell_size,
                           rms=input_data.rms_res,
                           save_name=input_data.filename.replace('.fits', ''))
        call_kwargs.update(input_data.get_eimshow_flux_kwargs())
        call_kwargs.update(_parse_kwargs(args.kwargs))
        mlibs.eimshow(input_data.filename, **call_kwargs)
    # ---------------------------------------------------------------------

    parser = argparse.ArgumentParser(description='Morphen.')
    parser.add_argument('-filename',  '--filename',  required=False, \
                        help='Image data.')
    parser.add_argument('-residualname',  '--residualname',  required=False, \
                        help='Associated Residual from Image Data')
    parser.add_argument('-psfname',  '--psfname',  required=False, \
                        help='PSF name.')
    parser.add_argument('-mode',  '--mode',
                        required=False,  nargs='?',  
                        default = '',
                        # default='general_decomp',
                        const=True,
                        help=' (DEV) Mode of operation. 2D_fit (general_decomp or radio_decomp) or other (general image analysis).')
    parser.add_argument('-image_stats',  '--image_stats',
                        required=False,  nargs='?',  default=False,
                        const=True,
                        help='Compute basic image statistics.')

    parser.add_argument('-find_sources',  '--find_sources',
                        required=False,  nargs='?',  default=False,
                        const=True,
                        # action='store_true',
                        help='Perform source extraction and basic photometry.')

    parser.add_argument('--sigma', type=float, default=10,
                        help='Sigma value for source extraction.')
    parser.add_argument('--ell_size_factor', type=float, default=2,
                        help='Factor of ellipse size for source plot artistics/statistics.')

    parser.add_argument('-sersic_radio',  '--sersic_radio',
                        required=False,  nargs='?',  default=False,
                        const=True,
                        help='Perform Sersic Image Fitting for radio data.')



    parser.add_argument('--solver2', type=str, default='least_squares',
                        help='2nd run solver method (nelder or least_squares).')

    parser.add_argument('-SFR-do',  '--SFR-do',
                        required=False,  nargs='?',  default=False,
                        const=True,
                        help='Compute SFR Estimates.')

    parser.add_argument('--redshift', type=float, default=0.01,
                        help='Redshift of the source.')

    parser.add_argument('-general_fit',  '--general_fit',
                        required=False,  nargs='?',  default=False,
                        const=True,
                        help='Perform Sersic Image Fitting.')

    parser.add_argument('-obs_type',  '--obs_type',
                        required=False,  nargs='?',  default='radio',
                        const=True,
                        help='Which kind of observations (radio or other?)')

    # parser.add_argument('-sersic_optical',  '--sersic_optical',
    #                     required=False,  nargs='?',  default=False,
    #                     const=True,
    #                     help='Perform Sersoc Image Fitting for optical data.')

    parser.add_argument('-noshow',  '--noshow', required=False,  nargs='?',
                        default=False,  const=True,
                        help='Do not show plots.')

    # --- (DEV) New-style run-mode flags -----------------------------------
    parser.add_argument('-plot_image', '--plot-image',
                        dest='plot_image', required=False,  nargs='?',
                        default=False,  const=True,
                        help='(DEV) Load --filename via read_data and plot it with eimshow.')
    parser.add_argument('--kwargs', nargs='*', default=[],
                        help="(DEV) Extra key=value pairs forwarded as kwargs to the "
                             "current mode's underlying function. Values are parsed as "
                             "Python literals first, then as Python expressions with "
                             "`mlibs` in scope (so mlibs.u, mlibs.np, mlibs.find_z_NED(...) "
                             "work). Quote each key=value pair as ONE argument if it "
                             "contains spaces, e.g.: "
                             "--kwargs vmax_factor=0.3 CM=viridis "
                             "\"scalebar_length=2 * mlibs.u.kpc\"")
    parser.add_argument('-show',  '--show',
                        dest='show', required=False,  nargs='?',
                        default=False,  const=True,
                        help='(DEV) Display figures on screen (default: no figures are shown, '
                             'e.g. for headless/batch runs).')
    # -----------------------------------------------------------------------
    # parser.add_argument('-filename',  '--filename',  required=False, \
    #                     help='The input array of the light profile of the galaxy.')

    args = parser.parse_args()
    # from matplotlib import use as mpluse
    # mpluse("Agg")

    # Default is headless (no figures shown), matching the current behaviour
    # of every mode below. `--show` opts in. Registered via atexit rather than
    # a plain trailing `plt.show()` call because several modes below (both
    # legacy and new) call `sys.exit()` partway through -- atexit still fires
    # on that path, so this one flag covers every mode without editing them.
    if args.show:
        atexit.register(mlibs.plt.show)

    # New-style modes take priority and exit right after running, leaving the
    # legacy `-mode`/`-image_stats`/`-find_sources`/... chain below untouched.
    for _flag_name, _handler in MODES.items():
        if getattr(args, _flag_name):
            _handler(args)
            sys.exit()

    if args.filename != None:
        if args.mode == "general_decomp":
            print("2D FIT MODE")

            if args.obs_type == "other":
                if args.psfname != None:
                    psf_name = args.psfname
                else:
                    #should return error and stop the code
                    psf_name = None
                    sys.exit("Error: Please, provide a psf image.")

                input_data=read_data(filename=args.filename,
                                     psfname=psf_name)
                
                _, mask_region = mlibs.t_mask_dilation(input_data.filename, 
                                            sigma=6, 
                                            dilation_size=3,
                                            iterations=3, show_figure=False,PLOT=False,
                                            use_distance_filter=False,
                                            rms=input_data.rms_img)

                _, mask_for_fit = mlibs.t_mask_dilation(input_data.filename, 
                                            sigma=8, 
                                            dilation_size=10,
                                            iterations=10, show_figure=False,PLOT=False,
                                            do_filtering=True,
                                            rms=input_data.rms_img)
                #this is a set of parameters that may work OKAY
                bwf, bhf = 0.5, 0.5 # standard value that works for a wide range of images. But, additional check is required!!!!
                fwf, fhf = 0.5, 0.5 # standard value that works for a wide range of images. But, additional check is required!!!!
                clean_param = 0.1
                deblend_cont = 1e-1
                deblend_nthresh = 15
                sigma_mask = 12
                sigma = 12
                ell_size_factor = 1
                minarea_factor = 4
                SE = source_extraction(input_data, 
                                        ell_size_factor = ell_size_factor, sigma = sigma,
                                        bwf = bwf,bhf = bhf,fwf = fwf, fhf = fhf,
                                        clean_param = clean_param, 
                                        deblend_cont = deblend_cont, 
                                        deblend_nthresh=deblend_nthresh,minarea_factor=minarea_factor,
                                        apply_mask=False,sigma_mask=sigma_mask,dilation_size=3,
                                        mask=mask_region,first_ID_only=True, sort_by = 'distance',
                                        show_petro_plots=False,algorithm='PF',show_bkg_map=False,
                                        dry_run=True, obs_type = 'other')
                # ids_to_add = ['1']
                ids_to_add = []
                comp_ids = ['1']
                SE = source_extraction(input_data, ids_to_add=ids_to_add,
                                        ell_size_factor = ell_size_factor, sigma = sigma,
                                        bwf=bwf, bhf=bhf, fwf=fwf, fhf=fhf,
                                        clean_param=clean_param, 
                                        deblend_cont=deblend_cont, minarea_factor=minarea_factor,
                                        deblend_nthresh=deblend_nthresh,
                                        apply_mask=False,sigma_mask=sigma_mask,
                                        mask = mask_region,
                                        show_petro_plots=False,obs_type = 'other',algorithm='PF', 
                                        SE_ref = SE,
                                        show_bkg_map=False,force_circular=False,
                                        dry_run=False)
                

                psf_fwhm = mlibs.psf_params(mlibs.load_fits_data(input_data.psfname))

                final_bkg, final_rms, diagnostics = mlibs.multiscale_segmentation_background(input_data.image_data_2D,
                                                                                    combine_method='finest_valid',
                                                                                    sigma_sigma_clip=2.0, maxiters=10,
                                                                                    beam_size_px=int(psf_fwhm+1), 
                                                                                    # sigma_sigma_clip=3.0, maxiters=15,
                                                                                    # beam_size_px=3, n_scales=8,
                                                                                    # do_plot=False,
                                                                                    # profile_type='azimuthal'
                                                                                    )
                rms_map = mlibs.estimate_RMS_map(input_data.image_data_2D, 
                                                filter_size=7, sigma=12.0,
                                                exclude_percentile=50,
                                                box_size=(7,7))
                fix_geometry = [True] * SE.n_IDs + [True,True]
                dr_fix = [5] * SE.n_IDs + [5,5]
                fix_x0_y0 = [True] * SE.n_IDs + [True,True]
                fix_value_n = [1.0] * SE.n_IDs + [1.0,1.0]
                fix_n = [False] * SE.n_IDs + [False,False]
                # reload_libs()
                smfg = sersic_multifit_general(input_data,
                                                SE, #source extraction object, from previous step
                                                convolution_mode='GPU',
                                                which_residual = 'user',
                                                mask=mask_region,
                                                mask_for_fit=mask_for_fit,
                                                Npsf=0.5,
                                                use_mask_for_fit=True,loss='cauchy',
                                                # bkg_map = mlibs.shuffle_2D(final_bkg), is_bkg_map_conv = False, <<UNSTABLE - DO NOT USE>>
                                                bkg_map = final_bkg, is_bkg_map_conv = True, 
                                                # bkg_map = np.ones_like(final_bkg)*1e-3, is_bkg_map_conv = True, ##TESTING <DO NOT USE>>
                                                # use_weights=True, 
                                                # rms_map = rms_map,
                                                tr_solver='exact',
                                                fix_geometry=fix_geometry, #for stability purposes, keep True for now. 
                                                comp_ids=comp_ids,# which component label is compact/bulge?
                                                dr_fix=dr_fix,#for each component, radial element size to fix (x0,y0) positions
                                                fix_x0_y0=fix_x0_y0,#for each component, fix or not the (x0,y0) positions.
                                                fix_value_n=fix_value_n,#for each component, the Sersic index value to be fixed. 
                                                fix_n=fix_n,#for each component, fix or not the Sersic index. 
                                                z = 0.1, # redshift; just an arbitrary value for now.
                                                )
                # mlibs.lmfit.report_fit(smfg.mini.result.params)
                sys.exit()
            else:
                print("<< In development >>")
                print("Not implemented yet.")
                sys.exit("Error: Only 'other' obs_type is allowed in general_decomp mode at the momment. ")


        else:
            input_data = read_data(filename=args.filename,
                                residualname=args.residualname,)
            if "--image_stats" in sys.argv:
                radio_image_analysis(input_data, z=args.redshift)
            if "--find_sources" in sys.argv:
                SE = source_extraction(input_data,sigma=args.sigma,
                                ell_size_factor=args.ell_size_factor,obs_type=args.obs_type)
                if "--sersic_radio" in sys.argv:
                    if args.residualname != None:
                        SMFR = sersic_multifit_radio(input_data, SE,
                                                    method2 = args.solver2)
                        # SMFR_decomp = decompose_emission(input_data,
                        #           SMFR,#seric multifit results object
                        #           z=SMFR.z)
                        if '--SFR-do' in sys.argv:
                            SFR = radio_star_formation(input_data,
                                                       SMFR, z=args.redshift)
                    else:
                        print("Error: Please, provide a residual image (e.g. the one "
                            "generate during interferometric deconvolution).")

                if "--general_fit" in sys.argv:
                    if args.psfname != None:
                        SMFR = sersic_multifit_general(input_data, SE,
                                                    method2=args.solver2)

                    else:
                        print("Error: Please, provide a psf image.")
