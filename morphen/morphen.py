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
__dates__ =  ('2024 03 25','2024 11 13')
__version__ = '0.5.0alpha-1'
__codename__ = 'Saurinho'
__author__ = 'Geferson Lucatelli'
__coauthors__ = ('Javier Moldon, Rob Beswick, '
                  'Fabricio Ferrari, Leonardo Ferreira')
__email__ = 'geferson.lucatelli@postgrad.manchester.ac.uk'
__date__ = '2024 12 18'
# print(__doc__)


import argparse
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

class config():
    """
    Configuration Class to specify basic parameters.
    """

    def reset_rc_params():
        """
        Global configuration for matplotlib.pyplot
        """
        mpl.rcParams.update({'font.size': 16,
                             'text.usetex': False,  #
                             'font.family': 'sans-serif',
                             'mathtext.fontset': 'stix',
                             'font.family': 'sans',
                             'font.weight': 'medium',  
                             'font.family': 'STIXGeneral',
                             'xtick.labelsize': 16,
                             'figure.figsize': (6, 4),
                             'ytick.labelsize': 16,
                             'axes.labelsize': 16,
                             'xtick.major.width': 1,
                             'ytick.major.width': 1,
                             'axes.linewidth': 1.5,
                             'axes.edgecolor':'orange',
                             'lines.linewidth': 2,
                             'legend.fontsize': 14,
                             'grid.linestyle': '--',
                             # 'grid.color':'black',
                             'axes.grid.which': 'major',  
                             'axes.grid.axis': 'both', 
                             'axes.spines.right': True,
                             'axes.grid': True,
                             'axes.titlesize' : 16,
                             'legend.framealpha': 1.0
                             })
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
                 imagelist = [],residuallist=[],):

        self.filename = filename
        self.residualname = residualname
        self.psfname = psfname
        self.print_names()
        self.get_data()
        try:
            self.get_info()
        except:
            pass

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

    def get_data(self):
        self.image_data_2D = None
        self.residual_data_2D = None
        self.psf_data_2D = None
        self.rms_img = None
        self.rms_res = None
        if self.filename != None:
            self.image_data_2D = mlibs.load_fits_data(self.filename)
            self.rms_img = mlibs.mad_std(self.image_data_2D)
        if self.residualname != None:
            self.residual_data_2D = mlibs.load_fits_data(self.residualname)
            self.rms_res = mlibs.mad_std(self.residual_data_2D)
        if self.psfname != None:
            self.psf_data_2D = mlibs.load_fits_data(self.psfname)
    
    def get_info(self):
        self.cell_size = mlibs.get_cell_size(self.filename)
        self.beam_area_px = mlibs.beam_area2(self.filename)
        


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
        # _logging_.logger.info("Computing image properties.")
        # self.levels, self.fluxes, self.agrow, self.plt_image, \
        #     self.omask, self.mask, self.results_im_props = \
        #     mlibs.compute_image_properties(img=self.input_data.filename,
        #                                    cell_size=self.cell_size,
        #                                    residual=self.input_data.residualname,
        #                                    sigma_mask=self.sigma_mask,
        #                                    dilation_size=self.dilation_size,
        #                                    crop=self.crop,
        #                                    iterations=config.mask_iterations,
        #                                    box_size=self.box_size,
        #                                    last_level=self.last_level,
        #                                    mask=self.mask,
        #                                    apply_mask=self.apply_mask,
        #                                    vmin_factor=self.vmin_factor,
        #                                    results=self.results,
        #                                    show_figure=self.show_figure,
        #                                    logger=_logging_.logger)
        #
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

        self.image_measures, self.mask, self.omask = \
            mlibs.measures(imagename=self.input_data.filename,
                           residualname=self.input_data.residualname,
                           z=self.z,
                           mask_component=self.mask_component,
                           sigma_mask=self.sigma_mask,
                           last_level=self.last_level,
                           vmin_factor=self.vmin_factor,
                           plot_catalog=True, 
                           data_2D=self.input_data.image_data_2D,
                           npixels=self.npixels,fwhm=self.fwhm,
                           kernel_size=self.kernel_size,
                           dilation_size=self.dilation_size,
                           main_feature_index=0,
                           results_final={},
                           crop=self.crop, box_size=self.box_size,
                           iterations=config.mask_iterations,
                           fracX=0.15, fracY=0.15,
                           deblend=False, bkg_sub=False,
                           bkg_to_sub=None, rms=None,
                           do_petro=self.do_petro,
                           apply_mask=self.apply_mask,
                           do_PLOT=True, SAVE=self.SAVE,
                           show_figure=self.show_figure,
                           mask=self.mask,
                           do_measurements='all',
                           compute_A=True,
                           add_save_name='',logger=_logging_.logger)


class source_extraction():
    """
    Source extraction class, responsible to find relevant regions of emission.

    For now, only SEP is implemented. Soon, other algorithms will be added
    in this class as alternatives for source extraction. These are:
        - PyBDSF
        - AstroDendro
        - Photutils
    """
    def __init__(self, input_data,z=0.05,ids_to_add=[1],
                 crop=False, box_size=256,
                 apply_mask=False, mask=None, dilation_size = None,
                 sigma_level=6, sigma_mask=6, vmin_factor=3, mask_component=None,
                 bwf=1, bhf=1, fwf=1, fhf=1,
                 segmentation_map = True, filter_type='conv',
                 deblend_nthresh=25, deblend_cont=1e-8,
                 clean_param=0.5, clean=True,
                 minarea_factor = 1.0,npixels=None,
                 sort_by='flux',  # sort detected source by distance
                 sigma=6,  # min rms to search for sources
                 mask_grow_iterations = 1,
                 ell_size_factor=None,  # unstable, please inspect!
                 obs_type = 'radio', algorithm='SEP',threshold_mode='sigma',
                 force_circular=False,
                 show_detection=False,show_petro_plots=False,
                 show_bkg_map=False,
                 SAVE=True, show_figure=True,dry_run = False,SE_ref=None):
        """
        Parameters
        ----------
        input_data : str
            Path to the image.
        z : float
            Redshift of the source.
        ids_to_add : list
            List of ids to be added to the catalogue.
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
        """

        self.input_data = input_data
        self.z = z
        self.ids_to_add = ids_to_add
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
        self.SAVE = SAVE
        self.show_figure = show_figure
        self.show_detection = show_detection
        self.show_bkg_map = show_bkg_map
        self.show_petro_plots = show_petro_plots
        self.obs_type = obs_type
        self.force_circular = force_circular
        self.mask_grow_iterations = mask_grow_iterations
        self.threshold_mode = threshold_mode
        
        if dry_run is True:
            self.show_detection = True
            self.get_sources()

        if dry_run is not True:
            self.contruct_source_properties()

    def get_sources(self):
        if self.algorithm == 'SEP':
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
            (self.masks, self.indices, self.bkg,
             self.seg_maps, self.objects) = \
                mlibs.phot_source_ext(self.input_data.filename,
                                      residual=self.input_data.residual_data_2D,
                               bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                               # filtering options for source detection
                               segmentation_map=self.segmentation_map,
                               filter_type=self.filter_type, mask=self.mask,
                               deblend_nthresh=self.deblend_nthresh,
                               deblend_cont=self.deblend_cont,
                               clean_param=self.clean_param,
                               clean=self.clean,
                               sort_by=self.sort_by,
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

    def contruct_source_properties(self):
        (self.sources_photometries, self.n_components, self.n_IDs, 
         self.masks, self.indices, self.objects,
         self.psf_name, self.mask, self.bkg) = \
            mlibs.prepare_fit(self.input_data.filename,
                              self.input_data.residualname,
                              z=self.z,ids_to_add = self.ids_to_add,
                              bw=self.bw, bh=self.bh, fw=self.fw, fh=self.fh,
                              sigma=self.sigma, sigma_mask=self.sigma_mask,
                              apply_mask=self.apply_mask,mask=self.mask,
                              mask_grow_iterations = self.mask_grow_iterations,
                              deblend_nthresh=self.deblend_nthresh,
                              ell_size_factor=self.ell_size_factor,
                              minarea=self.minarea,npixels = self.npixels,
                              minarea_factor=self.minarea_factor,
                              deblend_cont=self.deblend_cont,
                              clean_param=self.clean_param,
                              obs_type=self.obs_type,algorithm=self.algorithm,
                              sort_by=self.sort_by,
                              force_circular = self.force_circular,
                              show_petro_plots=self.show_petro_plots)

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
                 fix_geometry = True,
                 comp_ids = [],
                 fix_n = None, 
                 fix_value_n = None, 
                 fix_max_value_Rn = None,
                 fix_min_value_Rn = False,
                 dr_fix = None,fix_x0_y0=None,
                 sigma=6.0, use_mask_for_fit=False,mask_for_fit=None,
                 mask=None,
                 tr_solver = "exact",loss = 'cauchy',
                 convolution_mode='GPU',method1='least_squares',
                 self_bkg = False, bkg_rms_map = None,
                 is_rms_map_conv = False,
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
        fix_geometry : bool
            If True, fix the geometry of the components.
        comp_ids : list
            List of component IDs to be fitted.
        fix_n : list
            List of booleans to fix the Sersic index.
        fix_value_n : list
            List of values to fix the Sersic index.
        dr_fix : list
            List of values to fix the dr parameter.
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

            Examples: 
                fix_n = [True, True, True,True, True, True,True, True, True]
                fix_value_n=[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
                dr_fix = [3, 3, 100, 100, 10, 10, 10, 10, 10]
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
        self.comp_ids = comp_ids
        self.fix_geometry = fix_geometry
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
        self.bkg_rms_map = bkg_rms_map
        self.is_rms_map_conv = is_rms_map_conv
        self.verbose = verbose

        self.fix_min_value_Rn = fix_min_value_Rn
        
        if fix_n is None:
            self.fix_n = [True] * self.SE.n_components
        else:
            self.fix_n = fix_n
            
        if fix_value_n is None:
            self.fix_value_n = [0.5] * self.SE.n_components
        else:
            self.fix_value_n = fix_value_n
            
        if fix_max_value_Rn is None:
            self.fix_max_value_Rn = [False] * self.SE.n_components
        else:
            self.fix_max_value_Rn = fix_max_value_Rn
            
        if dr_fix is None:
            self.dr_fix = [10] * self.SE.n_components
        else:
            self.dr_fix = dr_fix
        if fix_x0_y0 is None:
            self.fix_x0_y0 = [True] * self.SE.n_components

        if fix_geometry is None:
            self.fix_geometry = [True] * self.SE.n_components
        else:
            self.fix_geometry = fix_geometry
        
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
         self.class_resuts, self.compact_model) = \
            mlibs.run_image_fitting(imagelist=[self.input_data.filename],
                                    residuallist=[self.input_data.residualname],
                                    aspect=self.aspect,
                                    which_residual=self.which_residual,
                                    comp_ids=self.comp_ids,# which IDs refers to compact components?
                                    sources_photometries=self.SE.sources_photometries,
                                    n_components=self.SE.n_components,
                                    z=self.z,
                                    convolution_mode=self.convolution_mode,
                                    method1=self.method1,
                                    method2=self.method2,
                                    mask=self.mask,
                                    use_mask_for_fit=self.use_mask_for_fit,
                                    mask_for_fit=self.mask_for_fit,
                                    bkg_rms_map=self.bkg_rms_map,
                                    self_bkg=self.self_bkg,
                                    is_rms_map_conv = self.is_rms_map_conv,
                                    save_name_append='',
                                    fix_n=self.fix_n,
                                    loss=self.loss,
                                    tr_solver = self.tr_solver,
                                    fix_value_n=self.fix_value_n,
                                    fix_max_value_Rn = self.fix_max_value_Rn,
                                    fix_min_value_Rn = self.fix_min_value_Rn,
                                    fix_geometry=self.fix_geometry,  # unstable if  False
                                    dr_fix=self.dr_fix,
                                    fix_x0_y0 = self.fix_x0_y0,
                                    parameters_mini_init = self.parameters_mini_init,
                                    sigma=self.sigma,
                                    logger=_logging_.logger,verbose=self.verbose)
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
                 fix_geometry = True,
                 comp_ids = ['1'],
                 fix_n = None,
                 fix_value_n = None, 
                 fix_max_value_Rn = None,
                 dr_fix = None,
                 constrained=True, self_bkg=False,
                 sigma=6.0, 
                 use_mask_for_fit=False,mask_for_fit=None,mask=None,
                 bkg_rms_map = None,
                 loss='cauchy', tr_solver='exact',
                 regularize=True, f_scale=1.0, ftol=1e-10,
                 xtol=1e-10, gtol=1e-10,
                 init_params=0.2, final_params=5.0,
                 which_residual = 'user',
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
            
        self.comp_ids = comp_ids
        self.fix_geometry = fix_geometry
        self.convolution_mode = convolution_mode
        self.constrained = constrained
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
        self.bkg_rms_map = bkg_rms_map

        if self.bkg_rms_map is not None:
            self.rms_map = self.bkg_rms_map
        else:
            if self.self_bkg == True:
                self.rms_map = self.SE.bkg
            else:
                self.rms_map = None

        if fix_n is None:
            self.fix_n = [False] * self.SE.n_components
        else:
            self.fix_n = fix_n

        if fix_value_n is None:
            self.fix_value_n = [1.5] * self.SE.n_components
        else:
            self.fix_value_n = fix_value_n
            
        if fix_max_value_Rn is None:
            self.fix_max_value_Rn = [False] * self.SE.n_components
        else:
            self.fix_max_value_Rn = fix_max_value_Rn

        if dr_fix is None:
            self.dr_fix = [10] * self.SE.n_components
        else:
            self.dr_fix = dr_fix

        if fix_geometry is None:
            self.fix_geometry = [True] * self.SE.n_components
        else:
            self.fix_geometry = fix_geometry

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
                           rms_map = self.rms_map,
                           which_residual=self.which_residual,
                           observation_type = self.SE.obs_type,
                           mask_region = self.mask_for_fit,
                           # rms_map=self.SE.bkg,
                           # rms_map=None,
                           fix_n=self.fix_n,
                           fix_value_n=self.fix_value_n,
                           fix_max_value_Rn=self.fix_max_value_Rn,
                           dr_fix=self.dr_fix,
                           convolution_mode=self.convolution_mode,
                           fix_geometry=self.fix_geometry, workers=-1,
                           method1=self.method1,
                           method2=self.method2,
                           init_params=self.init_params,
                           final_params=self.final_params,
                           loss=self.loss, tr_solver=self.tr_solver,
                           regularize=self.regularize, f_scale=self.f_scale,
                           ftol=self.ftol,
                           xtol=self.xtol, gtol=self.gtol,
                           save_name_append=self.save_name_append,
                           logger=_logging_.logger)

        all_comps_ids = np.arange(1, self.SE.n_components + 1).astype('str')
        mask_compact_ids = np.isin(all_comps_ids, np.asarray(self.comp_ids))
        ext_ids = list(all_comps_ids[~mask_compact_ids])

        special_name = ''
        compact_model = 0
        extended_model = 0
        compact_model_deconv = 0
        extended_model_deconv = 0

        if self.input_data.rms_res == None:
            rms_std_res = self.input_data.rms_img
        else:
            rms_std_res = self.input_data.rms_res

        for lc in self.comp_ids:
            compact_model = (compact_model +
                                self.model_dict['model_c' + lc + '_conv'])
            compact_model_deconv = (compact_model_deconv +
                                    self.model_dict['model_c' + lc])
        # if ext_ids is not None:
        if ext_ids == []:
            extended_model = 0
            extended_model_deconv = 0
            nfunctions = 1
        else:
            for le in ext_ids:
                extended_model = (extended_model +
                                    self.model_dict['model_c' + le + '_conv'])
                extended_model_deconv = (extended_model_deconv +
                                            self.model_dict['model_c' + le])
                nfunctions = None

        self.decomposition_results = mlibs.plot_decomp_results(imagename=self.input_data.filename,
                                                   compact=compact_model,
                                                   extended_model=extended_model,
                                                   rms=rms_std_res,
                                                   nfunctions=nfunctions,
                                                   obs_type=self.SE.obs_type,
                                                   special_name=special_name)

        mlibs.plot_fit_results(self.input_data.filename, 
                               self.model_dict,
                               self.image_results_conv,
                               self.SE.sources_photometries,
                               bkg_image=self.bkg_images[-1],
                               crop=False, box_size=200,
                               mask=self.mask_for_fit,
                               obs_type=self.SE.obs_type,
                               vmax_factor=0.3, vmin_factor=1.0)

        mlibs.plot_slices(data_2D=self.input_data.image_data_2D,
                          model_dict=self.model_dict,
                          image_results_conv=self.image_results_conv[-2],
                          Rp_props=self.SE.sources_photometries,
                          residual_2D=None)
        self.parameter_results = self.result_mini.params.valuesdict().copy()
        try:
            for param in self.result_mini.params.valuesdict().keys():
                self.parameter_results[param+'_err'] = self.result_mini.params[param].stderr
        except:
            pass
        self.parameter_results['#imagename'] = os.path.basename(self.input_data.filename)
        self.results_fit = {**self.parameter_results, **self.decomposition_results}
        # self.results_fit.append(all_results)
        pass

class morphometry():
    """
    Core functionalities from Morfometryka. 
    Morfometryka is not publically available yet, 
    so these functions will be added in a later stage. 
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
                                f" +/- {1000*self.SMFR.data_properties['flux_error_res_3'].values[0]:.2f} mJy")
        
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
            mlibs.print_logger_header(title="Decomposition Results",
                    logger=_logging_.logger)
            

            
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            _logging_.logger.info(f'Compact component ID {component_id}.')
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            
            # results from sersic fitting
            compact_dec_summary[f'f{component_id}_Rn'] = \
                self.SMFR.result_mini.params[f'f{component_id}_Rn'].value * self.pix_to_pc
            compact_dec_summary[f'f{component_id}_Rn_err'] = \
                self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr * self.pix_to_pc
            A_Rn,A_Rn_err = mlibs.radii_to_area(self.SMFR.result_mini.params[f'f{component_id}_Rn'].value,
                                                self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr)
            compact_dec_summary[f'f{component_id}_A_Rn'] = \
                mlibs.pix_area_to_kpc_area(A_Rn,self.pix_to_pc)
            compact_dec_summary[f'f{component_id}_A_Rn_err'] = \
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
            
            
            peak_snr = self.SMFR.components_conv_props['peak_of_flux'].loc[f'model_c_conv_{component_id}_props'] \
                / self.SMFR.components_conv_props['rms_residual'].loc[f'model_c_conv_{component_id}_props']
                
            # theta_Rn_maj_asec_err = self.SMFR.components_conv_props['bmajor'].iloc[0]/(2*peak_snr*theta_Rn_maj_asec)
            # theta_Rn_min_asec_err = self.SMFR.components_conv_props['bminor'].iloc[0]/(2*peak_snr*theta_Rn_min_asec)
            theta_Rn_maj_asec_err = self.SMFR.components_conv_props['bmajor'].iloc[0]/(2*peak_snr)
            theta_Rn_min_asec_err = self.SMFR.components_conv_props['bminor'].iloc[0]/(2*peak_snr)
            theta_Rn_asec_err = np.sqrt(theta_Rn_maj_asec_err**2 + theta_Rn_min_asec_err**2)

            
            Snu = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['total_flux_mask']
            Snu_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['flux_error_res_3']
            compact_dec_summary[f'comp_{component_id}_Snu'] = Snu
            compact_dec_summary[f'comp_{component_id}_Snu_err'] = Snu_err
            
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

            
            compact_dec_summary[f'f{component_id}_theta_Rn_maj_pc'] = theta_Rn_maj_pc
            compact_dec_summary[f'f{component_id}_theta_Rn_min_pc'] = theta_Rn_min_pc
            compact_dec_summary[f'f{component_id}_theta_Rn_pc'] = theta_Rn_pc
            compact_dec_summary[f'f{component_id}_theta_Rn_maj_pc_err'] = theta_Rn_maj_pc_err 
            compact_dec_summary[f'f{component_id}_theta_Rn_min_pc_err'] = theta_Rn_min_pc_err 
            compact_dec_summary[f'f{component_id}_theta_Rn_pc_err'] = theta_Rn_pc_err

            compact_dec_summary[f'f{component_id}_theta_Rn_maj_asec'] = theta_Rn_maj_asec
            compact_dec_summary[f'f{component_id}_theta_Rn_min_asec'] = theta_Rn_min_asec
            compact_dec_summary[f'f{component_id}_theta_Rn_asec'] = theta_Rn_asec
            compact_dec_summary[f'f{component_id}_theta_Rn_asec_err'] = theta_Rn_asec_err
            compact_dec_summary[f'f{component_id}_theta_Rn_maj_asec_err'] = theta_Rn_maj_asec_err
            compact_dec_summary[f'f{component_id}_theta_Rn_min_asec_err'] = theta_Rn_min_asec_err
            
            compact_dec_summary[f'f{component_id}_theta_R50_pc'] = theta_R50_pc
            compact_dec_summary[f'f{component_id}_theta_R50_maj_pc'] = theta_R50_maj_pc
            compact_dec_summary[f'f{component_id}_theta_R50_min_pc'] = theta_R50_min_pc
            compact_dec_summary[f'f{component_id}_theta_R50_pc_err'] = theta_R50_pc_err
            compact_dec_summary[f'f{component_id}_theta_R50_maj_pc_err'] = theta_R50_maj_pc_err
            compact_dec_summary[f'f{component_id}_theta_R50_min_pc_err'] = theta_R50_min_pc_err
            
            compact_dec_summary[f'f{component_id}_theta_R50_asec'] = theta_R50_asec
            compact_dec_summary[f'f{component_id}_theta_R50_maj_asec'] = theta_R50_maj_asec
            compact_dec_summary[f'f{component_id}_theta_R50_min_asec'] = theta_R50_min_asec
            compact_dec_summary[f'f{component_id}_theta_R50_asec_err'] = theta_R50_asec_err
            compact_dec_summary[f'f{component_id}_theta_R50_maj_asec_err'] = theta_R50_maj_asec_err
            compact_dec_summary[f'f{component_id}_theta_R50_min_asec_err'] = theta_R50_min_asec_err
            
            
            
            # compact_dec_summary[f'f{component_id}_theta_Rn_maj_pc_err'] = theta_Rn_maj_asec_err * self.pix_to_pc/self.cell_size
            # compact_dec_summary[f'f{component_id}_theta_Rn_min_pc_err'] = theta_Rn_min_asec_err * self.pix_to_pc/self.cell_size

            compact_dec_summary[f'f{component_id}_Tb_Rn'] = TB_Rn
            compact_dec_summary[f'f{component_id}_Tb_Rn_err'] = TB_Rn_err
            compact_dec_summary[f'f{component_id}_Tb_R50'] = TB_R50
            compact_dec_summary[f'f{component_id}_Tb_R50_err'] = TB_R50_err
            
        
            #properties from model images
            compact_dec_summary[f'comp_{component_id}_R50_d'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii'] * self.pix_to_pc
            compact_dec_summary[f'comp_{component_id}_R50_d_err'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii_err'] * self.pix_to_pc
            compact_dec_summary[f'comp_{component_id}_R95_d'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C95radii'] * self.pix_to_pc
            compact_dec_summary[f'comp_{component_id}_R95_d_err'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C95radii_err'] * self.pix_to_pc
            


            A_R50_d = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix50']
            A_R50_d_err = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix50_err']
            
            compact_dec_summary[f'comp_{component_id}_A_R50_d'] = mlibs.pix_area_to_kpc_area(A_R50_d,self.pix_to_pc)
            compact_dec_summary[f'comp_{component_id}_A_R50_d_err'] = mlibs.pix_area_to_kpc_area(A_R50_d_err,self.pix_to_pc)
                
            A_R95_d = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix95']
            A_R95_d_err = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix95_err']

            compact_dec_summary[f'comp_{component_id}_A_R95_d'] = mlibs.pix_area_to_kpc_area(A_R95_d,self.pix_to_pc)
            compact_dec_summary[f'comp_{component_id}_A_R95_d_err'] = mlibs.pix_area_to_kpc_area(A_R95_d_err,self.pix_to_pc)


            compact_dec_summary[f'comp_{component_id}_s_Snu_A50_d'] = 0.5*Snu / compact_dec_summary[f'comp_{component_id}_A_R50_d']
            compact_dec_summary[f'comp_{component_id}_s_Snu_A50_d_err'] = 0.5*Snu_err / compact_dec_summary[f'comp_{component_id}_A_R50_d']
            compact_dec_summary[f'comp_{component_id}_s_Snu_A95_d'] = Snu / compact_dec_summary[f'comp_{component_id}_A_R95_d']
            compact_dec_summary[f'comp_{component_id}_s_Snu_A95_d_err'] = Snu_err / compact_dec_summary[f'comp_{component_id}_A_R95_d']
            compact_dec_summary[f'comp_{component_id}_s_Snu_A_Rn'] = 0.5*Snu / compact_dec_summary[f'f{component_id}_A_Rn']
            compact_dec_summary[f'comp_{component_id}_s_Snu_A_Rn_err'] = 0.5*Snu_err / compact_dec_summary[f'f{component_id}_A_Rn']

            #CONVOLVED QUANTITIES
            compact_conv_summary[f'comp_{component_id}_R50'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C50radii'] * self.pix_to_pc
            compact_conv_summary[f'comp_{component_id}_R50_err'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C50radii_err'] * self.pix_to_pc
            compact_conv_summary[f'comp_{component_id}_R95'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C95radii'] * self.pix_to_pc
            compact_conv_summary[f'comp_{component_id}_R95_err'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C95radii_err'] * self.pix_to_pc



            A_R50 = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix50']
            A_R50_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix50_err']
            
            compact_conv_summary[f'comp_{component_id}_A_R50'] = mlibs.pix_area_to_kpc_area(A_R50,self.pix_to_pc)
            compact_conv_summary[f'comp_{component_id}_A_R50_err'] = mlibs.pix_area_to_kpc_area(A_R50_err,self.pix_to_pc)
                
            A_R95 = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix95']
            A_R95_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix95_err']

            compact_conv_summary[f'comp_{component_id}_A_R95'] = mlibs.pix_area_to_kpc_area(A_R95,self.pix_to_pc)
            compact_conv_summary[f'comp_{component_id}_A_R95_err'] = mlibs.pix_area_to_kpc_area(A_R95_err,self.pix_to_pc)
            compact_conv_summary[f'comp_{component_id}_s_Snu_A95'] = Snu / compact_conv_summary[f'comp_{component_id}_A_R95']
            compact_conv_summary[f'comp_{component_id}_s_Snu_A95_err'] = Snu_err / compact_conv_summary[f'comp_{component_id}_A_R95']
            
            _logging_.logger.info(f"Snu = {1000*compact_dec_summary[f'comp_{component_id}_Snu']:.2f}" 
                                  f" +/- {1000*compact_dec_summary[f'comp_{component_id}_Snu_err']:.2f} mJy")
            _logging_.logger.info(f"Rn = {compact_dec_summary[f'f{component_id}_Rn']:.2f}" 
                                  f" +/- {compact_dec_summary[f'f{component_id}_Rn_err']:.2f} pc")
            _logging_.logger.info(f"R50_d = {compact_dec_summary[f'comp_{component_id}_R50_d']:.2f}"
                                  f" +/- {compact_dec_summary[f'comp_{component_id}_R50_d_err']:.2f} pc")
            # _logging_.logger.info(f"R95_d = {compact_dec_summary[f'comp_{component_id}_R95_d']:.2f} pc")
            _logging_.logger.info(f"R95_d = {compact_dec_summary[f'comp_{component_id}_R95_d']:.2f}"
                                  f" +/- {compact_dec_summary[f'comp_{component_id}_R95_d_err']:.2f} pc")
            
            _logging_.logger.info(f"R50 = {compact_conv_summary[f'comp_{component_id}_R50']:.2f}"
                                  f" +/- {compact_conv_summary[f'comp_{component_id}_R50_err']:.2f} pc")
            _logging_.logger.info(f"R95 = {compact_conv_summary[f'comp_{component_id}_R95']:.2f}"
                                  f" +/- {compact_conv_summary[f'comp_{component_id}_R95_err']:.2f} pc")            

            _logging_.logger.info(f"Theta_Rn = {compact_dec_summary[f'f{component_id}_theta_Rn_pc']:.2f}" 
                                  f" +/- {compact_dec_summary[f'f{component_id}_theta_Rn_pc_err']:.2f} pc")
            _logging_.logger.info(f"Theta_R50 = {compact_dec_summary[f'f{component_id}_theta_R50_pc']:.2f}"
                                    f" +/- {compact_dec_summary[f'f{component_id}_theta_R50_pc_err']:.2f} pc")
            _logging_.logger.info(f"Tb_Rn = ({compact_dec_summary[f'f{component_id}_Tb_Rn']:.2f}" 
                                  f" +/- {compact_dec_summary[f'f{component_id}_Tb_Rn_err']:.2f}) x 10^5 K")
            _logging_.logger.info(f"Tb_R50 = ({compact_dec_summary[f'f{component_id}_Tb_R50']:.2f}"
                                    f" +/- {compact_dec_summary[f'f{component_id}_Tb_R50_err']:.2f}) x 10^5 K")


            _logging_.logger.info(f"A_Rn = ({1e3*compact_dec_summary[f'f{component_id}_A_Rn']:.3f}" 
                                  f" +/- {1e3*compact_dec_summary[f'f{component_id}_A_Rn_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"A50_d = ({1e3*compact_dec_summary[f'comp_{component_id}_A_R50_d']:.3f}"
                                  f" +/- {1e3*compact_dec_summary[f'comp_{component_id}_A_R50_d_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"A95_d = ({1e3*compact_dec_summary[f'comp_{component_id}_A_R95_d']:.3f}"
                                  f" +/- {1e3*compact_dec_summary[f'comp_{component_id}_A_R95_d_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"A50 = ({1e3*compact_conv_summary[f'comp_{component_id}_A_R50']:.3f}"
                                  f" +/- {1e3*compact_conv_summary[f'comp_{component_id}_A_R50_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"A95 = ({1e3*compact_conv_summary[f'comp_{component_id}_A_R95']:.3f}"
                                  f" +/- {1e3*compact_conv_summary[f'comp_{component_id}_A_R95_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"sSnu_A95_d = {compact_dec_summary[f'comp_{component_id}_s_Snu_A95_d']:.3f}"
                                  f" +/- {compact_dec_summary[f'comp_{component_id}_s_Snu_A95_d_err']:.3f} Jy/kpc^2")
            _logging_.logger.info(f"sSnu_A95 = {compact_conv_summary[f'comp_{component_id}_s_Snu_A95']:.3f}"
                                  f" +/- {compact_conv_summary[f'comp_{component_id}_s_Snu_A95_err']:.3f} Jy/kpc^2")

            # _logging_.logger.info(f"Rn =  {compact_dec_summary[f'f{component_id}_Rn']}")




        compact_conv_summary[f'Snu_comp_total'] = self.SMFR.results_compact_conv_morpho['total_flux_mask'].values[0]
        compact_conv_summary[f'Snu_comp_total_err'] = self.SMFR.results_compact_conv_morpho['flux_error_res_3'].values[0]
        _logging_.logger.info(f"Snu Compact Total = {1000*compact_conv_summary[f'Snu_comp_total']:.2f}" 
                                f" +/- {1000*compact_conv_summary[f'Snu_comp_total_err']:.2f} mJy")
         
        self.compact_dec_summary = compact_dec_summary
        self.compact_conv_summary = compact_conv_summary
        
        pass
    
    def compute_extended_properties(self):
        all_comps_ids = np.arange(1, self.SMFR.SE.n_components + 1).astype('str')
        mask_compact_ids = np.isin(all_comps_ids, np.asarray(self.SMFR.comp_ids))
        self.SMFR.ext_ids = list(all_comps_ids[~mask_compact_ids])
        
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
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            _logging_.logger.info(f'Extended component ID {component_id}.')
            _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
            # results from sersic fitting
            ext_dec_summary[f'f{component_id}_Rn'] = \
                self.SMFR.result_mini.params[f'f{component_id}_Rn'].value * self.pix_to_pc
            ext_dec_summary[f'f{component_id}_Rn_err'] = \
                self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr * self.pix_to_pc
            A_Rn,A_Rn_err = mlibs.radii_to_area(self.SMFR.result_mini.params[f'f{component_id}_Rn'].value,
                                                self.SMFR.result_mini.params[f'f{component_id}_Rn'].stderr)
            ext_dec_summary[f'f{component_id}_A_Rn'] = \
                mlibs.pix_area_to_kpc_area(A_Rn,self.pix_to_pc)
            ext_dec_summary[f'f{component_id}_A_Rn_err'] = \
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
            
            peak_snr = self.SMFR.components_conv_props['peak_of_flux'].loc[f'model_c_conv_{component_id}_props'] \
                / self.SMFR.components_conv_props['rms_residual'].loc[f'model_c_conv_{component_id}_props']
                
            # theta_Rn_maj_asec_err = self.SMFR.components_conv_props['bmajor'].iloc[0]/(2*peak_snr*theta_Rn_maj_asec)
            # theta_Rn_min_asec_err = self.SMFR.components_conv_props['bminor'].iloc[0]/(2*peak_snr*theta_Rn_min_asec)
            theta_Rn_maj_asec_err = self.SMFR.components_conv_props['bmajor'].iloc[0]/(2*peak_snr)
            theta_Rn_min_asec_err = self.SMFR.components_conv_props['bminor'].iloc[0]/(2*peak_snr)
            theta_Rn_asec_err = np.sqrt(theta_Rn_maj_asec_err**2 + theta_Rn_min_asec_err**2)
            
            Snu = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['total_flux_mask']
            Snu_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['flux_error_res_3']
            ext_dec_summary[f'comp_{component_id}_Snu'] = Snu
            ext_dec_summary[f'comp_{component_id}_Snu_err'] = Snu_err
            
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
            
            ext_dec_summary[f'f{component_id}_theta_Rn_maj_pc'] = theta_Rn_maj_pc
            ext_dec_summary[f'f{component_id}_theta_Rn_min_pc'] = theta_Rn_min_pc
            ext_dec_summary[f'f{component_id}_theta_Rn_pc'] = theta_Rn_pc
            ext_dec_summary[f'f{component_id}_theta_Rn_maj_pc_err'] = theta_Rn_maj_pc_err 
            ext_dec_summary[f'f{component_id}_theta_Rn_min_pc_err'] = theta_Rn_min_pc_err 
            ext_dec_summary[f'f{component_id}_theta_Rn_pc_err'] = theta_Rn_pc_err
            
            ext_dec_summary[f'f{component_id}_theta_Rn_maj_asec'] = theta_Rn_maj_asec
            ext_dec_summary[f'f{component_id}_theta_Rn_min_asec'] = theta_Rn_min_asec
            ext_dec_summary[f'f{component_id}_theta_Rn_asec'] = theta_Rn_asec
            ext_dec_summary[f'f{component_id}_theta_Rn_asec_err'] = theta_Rn_asec_err
            ext_dec_summary[f'f{component_id}_theta_Rn_maj_asec_err'] = theta_Rn_maj_asec_err
            ext_dec_summary[f'f{component_id}_theta_Rn_min_asec_err'] = theta_Rn_min_asec_err
            
            ext_dec_summary[f'f{component_id}_theta_R50_pc'] = theta_R50_pc
            ext_dec_summary[f'f{component_id}_theta_R50_maj_pc'] = theta_R50_maj_pc
            ext_dec_summary[f'f{component_id}_theta_R50_min_pc'] = theta_R50_min_pc
            ext_dec_summary[f'f{component_id}_theta_R50_pc_err'] = theta_R50_pc_err
            ext_dec_summary[f'f{component_id}_theta_R50_maj_pc_err'] = theta_R50_maj_pc_err
            ext_dec_summary[f'f{component_id}_theta_R50_min_pc_err'] = theta_R50_min_pc_err
            
            ext_dec_summary[f'f{component_id}_theta_R50_asec'] = theta_R50_asec
            ext_dec_summary[f'f{component_id}_theta_R50_maj_asec'] = theta_R50_maj_asec
            ext_dec_summary[f'f{component_id}_theta_R50_min_asec'] = theta_R50_min_asec
            ext_dec_summary[f'f{component_id}_theta_R50_asec_err'] = theta_R50_asec_err
            ext_dec_summary[f'f{component_id}_theta_R50_maj_asec_err'] = theta_R50_maj_asec_err
            ext_dec_summary[f'f{component_id}_theta_R50_min_asec_err'] = theta_R50_min_asec_err

            ext_dec_summary[f'f{component_id}_Tb_Rn'] = TB_Rn
            ext_dec_summary[f'f{component_id}_Tb_Rn_err'] = TB_Rn_err
            ext_dec_summary[f'f{component_id}_Tb_R50'] = TB_R50
            ext_dec_summary[f'f{component_id}_Tb_R50_err'] = TB_R50_err 
        
            #properties from model images
            ext_dec_summary[f'comp_{component_id}_R50_d'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii'] * self.pix_to_pc
            ext_dec_summary[f'comp_{component_id}_R50_d_err'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C50radii_err'] * self.pix_to_pc
            ext_dec_summary[f'comp_{component_id}_R95_d'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C95radii'] * self.pix_to_pc
            ext_dec_summary[f'comp_{component_id}_R95_d_err'] = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['C95radii_err'] * self.pix_to_pc


            A_R50_d = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix50']
            A_R50_d_err = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix50_err']
            ext_dec_summary[f'comp_{component_id}_A_R50_d'] = mlibs.pix_area_to_kpc_area(A_R50_d,self.pix_to_pc)
            ext_dec_summary[f'comp_{component_id}_A_R50_d_err'] = mlibs.pix_area_to_kpc_area(A_R50_d_err,self.pix_to_pc)

            
            A_R95_d = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix95']
            A_R95_d_err = self.SMFR.components_deconv_props.loc[f'model_c_deconv_{component_id}_props']['npix95_err']
            

            
            ext_dec_summary[f'comp_{component_id}_A_R95_d'] = mlibs.pix_area_to_kpc_area(A_R95_d,self.pix_to_pc)
            ext_dec_summary[f'comp_{component_id}_A_R95_d_err'] = mlibs.pix_area_to_kpc_area(A_R95_d_err,self.pix_to_pc)



            ext_dec_summary[f'comp_{component_id}_s_Snu_A50_d'] = 0.5*Snu / ext_dec_summary[f'comp_{component_id}_A_R50_d']
            ext_dec_summary[f'comp_{component_id}_s_Snu_A50_d_err'] = 0.5*Snu_err / ext_dec_summary[f'comp_{component_id}_A_R50_d']
            ext_dec_summary[f'comp_{component_id}_s_Snu_A95_d'] = Snu / ext_dec_summary[f'comp_{component_id}_A_R95_d']
            ext_dec_summary[f'comp_{component_id}_s_Snu_A95_d_err'] = Snu_err / ext_dec_summary[f'comp_{component_id}_A_R95_d']
            ext_dec_summary[f'comp_{component_id}_s_Snu_A_Rn'] = 0.5*Snu / ext_dec_summary[f'f{component_id}_A_Rn']
            ext_dec_summary[f'comp_{component_id}_s_Snu_A_Rn_err'] = 0.5*Snu_err / ext_dec_summary[f'f{component_id}_A_Rn']

            #CONVOLVED QUANTITIES
            ext_conv_summary[f'comp_{component_id}_R50'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C50radii'] * self.pix_to_pc
            ext_conv_summary[f'comp_{component_id}_R50_err'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C50radii_err'] * self.pix_to_pc
            ext_conv_summary[f'comp_{component_id}_R95'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C95radii'] * self.pix_to_pc
            ext_conv_summary[f'comp_{component_id}_R95_err'] = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['C95radii_err'] * self.pix_to_pc
            
            A_R50 = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix50']
            A_R50_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix50_err']
            
            
            ext_conv_summary[f'comp_{component_id}_A_R50'] = mlibs.pix_area_to_kpc_area(A_R50,self.pix_to_pc)
            ext_conv_summary[f'comp_{component_id}_A_R50_err'] = mlibs.pix_area_to_kpc_area(A_R50_err,self.pix_to_pc)


                
            A_R95 = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix95']
            A_R95_err = self.SMFR.components_conv_props.loc[f'model_c_conv_{component_id}_props']['npix95_err']    
            

            ext_conv_summary[f'comp_{component_id}_A_R95'] = mlibs.pix_area_to_kpc_area(A_R95,self.pix_to_pc)
            ext_conv_summary[f'comp_{component_id}_A_R95_err'] = mlibs.pix_area_to_kpc_area(A_R95_err,self.pix_to_pc)




            ext_conv_summary[f'comp_{component_id}_s_Snu_A50'] = 0.5*Snu / ext_conv_summary[f'comp_{component_id}_A_R50']
            ext_conv_summary[f'comp_{component_id}_s_Snu_A50_err'] = 0.5*Snu_err / ext_conv_summary[f'comp_{component_id}_A_R50']
            ext_conv_summary[f'comp_{component_id}_s_Snu_A95'] = Snu / ext_conv_summary[f'comp_{component_id}_A_R95']
            ext_conv_summary[f'comp_{component_id}_s_Snu_A95_err'] = Snu_err / ext_conv_summary[f'comp_{component_id}_A_R95']
        
            _logging_.logger.info(f"Snu = {1000*ext_dec_summary[f'comp_{component_id}_Snu']:.2f}" 
                                  f" +/- {1000*ext_dec_summary[f'comp_{component_id}_Snu_err']:.2f} mJy")
            _logging_.logger.info(f"Rn = {ext_dec_summary[f'f{component_id}_Rn']:.2f}" 
                                  f" +/- {ext_dec_summary[f'f{component_id}_Rn_err']:.2f} pc")
            _logging_.logger.info(f"R50_d = {ext_dec_summary[f'comp_{component_id}_R50_d']:.2f} pc"
                                    f" +/- {ext_dec_summary[f'comp_{component_id}_R50_d_err']:.2f} pc")
            _logging_.logger.info(f"R95_d = {ext_dec_summary[f'comp_{component_id}_R95_d']:.2f}"
                                  f" +/- {ext_dec_summary[f'comp_{component_id}_R95_d_err']:.2f} pc")
            
            _logging_.logger.info(f"R50 = {ext_conv_summary[f'comp_{component_id}_R50']:.2f}"
                                  f" +/- {ext_conv_summary[f'comp_{component_id}_R50_err']:.2f} pc")
            _logging_.logger.info(f"R95 = {ext_conv_summary[f'comp_{component_id}_R95']:.2f}"
                                  f" +/- {ext_conv_summary[f'comp_{component_id}_R95_err']:.2f} pc")

            _logging_.logger.info(f"Theta_Rn = {ext_dec_summary[f'f{component_id}_theta_Rn_pc']:.2f}" 
                                  f" +/- {ext_dec_summary[f'f{component_id}_theta_Rn_pc_err']:.2f} pc")
            _logging_.logger.info(f"Theta_R50 = {ext_dec_summary[f'f{component_id}_theta_R50_pc']:.2f}"
                                    f" +/- {ext_dec_summary[f'f{component_id}_theta_R50_pc_err']:.2f} pc")
            _logging_.logger.info(f"Tb_Rn = ({ext_dec_summary[f'f{component_id}_Tb_Rn']:.2f}" 
                                  f" +/- {ext_dec_summary[f'f{component_id}_Tb_Rn_err']:.2f}) x 10^5 K")
            _logging_.logger.info(f"Tb_R50 = ({ext_dec_summary[f'f{component_id}_Tb_R50']:.2f}"
                                    f" +/- {ext_dec_summary[f'f{component_id}_Tb_R50_err']:.2f}) x 10^5 K")

            _logging_.logger.info(f"A_Rn = ({1e3*ext_dec_summary[f'f{component_id}_A_Rn']:.3f}" 
                                  f" +/- {1e3*ext_dec_summary[f'f{component_id}_A_Rn_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"A50_d = ({1e3*ext_dec_summary[f'comp_{component_id}_A_R50_d']:.3f}"
                                  f" +/- {1e3*ext_dec_summary[f'comp_{component_id}_A_R50_d_err']:.3f}) x 1e-3 kpc^2")            
            _logging_.logger.info(f"A95_d = ({1e3*ext_dec_summary[f'comp_{component_id}_A_R95_d']:.3f}"
                                  f" +/- {1e3*ext_dec_summary[f'comp_{component_id}_A_R95_d_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"A50 = ({1e3*ext_conv_summary[f'comp_{component_id}_A_R50']:.3f}"
                                  f" +/- {1e3*ext_conv_summary[f'comp_{component_id}_A_R50_err']:.3f}) x 1e-3 kpc^2")
            _logging_.logger.info(f"A95 = ({1e3*ext_conv_summary[f'comp_{component_id}_A_R95']:.3f}"
                                  f" +/- {1e3*ext_conv_summary[f'comp_{component_id}_A_R95_err']:.3f}) x 1e-3 kpc^2")        
            _logging_.logger.info(f"sSnu_A95_d = {ext_dec_summary[f'comp_{component_id}_s_Snu_A95_d']:.3f}"
                                  f" +/- {ext_dec_summary[f'comp_{component_id}_s_Snu_A95_d_err']:.3f} Jy/kpc^2")
            _logging_.logger.info(f"sSnu_A95 = {ext_conv_summary[f'comp_{component_id}_s_Snu_A95']:.3f}"
                                  f" +/- {ext_conv_summary[f'comp_{component_id}_s_Snu_A95_err']:.3f} Jy/kpc^2")
        
        # ext_conv_summary[f'Snu_ext_total'] = self.SMFR.results_ext_conv_morpho['total_flux_mask'].values[0]
        # ext_conv_summary[f'Snu_ext_total_err'] = self.SMFR.results_ext_conv_morpho['flux_error_res_3'].values[0]
        ext_conv_summary[f'Snu_ext_total'] = self.SMFR.data_properties['total_flux_mask'].values[0]-self.SMFR.results_compact_conv_morpho['total_flux_mask'].values[0]
        ext_conv_summary[f'Snu_ext_total_2'] = self.SMFR.results_fit['flux_density_ext']
        ext_conv_summary[f'Snu_ext_total_err'] = self.SMFR.data_properties['flux_error_res_3'].values[0]
        _logging_.logger.info(f"Snu Extended Total = {1000*ext_conv_summary[f'Snu_ext_total']:.2f}" 
                                f" +/- {1000*ext_conv_summary[f'Snu_ext_total_err']:.2f} mJy")
        
        _logging_.logger.info(f'-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+')
        
        A50_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix50'].values[0],self.pix_to_pc)
        A50_err_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix50_err'].values[0],self.pix_to_pc)
        A95_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix95'].values[0],self.pix_to_pc)
        A95_err_kpc_ext_conv = mlibs.pix_area_to_kpc_area(self.SMFR.results_ext_conv_morpho['npix95_err'].values[0],self.pix_to_pc)
        
        _logging_.logger.info(f"Snu Extended Data = {1000*self.SMFR.results_ext_conv_morpho[f'total_flux_mask'].values[0]:.2f}" 
                                f" +/- {1000*self.SMFR.results_ext_conv_morpho[f'flux_error_res_3'].values[0]:.2f} mJy")

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
    parser = argparse.ArgumentParser(description='Morphen.')
    parser.add_argument('-filename',  '--filename',  required=False, \
                        help='Image data.')
    parser.add_argument('-residualname',  '--residualname',  required=False, \
                        help='Associated Residual from Image Data')
    parser.add_argument('-psfname',  '--psfname',  required=False, \
                        help='PSF name.')

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
    # parser.add_argument('-filename',  '--filename',  required=False, \
    #                     help='The input array of the light profile of the galaxy.')

    args = parser.parse_args()

    if args.filename != None:
        input_data = read_data(filename=args.filename,
                               residualname=args.residualname,)
        if "--image_stats" in sys.argv:
            radio_image_analysis(input_data)
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
                                                   SMFR,
                                                   z=args.redshift)
                else:
                    print("Error: Please, provide a residual image (e.g. the one "
                          "generate during interferometric deconvolution).")

            if "--general_fit" in sys.argv:
                if args.psfname != None:
                    SMFR = sersic_multifit_general(input_data, SE,
                                                 method2=args.solver2)

                else:
                    print("Error: Please, provide a psf image.")
