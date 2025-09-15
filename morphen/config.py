#setting the GPU memory fraction to be used of 25% should be fine!
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '0.25'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
# os.environ["XLA_FLAGS"] = '--xla_force_host_platform_device_count=6'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'



# os.environ["NUM_CPUS"] = "12"  # Set the desired number of CPU threads here
# # os.environ["XLA_FLAGS"] = "--xla_cpu_threads=2"
os.environ["TF_XLA_FLAGS"] = "--xla_cpu_threads=6"  # Set the desired number of CPU
# # threads here
#
os.environ['MKL_NUM_THREADS']='6'
os.environ['OPENBLAS_NUM_THREADS']='6'
os.environ["NUM_INTER_THREADS"]="6"
# os.environ["NUM_INTRA_THREADS"]="12"
#
os.environ["XLA_FLAGS"] = ("--xla_cpu_multi_thread_eigen=false "
                           "intra_op_parallelism_threads=6")



"""
#Config
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
                            # 'text.usetex' : True,
                            # 'font.family' : 'serif',
                            # 'font.serif' : ['Garamond Libre', 'EB Garamond', 'Cormorant Garamond', 'serif'],
                            #  'text.latex.preamble': r'''
                            #     \usepackage{ebgaramond-maths}
                            #     \usepackage{garamondlibre}
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
    # mpl.rcParams.update({'font.size': 16,
    #                      "text.usetex": False,  #
    #                      "font.family": "sans-serif",
    #                      'mathtext.fontset': 'stix',
    #                      "font.family": "sans",
    #                      'font.weight': 'medium',
    #                      'font.family': 'STIXGeneral',
    #                      'xtick.labelsize': 16,
    #                      'figure.figsize': (6, 4),
    #                      'ytick.labelsize': 16,
    #                      'axes.labelsize': 16,
    #                      'xtick.major.width': 1,
    #                      'ytick.major.width': 1,
    #                      'axes.linewidth': 1.5,
    #                      'axes.edgecolor': 'orange',
    #                      'lines.linewidth': 2,
    #                      'legend.fontsize': 14,
    #                      'grid.linestyle': '--',
    #                      # 'grid.color':'black',
    #                      'axes.grid.which': 'major',
    #                      'axes.grid.axis': 'both',
    #                      'axes.spines.right': True,
    #                      'axes.grid': True,
    #                      'axes.titlesize': 16,
    #                      'legend.framealpha': 1.0,
    #                      })
    pass


def reset_rc_params_simple():
    """
    Global configuration for matplotlib.pyplot
    """
    global_font_size = 14
    mpl.rcParams.update({
        'font.size': global_font_size,
        'text.usetex': False, 
        'font.family': 'sans-serif',
        'mathtext.fontset': 'stix',
        'font.family': 'sans',
        'font.weight': 'medium',  
        'font.family': 'STIXGeneral',
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
        'axes.grid.which': 'major',  
        'axes.grid.axis': 'both', 
        'axes.spines.right': True,
        'axes.grid': True,
        'axes.titlesize' : global_font_size,
        'legend.framealpha': 1.0
    })
    pass


def set_consistent_style(figsize,factor = 3.0):
    """Set matplotlib style for consistent appearance"""
    base_size = (figsize[0] * figsize[1]) ** 0.5 * factor
    
    plt.rcParams.update({
        'font.size': base_size*0.9,
        'axes.titlesize': base_size,
        'axes.labelsize': base_size,
        'xtick.labelsize': base_size,
        'ytick.labelsize': base_size,
        'legend.fontsize': base_size-1,
    })


cb_colors = {
    # Colourblind friendly attempt.
    'blue':         '#0077BB',    # Blue
    'orange':       '#E69F00',    # Orange
    'darker_yellow':'#E6C700',    # Darker yellow
    'light_yellow': '#F9E949',    # Light yellow
    'lime_green':   '#88BB44',    # Lime-green
    'purple':       '#994F9F',    # Purple
    'bamboo':       '#D5A87F',    # Light brown/tan
    'turquoise':    '#5DC1C0',    # Turquoise
    'grey':         '#777777',    # Medium grey
    'red':          '#CC3311',    # Red
    'magenta':      '#EE3377',    # Magenta/Pink
    'dark_blue':    '#33557B',    # Dark blue
    'mint':         '#56C667',    # Mint green
    'orange':         'orange',     # Light brown/tan
    'dark_teal':    '#006B5C',    # Dark teal
    'dark_orange':  '#FF8425',    # Light yellow
}