#setting the GPU memory fraction to be used of 25% should be fine!
os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '0.25'
os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'true'
# os.environ["XLA_FLAGS"] = '--xla_force_host_platform_device_count=6'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '0'


# os.environ["NUM_CPUS"] = "12"  # Set the desired number of CPU threads here
# # os.environ["XLA_FLAGS"] = "--xla_cpu_threads=2"
os.environ["TF_XLA_FLAGS"] = "--xla_cpu_threads=8"  # Set the desired number of CPU
# # threads here
#
os.environ['MKL_NUM_THREADS']='8'
os.environ['OPENBLAS_NUM_THREADS']='8'
os.environ["NUM_INTER_THREADS"]="8"
# os.environ["NUM_INTRA_THREADS"]="12"
# #
os.environ["XLA_FLAGS"] = ("--xla_cpu_multi_thread_eigen=false "
                           "intra_op_parallelism_threads=8")



"""
#Config
"""
def reset_rc_params():
    """
    Global configuration for matplotlib.pyplot
    """
    mpl.rcParams.update({'font.size': 16,
                         "text.usetex": False,  #
                         "font.family": "sans-serif",
                         'mathtext.fontset': 'stix',
                         "font.family": "sans",
                         'font.weight': 'medium',
                         'font.family': 'STIXGeneral',
                         'xtick.labelsize': 16,
                         'figure.figsize': (6, 4),
                         'ytick.labelsize': 16,
                         'axes.labelsize': 16,
                         'xtick.major.width': 1,
                         'ytick.major.width': 1,
                         'axes.linewidth': 1.5,
                         'axes.edgecolor': 'orange',
                         'lines.linewidth': 2,
                         'legend.fontsize': 14,
                         'grid.linestyle': '--',
                         # 'grid.color':'black',
                         'axes.grid.which': 'major',
                         'axes.grid.axis': 'both',
                         'axes.spines.right': True,
                         'axes.grid': True,
                         'axes.titlesize': 16,
                         'legend.framealpha': 1.0
                         })
    pass
