import pandas as pd
import numpy as np
import argparse
import sys
import os
import subprocess
import shlex
import glob
import time
import uuid
import shutil
from pathlib import Path

sys.path.append(os.path.dirname(os.path.abspath(__file__)))
from wsclean_container import resolve_container, DEFAULT_CONTAINER
import channel_division

import resource
soft, hard = resource.getrlimit(resource.RLIMIT_NPROC)
new_soft_limit = hard
resource.setrlimit(resource.RLIMIT_NPROC, (new_soft_limit, hard))
soft, hard = resource.getrlimit(resource.RLIMIT_NPROC)


def imaging(g_name, field, uvtaper, robust, base_name='clean_image',
            continue_clean='False',nc=4):
    g_vis = g_name + '.ms'
    """
    # uvtaper_mode+uvtaper_args+'.'+uvtaper_addmode+uvtaper_addargs+
    """
    print(uvtaper_addmode, uvtaper_addargs, robust)
    if uvtaper != '':
        taper = 'taper_'
    else:
        taper = ''

    if continue_clean == 'True':
        print('*************************************')
        image_to_continue = glob.glob(f"{root_dir_sys}/*MFS-image.fits")
        image_to_continue.sort(key=os.path.getmtime, reverse=False)
        image_to_continue = os.path.basename(image_to_continue[-1])
        image_deepclean_name = image_to_continue.replace('-MFS-image.fits','')
        print('Using prefix from previous image: ', image_deepclean_name)


    if continue_clean == 'False':
        image_deepclean_name = (base_name + '_' + g_name + '_' +
                                imsizex + 'x' + imsizey + '_' + \
                                'am_' + str(args.nsigma_automask) + '_at_' + str(args.nsigma_autothreshold) + '_' + \
                                cell + '_' + \
                                deconvolver[1:] + '_' + taper + \
                                uvtaper + '_r' + str(robust))

    ext = ''
    if '-join-channels' in deconvolver_args:
        print('Using mtmfs method.')
        ext = ext + '-MFS'
    ext = ext + '-image.fits'

    print(image_deepclean_name)

    if not os.path.exists(root_dir_sys + image_deepclean_name + ext) or continue_clean == 'True':
        # if nc < 4:
        #     _nc = 4
        # elif nc > 6:
        #     _nc = nc
        # else:
        #     _nc = 6
        if nc > 16:
            _nc = 16
        else:
            # _nc = 8
            _nc = nc
        if running_container == 'native':
            # 'mpirun -np 4 wsclean-mp'

            command_exec = (
                'mpirun -np '+str(_nc)+' wsclean-mp -name ' + root_dir + image_deepclean_name +
                ' -size ' + imsizex + ' ' + imsizey + ' -scale ' + cell +
                ' ' + gain_args + ' -niter ' + niter + ' -weight ' + weighting +
                ' ' + robust + ' ' + auto_mask + ' ' + auto_threshold + mask_file +
                ' ' + deconvolver + ' ' + deconvolver_options +
                ' ' + deconvolver_args + ' ' + taper_mode + uvtaper +
                ' ' + opt_args + ' ' + data_column + ' ' + root_dir + g_vis)
            print(' ++==>> Command to be executed by WSClean: ')
            print(command_exec)
            subprocess.run(shlex.split(command_exec), check=True)

        # if running_container == 'singularity':
        #     command_exec = (
        #     'singularity exec --nv '
        #     '--bind /usr/local/cuda-13.0:/usr/local/cuda-13.0 ' \
        #     '--bind /usr/local/cuda-13.0/lib64:/usr/local/cuda/lib64 ' \
        #     '--bind ' + mount_dir + ' ' + wsclean_dir +
        #     ' ' + 'wsclean -name ' + root_dir +
        #     # ' ' + 'mpirun -np ' + str(_nc) + ' wsclean-mp -name ' + root_dir +
        #     # ' ' + 'mpirun --use-hwthread-cpus wsclean-mp -name ' + root_dir +
        #     image_deepclean_name +
        #     ' -size ' + imsizex + ' ' + imsizey + ' -scale ' + cell +
        #     ' ' + gain_args + ' -niter ' + niter + ' -weight ' + weighting +
        #     ' ' + robust + ' ' + auto_mask + ' ' + auto_threshold + mask_file +
        #     ' ' + deconvolver + ' ' + deconvolver_options +
        #     ' ' + deconvolver_args + ' ' + taper_mode + uvtaper +
        #     ' ' + opt_args + ' ' + data_column + ' ' + root_dir + g_vis)

        #     print(' ++==>> Command to be executed by Singularity > WSClean: ')
        #     print(command_exec)
        #     os.system(command_exec)

        # Define a scratch dir on the host
        # scratch_dir = '/media/sagauga/starbyte/wsclean_temp/'
        # scratch_dir = f'/media/sagauga/starbyte/wsclean_temp/{uuid.uuid4().hex}/'
        # scratch_dir = f'/home/lucatelli/astronomical_data/wsclean_temp/{uuid.uuid4().hex}/'
        if os.path.exists('/media/sagauga/temp_wsclean/'):
            scratch_dir = f'/media/sagauga/temp_wsclean/wsclean_temp/{uuid.uuid4().hex}/'
        # if os.path.exists('/media/sagauga/void/'):
        #     scratch_dir = f'/media/sagauga/void/wsclean_temp/{uuid.uuid4().hex}/'
        # elif os.path.exists('/media/sagauga/galnet/'):
        #     scratch_dir = f'/media/sagauga/galnet/wsclean_temp/{uuid.uuid4().hex}/'
            # scratch_dir = f'/media/sagauga/starbyte/wsclean_temp/{uuid.uuid4().hex}/'
        # if os.path.exists('/media/sagauga/starbyte/'):
        #     scratch_dir = f'/media/sagauga/starbyte/wsclean_temp/{uuid.uuid4().hex}/'
        #     # scratch_dir = f'/media/sagauga/galnet/wsclean_temp/{uuid.uuid4().hex}/'
        elif os.path.exists('/mnt/scratch/lucatelli/temp_wsclean/'):
            scratch_dir = f'/mnt/scratch/lucatelli/temp_wsclean/wsclean_temp/{uuid.uuid4().hex}/'
        elif os.path.exists('/home/lucatelli/astronomical_data/'):
            scratch_dir = f'/home/lucatelli/astronomical_data/wsclean_temp/{uuid.uuid4().hex}/'
        else:
            # scratch_dir = str(Path.home() / f'wsclean_temp/{uuid.uuid4().hex}/')
            scratch_dir = str(f'./wsclean_temp/{uuid.uuid4().hex}/')
        os.makedirs(scratch_dir, exist_ok=True)
        print(' >> Using scratch directory for WSClean temporary files: ', scratch_dir)


        if running_container == 'singularity':
            command_exec = (
            'singularity exec --nv '
            # '--bind /usr/local/cuda-13.0:/usr/local/cuda-13.0 ' \
            # '--bind /usr/local/cuda-13.0/lib64:/usr/local/cuda/lib64 ' \
            '--bind ' + scratch_dir + ':/wsclean_tmp ' +  # <-- bind scratch
            '--bind ' + mount_dir + ' ' + wsclean_dir +
            ' ' + 'wsclean -name ' + root_dir + image_deepclean_name +
            # ' ' + 'mpirun -np ' + str(_nc) + ' wsclean-mp -name ' + root_dir + image_deepclean_name +
            ' -size ' + imsizex + ' ' + imsizey + ' -scale ' + cell +
            ' ' + gain_args + ' -niter ' + niter + ' -weight ' + weighting +
            ' ' + robust + ' ' + auto_mask + ' ' + auto_threshold + mask_file +
            ' ' + deconvolver + ' ' + deconvolver_options +
            ' ' + deconvolver_args + ' ' + taper_mode + uvtaper +
            ' -temp-dir /wsclean_tmp ' +                  # <-- use it
            ' ' + opt_args + ' ' + data_column + ' ' + root_dir + g_vis)
            print(' ++==>> Command to be executed by Singularity > WSClean: ')
            print(command_exec)
            subprocess.run(shlex.split(command_exec), check=True)

        shutil.rmtree(scratch_dir, ignore_errors=True)
        image_stats = {
            "#basename": image_deepclean_name + ext}  # get_image_statistics(image_deep_selfcal  + ext)
        image_stats['imagename'] = image_deepclean_name + ext
        '''
        save dictionary to file
        '''
        return (image_stats)
    else:
        print('Skipping imaging; already done.')
        return (None)

    # pass


def parse_float_list(str_values):
    return [float(val.strip()) for val in str_values.strip('[]').split(',')]


def parse_str_list(str_values):
    return [s.strip() for s in str_values.strip('[]').split(',')]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Helper for wsclean imaging.")
    parser.add_argument("--p", type=str, help="The path to the MS file.")
    parser.add_argument("--f", nargs='?', default=False,
                        const=True, help="The name of the ms file")
    parser.add_argument("--r",
                        type=parse_float_list, nargs='?',
                        const=True, default=[0.5],
                        help="List of robust values")
    parser.add_argument("--t", type=parse_str_list, nargs='?',
                        const=True, default=[''],
                        help="List of sky-tapers values")

    parser.add_argument("--data", type=str, nargs='?', default='DATA',  # 'CORRECTED_DATA'
                        help="Which data column to use")

    parser.add_argument("--wsclean_install", type=str, nargs='?', default='singularity',
                        help="How wsclean was installed (singularity or native)?")

    parser.add_argument("--wsclean_sif", type=str, nargs='?', default=None,
                        help="Path to the WSClean singularity/apptainer image. "
                             "If omitted, the image is looked up in $PH4SER_WSCLEAN_SIF, "
                             "$PH4SER_CONTAINER_DIR, ph4ser_config, the ph4ser module "
                             "directory and ~/.ph4ser/containers, in that order. "
                             "Run 'python wsclean_container.py' to download it.")

    parser.add_argument("--update_model", type=str, nargs='?', default='False',
                        help="Update model after cleaning?")

    parser.add_argument("--deconvolution_mode", type=str, nargs='?', default='good',
                        help="This is not a wsclean parameter. "
                             "If 'good' will use a proper set of wsclean arguments to perform a "
                             "good deconvolution, but slower. If 'fast' will not use complex "
                             "deconvolution, "
                             "e.g. no multi-scale, no MFS, etc. This is not intended for final "
                             "science data. It is intended for quick-look data.")


    parser.add_argument("--with_multiscale", type=str, nargs='?', default='False',
                        help="Use multiscale deconvolver?")

    parser.add_argument("--shift", type=str, nargs='?', default=None,
                        help="New phase center to shift for imaging."
                             "Eg. --shift 13:15:30.68 +62.07.45.357")

    parser.add_argument("--scales", type=str, nargs='?', default="None",
                        help="Scales to be used with the multiscale deconvolver in WSClean. "
                             "If None, scales will be determined automatically by WSClean.")
    parser.add_argument("--maxmscales", type=str, nargs='?', default="8",
                        help="Maximum number of scales to be used with the multiscale deconvolver in WSClean. ")

    parser.add_argument("--sx", type=str, nargs='?', default='2048',
                        help="Image Size x-axis")
    parser.add_argument("--sy", type=str, nargs='?', default='2048',
                        help="Image Size y-axis")
    parser.add_argument("--cellsize", type=str, nargs='?', default='0.05asec',
                        help="Cell size")
    parser.add_argument("--niter", type=str, nargs='?', default='50000',
                        help="Number of iterations during cleaning.")

    parser.add_argument("--maxuv_l", type=str, nargs='?', default=None,
                        help="Max uv distance in lambda.")

    parser.add_argument("--minuv_l", type=str, nargs='?', default=None,
                        help="Min uv distance in lambda.")

    parser.add_argument("--nsigma_automask", type=str, nargs='?', default='4.0',
                        help="Sigma level for automasking in wsclean.")

    parser.add_argument("--nsigma_autothreshold", type=str, nargs='?', default='2.0',
                        help="Sigma level for autothreshold in WSClean.")

    parser.add_argument("--mask", nargs='?', default=None,
                        const=True, help="A fits-file mask to be used.")

    parser.add_argument("--nc", type=int, nargs='?', default=4,
                        help="Number of channels division to be used in "
                             "the MFS deconvolution.")

    parser.add_argument("--channel_division", type=str, nargs='?', default='auto',
                        help="How the bandwidth is divided into --nc sub-bands. "
                             "'auto': split frequencies computed from the MS "
                             "(gaps + region where all instruments overlap; "
                             "see channel_division.py). 'default': WSClean's "
                             "own division. 'gap': -gap-channel-division. "
                             "Anything else is passed as a comma-separated "
                             "list of split frequencies in Hz to "
                             "-channel-division-frequencies.")

    parser.add_argument("--channel_division_mode", type=str, nargs='?', default='bandwidth',
                        help="For --channel_division auto: split each frequency "
                             "block into equal 'bandwidth' or equal 'weight'.")

    parser.add_argument("--negative_arg", type=str, nargs='?', default='negative',
                        help="Enable/disable negative clean components during cleaning.")

    parser.add_argument("--quiet", type=str, nargs='?', default='False',
                        help="Enable/disable WSClean outputs.")


    parser.add_argument("--continue_clean", type=str, nargs='?', default='False',
                        help="Continue cleaning?")


    parser.add_argument("--opt_args", type=str, nargs='?', default='',
                        help="Optional/additional arguments to be passed to WSClean. "
                             "Warning: Do not repeat previously defined arguments."
                             "Example: ' -apply-facet-beam -dd-psf-grid 6 6 -facet-beam-update 60 '")


    parser.add_argument("--save_basename", type=str, nargs='?', default='image',
                        help="optional basename for saving image files.")

    args = parser.parse_args()

    if args.update_model == 'True':
        update_model_option = ' -update-model-required '
    else:
        update_model_option = ' -no-update-model-required '

    running_container = args.wsclean_install

    if running_container == 'native':
        os.environ['OPENBLAS_NUM_THREADS'] = '1'

    threads = os.cpu_count()
    # threads_proc = int(threads - 8) if threads > 16 else int(threads)
    # threads_proc = int(16)
    threads_proc = int(threads)
    
    print(f' >> Using {threads_proc} threads for processing. Total threads available: {threads}')

    # for i in range(len(image_list)):
    field = os.path.basename(args.f).replace('.ms', '')
    g_name = field
    root_dir_sys = os.path.dirname(args.f) + '/'
    robusts = args.r
    tapers = args.t

    if running_container == 'singularity':
        mount_dir = root_dir_sys + ':/mnt'
        root_dir = '/mnt/'
        # The image is not shipped with the repository: resolve it from the
        # command line, the environment, the config, the module directory or
        # the shared cache. See wsclean_container.py.
        wsclean_dir = resolve_container(DEFAULT_CONTAINER,
                                        explicit=args.wsclean_sif)
        print('Using WSClean from singularity image: ', wsclean_dir)
        print('WSClean version from singularity image: ')
        subprocess.run([
            "singularity", "exec", "--nv",
            # "--bind", "/usr/local/cuda-13.0:/usr/local/cuda-13.0",
            # "--bind", "/usr/local/cuda-13.0/lib64:/usr/local/cuda/lib64",
            "--bind", mount_dir,
            wsclean_dir, "wsclean", "--version"
        ], check=True)
    if running_container == 'native':
        mount_dir = ''
        root_dir = root_dir_sys
        os.environ['OPENBLAS_NUM_THREADS'] = '1'

    base_name = args.save_basename

    ## Setting image and deconvolution noptions.
    ### Cleaning arguments
    auto_mask = ' -auto-mask ' + args.nsigma_automask
    # auto_mask = ' '
    auto_threshold = ' -auto-threshold ' + args.nsigma_autothreshold
    if args.mask == 'None' or args.mask == None:
        mask_file = ' '
    else:
        # if args.mask != 'None' or args.mask != None:
        if running_container == 'native':
            mask_file = ' -fits-mask ' + args.mask + ' '
        if running_container == 'singularity':
            mask_file = ' -fits-mask ' + root_dir+os.path.basename(args.mask) + ' '

    # data to run deconvolution
    data_column = ' -data-column ' + args.data
    with_multiscale = args.with_multiscale
    ### Selecting the deconvolver
    # deconvolution_mode = 'good'
    if args.deconvolution_mode == 'good':
        if with_multiscale == True or with_multiscale == 'True':
            deconvolver = '-multiscale'
            deconvolver_options = ' -multiscale-scale-bias 0.65 -multiscale-gain 0.05 -multiscale-convolution-padding 1.3 -padding 1.3 -clean-border 1 '
            if (args.scales is None) or (args.scales == 'None'):
                deconvolver_options = deconvolver_options + ' -multiscale-max-scales ' + args.maxmscales + ' '
            else:
                deconvolver_options = (deconvolver_options + ' -multiscale-scales ' + args.scales + ' ')
        else:
            deconvolver = ' '
            deconvolver_options = (' ')
        nc = args.nc
        # nc = 8 #debug
        negative_arg = '-'+args.negative_arg
        if nc > 3:
            dec_chan = int(np.nanmax([4, nc//2]))
        else:
            dec_chan = nc
        dec_chan = nc #debug
        deconvolver_args = (' '
                            # Setting some important arguments to custom values
                            '-deconvolution-threads ' + str(threads_proc) + ' -j ' + str(threads_proc) + ' '
                            '-parallel-reordering ' + str(threads_proc) + ' '
                            '-parallel-deconvolution 1024 '
                            # '-parallel-deconvolution 1024 '
                            '-channels-out '+str(nc)+' -join-channels ' + negative_arg + ' '
                            '-no-mf-weighting ' # this must be set to generate science images 
                            '-fit-spectral-pol  ' +str(3)+' '+' -deconvolution-channels ' +str(dec_chan)+' '
                            # '-fit-spectral-log-pol ' +str(2)+' '+' -deconvolution-channels ' +str(dec_chan)+' '
                            '-gridder wgridder -wgridder-accuracy 1e-5 -parallel-gridding 8 '
                            '-apply-primary-beam '
                            #updates 2026
                            '-weighting-rank-filter 3 -weighting-rank-filter-size 128 '
                            # '-gridder idg -idg-mode hybrid -grid-with-beam -circular-beam -save-psf-pb -save-uv '
                            #cleanup 2026
                            # '-wstack-nwlayers-factor 3 ' +' -deconvolution-channels ' +str(nc)
                            # '-fit-spectral-pol  ' +str(4)+' '
                            # '-save-psf-pb -save-weights '
                            # Other arguments that may be useful
                            # '-store-imaging-weights -save-weights '
                            # '-no-negative -abs-threshold 3.5e-6 '
                            # Some testing parameters
                            # '-save-weights -local-rms -local-rms-window 50 '
                            # '-beam-fitting-size 0.1 '
                            # ' -circular-beam -beam-size 0.1arcsec -beam-fitting-size = 0.7 ' 
                            # '-facet-beam-update 60 -save-aterms -diagonal-solutions '
                            # '-apply-facet-solutions ' 
                            # '-gridder idg -idg-mode hybrid -apply-primary-beam 
                            # '-apply-facet-beam -facet-beam-update 120 -save-aterms '
                            # -diagonal-solutions ' 
                            # '-save-source-list '
                            )
    if args.deconvolution_mode == 'fast':
        deconvolver = ' '
        deconvolver_options = (' ')
        deconvolver_args = (' -save-source-list '
                            '-deconvolution-threads 16 -j 16 -parallel-reordering 16 '
                            '-parallel-deconvolution 2048') #-parallel-gridding 24

    # image parameters
    weighting = 'briggs'
    imsizex = args.sx
    imsizey = args.sy
    cell = args.cellsize
    niter = args.niter

    """
    # taper options (this is a to-do)
    uvtaper_mode = '-taper-tukey'
    uvtaper_args = '900000'
    uvtaper_addmode = '-maxuv-l'
    uvtaper_addargs = '800000'
    taper_mode='-taper-gaussian '
    uvtaper_mode = '-taper-gaussian'
    uvtaper_args = '0.05asec'
    """
    uvtaper_addmode = ''
    uvtaper_addargs = ''
    # uvtaper = uvtaper_mode + ' '+ uvtaper_args + ' ' +uvtaper_addmode + ' ' +
    # uvtaper_addargs
    uvselection = ''
    if args.maxuv_l is not None:
        uvselection = ' -maxuv-l ' + args.maxuv_l + ' '
    if args.minuv_l is not None:
        uvselection = uvselection + ' -minuv-l ' + args.minuv_l + ' '
    # general arguments
    gain_args = ' -mgain 0.5 -gain 0.05 -nmiter 20 ' # -super-weight 9.0

    if args.shift == 'None' or args.shift == None:
        # if args.shift != ' ':
        shift_options = ' '
    else:
        shift_options = ' -shift ' + args.shift + ' '
    # shift_options = ' '  # -shift 13:15:30.68  +62.07.45.357 '#' -shift 13:15:28.903
    # +62.07.11.886 '
    if args.quiet == 'True':
        quiet = ' -quiet '
        # quiet = ' ' #debug
    else:
        quiet = ' -v '

    # quiet = ' '

    if args.continue_clean == 'True':
        """
        Not HANDLED PROPERLY YET
        """
        continue_clean = ' --continue '
    else:
        continue_clean = ' '
    opt_args = (
                # ' -mem 80 -abs-mem 35 '
                # ' -abs-mem 80 '
                # '-pol RL,LR -no-negative -circular-beam -no-reorder '
                # ' -save-first-residual -save-weights -save-uv '-maxuv-l 3150000
                ' '+uvselection+continue_clean+args.opt_args+' '
                ' -log-time -field 0 ' + quiet + update_model_option + ' ')
    opt_args = opt_args + shift_options

    # sub-band division for -channels-out
    if args.deconvolution_mode == 'good' and nc > 1:
        if args.channel_division == 'auto':
            try:
                division_args = channel_division.auto_division_args(
                    args.f, nc, mode=args.channel_division_mode,
                    spws=channel_division.parse_spws(opt_args),
                    field=channel_division.parse_field(opt_args),
                    csv_file=(root_dir_sys + base_name + '_' + g_name +
                              '_nc' + str(nc) + '_channel_division.csv'))
            except Exception as e:
                print(' !! channel_division failed, using WSClean default '
                      'division: ', e)
                division_args = ' '
        elif args.channel_division == 'gap':
            division_args = ' -gap-channel-division '
        elif args.channel_division in ('default', 'None', ''):
            division_args = ' '
        else:
            division_args = (' -channel-division-frequencies ' +
                             args.channel_division + ' ')
        deconvolver_args = deconvolver_args + division_args

    for robust in robusts:
        for uvtaper in tapers:
            if uvtaper == '':
                taper_mode = ''
            else:
                taper_mode = '-taper-gaussian '

            startTime = time.time()
            image_statistics = imaging(g_name=g_name,
                                       # base_name='2_selfcal_update_model',
                                       # base_name='image',
                                       base_name=base_name,
                                       field=field, robust=str(robust),
                                       uvtaper=uvtaper,
                                       continue_clean=args.continue_clean,
                                       nc=int(1*nc))

            exec_time = time.time() - startTime
            print(f" ++==>> Exec time cleaning = {exec_time:.1f} s")

            if image_statistics is not None:
                image_statistics['robust'] = robust
                # image_statistics['vwt'] = vwt
                image_statistics['uvtaper'] = uvtaper
                df = pd.DataFrame.from_dict(image_statistics, orient='index').T
                df.to_csv(root_dir_sys + image_statistics['imagename'].replace('.fits',
                                                                               '_data.csv'),
                          header=True,
                          index=False)

            else:
                pass