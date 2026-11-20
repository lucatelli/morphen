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


def resolve_ms_path(path):
    """
    Full path to the measurement set named by --f.

    Pipeline callers usually pass --f *without* the .ms extension (the script
    appends '.ms' itself when it builds the WSClean command), while a manual
    run usually includes it. Anything that has to open the MS directly -- the
    channel division, for one -- needs the real directory, so accept both
    spellings.

    Args:
        path: The value of --f, with or without the .ms extension.

    Returns:
        str: The first candidate that exists on disk; otherwise `path` with
        '.ms' appended, which keeps the error message informative.
    """
    candidates = [path] if path.endswith('.ms') else [path + '.ms', path]
    for candidate in candidates:
        if os.path.isdir(candidate):
            return candidate
    return candidates[0]


def parse_ms_list(value):
    """
    Split the value of --f into a list of measurement sets.

    A single path (no comma) is returned untouched, so the single-MS code path
    sees exactly the string it always did. A comma-separated list, with or
    without brackets (same style as --r/--t), is split and stripped; exact
    duplicates are dropped, keeping the first occurrence.

    Args:
        value: The value of --f.

    Returns:
        list of str: One entry per MS, as given on the command line.
    """
    if ',' not in value:
        return [value]
    items = [s for s in parse_str_list(value) if s != '']
    return list(dict.fromkeys(items))


def ms_stem(path):
    """Basename of an MS path without the trailing '.ms'."""
    name = os.path.basename(os.path.normpath(path))
    return name[:-3] if name.endswith('.ms') else name


def default_multi_ms_name(ms_paths, max_len=80):
    """
    Tag used in the image name when several MS are imaged together.

    The MS stems joined with '+' (e.g. 'obs1+obs2'); if that gets longer than
    `max_len`, '<first_stem>+<N-1>more' instead.
    """
    stems = [ms_stem(p) for p in ms_paths]
    name = '+'.join(stems)
    if len(name) > max_len or len(set(stems)) < len(stems):
        name = stems[0] + '+' + str(len(stems) - 1) + 'more'
    return name


def build_container_mounts(out_dir, ms_paths, mask=None):
    """
    Bind mounts needed to see the output dir, every MS and the mask inside the
    container.

    The output directory is mounted on /mnt, as in the single-MS case, so the
    image name, the mask and --continue_clean keep resolving there. Each other
    directory that holds an MS gets its own top-level target (/mnt_ms0,
    /mnt_ms1, ...), so MS with the same basename in different directories do
    not collide. A mask outside every mounted directory goes on /mnt_mask.

    MS that share a basename (e.g. several '1344+5553.calibrated.ms' from
    different epochs) must also get distinct names inside the container:
    with -temp-dir, WSClean names its reordered files after the MS basename
    only, so same-named MS overwrite each other's temp files and WSClean
    crashes (SIGSEGV) when reading them back. Each repeated one is therefore
    bound by itself under a new name, '<stem>_ms<i>.ms'.

    Args:
        out_dir: Host directory where WSClean writes the images.
        ms_paths: Absolute host paths of the measurement sets.
        mask: Host path of the FITS mask, or None.

    Returns:
        tuple: (binds, container_ms_paths, container_mask), where `binds` is a
        list of 'host:target' strings for singularity --bind and
        `container_mask` is None when no mask is given.
    """
    out_dir = os.path.abspath(out_dir)
    targets = {out_dir: '/mnt'}
    binds = [out_dir + ':/mnt']
    container_ms_paths = []
    seen_names = set()
    for i, ms in enumerate(ms_paths):
        if os.path.basename(ms) in seen_names:
            target = ('/mnt_ms' + str(len(targets) - 1) + '/' + ms_stem(ms) +
                      '_ms' + str(i) + '.ms')
            targets[ms] = target
            binds.append(ms + ':' + target)
            container_ms_paths.append(target)
            continue
        seen_names.add(os.path.basename(ms))
        host_dir = os.path.dirname(ms)
        if host_dir not in targets:
            targets[host_dir] = '/mnt_ms' + str(len(targets) - 1)
            binds.append(host_dir + ':' + targets[host_dir])
        container_ms_paths.append(targets[host_dir] + '/' + os.path.basename(ms))

    container_mask = None
    if mask is not None:
        mask_dir = os.path.dirname(os.path.abspath(mask))
        if mask_dir not in targets:
            targets[mask_dir] = '/mnt_mask'
            binds.append(mask_dir + ':/mnt_mask')
        container_mask = targets[mask_dir] + '/' + os.path.basename(mask)
    return binds, container_ms_paths, container_mask


def imaging(g_name, field, uvtaper, robust, base_name='clean_image',
            continue_clean='False',nc=4, vis_list=None):
    g_vis = g_name + '.ms'
    # vis_list: MS paths as WSClean must see them (container paths when
    # running in singularity). None keeps the single-MS behaviour.
    if vis_list is None:
        vis_args = root_dir + g_vis
    else:
        vis_args = ' '.join(vis_list)
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
                ' ' + opt_args + ' ' + data_column + ' ' + vis_args)
            print(' ++==>> Command to be executed by WSClean: ')
            print(command_exec)
            if args.dry_run == 'True':
                print(' >> --dry_run True: WSClean not executed.')
                return (None)
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
        #     ' ' + opt_args + ' ' + data_column + ' ' + vis_args)

        #     print(' ++==>> Command to be executed by Singularity > WSClean: ')
        #     print(command_exec)
        #     os.system(command_exec)

        # Define a scratch dir on the host
        # scratch_dir = '/media/sagauga/starbyte/wsclean_temp/'
        # scratch_dir = f'/media/sagauga/starbyte/wsclean_temp/{uuid.uuid4().hex}/'
        # scratch_dir = f'/home/lucatelli/astronomical_data/wsclean_temp/{uuid.uuid4().hex}/'
        if os.path.exists('/media/sagauga/temp_wsclean/'):
            scratch_dir = f'/media/sagauga/temp_wsclean/wsclean_temp/{uuid.uuid4().hex}/'
        elif os.path.exists('/dev/shm/wsclean_temp/'):
            scratch_dir = f'/dev/shm/wsclean_temp/{uuid.uuid4().hex}/'
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
        # scratch_dir = f'/dev/shm/wsclean_temp/{uuid.uuid4().hex}/'
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
            ' ' + opt_args + ' ' + data_column + ' ' + vis_args)
            print(' ++==>> Command to be executed by Singularity > WSClean: ')
            print(command_exec)
            if args.dry_run == 'True':
                print(' >> --dry_run True: WSClean not executed.')
                shutil.rmtree(scratch_dir, ignore_errors=True)
                return (None)
            try:
                subprocess.run(shlex.split(command_exec), check=True)
            finally:
                # also when WSClean fails: its temp files can be several GB
                shutil.rmtree(scratch_dir, ignore_errors=True)

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


def _predict_scratch_dir():
    """Host scratch directory for WSClean temporary files (as in imaging())."""
    for base in ('/media/sagauga/temp_wsclean/', '/dev/shm/wsclean_temp/',
                 '/mnt/scratch/lucatelli/temp_wsclean/',
                 '/home/lucatelli/astronomical_data/'):
        if os.path.exists(base):
            return f'{base}wsclean_temp/{uuid.uuid4().hex}/'
    return f'./wsclean_temp/{uuid.uuid4().hex}/'


def predict(prefix, g_name, nc=4, division_args=' ', vis_list=None):
    """
    Predict visibilities from existing model images (wsclean -predict).

    The model images must be in the output directory and be named as WSClean
    names them in a normal imaging run: ``<prefix>-NNNN-model.fits`` for
    ``nc`` > 1 (``<prefix>-model.fits`` for ``nc`` = 1). No imaging or
    cleaning is done; the predicted visibilities are written to MODEL_DATA.

    ``nc`` and ``division_args`` must be the ones of the run that made the
    model images, so that each model image is applied to the same channels.
    """
    if vis_list is None:
        vis_args = root_dir + g_name + '.ms'
    else:
        vis_args = ' '.join(vis_list)
    if nc > 16:
        nc = 16
    predict_args = (' -predict -name ' + root_dir + prefix +
                    ' -channels-out ' + str(nc) + ' ' + division_args +
                    ' -gridder wgridder -wgridder-accuracy 1e-5 '
                    '-parallel-gridding 8 -padding 1.3 '
                    '-j ' + str(threads_proc) + ' -log-time -field 0 ' +
                    quiet + ' ')
    scratch_dir = _predict_scratch_dir()
    os.makedirs(scratch_dir, exist_ok=True)
    if running_container == 'native':
        command_exec = ('wsclean' + predict_args + '-temp-dir ' + scratch_dir +
                        ' ' + vis_args)
    else:
        command_exec = ('singularity exec --nv '
                        '--bind ' + scratch_dir + ':/wsclean_tmp ' +
                        '--bind ' + mount_dir + ' ' + wsclean_dir +
                        ' wsclean' + predict_args +
                        '-temp-dir /wsclean_tmp ' + vis_args)
    print(' ++==>> Command to be executed by WSClean (predict): ')
    print(command_exec)
    try:
        if args.dry_run == 'True':
            print(' >> --dry_run True: WSClean not executed.')
            return
        subprocess.run(shlex.split(command_exec), check=True)
    finally:
        shutil.rmtree(scratch_dir, ignore_errors=True)


def parse_float_list(str_values):
    return [float(val.strip()) for val in str_values.strip('[]').split(',')]


def parse_str_list(str_values):
    return [s.strip() for s in str_values.strip('[]').split(',')]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Helper for wsclean imaging.")
    parser.add_argument("--p", type=str, help="The path to the MS file.")
    parser.add_argument("--f", nargs='?', default=False,
                        const=True,
                        help="The name of the ms file (the .ms extension is "
                             "optional). Several MS can be imaged together as a "
                             "comma-separated list, e.g. "
                             "--f /data/a/obs1.ms,/data/b/obs2.ms (brackets are "
                             "allowed). All of them must have the --data column, "
                             "and -field 0 is applied to each one. With more than "
                             "one MS, --channel_division auto is not available "
                             "and falls back to WSClean's default division.")

    parser.add_argument("--outdir", type=str, nargs='?', default=None,
                        help="Directory for the images and the CSV file. "
                             "Default: the directory of the (first) MS.")

    parser.add_argument("--name", type=str, nargs='?', default=None,
                        help="Tag used in the image names in place of the MS "
                             "name. Default with several MS: the MS names "
                             "joined with '+' (or '<first>+<N-1>more' if "
                             "that is too long).")

    parser.add_argument("--dry_run", type=str, nargs='?', default='False',
                        help="If 'True', print the WSClean command (including "
                             "the container mounts) without running it.")
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
    parser.add_argument("--maxmscales", type=str, nargs='?', default="6",
                        help="Maximum number of scales to be used with the multiscale deconvolver in WSClean. ")

    parser.add_argument("--sx", type=str, nargs='?', default='2048',
                        help="Image Size x-axis")
    parser.add_argument("--sy", type=str, nargs='?', default='2048',
                        help="Image Size y-axis")
    parser.add_argument("--cellsize", type=str, nargs='?', default='0.05asec',
                        help="Cell size")
    parser.add_argument("--niter", type=str, nargs='?', default='200000',
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

    parser.add_argument("--channel_division", type=str, nargs='?', default='default',
                        help="How the bandwidth is divided into --nc sub-bands. "
                             "'auto': split frequencies computed from the MS "
                             "(gaps + region where all instruments overlap; "
                             "see channel_division.py). 'default': WSClean's "
                             "own division. 'gap': -gap-channel-division. "
                             "Anything else is passed as a comma-separated "
                             "list of split frequencies in Hz to "
                             "-channel-division-frequencies.")

    parser.add_argument("--channel_division_mode", type=str, nargs='?', default='weight',
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


    parser.add_argument("--predict", type=str, nargs='?', default=None,
                        help="Predict only (wsclean -predict): write to "
                             "MODEL_DATA the visibilities of the model images "
                             "<outdir>/<PREDICT>-NNNN-model.fits "
                             "(<PREDICT>-model.fits with --nc 1). No imaging "
                             "is done. --nc and --channel_division must be "
                             "the ones used to make the model images.")

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

    ms_inputs = parse_ms_list(args.f)
    # A single MS with no --outdir/--name keeps the original code path, so
    # existing callers (e.g. radio_utils.run_wsclean) see no change at all.
    legacy_single = (len(ms_inputs) == 1 and args.outdir is None
                     and args.name is None)
    mask_given = not (args.mask == 'None' or args.mask is None)

    if legacy_single:
        f_path = ms_inputs[0]
        field = os.path.basename(f_path).replace('.ms', '')
        g_name = field
        root_dir_sys = os.path.dirname(f_path) + '/'
        vis_path = resolve_ms_path(f_path)   # --f may omit the .ms extension
        ms_paths = [vis_path]
        vis_list = None
    else:
        ms_paths = [os.path.abspath(resolve_ms_path(p)) for p in ms_inputs]
        ms_paths = list(dict.fromkeys(ms_paths))
        missing = [p for p in ms_paths if not os.path.isdir(p)]
        if missing:
            sys.exit(' !! Measurement set(s) not found: ' + ', '.join(missing))
        if mask_given and not os.path.isfile(args.mask):
            sys.exit(' !! Mask not found: ' + args.mask)
        if args.name is not None:
            g_name = args.name
        elif len(ms_paths) == 1:
            g_name = ms_stem(ms_paths[0])
        else:
            g_name = default_multi_ms_name(ms_paths)
        field = g_name
        out_dir = (os.path.abspath(args.outdir) if args.outdir is not None
                   else os.path.dirname(ms_paths[0]))
        os.makedirs(out_dir, exist_ok=True)
        root_dir_sys = out_dir + '/'
        vis_path = ms_paths[0]
    multi_ms = len(ms_paths) > 1
    robusts = args.r
    tapers = args.t

    if running_container == 'singularity':
        if legacy_single:
            bind_list = [root_dir_sys + ':/mnt']
        else:
            bind_list, vis_list, mask_in_container = build_container_mounts(
                root_dir_sys, ms_paths,
                mask=args.mask if mask_given else None)
        # Joined so that the '--bind ' + mount_dir in imaging() expands to one
        # --bind flag per directory.
        mount_dir = ' --bind '.join(bind_list)
        root_dir = '/mnt/'
        # The image is not shipped with the repository: resolve it from the
        # command line, the environment, the config, the module directory or
        # the shared cache. See wsclean_container.py.
        wsclean_dir = resolve_container(DEFAULT_CONTAINER,
                                        explicit=args.wsclean_sif)
        print('Using WSClean from singularity image: ', wsclean_dir)
        print('WSClean version from singularity image: ')
        bind_flags = [x for b in bind_list for x in ("--bind", b)]
        subprocess.run([
            "singularity", "exec", "--nv",
            # "--bind", "/usr/local/cuda-13.0:/usr/local/cuda-13.0",
            # "--bind", "/usr/local/cuda-13.0/lib64:/usr/local/cuda/lib64",
            *bind_flags,
            wsclean_dir, "wsclean", "--version"
        ], check=True)
    if running_container == 'native':
        mount_dir = ''
        root_dir = root_dir_sys
        os.environ['OPENBLAS_NUM_THREADS'] = '1'
        if not legacy_single:
            vis_list = ms_paths
            mask_in_container = args.mask if mask_given else None

    if not legacy_single:
        print(' >> Imaging ' + str(len(ms_paths)) + ' measurement set(s):')
        for host_ms, wsclean_ms in zip(ms_paths, vis_list):
            print('      ' + host_ms + '  ->  ' + wsclean_ms)
        print(' >> Output directory: ' + root_dir_sys)
        print(' >> Image name tag: ' + g_name)

    base_name = args.save_basename

    ## Setting image and deconvolution noptions.
    ### Cleaning arguments
    auto_mask = ' -auto-mask ' + args.nsigma_automask
    # auto_mask = ' '
    auto_threshold = ' -auto-threshold ' + args.nsigma_autothreshold
    if args.mask == 'None' or args.mask == None:
        mask_file = ' '
    elif not legacy_single:
        mask_file = ' -fits-mask ' + mask_in_container + ' '
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
            deconvolver_options = ' -multiscale-scale-bias 0.75 -multiscale-gain 0.1 -multiscale-convolution-padding 1.3 -padding 1.3 -clean-border 1 '
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
                            '-parallel-deconvolution 2048 '
                            # '-parallel-deconvolution 1024 '
                            '-channels-out '+str(nc)+' -join-channels ' + negative_arg + ' '
                            '-no-mf-weighting ' # this must be set to generate science images 
                            '-fit-spectral-pol  ' +str(2)+' '+' -deconvolution-channels ' +str(dec_chan)+' '
                            # '-fit-spectral-log-pol ' +str(2)+' '+' -deconvolution-channels ' +str(dec_chan)+' '
                            '-gridder wgridder -wgridder-accuracy 1e-5 -parallel-gridding 8 '
                            # '-apply-primary-beam -circular-beam '
                            #updates 2026
                            '-weighting-rank-filter 3 -weighting-rank-filter-size 128 '
                            # '-weighting-rank-filter 3 -weighting-rank-filter-size 128 -local-rms -local-rms-window 25 -local-rms-strength 0.5 '
                            # '-gridder idg -idg-mode hybrid -grid-with-beam -circular-beam -save-psf-pb -save-uv '
                            #cleanup 2026
                            # '-wstack-nwlayers-factor 3 ' +' -deconvolution-channels ' +str(nc)
                            # '-fit-spectral-pol  ' +str(4)+' '
                            # '-save-psf-pb -save-weights -use-weights-as-taper '
                            # Other arguments that may be useful
                            # '-store-imaging-weights -save-weights '
                            # '-no-negative -abs-threshold 3.5e-6 '
                            # Some testing parameters
                            # '-save-weights -local-rms -local-rms-window 100 -local-rms-strength 0.5 '
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
    gain_args = ' -mgain 0.60 -gain 0.1 -nmiter 30 ' # -super-weight 9.0

    if args.shift == 'None' or args.shift == None:
        # if args.shift != ' ':
        shift_options = ' '
    else:
        shift_options = ' -shift ' + args.shift + ' '
    # shift_options = ' '  # -shift 13:15:30.68  +62.07.45.357 '#' -shift 13:15:28.903
    # +62.07.11.886 '
    if args.quiet == 'True':
        quiet = ' -quiet '
    else:
        quiet = ' -v '

    # quiet = ' -v ' #debug override

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
    division_args = ' '
    if args.deconvolution_mode == 'good' and nc > 1:
        if args.channel_division == 'auto' and multi_ms:
            # channel_division reads a single MS; not supported for a list yet.
            print(' !! --channel_division auto is not available with several '
                  'MS, using WSClean default division.')
            division_args = ' '
        elif args.channel_division == 'auto':
            try:
                division_args = channel_division.auto_division_args(
                    vis_path, nc, mode=args.channel_division_mode,
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

    if args.predict not in (None, 'None'):
        # Predict only: no imaging loop.
        startTime = time.time()
        predict(args.predict, g_name, nc=int(nc), division_args=division_args,
                vis_list=vis_list)
        print(f" ++==>> Exec time predict = {time.time() - startTime:.1f} s")
        sys.exit(0)

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
                                       nc=int(1*nc),
                                       vis_list=vis_list)

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