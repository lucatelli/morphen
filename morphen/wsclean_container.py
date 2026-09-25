#!/usr/bin/env python3
"""
Locate (and, if needed, download) the WSClean Singularity/Apptainer container
used by ph4ser.

Container images are ~1.5 GB, far above the 100 MB per-file limit GitHub
imposes on repositories, so they are not tracked by git. They are published as
release assets and fetched on demand into the clone, or into a shared cache
directory so that several clones share one copy.

Use as a library
----------------
    from wsclean_container import resolve_container
    sif = resolve_container()                      # default container
    sif = resolve_container(explicit=args.wsclean_sif)

Use from the command line
-------------------------
    python wsclean_container.py                    # download the default container
    python wsclean_container.py --list             # show known containers
    python wsclean_container.py --where            # print the resolved path and exit
    python wsclean_container.py --name wsclean37_cpu_native.sif
    python wsclean_container.py --dest ~/.ph4ser/containers
    python wsclean_container.py --url https://.../custom.sif --no-verify
"""

import argparse
import hashlib
import os
import shutil
import subprocess
import sys


# ---------------------------------------------------------------------------
# Known containers
# ---------------------------------------------------------------------------
# `url`     : direct download URL; leave '' until the image has been published.
# `sha256`  : checksum used to verify a download; None disables verification.
# `size`    : expected size in bytes, used only to report download progress.
#
# Containers are published as GitHub release assets (free, no bandwidth cap on
# public repos, 2 GiB per asset). To add one, upload it to a release -- either
# by drag-and-drop on the releases page, or:
#     gh release create v0.13 --title 'ph4ser v0.13' --notes '...'
#     gh release upload v0.13 ph4ser/<image>.sif
# the URL is then
#     https://github.com/lucatelli/ph4ser/releases/download/<tag>/<image>.sif
# and its checksum comes from `sha256sum <image>.sif`.
# ---------------------------------------------------------------------------

DEFAULT_CONTAINER = 'wsclean36_cpu_portable_emerlin.sif'

CONTAINERS = {
    'wsclean36_cpu_portable_emerlin.sif': {
        'url': 'https://github.com/lucatelli/ph4ser/releases/download/'
               'v0.13/wsclean36_cpu_portable_emerlin.sif',
        'sha256': '54268390b2de7e11072529fe04cbab0fe30d7def90d334efee4c551c289b8373',
        'size': 1631297536,
        'description': 'WSClean 3.6, CPU, portable build (IDG + EveryBeam, e-MERLIN ready).',
    },
    'wsclean36_gpu_portable_emerlin.sif': {
        'url': '',
        'sha256': None,
        'size': None,
        'description': 'WSClean 3.6, GPU, portable build (IDG + EveryBeam, e-MERLIN ready).',
    },
    'wsclean37_cpu_native.sif': {
        'url': '',
        'sha256': None,
        'size': None,
        'description': 'WSClean 3.7, CPU, native build (best performance, host-specific).',
    },
}

# Shared cache, so that several clones of the repository do not each hold their
# own 1.5 GB copy of the same image.
CACHE_DIR = os.path.join(os.path.expanduser('~'), '.ph4ser', 'containers')

# Directory holding this file (and, by default, the container next to it).
MODULE_DIR = os.path.dirname(os.path.abspath(__file__))


# ---------------------------------------------------------------------------
# Resolution
# ---------------------------------------------------------------------------

def container_search_path(name=DEFAULT_CONTAINER):
    """
    Directories/files searched for a container, in order of precedence.

    Args:
        name: Container file name (not a path).

    Returns:
        list: List of ``(path, origin)`` tuples, most specific first, where
        `origin` is a short human-readable description of where the candidate
        came from (used in error messages).
    """
    candidates = []

    env_sif = os.environ.get('PH4SER_WSCLEAN_SIF')
    if env_sif:
        candidates.append((os.path.abspath(os.path.expanduser(env_sif)),
                           'environment variable PH4SER_WSCLEAN_SIF'))

    env_dir = os.environ.get('PH4SER_CONTAINER_DIR')
    if env_dir:
        candidates.append((os.path.join(os.path.abspath(os.path.expanduser(env_dir)), name),
                           'environment variable PH4SER_CONTAINER_DIR'))

    # ph4ser_config.py may define general_settings['wsclean_sif']
    try:
        import ph4ser_config as _cf
        cfg_sif = getattr(_cf, 'general_settings', {}).get('wsclean_sif')
        if cfg_sif:
            candidates.append((os.path.abspath(os.path.expanduser(cfg_sif)),
                               "ph4ser_config general_settings['wsclean_sif']"))
    except Exception:
        pass

    # next to this file -- the layout used during development
    candidates.append((os.path.join(MODULE_DIR, name), 'ph4ser module directory'))
    # shared cache
    candidates.append((os.path.join(CACHE_DIR, name), 'shared cache ~/.ph4ser/containers'))

    return candidates


def auto_fetch_enabled():
    """
    Whether a missing container may be downloaded automatically.

    Controlled by ``$PH4SER_AUTO_FETCH``: set it to 0/false/no/off to require
    that containers be installed by hand.

    Returns:
        bool: True unless auto-fetching has been switched off.
    """
    return os.environ.get('PH4SER_AUTO_FETCH', '1').strip().lower() \
        not in ('0', 'false', 'no', 'off')


def resolve_container(name=DEFAULT_CONTAINER, explicit=None, required=True,
                      verbose=True, auto_fetch=True):
    """
    Return the path to a usable WSClean container, downloading it if needed.

    The first existing candidate wins, searched in this order:

      1. `explicit` (e.g. the --wsclean_sif command-line argument)
      2. ``$PH4SER_WSCLEAN_SIF``          (full path to an image)
      3. ``$PH4SER_CONTAINER_DIR``        (directory holding images)
      4. ``ph4ser_config.general_settings['wsclean_sif']``
      5. the ph4ser module directory      (the development layout)
      6. ``~/.ph4ser/containers``         (shared cache)

    If none of those exist and a download URL is registered for `name`, the
    image is fetched automatically, so that a fresh clone images without any
    manual setup step. The download is resumable and checksum-verified; set
    ``$PH4SER_AUTO_FETCH=0`` to turn it off and get the old "not found" error
    with installation instructions instead.

    Args:
        name: Container file name to look for.
        explicit: Path given explicitly by the caller; takes precedence over
            everything else and is an error if it does not exist.
        required: If True, raise when nothing is found and nothing could be
            downloaded; if False, return None.
        verbose: If True, report which container was selected.
        auto_fetch: If False, never download (used by the informational
            command-line options).

    Returns:
        str: Absolute path to the container, or None if `required` is False and
        nothing was found.

    Raises:
        FileNotFoundError: If nothing is found or downloaded and `required` is
            True, or if `explicit` was given but does not exist.
    """
    if explicit:
        path = os.path.abspath(os.path.expanduser(explicit))
        if not os.path.isfile(path):
            raise FileNotFoundError(
                f'WSClean container not found at the requested path: {path}')
        if verbose:
            print(f' >> Using WSClean container: {path} (explicitly requested)')
        return path

    for path, origin in container_search_path(name):
        if os.path.isfile(path):
            if verbose:
                print(f' >> Using WSClean container: {path} (from {origin})')
            return path

    # nothing installed: fetch it, unless the caller or the user said not to
    url = CONTAINERS.get(name, {}).get('url')
    blocked = None
    if not (auto_fetch and required):
        blocked = 'not requested by the caller'
    elif not auto_fetch_enabled():
        blocked = 'disabled by PH4SER_AUTO_FETCH'
    elif not url:
        blocked = f'no download URL is registered for "{name}"'

    if blocked is None:
        size = CONTAINERS.get(name, {}).get('size')
        print(f'++==> WSClean container "{name}" is not installed yet; '
              f'downloading it now'
              + (f' ({size / 1024 ** 3:.2f} GiB).' if size else '.'))
        print('      This happens once. Set PH4SER_AUTO_FETCH=0 to disable, or '
              'point PH4SER_WSCLEAN_SIF at an existing image.')
        try:
            return fetch_container(name=name)
        except Exception as exc:
            print(f'     !!==> Automatic download failed: {exc}')
            blocked = 'the automatic download failed'

    if not required:
        return None

    searched = '\n'.join(f'      - {path}   [{origin}]'
                         for path, origin in container_search_path(name))
    raise FileNotFoundError(
        f'\nWSClean container "{name}" not found, and it was not downloaded '
        f'({blocked}).\n'
        f'   Searched:\n{searched}\n\n'
        f'   Container images are not distributed through git (they exceed the\n'
        f'   100 MB GitHub file limit). Fetch the image with:\n\n'
        f'       python {os.path.join(MODULE_DIR, "wsclean_container.py")}\n\n'
        f'   or point ph4ser at an existing image:\n\n'
        f'       export PH4SER_WSCLEAN_SIF=/path/to/{name}\n')


# ---------------------------------------------------------------------------
# Download
# ---------------------------------------------------------------------------

def sha256sum(path, chunk_size=1 << 20, verbose=True):
    """
    SHA-256 checksum of a file, read in chunks so that large images do not
    have to be held in memory.

    Args:
        path: File to hash.
        chunk_size: Read size in bytes.
        verbose: If True, print a short progress note.

    Returns:
        str: Hexadecimal digest.
    """
    if verbose:
        print(f' >> Verifying checksum of {os.path.basename(path)} ...')
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(chunk_size), b''):
            digest.update(chunk)
    return digest.hexdigest()


def fetch_container(name=DEFAULT_CONTAINER, dest=None, url=None, verify=True,
                    force=False):
    """
    Download a container image and verify it.

    The download is resumable: an interrupted transfer is continued rather than
    restarted (``curl -C -`` / ``wget -c``). The file is written to a
    ``.part`` temporary name and only moved into place once the checksum
    matches, so an aborted or corrupted download can never masquerade as a
    usable container.

    Args:
        name: Container file name (a key of CONTAINERS, or any name if `url`
            is given explicitly).
        dest: Destination directory; defaults to the ph4ser module directory,
            falling back to the shared cache if that is not writable.
        url: Download URL; defaults to the registered URL for `name`.
        verify: If True, check the SHA-256 against the registered value.
        force: If True, re-download even if the container is already present.

    Returns:
        str: Absolute path to the downloaded container.

    Raises:
        ValueError: If no URL is known for `name`.
        RuntimeError: If the download fails or the checksum does not match.
    """
    entry = CONTAINERS.get(name, {})
    url = url or entry.get('url')

    if not url:
        raise ValueError(
            f'\nNo download URL is registered for "{name}".\n'
            f'   Edit CONTAINERS in {os.path.abspath(__file__)} and set its "url",\n'
            f'   or pass one explicitly:\n\n'
            f'       python wsclean_container.py --name {name} --url <URL>\n')

    if dest is None:
        dest = MODULE_DIR if os.access(MODULE_DIR, os.W_OK) else CACHE_DIR
    dest = os.path.abspath(os.path.expanduser(dest))
    os.makedirs(dest, exist_ok=True)

    target = os.path.join(dest, name)
    if os.path.isfile(target) and not force:
        print(f' >> Container already present: {target}')
        print('    (use --force to download it again)')
        return target

    partial = target + '.part'
    expected = entry.get('sha256') if verify else None
    size = entry.get('size')

    print(f' >> Downloading {name}')
    print(f'    from {url}')
    print(f'    to   {target}')
    if size:
        note = ' -- this will take a while' if size > 1024 ** 3 else ''
        print(f'    ({size / 1024 ** 3:.2f} GiB{note})')

    if shutil.which('curl'):
        command = ['curl', '-L', '--fail', '-C', '-', '-o', partial, url]
    elif shutil.which('wget'):
        command = ['wget', '-c', '-O', partial, url]
    else:
        raise RuntimeError('Neither curl nor wget is available to download the container.')

    try:
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(
            f'Download failed ({exc}). The partial file was kept at {partial}; '
            f'run the command again to resume it.')

    if expected:
        found = sha256sum(partial)
        if found != expected:
            raise RuntimeError(
                f'\nChecksum mismatch for {name}:\n'
                f'    expected {expected}\n'
                f'    found    {found}\n'
                f'   The partial file was kept at {partial} for inspection; '
                f'delete it and retry.')
        print('    Checksum OK.')
    else:
        print('    No checksum registered for this container; skipping verification.')

    os.replace(partial, target)
    print(f' >> Container ready: {target}')
    return target


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Locate or download the WSClean container used by ph4ser.',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--name', type=str, default=DEFAULT_CONTAINER,
                        help=f'Container file name (default: {DEFAULT_CONTAINER}).')
    parser.add_argument('--dest', type=str, default=None,
                        help='Destination directory (default: the ph4ser module '
                             'directory, or ~/.ph4ser/containers if it is not writable).')
    parser.add_argument('--url', type=str, default=None,
                        help='Download URL, overriding the registered one.')
    parser.add_argument('--no-verify', action='store_true',
                        help='Skip SHA-256 verification of the download.')
    parser.add_argument('--force', action='store_true',
                        help='Download again even if the container is already present.')
    parser.add_argument('--where', action='store_true',
                        help='Print the path of the container that would be used, then exit.')
    parser.add_argument('--list', action='store_true',
                        help='List the known containers, then exit.')
    args = parser.parse_args()

    if args.list:
        print('Known containers:\n')
        for name, entry in CONTAINERS.items():
            marker = ' (default)' if name == DEFAULT_CONTAINER else ''
            local = resolve_container(name, required=False, verbose=False)
            print(f'  {name}{marker}')
            print(f'      {entry["description"]}')
            print(f'      url    : {entry["url"] or "(not published yet)"}')
            print(f'      local  : {local or "(not present)"}\n')
        return 0

    if args.where:
        try:
            print(resolve_container(args.name, verbose=False, auto_fetch=False))
            return 0
        except FileNotFoundError as exc:
            print(exc, file=sys.stderr)
            return 1

    existing = resolve_container(args.name, required=False, verbose=False)
    if existing and not args.force:
        print(f' >> Container already available: {existing}')
        print('    (use --force to download it again)')
        return 0

    try:
        fetch_container(name=args.name, dest=args.dest, url=args.url,
                        verify=not args.no_verify, force=args.force)
    except (ValueError, RuntimeError) as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
