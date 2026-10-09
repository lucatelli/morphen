# CASA site configuration for the morphen Singularity images.
#
# Repo name is casasiteconfig_morphen.py (installed as /opt/casa/casasiteconfig.py)
# because casaconfig also imports a module named `casasiteconfig` from sys.path:
# under its real name, any python started in this directory would load it.
#
# casaconfig loads /opt/casa/casasiteconfig.py automatically, before the user's
# ~/.casa/config.py; the user file then overrides only the keys it sets.
# This file is plain Python run at `import casatools`, so it can pick a data
# location at runtime instead of hard-coding one at build time.
#
# Two modes:
#
#   1. Default (read-only): casarundata + measures were pulled into
#      /opt/casa/data when the image was built. CASA reads them from there and
#      never tries to write, so it works with a read-only, empty, or fake $HOME
#      (--containall, --no-home, HPC). The IERS/ephemeris tables are then as old
#      as the image.
#
#   2. Updating measures: export MORPHEN_CASADATA=/some/writable/dir (usually a
#      bind mount, e.g. -B /scratch/me/casadata:/casadata). The dir becomes
#      measurespath with measures auto-updates on. If it is empty it is seeded
#      once (~30 MB, under a second) with a *measures overlay*: real, writable
#      copies of the tables the daily measures tarball replaces (geodetic/ and
#      ephemerides/DE200, DE405), and symlinks into /opt/casa/data for
#      everything else (alma, catalogs, nrao, Lines, Sources, JPL-Horizons...).
#      casaconfig then sees a normal casarundata tree and only refreshes the
#      measures (~27 MB download, at most once a day). casarundata itself stays
#      at the image's version (data_auto_update = False): it is pinned with
#      the CASA release anyway, and an update would try to write through the
#      symlinks into the read-only image.
#      Concurrent jobs sharing the dir are serialised by a seed lock here and
#      by casaconfig's own data_update.lock for the download. The dir must be
#      owned by the user running CASA (a casaconfig rule for auto-updates);
#      otherwise this falls back to mode 1 with a message.
#      Pointing MORPHEN_CASADATA at an existing full casarundata tree that
#      casaconfig already maintains (e.g. a host ~/.casa/data) also works: it
#      is used as-is, with measures updates only.
#
# A user ~/.casa/config.py that sets measurespath or *_auto_update still wins.

import os as _os
import sys as _sys
import tempfile as _tempfile

_BAKED = '/opt/casa/data'
# Tables shipped in the NRAO/ASTRON measures tarball (what measures_update
# rewrites); everything else in casarundata is static between CASA releases.
_MEASURES_EPHEM = ('DE200', 'DE405')
_SEED_LOCK = '.morphen-seed.lock'


def _say(msg):
    print('morphen: ' + msg, file=_sys.stderr)


def _user():
    # --containall/--cleanenv drop $USER; fall back to the passwd entry or uid.
    try:
        import getpass as _getpass
        return _getpass.getuser()
    except Exception:
        return str(_os.getuid())


def _writable_dir(path):
    try:
        _os.makedirs(path, exist_ok=True)
    except OSError:
        return False
    return _os.path.isdir(path) and _os.access(path, _os.W_OK | _os.X_OK)


def _scratch(sub):
    base = _os.environ.get('TMPDIR') or _tempfile.gettempdir()
    path = _os.path.join(base, '%s-morphen' % _user(), sub)
    _os.makedirs(path, exist_ok=True)
    return path


def _seed_overlay(rw):
    # Called with the seed lock held, on an empty rw or one holding only the
    # leftovers of an interrupted seed (copies/links are redone idempotently).
    import shutil as _shutil
    _shutil.copytree(_os.path.join(_BAKED, 'geodetic'), _os.path.join(rw, 'geodetic'),
                     symlinks=True, dirs_exist_ok=True)
    eph_rw = _os.path.join(rw, 'ephemerides')
    _os.makedirs(eph_rw, exist_ok=True)
    for name in sorted(_os.listdir(_os.path.join(_BAKED, 'ephemerides'))):
        src = _os.path.join(_BAKED, 'ephemerides', name)
        dst = _os.path.join(eph_rw, name)
        if name in _MEASURES_EPHEM:
            _shutil.copytree(src, dst, symlinks=True, dirs_exist_ok=True)
        elif not _os.path.lexists(dst):
            _os.symlink(src, dst)
    for name in sorted(_os.listdir(_BAKED)):
        if name in ('geodetic', 'ephemerides', 'readme.txt', 'data_update.lock'):
            continue
        dst = _os.path.join(rw, name)
        if not _os.path.lexists(dst):
            _os.symlink(_os.path.join(_BAKED, name), dst)
    # The casarundata readme goes in last: it is what marks the dir as seeded
    # (and what casaconfig checks), so an interrupted seed is redone next time.
    _shutil.copy2(_os.path.join(_BAKED, 'readme.txt'), _os.path.join(rw, 'readme.txt'))


def _prepare_rw(rw):
    # Returns True when rw can be used as an auto-updating measurespath.
    if _os.stat(rw).st_uid != _os.getuid():
        _say('MORPHEN_CASADATA=%s is not owned by %s, and casaconfig only '
             'auto-updates a measurespath owned by the user; using the read-only '
             'data baked into the image' % (rw, _user()))
        return False
    import fcntl as _fcntl
    with open(_os.path.join(rw, _SEED_LOCK), 'a') as lock:
        _fcntl.flock(lock, _fcntl.LOCK_EX)
        if _os.path.exists(_os.path.join(rw, 'readme.txt')):
            return True  # seeded earlier, or a casaconfig-maintained casarundata
        others = [f for f in _os.listdir(rw) if f not in (_SEED_LOCK, 'data_update.lock')]
        if not _os.path.isdir(_BAKED):
            _say('no %s in this image to seed MORPHEN_CASADATA from' % _BAKED)
            return False
        if others:
            # A leftover of an interrupted seed is ours to finish; anything
            # else is not casaconfig data and must not be touched.
            ours = set(others) <= set(_os.listdir(_BAKED))
            if not ours:
                _say('MORPHEN_CASADATA=%s is not empty and holds no casaconfig '
                     'readme.txt; using the read-only data baked into the image' % rw)
                return False
        _say('seeding %s with a measures overlay of %s (one-off, ~30 MB)' % (rw, _BAKED))
        _seed_overlay(rw)
    return True


_rw = _os.environ.get('MORPHEN_CASADATA', '').strip()
_rw_ok = False
if _rw:
    _rw = _os.path.abspath(_os.path.expanduser(_rw))
    if not _writable_dir(_rw):
        _say('MORPHEN_CASADATA=%s is not writable, using the read-only data '
             'baked into the image' % _rw)
    else:
        try:
            _rw_ok = _prepare_rw(_rw)
        except Exception as _exc:
            _say('could not prepare MORPHEN_CASADATA=%s (%s); using the read-only '
                 'data baked into the image' % (_rw, _exc))

# datapath is set explicitly (casaconfig would default it to [measurespath])
# so that the baked tree stays a fallback: if a user config points
# measurespath somewhere empty or missing, casatools still finds the IERS
# tables via datapath and imports, instead of raising ImportError.
if _rw_ok:
    measurespath = _rw
    measures_auto_update = True
    data_auto_update = False
    datapath = [_rw, _BAKED]
elif _os.path.isdir(_BAKED):
    measurespath = _BAKED
    measures_auto_update = False
    data_auto_update = False
    datapath = [_BAKED]
# else: no baked data in this image, keep casaconfig's defaults (~/.casa/data).

# cachedir (~/.casa) holds rc/telemetry files; casa logs default to the cwd.
# Redirect both when they are not writable. A read-only cwd does not crash
# casatasks (the log is just silently not written), but the log is lost.
# Guarded: an exception anywhere in this file makes casaconfig discard *all*
# of it (back to ~/.casa/data defaults), so a bad TMPDIR must not do that.
try:
    _cache = _os.path.expanduser('~/.casa')
    if not _writable_dir(_cache):
        cachedir = _scratch('casa')
    if not _os.access(_os.getcwd(), _os.W_OK):
        import time as _time
        logfile = _os.path.join(_scratch('casa'),
                                'casa-%s.log' % _time.strftime('%Y%m%d-%H%M%S', _time.gmtime()))
except OSError as _exc:
    _say('no writable scratch dir for CASA cache/logs (%s)' % _exc)
