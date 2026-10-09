"""Smoke test for the morphen Singularity images.

Run inside the container, from a writable working dir:

    singularity exec morphen.sif python /opt/morphen/tests/container_smoke.py

Each check prints PASS/FAIL; exit status is the number of failures. Flags:
    --quick     imports only (used by the images' %test section)
    --no-gui    skip the plotms/imview headless launch checks
    --gpu       require JAX to run on a GPU (GPU image, run with --nv)
"""
import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import time
import traceback

RESULTS = []


def check(name):
    def deco(fn):
        def run(*a, **k):
            t0 = time.time()
            try:
                info = fn(*a, **k)
                RESULTS.append((name, True, info))
                print('PASS  %-34s %6.1fs  %s' % (name, time.time() - t0, info or ''), flush=True)
            except Exception as exc:  # noqa: BLE001 - report every failure kind
                RESULTS.append((name, False, repr(exc)))
                print('FAIL  %-34s %6.1fs  %r' % (name, time.time() - t0, exc), flush=True)
                traceback.print_exc()
        return run
    return deco


@check('python / platform')
def t_python():
    return '%s %s' % (sys.version.split()[0], sys.executable)


@check('scientific stack imports')
def t_stack():
    import numpy, scipy, astropy, pandas, matplotlib, skimage, sklearn  # noqa: F401
    import photutils, sep, fitsio, lmfit, emcee, dynesty, petrofit, corner  # noqa: F401
    import bilby, arviz, h5py, datashader, image_registration, astroquery  # noqa: F401
    return 'numpy %s astropy %s photutils %s bilby %s' % (
        numpy.__version__, astropy.__version__, photutils.__version__, bilby.__version__)


@check('jax jit (default device)')
def t_jax():
    import jax
    import jax.numpy as jnp

    @jax.jit
    def sersic(r, n, re):
        bn = 2.0 * n - 1.0 / 3.0
        return jnp.exp(-bn * ((r / re) ** (1.0 / n) - 1.0))

    out = sersic(jnp.linspace(0.1, 5, 64), 2.0, 1.5).block_until_ready()
    assert out.shape == (64,)
    return 'jax %s devices=%s' % (jax.__version__, jax.devices())


@check('jax gpu (matmul + fft conv)')
def t_jax_gpu():
    import jax
    import jax.numpy as jnp
    assert jax.default_backend() == 'gpu', 'backend is %s (run with --nv?)' % jax.default_backend()
    k = jax.random.PRNGKey(0)
    a = jax.random.normal(k, (2048, 2048))
    conv = jax.jit(lambda x: jnp.real(jnp.fft.ifft2(jnp.fft.fft2(x) * jnp.fft.fft2(x.T))))
    m = jax.jit(lambda x: (x @ x.T).sum())
    conv(a).block_until_ready()
    m(a).block_until_ready()
    t0 = time.time()
    conv(a).block_until_ready()
    m(a).block_until_ready()
    return '%s, %.1f ms' % (jax.devices()[0].device_kind, (time.time() - t0) * 1e3)


@check('casatools import (no writes)')
def t_casatools():
    import casaconfig.config as cfg
    import casatools
    return 'measurespath=%s auto=%s/%s loaded=%s' % (
        cfg.measurespath, cfg.measures_auto_update, cfg.data_auto_update,
        [os.path.basename(f) for f in cfg.load_success()])


@check('casa measures (IERS frame conv.)')
def t_measures():
    import casatools
    me = casatools.measures()
    qa = casatools.quanta()
    me.doframe(me.observatory('VLA'))
    me.doframe(me.epoch('utc', qa.quantity('2024-01-01T00:00:00')))
    d = me.measure(me.direction('J2000', '12h00m00', '45d00m00'), 'AZEL')
    ut1 = me.measure(me.epoch('utc', 'today'), 'ut1')
    return 'az=%.4f rad, ut1 ok=%s' % (d['m0']['value'], bool(ut1))


@check('casatasks imstat/imhead/export')
def t_casatasks(tmp):
    import numpy as np
    from astropy.io import fits
    import casatasks
    fitsfile = os.path.join(tmp, 'img.fits')
    data = np.random.default_rng(0).normal(0, 1e-4, (128, 128)).astype('f4')
    data[60:68, 60:68] += 1e-2
    hdr = fits.Header()
    hdr.update(CTYPE1='RA---SIN', CTYPE2='DEC--SIN', CRVAL1=180.0, CRVAL2=45.0,
               CDELT1=-1e-4, CDELT2=1e-4, CRPIX1=64, CRPIX2=64, BUNIT='Jy/beam',
               BMAJ=3e-4, BMIN=3e-4, BPA=0.0, CUNIT1='deg', CUNIT2='deg')
    fits.writeto(fitsfile, data, hdr, overwrite=True)
    im = os.path.join(tmp, 'img.image')
    casatasks.importfits(fitsimage=fitsfile, imagename=im, overwrite=True)
    st = casatasks.imstat(imagename=im)
    casatasks.imhead(imagename=im, mode='list')
    casatasks.exportfits(imagename=im, fitsimage=os.path.join(tmp, 'rt.fits'), overwrite=True)
    return 'max=%.3g' % st['max'][0]


@check('morphen: import mlibs + mp')
def t_morphen():
    import mlibs
    import morphen as mp
    assert hasattr(mlibs, 'compute_image_properties')
    assert hasattr(mp, 'source_extraction')
    return 'mlibs from %s' % os.path.dirname(mlibs.__file__)


@check('astropy IERS / cache dir')
def t_astropy():
    from astropy.config import get_cache_dir
    from astropy.time import Time
    t = Time('2024-06-01')
    return 'ut1-utc=%s cache=%s' % (t.delta_ut1_utc, get_cache_dir())


def _launch_gui(app_path, args, seconds=8):
    # xvfb-run when available; otherwise the morphen AppImage shim renders
    # offscreen by itself (QT_QPA_PLATFORM=offscreen when there is no X).
    cmd = ([app_path] if shutil.which('xvfb-run') is None
           else ['xvfb-run', '-a', app_path])
    p = subprocess.Popen(cmd + args, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    time.sleep(seconds)
    alive = p.poll() is None
    if alive:
        p.kill()
    out = p.communicate()[0].decode(errors='replace')[-400:]
    if not alive and p.returncode != 0:
        raise RuntimeError('exited %s: %s' % (p.returncode, out))
    return 'alive after %ss' % seconds if alive else 'exit 0'


@check('casaplotms binary (headless)')
def t_plotms():
    import casaplotms
    app = os.path.join(os.path.dirname(casaplotms.__file__), '__bin__', 'casaplotms-x86_64.AppImage')
    return _launch_gui(app, ['--nopopups'])


@check('casaviewer binary (headless)')
def t_viewer():
    import casaviewer
    app = os.path.join(os.path.dirname(casaviewer.__file__), '__bin__', 'casaviewer-x86_64.AppImage')
    return _launch_gui(app, [])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--no-gui', action='store_true')
    ap.add_argument('--gpu', action='store_true')
    a = ap.parse_args()
    t_python()
    t_stack()
    t_jax()
    if a.gpu:
        t_jax_gpu()
    t_casatools()
    if not a.quick:
        tmp = tempfile.mkdtemp(prefix='morphen-smoke-')
        cwd0 = os.getcwd()
        os.chdir(tmp)
        t_measures()
        t_casatasks(tmp)
        t_morphen()
        t_astropy()
        if not a.no_gui:
            t_plotms()
            t_viewer()
        # Leave tmp before deleting it: casacore's exit handler calls
        # getcwd() and aborts (rc 134) if the cwd no longer exists.
        os.chdir(cwd0)
        shutil.rmtree(tmp, ignore_errors=True)
    nfail = sum(not ok for _, ok, _ in RESULTS)
    print('\n%d/%d checks passed' % (len(RESULTS) - nfail, len(RESULTS)))
    sys.exit(nfail)


if __name__ == '__main__':
    main()
