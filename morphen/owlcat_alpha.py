import os
import sys
import numpy
#import pylab
import glob
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel
from astropy.convolution import convolve
from scipy.optimize import curve_fit
import os
import re
import glob
import numpy
import matplotlib.pyplot as plt
import os
import glob
import pylab
import sys
import scipy.signal
from shutil import copyfile
from astropy.convolution import Gaussian2DKernel
from astropy.io import fits
from astropy.wcs import WCS
from astropy.convolution import convolve
from optparse import OptionParser
from scipy import ndimage, stats
from scipy.optimize import curve_fit
import lmfit

# overwrite = True


def linear_function(x, a, b):
    return a + b*x #+ c*x*x
def do_fit_spec(data,freqs):
    def linear_function(x, a, b):
        return a + b*x #+ c*x*x
    
    def min_func(params):
        a = params['a']
        b = params['b']
        res = data - linear_function(freqs, a, b)
        return res.copy()
        
    
    fit_params = lmfit.Parameters()
    fit_params.add("a", value=0.0, min=-100, max=100)
    fit_params.add("b", value=-1.0, min=-2.5, max=1.5)
    
    mini = lmfit.Minimizer(min_func, fit_params, max_nfev=5000,
                        nan_policy='omit', reduce_fcn='neglogcauchy')
    
    result_1 = mini.minimize(method='least_squares',
                           max_nfev=200000, #f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12, 
                           verbose=0
                           )
    second_run_params = result_1.params
    
    result = mini.minimize(method='least_squares',
                           params=second_run_params,
                           max_nfev=200000, #f_scale = 1.0,
                           loss="cauchy", tr_solver="exact",
                           ftol=1e-12, xtol=1e-12, gtol=1e-12, 
                           verbose=0
                        )
    
    return result

def get_image(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
        image = numpy.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        image = numpy.array(input_hdu.data[0, :, :])
    elif len(input_hdu.data.shape) == 4:
        image = numpy.array(input_hdu.data[0, 0, :, :])
    else:
        image = numpy.array(input_hdu.data[0, 0, 0, :, :])
    return image


def deg2rad(xx):
    return numpy.pi * xx / 180.0


def get_beam(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    hdr = input_hdu.header
    bmaj = hdr.get('BMAJ')
    bmin = hdr.get('BMIN')
    bpa = hdr.get('BPA')
    pixscale = hdr.get('CDELT2')
    return bmaj, bmin, bpa, pixscale


def flush_fits(newimage, fitsfile):
    f = fits.open(fitsfile, mode='update')
    input_hdu = f[0]
    if len(input_hdu.data.shape) == 2:
        input_hdu.data[:, :] = newimage
    elif len(input_hdu.data.shape) == 3:
        input_hdu.data[0, :, :] = newimage
    elif len(input_hdu.data.shape) == 4:
        input_hdu.data[0, 0, :, :] = newimage
    else:
        input_hdu.data[0, 0, 0, :, :] = newimage
    f.flush()


def bb(xx):
    return str(round((xx * 3600.0), 3))


def get_target_beam(residuals,mode='max'):
    majors = []
    minors = []
    pas = []
    for residual in residuals:
        bmaj, bmin, bpa, pixscale = get_beam(residual)
        majors.append(bmaj)
        minors.append(bmin)
        pas.append(bpa)
    if mode == 'max':
        target_bmaj = numpy.max(majors)
        target_bmin = numpy.max(minors)
        target_bpa = numpy.median(numpy.array(pas))
    if mode == 'mean':
        target_bmaj = numpy.mean(majors)
        target_bmin = numpy.mean(minors)
        target_bpa = numpy.median(numpy.array(pas))
    if mode == 'median':
        target_bmaj = numpy.median(majors)
        target_bmin = numpy.median(minors)
        target_bpa = numpy.median(numpy.array(pas))
    
    print(
        '       | Target beam is: ' + bb(target_bmaj) + ',' + bb(target_bmin) + ',' + str(
            round(target_bpa, 3)))
    return target_bmaj, target_bmin, target_bpa


def drop_deg_axes(fitsfile):
    input_hdu = fits.open(fitsfile, mode='update')
    data = input_hdu[0].data[0, 0, :, :]
    input_hdu[0].data = data
    input_hdu.flush()
    hdr = input_hdu[0].header


def fix_beam_header(fitsfile, bmaj, bmin, bpa):
    hdu = fits.open(fitsfile, mode='update')
    hdr = hdu[0].header
    hdr.set('BMAJ', bmaj)
    hdr.set('BMIN', bmin)
    hdr.set('BPA', bpa)
    hdu.flush()


def restore_convolve(residuals,kernel_size,ref_target_beam = None,mode='max',target_beam=None,
                     do_debug=True,overwrite=True):
    if ref_target_beam is not None:
        target_bmaj, target_bmin, target_bpa = get_target_beam([ref_target_beam])
    elif target_beam is not None:
        target_bmaj, target_bmin, target_bpa = target_beam['target_bmaj'], target_beam['target_bmin'], target_beam['target_bpa']
    else:
        target_bmaj, target_bmin, target_bpa = get_target_beam(residuals,mode)
    for residual_fits in residuals:
        print('       | ' + residual_fits)
        # Get restoring beam
        bmaj, bmin, bpa, pixscale = get_beam(residual_fits)

        # Setup FITS files for beams and kernels
        print('++++++++++++++++++++++++++++')
        print(residual_fits)
        print(residual_fits.replace('.fits', '_template.fits'))
        template_fits = residual_fits.replace('.fits', '_template.fits')
        target_beam_fits = template_fits.replace('.fits', '_target_beam.fits')
        restoring_beam_fits = template_fits.replace('.fits', '_restoring_beam.fits')
        kernel_fits = template_fits.replace('.fits', '_kernel.fits')

        template_kernel_fits = template_fits.replace('.fits', '_kernel.fits')
        
        if os.path.exists(template_kernel_fits) and overwrite:
            os.system(f"rm {template_kernel_fits}")
        
        restored_fits = residual_fits.replace('residual', 'image-conv')
        pbcor_fits = restored_fits.replace('.fits', '_pbcor.fits')
        # breakpoint()
        print(' ----> | Creating kernel template image')
        os.system('fitstool.py -f -z ' + str(kernel_size) + ' -o ' + template_fits + ' ' +
                  residuals[0])
        drop_deg_axes(template_fits)

        # Target beam image
        print(' ----> | Creating target beam image')
        copyfile(template_fits, target_beam_fits)
        target_xstd = target_bmin / (2.3548 * pixscale)
        target_ystd = target_bmaj / (2.3548 * pixscale)
        target_theta = deg2rad(target_bpa)
        target_gaussian = Gaussian2DKernel(x_stddev=target_xstd, y_stddev=target_ystd,
                                           theta=target_theta, x_size=kernel_size,
                                           y_size=kernel_size, mode='center')
        target_image = target_gaussian.array
        target_image = target_image / numpy.max(target_image)
        flush_fits(target_image, target_beam_fits)

        # Restoring beam image
        print(' ----> | Creating restoring beam image')
        copyfile(template_fits, restoring_beam_fits)
        xstd = bmin / (2.3548 * pixscale)
        ystd = bmaj / (2.3548 * pixscale)
        theta = deg2rad(bpa)
        restoring_gaussian = Gaussian2DKernel(x_stddev=xstd, y_stddev=ystd, theta=theta,
                                              x_size=kernel_size, y_size=kernel_size,
                                              mode='center')
        restoring_image = restoring_gaussian.array
        restoring_image = restoring_image / numpy.max(restoring_image)
        flush_fits(restoring_image, restoring_beam_fits)
        # Run PyPHER to get homogenising kernel
        os.system(
            'pypher ' + restoring_beam_fits + ' ' + target_beam_fits + ' ' + kernel_fits)

        print(' <---- | Reading kernel image')
        kernel_image = get_image(kernel_fits)
        print(' <---- | Reading residual image')
        residual_image = get_image(residual_fits)
        print('  (X)  | Convolving residual with homogenising kernel')
        residual_conv_image = scipy.signal.fftconvolve(residual_image, kernel_image,
                                                       mode='same')

        print(' <---- | Reading model image')
        model_image = get_image(residual_fits.replace('residual', 'model'))
        print('  (X)  | Convolving model with target beam')
        model_conv_image = scipy.signal.fftconvolve(model_image, target_image,
                                                    mode='same')

        print('   +   | Restoring model to residual')
        restored_image = residual_conv_image + model_conv_image
        print(' ----> | Creating restored image')
        copyfile(residual_fits, restored_fits)
        flush_fits(restored_image, restored_fits)
        fix_beam_header(restored_fits, target_bmaj, target_bmin, target_bpa)

        #    print('   /   | Applying primary beam correction')
        #    pbcor_image = pbcor(pbdir,restored_fits)
        #    print(' ----> | Creating primary beam corrected image')
        #    copyfile(restored_fits,pbcor_fits)
        #    flush_fits(pbcor_image,pbcor_fits)

        # Clean up templates
        if not do_debug:
            os.system(f"rm {template_fits}")
            os.system(f"rm {target_beam_fits}")
            os.system(f"rm {restoring_beam_fits}")
            os.system(f"rm {kernel_fits}")
            os.system(f"rm {kernel_fits.replace('.fits', '.log')}")
        print('       | Done')
        print('::::::::::::::::::::::::::::::::::::::::::::::::::::::::::')


def fitfunc(x, A, B):
    return A * x + B


def getfreqs(fitslist):
    freqs = []
    for fitsfile in fitslist:
        hdu = fits.open(fitsfile)
        hdr = hdu[0].header
        freq = hdr['CRVAL3']
        freqs.append(freq)
    freqs = numpy.array(freqs)
    return freqs


def getImage(fitsfile):
    input_hdu = fits.open(fitsfile)[0]
    if len(input_hdu.data.shape) == 2:
        image = numpy.array(input_hdu.data[:, :])
    elif len(input_hdu.data.shape) == 3:
        image = numpy.array(input_hdu.data[0, :, :])
    else:
        image = numpy.array(input_hdu.data[0, 0, :, :])
    return image


def makecube(fitslist):
    temp = []
    for fitsfile in fitslist:
        img = getImage(fitsfile)
        temp.append(img)
    cube = numpy.dstack(temp)
    return cube


def make_alpha_maps(fitslist, reslist, prefix,threshold=None,mask=None):
    # alphafits = prefix + '_alpha_' + str(threshold).replace('.', 'p') + '.fits'
    # alphaerrorfits = prefix + '_alphaerror_' + str(threshold).replace('.', 'p') + '.fits'
    # f0fits = prefix + '_f0_' + str(threshold).replace('.', 'p') + '.fits'
    # nfluxfits = prefix + '_nflux_' + str(threshold).replace('.', 'p') + '.fits'
    alphafits = prefix + '_alpha.fits'
    alphaerrorfits = prefix + '_alphaerror.fits'
    f0fits = prefix + '_f0.fits'
    nfluxfits = prefix + '_nflux.fits'

    print(fitslist)
    os.system('cp ' + fitslist[0] + ' ' + alphafits)
    os.system('cp ' + fitslist[0] + ' ' + alphaerrorfits)
    # os.system('cp '+fitslist[0]+' '+f0fits)
    # os.system('cp '+fitslist[0]+' '+nfluxfits)

    freqs = getfreqs(fitslist)
    logf = numpy.log10(freqs)
    print(freqs)
    print(logf)
    cube = makecube(fitslist)
    cube_res = makecube(reslist)
    
    if (threshold is None) and (mask is None):
        raise ValueError("At least a mask or threshold must be provided.")
    
    if (threshold is not None) and (mask is None):
        cube = makecube(fitslist)
        mask = cube > threshold
        cube[~mask] = numpy.nan
        flattened_mask = numpy.mean(mask, axis=2)
        skymask = flattened_mask > 0.5
        idx = numpy.column_stack(numpy.where(skymask == True))
    else:
        idx = numpy.column_stack(numpy.where(mask==True))
        mask_3d = mask[:, :, numpy.newaxis]
        inv_mask = ~mask
        inv_mask_3d = inv_mask[:, :, numpy.newaxis]
        cube = numpy.where(inv_mask_3d, numpy.nan, cube)
        mask = inv_mask_3d.copy()
        # masked_cube_res = numpy.where(inv_mask_3d, numpy.nan, conv_cube_res)
        
    


    alphaimage = numpy.empty((cube.shape[0], cube.shape[1]))
    alphaerrorimage = numpy.empty((cube.shape[0], cube.shape[1]))
    f0image = numpy.empty((cube.shape[0], cube.shape[1]))
    nfluximage = numpy.empty((cube.shape[0], cube.shape[1]))

    alphaimage[:] = numpy.nan
    alphaerrorimage[:] = numpy.nan
    f0image[:] = numpy.nan
    nfluximage[:] = numpy.nan

    nspec = len(idx)
    tenpc = int(nspec / 10.0)
    count = 0
    pcount = 0
    for i, j in idx:        
        spec = numpy.array(cube[i, j, :])
        spec_res = numpy.array(cube_res[i, j, :])
        if numpy.isnan(spec.sum()):
            specmask = ~numpy.isnan(spec)
        else:
            specmask = [True] * len(freqs)
        # print(spec[specmask])
        # print(freqs[specmask])
        num_points = len(spec[specmask])
        f0 = numpy.mean(freqs[specmask])
        logspec = numpy.log10(abs(spec[specmask]))
        logspec_res = numpy.log10(abs(spec_res[specmask]))
        # print(logspec)
        # results_fit = do_fit_spec(logspec,logf[specmask])
        # alpha_err = results_fit.params['a'].value
        # alphaimage[i, j] = results_fit.params['b'].value
        # alphaerrorimage[i, j] = alpha_err
        
        popt, pcov = curve_fit(fitfunc, logf[specmask], logspec,
                            #    sigma=logspec_res, 
                               sigma=abs(logspec*0.05),
                               absolute_sigma=True)
        A = popt[0]
        alpha_err = numpy.sqrt(numpy.diag(pcov))[0]
        alphaimage[i, j] = A
        alphaerrorimage[i, j] = alpha_err

        f0image[i, j] = f0
        nfluximage[i, j] = num_points

        if count == 0:
            print(str(pcount) + '%...')
            model_best = fitfunc(logf[specmask], popt[0], popt[1])
            plt.figure()
            plt.errorbar(logf[specmask], 
                        logspec, 
                        yerr=abs(logspec*0.05), 
                        # yerr=abs(logspec_res), 
                        fmt='o',
                        label='Data', color='k', ecolor='gray',
                        alpha=0.5)

            plt.plot(logf[specmask], model_best,
                    color='red', ls='-.', label='Best-fit model')
            plt.legend()
            plt.show()
        
        count += 1
        if count == tenpc:
            count = 0
            pcount += 10

    flush_fits(alphaimage, alphafits)
    flush_fits(alphaerrorimage, alphaerrorfits)
    # flush_fits(f0image,f0fits)
    # flush_fits(nfluximage,nfluxfits)
    return(alphaimage,alphaerrorimage)