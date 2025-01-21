# Radio Spectral Index Maps
## Introduction

The radio spectral index is a measure of the steepness of the radio spectrum of a source. It is defined as the slope of the log-log plot of the flux density of a source as a function of frequency. The spectral index is defined as:

$$
\alpha = \frac{\log(S_2/S_1)}{\log(\nu_2/\nu_1)}
$$

where $S_1$ and $S_2$ are the flux densities at frequencies $\nu_1$ and $\nu_2$ respectively. The spectral index is negative for a steep spectrum source and positive for a flat spectrum source. The spectral index is a useful probe to study the physical properties of a radio source. For example, a steep spectrum source is usually associated with synchrotron radiation from relativistic electrons, while a flat spectrum source is usually associated with thermal radiation from ionized gas.


## Radio Spectral Index Maps

With `morphen`, you can create radio spectral index maps from radio images at different frequencies. The radio spectral index map is created by calculating the spectral index at each pixel of the radio images. The radio spectral index map can be used to study the spectral properties of the source and to identify regions of different spectral indices.

As input, `morphen` uses two or more radio images at different frequencies. For multiple images, `morphen` uses sub-band images generated with `WSClean`. 

### Input of $uv$ and beam matched images

The ideal case scenario when creating radio spectral index maps is to have $uv$ and beam matched images at different frequencies. This ensures that the spectral index map is not affected by differences in the $uv$ coverage and beam size of the images. Providing $uv$ and beam matched images as input to `morphen` is recommended. 

As an example, the Jupyter notebook [morphen_spec_index_Mrk331_example.ipynb](morphen_spec_index_Mrk331_example.ipynb)
contains a step-by-step guide to create radio spectral index maps Mrk 331 using $uv$ and beam matched WSClean images, from 1.1 up to 35GHz.

### Input of images with different beam size
However, for cases where beam matched images are not available, `morphen` will create radio spectral index maps from images with different beam size. To transform all images into a common beam size, `morphen` uses the `alpha` code [https://github.com/IanHeywood/alpha](https://github.com/IanHeywood/alpha), with improvements made by Javier Moldon and Geferson Lucatelli.

### ** WORK IN PROGRESS ** 
### Power-law spectral index map.

### Non-thermal spectral index and thermal fraction maps.



<!-- ## Minimisation
This implementation uses a robust fit, with the `LMFIT` package.  -->