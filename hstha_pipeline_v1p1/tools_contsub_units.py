from astropy.io import fits
import astropy.wcs as wcs
import astropy.units as u
from astropy.convolution import convolve_fft
from astropy.stats import mad_std
import numpy as np
from radio_beam import Beam
from reproject import reproject_interp
from matplotlib import pyplot as plt
from scipy.ndimage import binary_dilation
from astropy.modeling import models, fitting
from astropy.table import Table, vstack 
from glob import glob 
from synphot import SpectralElement, units
from matplotlib.ticker import (MultipleLocator)
import os
import warnings 
warnings.filterwarnings('ignore')


def get_electrons_2_ergcm2sA(hdu, photflam=None):

    data = hdu.data.copy()

    if photflam == None: 
        # Get the necessary header keywords for scaling and conversion
        photflam = hdu.header['PHOTFLAM']
    
    # Scale the data using photflam and photbw
    data_conv = data * photflam

    hdu.data = np.array(data_conv, dtype=np.float32) *1e20
    hdu.header['BUNIT'] = ('erg/s/cm2/A/pixel', '1e-20 erg/s/cm2/A')

    return(hdu)


def get_Jy_2_ergcm2sA(hdu, photplam):

    data = hdu.data.copy()
    
    w = photplam * u.AA
    a = 1. * u.Jy
    b = a.to(u.erg / u.cm**2 / u.s / u.AA, u.spectral_density(w))
    data_conv = data * b.value

    hdu.data = np.array(data_conv, dtype=np.float32) *1e20
    hdu.header['BUNIT'] = ('erg/s/cm2/A/pixel', '1e-20 erg/s/cm2/A')

    return(hdu)


def get_ergcm2sA_2_ergcm2s(hdu, photbw):

    data = hdu.data.copy()
    data_conv = data * photbw
    hdu.data = np.array(data_conv, dtype=np.float32)
    hdu.header['BUNIT'] = ('erg/s/cm2/pixel', '1e-20 erg/s/cm2')

    return(hdu)