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
    hdu.header['BUNIT'] = ('erg/s/cm2/A/pixel', '1e-20 erg/s/cm2/A/pixel')

    return(hdu)


def get_Jy_2_ergcm2sA(hdu, photplam):

    data = hdu.data.copy()
    
    w = photplam * u.AA
    a = 1. * u.Jy
    b = a.to(u.erg / u.cm**2 / u.s / u.AA, u.spectral_density(w))
    data_conv = data * b.value

    hdu.data = np.array(data_conv, dtype=np.float32) *1e20
    hdu.header['BUNIT'] = ('erg/s/cm2/A/pixel', '1e-20 erg/s/cm2/A/pixel')

    return(hdu)


def get_ergcm2sA_2_ergcm2s(hdu, photbw):

    hdu_new = hdu.copy()

    data = hdu_new.data
    data_conv = data * photbw

    hdu_new.data = np.array(data_conv, dtype=np.float32)
    hdu_new.header['BUNIT'] = ('erg/s/cm2/pixel', '1e-20 erg/s/cm2/pixel')

    return(hdu_new)


def convert_perpix_to_perarcsec(hdu):

    hdu = hdu.copy()  

    keys = list(hdu.header.keys())
    if 'CD1_1' in keys:
        keyword = 'CD1_1'
    elif 'PC1_1' in keys:
        keyword = 'PC1_1'
    else: 
        print('Warning - using CDELT!')
        keyword = 'CDELT1'    

    pix_size = (hdu.header[keyword]*u.deg).to('arcsec')
    pix_area = pix_size**2

    hdu.data = hdu.data / pix_area.value
    hdu.header['BUNIT'] =  ('erg/s/cm2/arcsec2', '10^-20 erg/s/cm2/arcsec2') 

    return(hdu)