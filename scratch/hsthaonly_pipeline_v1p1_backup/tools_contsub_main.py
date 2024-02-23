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


def get_contsub(hdu_halpha, hdu_cont1, hdu_cont2, 
                photplam_halpha=None, photplam_cont1=None, photplam_cont2=None):

    if photplam_halpha == None:
        photplam_halpha = hdu_halpha.header['PHOTPLAM']
        photplam_cont1 = hdu_cont1.header['PHOTPLAM']
        photplam_cont2 = hdu_cont2.header['PHOTPLAM']

    weight_cont1 = abs(photplam_cont2 - photplam_halpha) / abs(photplam_cont1 - photplam_cont2)
    weight_cont2 = abs(photplam_cont1 - photplam_halpha) / abs(photplam_cont1 - photplam_cont2)

    coef_cont1 = weight_cont1
    coef_cont2 = weight_cont2

    data_cont1 = hdu_cont1.data
    data_cont2 = hdu_cont2.data

    data_cont1[data_cont1<=0] = np.nan
    data_cont2[data_cont2<=0] = np.nan

    data_cont1 = np.log10(data_cont1)
    data_cont2 = np.log10(data_cont2)

    hdu_cont1.data = data_cont1 * coef_cont1
    hdu_cont2.data = data_cont2 * coef_cont2

    data_cont = 10**(hdu_cont1.data + hdu_cont2.data)
    data_cont[np.isnan(data_cont)] = 0

    hdu_halpha_cont = hdu_halpha.copy()
    hdu_halpha_contsub = hdu_halpha.copy()

    hdu_halpha_cont.data = data_cont
    hdu_halpha_contsub.data = hdu_halpha.data - data_cont

    hdu_halpha_cont.data = np.array(hdu_halpha_cont.data, dtype=np.float32)
    hdu_halpha_contsub.data = np.array(hdu_halpha_contsub.data, dtype=np.float32)

    return(hdu_halpha_contsub, hdu_halpha_cont)