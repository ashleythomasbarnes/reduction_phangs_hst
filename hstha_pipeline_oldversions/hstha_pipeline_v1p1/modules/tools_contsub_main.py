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

    """
    Perform continuum subtraction on H-alpha image data using two continuum images.

    This function takes in three HDU (Header/Data Unit) objects: one for the H-alpha image 
    and two for the continuum images. It calculates the continuum-subtracted H-alpha image 
    by scaling and combining the continuum images based on their photometric wavelengths.

    Parameters:
    hdu_halpha (HDU): HDU object containing the H-alpha image data.
    hdu_cont1 (HDU): HDU object containing the first continuum image data.
    hdu_cont2 (HDU): HDU object containing the second continuum image data.
    photplam_halpha (float, optional): Photometric wavelength of the H-alpha image. 
                                        If None, it is extracted from the HDU header.
    photplam_cont1 (float, optional): Photometric wavelength of the first continuum image. 
                                        If None, it is extracted from the HDU header.
    photplam_cont2 (float, optional): Photometric wavelength of the second continuum image. 
                                        If None, it is extracted from the HDU header.

    Returns:
    tuple: A tuple containing two HDU objects:
        - hdu_halpha_contsub (HDU): HDU object with the continuum-subtracted H-alpha image data.
        - hdu_halpha_cont (HDU): HDU object with the combined continuum image data.
    """

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


def get_contsub_err(hdu_halpha, hdu_cont1, hdu_cont2, 
                hdu_halpah_err=None, hdu_cont1_err=None, hdu_cont2_err=None,
                photplam_halpha=None, photplam_cont1=None, photplam_cont2=None,
                relwt_cont1=1, relwt_cont2=1, relwt_cont=1):

    """
    Perform continuum subtraction on H-alpha image data using two continuum images WITH associated error.

    This function takes in H-alpha and continuum images along with their errors,
    and computes the continuum-subtracted H-alpha image and its error. The 
    continuum images are weighted based on their photometric wavelengths.

    Parameters:
    hdu_halpha (HDU): HDU object for the H-alpha image.
    hdu_cont1 (HDU): HDU object for the first continuum image.
    hdu_cont2 (HDU): HDU object for the second continuum image.
    hdu_halpah_err (HDU, optional): HDU object for the H-alpha error image. Default is None.
    hdu_cont1_err (HDU, optional): HDU object for the first continuum error image. Default is None.
    hdu_cont2_err (HDU, optional): HDU object for the second continuum error image. Default is None.
    photplam_halpha (float, optional): Photometric wavelength of the H-alpha image. 
                                        If None, it is extracted from the HDU header.
    photplam_cont1 (float, optional): Photometric wavelength of the first continuum image. 
                                        If None, it is extracted from the HDU header.
    photplam_cont2 (float, optional): Photometric wavelength of the second continuum image. 
                                        If None, it is extracted from the HDU header.
    relwt_cont1 (float, optional): Relative weight of the first continuum image. Default is 1.
    relwt_cont2 (float, optional): Relative weight of the second continuum image. Default is 1.
    relwt_cont (float, optional): Relative weight of the combined continuum image. Default is 1.
                                        

    Returns:
    tuple: A tuple containing:
    - (hdu_halpha_contsub, hdu_halpha_cont): HDU objects for the continuum-subtracted H-alpha image and the combined continuum image.
    - (hdu_halpha_contsub_err, hdu_halpha_cont_err): HDU objects for the error images of the continuum-subtracted H-alpha and the combined continuum.
    """

    if photplam_halpha == None:
        photplam_halpha = hdu_halpha.header['PHOTPLAM']
        photplam_cont1 = hdu_cont1.header['PHOTPLAM']
        photplam_cont2 = hdu_cont2.header['PHOTPLAM']

    wt_cont1 = abs(photplam_cont2 - photplam_halpha) / abs(photplam_cont1 - photplam_cont2)
    wt_cont2 = abs(photplam_cont1 - photplam_halpha) / abs(photplam_cont1 - photplam_cont2)

    data_halpha = hdu_halpha.data
    data_halphaerr = hdu_halpah_err.data
    data_cont1 = hdu_cont1.data
    data_cont2 = hdu_cont2.data
    data_cont1err = hdu_cont1_err.data
    data_cont2err = hdu_cont2_err.data
    data_cont1relerr = data_cont1err / data_cont1
    data_cont2relerr = data_cont2err / data_cont2

    data_cont1[data_cont1<=0] = np.nan
    data_cont2[data_cont2<=0] = np.nan

    data_cont1 = np.log10(data_cont1)
    data_cont2 = np.log10(data_cont2)
    data_cont1err = np.abs(data_cont1relerr / np.log(10))
    data_cont2err = np.abs(data_cont2relerr / np.log(10))

    data_cont1_wt = data_cont1 * wt_cont1 * relwt_cont1
    data_cont2_wt = data_cont2 * wt_cont2 * relwt_cont2
    data_cont1err_wt = data_cont1err * wt_cont1
    data_cont2err_wt = data_cont2err * wt_cont2

    data_cont = 10**(data_cont1_wt + data_cont2_wt)
    data_conterr_log = np.sqrt(data_cont1err_wt**2 + data_cont2err_wt**2)
    data_conterr = data_cont * np.log(10) * data_conterr_log

    data_cont[np.isnan(data_cont)] = 0
    data_conterr[np.isnan(data_conterr)] = 0

    data_contsub = data_halpha - (data_cont * relwt_cont)
    data_contsuberr = np.sqrt(data_halphaerr**2 + data_conterr**2)

    data_contsuberr[np.isnan(data_contsub)] = np.nan

    hdu_halpha_cont = fits.PrimaryHDU(data_cont, header=hdu_halpha.header)
    hdu_halpha_contsub = fits.PrimaryHDU(data_contsub, header=hdu_halpha.header)
    hdu_halpha_cont_err = fits.PrimaryHDU(data_conterr, header=hdu_halpha.header)
    hdu_halpha_contsub_err = fits.PrimaryHDU(data_contsuberr, header=hdu_halpha.header)

    hdu_halpha_cont.data = np.array(hdu_halpha_cont.data, dtype=np.float32)
    hdu_halpha_contsub.data = np.array(hdu_halpha_contsub.data, dtype=np.float32)
    hdu_halpha_cont_err.data = np.array(hdu_halpha_cont_err.data, dtype=np.float32)
    hdu_halpha_contsub_err.data = np.array(hdu_halpha_contsub_err.data, dtype=np.float32)

    # Remove inf values
    hdu_halpha_cont_err.data[hdu_halpha_cont_err.data==np.inf] = np.nan
    hdu_halpha_contsub_err.data[hdu_halpha_contsub_err.data==np.inf] = np.nan

    return((hdu_halpha_contsub, hdu_halpha_cont), (hdu_halpha_contsub_err, hdu_halpha_cont_err))