# Adapted by Ashley Barnes
# Thanks to Stephen Hannon - shann004@ucr.edu; UC Riverside/Caltech/IPAC
# Translated into Python from the original work of Janice Lee (some original code is included but commented out)

import numpy as np
from astropy.io import fits
import os

def perform_continuum_subtraction_onecont(hdu_halpha, hdu_cont1, halpha_outputfilename, cont_outputfilename):
    """
    Perform continuum subtraction on a single emission line image using one continuum image.

    Args:
        hdu_halpha (astropy.io.fits.ImageHDU): HDU containing the emission line image (e.g., H-alpha).
        hdu_cont1 (astropy.io.fits.ImageHDU): HDU containing the continuum image.
        photflam_cont1 (float): Photometric sensitivity of the continuum image.
        photplam_cont1 (float): Pivot wavelength of the continuum image.
        zpt_ab_halpha (float): AB zeropoint of the emission line image.
        halpha_outputfilename (str): Output filename for the continuum-subtracted emission line image.
        cont_outputfilename (str): Output filename for the scaled continuum image.

    Returns:
        None
    """

    # Retrieve photflam & photplam for halpha files
    photflam_halpha = float(hdu_halpha.header['PHOTFLAM']) 
    photplam_halpha = float(hdu_halpha.header['PHOTPLAM']) 

    # Retrieve photflam & photplam for continuum files
    photflam_cont1 = float(hdu_cont1.header['PHOTFLAM']) 
    photplam_cont1 = float(hdu_cont1.header['PHOTPLAM']) 

    # Calculate ZP for halpha
    zpt_ab_halpha = -2.5 * np.log10(photflam_halpha) - 5 * np.log10(photplam_halpha) - 2.408

    # Calculate ZP for continuum and scale
    zpt_ab_cont = -2.5 * np.log10(photflam_cont1) - 5 * np.log10(photplam_cont1) - 2.408
    scale_cont = 10 ** (-0.4 * (zpt_ab_halpha - zpt_ab_cont))
    coef_cont = 1 / scale_cont

    # Scale continuum data & save file
    hdu_cont1.data = hdu_cont1.data * coef_cont

    # Subtract scaled continuum from emission line data & save file
    hdu_halpha.data = hdu_halpha.data - hdu_cont1.data

    hdu_cont1.writeto(cont_outputfilename, overwrite=True)
    hdu_halpha.writeto(halpha_outputfilename, overwrite=True)

    return(None)

def perform_continuum_subtraction_twocont(hdu_halpha, hdu_cont1, hdu_cont2, halpha_outputfilename, cont_outputfilename):
    """
    Perform continuum subtraction on an emission line image using a single continuum image.

    Args:
        hdu_halpha (fits.PrimaryHDU or fits.ImageHDU): HDU object of the emission line image.
        hdu_cont1 (fits.PrimaryHDU or fits.ImageHDU): HDU object of the continuum image.
        hdu_cont2 (fits.PrimaryHDU or fits.ImageHDU): HDU object of the second continuum image.
        halpha_outputfilename (str): Output filename for the continuum-subtracted emission line image.
        cont_outputfilename (str): Output filename for the scaled continuum image.

    Returns:
        None
    """
    # Retrieve photflam & photplam for Halpha files
    photflam_halpha = float(hdu_halpha.header['PHOTFLAM'])
    photplam_halpha = float(hdu_halpha.header['PHOTPLAM'])

    # Retrieve photflam & photplam for continuum files
    photflam_cont1 = float(hdu_cont1.header['PHOTFLAM'])
    photplam_cont1 = float(hdu_cont1.header['PHOTPLAM'])
    photflam_cont2 = float(hdu_cont2.header['PHOTFLAM'])
    photplam_cont2 = float(hdu_cont2.header['PHOTPLAM'])

    # Calculate zero point (ZP) for Halpha
    zpt_ab_halpha = -2.5 * np.log10(photflam_halpha) - 5 * np.log10(photplam_halpha) - 2.408

    # Calculate ZP and scaling factor for continuum 1
    zpt_ab_cont1 = -2.5 * np.log10(photflam_cont1) - 5 * np.log10(photplam_cont1) - 2.408
    scale_cont1 = 10 ** (-0.4 * (zpt_ab_halpha - zpt_ab_cont1))

    # Calculate ZP and scaling factor for continuum 2
    zpt_ab_cont2 = -2.5 * np.log10(photflam_cont2) - 5 * np.log10(photplam_cont2) - 2.408
    scale_cont2 = 10 ** (-0.4 * (zpt_ab_halpha - zpt_ab_cont2))

    # Calculate weights for continuum 1 and continuum 2
    weight_cont1 = abs(photplam_cont2 - photplam_halpha) / abs(photplam_cont1 - photplam_cont2)
    weight_cont2 = abs(photplam_cont1 - photplam_halpha) / abs(photplam_cont1 - photplam_cont2)

    # Calculate coefficients for continuum scaling
    coef_conti = 1 / scale_cont1
    coef_cont1 = weight_cont1 / scale_cont1
    coef_cont2 = weight_cont2 / scale_cont2

    # Perform continuum scaling
    hdu_conti = hdu_cont1.copy()
    hdu_conti.data = hdu_conti.data * coef_conti
    hdu_cont1.data = hdu_cont1.data * coef_cont1
    hdu_cont2.data = hdu_cont2.data * coef_cont2

    # Create a mask to filter out pixels less than or equal to 0
    cont_ab_mask = hdu_cont1.data * hdu_cont2.data
    mask_cont = cont_ab_mask > 0

    # Prepare weighted continuum data
    cont_ab_wt = hdu_cont1.data + hdu_cont2.data

    # Prepare scaled continuum data
    cont_ab_scaled = hdu_conti.data

    # Perform continuum subtraction
    data_cont = np.where(mask_cont, cont_ab_wt, cont_ab_scaled)

    # Write final continuum file 
    hdu_cont = fits.PrimaryHDU(data_cont,hdu_cont1.header)
    hdu_cont.writeto(cont_outputfilename, overwrite=True)
    print('[INFO] Continuum file saved:', cont_outputfilename + '.fits')

    # Subtract scaled continuum from emission line data & save file
    hdu_halpha.data = hdu_halpha.data - hdu_cont.data
    hdu_halpha.writeto(halpha_outputfilename, overwrite=True)

    return(None)

def continuum_subtraction(halpha_inputfilename, cont1_inputfilename, halpha_outputfilename, cont_outputfilename, cont2_inputfilename=None):
    """
    Perform continuum subtraction on an emission line image using one or two continuum images.

    Args:
        halpha_inputfilename (str): Input filename for the emission line image.
        cont1_inputfilename (str): Input filename for the first continuum image.
        cont2_inputfilename (str): Input filename for the second continuum image (optional).
        halpha_outputfilename (str): Output filename for the continuum-subtracted emission line image.
        cont_outputfilename (str): Output filename for the scaled continuum image.

    Returns:
        None
    """
    print(f'[INFO] Loading Continuum 1: {cont1_inputfilename}')
    print(f'[INFO] Loading Emission Line Image: {halpha_inputfilename}')

    # Open halpha and continuum files
    with fits.open(halpha_inputfilename) as hdus_halpha, fits.open(cont1_inputfilename) as hdus_cont1:
        hdu_halpha = hdus_halpha[0].copy() if hdus_halpha[0].data is not None else hdus_halpha[1].copy()
        hdu_cont1 = hdus_cont1[0].copy() if hdus_cont1[0].data is not None else hdus_cont1[1].copy()

    # Open halpha and second continuum files
    if cont2_inputfilename is not None:
        print(f'[INFO] Loading Continuum 2: {cont2_inputfilename}')
        with fits.open(cont2_inputfilename) as hdus_cont2:
            hdu_cont2 = hdus_cont2[0].copy() if hdus_cont2[0].data is not None else hdus_cont2[1].copy()

    # Perform continuum subtraction
    if cont2_inputfilename is None:
        perform_continuum_subtraction_onecont(hdu_halpha, hdu_cont1, halpha_outputfilename, cont_outputfilename)
    else: 
        perform_continuum_subtraction_twocont(hdu_halpha, hdu_cont1, hdu_cont2, halpha_outputfilename, cont_outputfilename)
        
    return(None)