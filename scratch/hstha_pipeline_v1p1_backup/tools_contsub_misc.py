from astropy.io import fits
import numpy as np
from glob import glob 
from synphot import SpectralElement, units
import os
import warnings 
from scipy.ndimage import binary_closing
warnings.filterwarnings('ignore')


def get_hdu(rootdir, filename, hdu_id=0, return_filename=False):
    filename_full = glob(rootdir+filename)[0]
    if hdu_id == 'all':
        hdu = fits.open(filename_full)
    else:
        hdu = fits.open(filename_full)[hdu_id]
    print(filename_full)

    if return_filename: 
        return(hdu, filename_full)
    else:  
        return(hdu)


def write_hdu(hdu, rootdir, filename, appdir='hst_contsub/'):
    filename_full = rootdir+appdir+filename
    print(filename_full)
    hdu.writeto(filename_full, overwrite=True)


def make_paths(rootdir, appdir='hst_contsub/'):
    
    print('[Info] Outputing to the following: '+rootdir+appdir)

    if not os.path.isdir(rootdir+appdir):
        os.mkdir(rootdir+appdir)  
    if not os.path.isdir(rootdir+appdir+'figs'):
        os.mkdir(rootdir+appdir+'figs')
    os.system('rm -rf '+rootdir+appdir+'*.fits')


def remove_paths(rootdir, appdir='hst_contsub/'):
    os.system('rm -rf %s%s/figs/' %(rootdir, appdir))
    os.system('mkdir %s%s/figs/' %(rootdir, appdir))
    

def get_bandpassinfo(rootdir_bp):

    files = glob('%s*.dat' %rootdir_bp)
    files.sort()

    bp = {}
    for file in files:

        print(file)

        area = 45238.93416 * units.AREA  # HST
        bp_ = SpectralElement.from_file(file)
        name = file.split('/')[-1].split('.dat')[0].replace('HST_', '').replace('.F', '_F')
        name = name.replace('WFC_', '')
        name = name.replace('WFC3_', '')
        name = name.replace('UVIS1', 'UVIS')
        name = name.replace('UVIS2', 'UVIS')

        bp[name] = {'equivwidth': bp_.equivwidth().value, 
                    'integrate': bp_.integrate().value, 
                    'rmswidth': bp_.rmswidth().value, 
                    'photbw': bp_.photbw().value, 
                    'fwhm': bp_.fwhm().value, 
                    'rectwidth': bp_.rectwidth().value, 
                    'pivot': bp_.pivot().value, 
                    'unit_response': bp_.unit_response(area).value}  

    return(bp)


def get_nanzeros(hdu):
    hdu.data[hdu.data == 0] = np.nan
    return(hdu)


def remove_nan_padding(hdu):
    """
    Remove padding of NaN values from the edges of an HDU image.
    
    Parameters:
    - hdu (HDU): The input HDU.
    
    Returns:
    - HDU: A new HDU with padding removed.
    """
    
    print('[INFO] Remove NaN values around edge of image...')

    data = hdu.data
    
    # Check if the entire image is NaNs
    if np.all(np.isnan(data)):
        print("Warning: Entire image is NaNs. Returning the original image.")
        return hdu
    
    # Find valid data indices along each axis
    valid_x = np.where(np.nansum(data, axis=0)!=0)[0]
    valid_y = np.where(np.nansum(data, axis=1)!=0)[0]

    # In the rare case there's still no valid data
    if len(valid_x) == 0 or len(valid_y) == 0:
        print("Warning: No valid data found. Returning the original image.")
        return hdu
    
    # Crop the data array
    cropped_data = data[valid_y[0]:valid_y[-1]+1, valid_x[0]:valid_x[-1]+1]
    
    # Create a new header by copying the old one
    new_header = hdu.header.copy()
    
    # Update reference pixel (CRPIX) values in the header
    if 'CRPIX1' in new_header:
        new_header['CRPIX1'] -= valid_x[0]
    if 'CRPIX2' in new_header:
        new_header['CRPIX2'] -= valid_y[0]
    
    # Create a new HDU with cropped data and updated header
    new_hdu = fits.PrimaryHDU(cropped_data, new_header)
    
    return new_hdu

def get_covmask(hdu_hst_f555w, hdu_hst_f65Xn, hdu_hst_f814w):

    """Mask to same coverage..."""

    mask = ~(np.isnan(hdu_hst_f65Xn.data) | np.isnan(hdu_hst_f555w.data) | np.isnan(hdu_hst_f814w.data))

    mask = mask*1
    mask_close = binary_closing(mask, structure=np.ones((10,10)), iterations=5)

    hdu_hst_f555w.data[~mask_close] = np.nan
    hdu_hst_f65Xn.data[~mask_close] = np.nan
    hdu_hst_f814w.data[~mask_close] = np.nan

    return(hdu_hst_f555w, hdu_hst_f65Xn, hdu_hst_f814w)