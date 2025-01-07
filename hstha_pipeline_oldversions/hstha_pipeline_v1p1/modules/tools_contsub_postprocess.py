from astropy.io import fits
from astropy.stats import mad_std
import numpy as np
from scipy.ndimage import binary_dilation, binary_closing
from astropy.convolution import Gaussian2DKernel, convolve, interpolate_replace_nans
from deepCR import deepCR
import os

from tools_contsub_smoothregrid import * 

def get_mask(hdu):
    mask = ~np.isnan(hdu.data)*1
    mask_close = binary_closing(mask, structure=np.ones((10,10)), iterations=1)
    mask_hdu = fits.PrimaryHDU(np.int32(mask_close*1), hdu.header)
    return(mask_hdu) 

def get_mask_stars(hdu, hdu_mask, hdu_mask_map, sigma=0.5, factor=np.sqrt(2)):
    
    def smooth_hdu_gaussian(data, sigma_x=0.5, sigma_y=0.5):
        kernel = Gaussian2DKernel(sigma_x, sigma_y)
        smoothed_data = convolve(data, kernel, normalize_kernel=True)
        return smoothed_data

    # hdu_masked, mask = mask_hdu_with_ds9(hdu, region_filename)
    hdu_masked = hdu.copy()
    mask = hdu_mask.data>0

    std = mad_std(hdu.data, ignore_nan=True)
    mean = np.nanmean(hdu.data[hdu.data < std * 5])

    noise = np.random.normal(mean, std, hdu_masked.data.shape)
    noise = smooth_hdu_gaussian(noise * factor, sigma_x=sigma, sigma_y=sigma) 
    hdu_masked.data[mask] = noise[mask]

    mask = hdu_mask_map.data == 0
    hdu_masked.data[mask] = np.nan

    return(hdu_masked)

def get_interp_negs(hdu, hdu_mask, sigma_high=5, sigma_low=1):

    # Create a copy of the anchored HST data HDU
    hdu_i = hdu.copy()

    # Replace negative values with NaNs
    std = mad_std(hdu_i.data, ignore_nan=True)
    mask = hdu_i.data <= 3*std
    std = mad_std(hdu_i.data[mask], ignore_nan=True)

    # Replace negative values with NaNs
    std = mad_std(hdu_i.data[mask], ignore_nan=True)
    mask_high = hdu_i.data <= -std*sigma_high
    mask_low = hdu_i.data <= -std*sigma_low
    # mask_low = hdu_i.data <= 0 
    mask = binary_dilation(mask_high, iterations=-1, mask=mask_low)
    mask = binary_closing(mask, iterations=1)
    hdu_i.data[mask] = np.nan

    # Perform interpolation using Gaussian kernels
    kernel = Gaussian2DKernel(x_stddev=1)
    hdu_i.data = interpolate_replace_nans(hdu_i.data, kernel)
    hdu_i.data = interpolate_replace_nans(hdu_i.data, kernel)
    hdu_i.data = interpolate_replace_nans(hdu_i.data, kernel)

    mask = hdu_mask.data == 0
    hdu_i.data[mask] = np.nan

    # Save the processed anchored HST image to the output file
    print(f"[INFO] Negative values processed")

    return hdu_i

def get_cosmicrays(hdu, hdu_mask, dilation_iterations=5, threshold=0.25, patch=1024,
                          model_path='/Users/abarnes/opt/anaconda3/lib/python3.9/site-packages/learned_models/mask/ACS-WFC-F606W.pth'):
    
    # Load the FITS file and extract image data and header
    data = hdu.data.copy()

    # Look elsewhere for model 
    if not os.path.isfile(model_path):
        model_path = '/lustre/opsw/work/abarnes/applications/anaconda3/lib/python3.11/site-packages/learned_models/mask/ACS-WFC-F606W.pth'

    # Create an instance of deepCR with specified model configuration
    mdl = deepCR(mask=model_path, device="CPU")

    print('[INFO] [deepCR] Running deepCR...')
    # Apply the model to the input image to detect cosmic rays and generate a mask
    if patch==None: 
        mask = mdl.clean(data, threshold=threshold, inpaint=False)
    else: 
        print('[INFO] [deepCR] Running with patch=%i' %patch)
        mask = mdl.clean(data, threshold=threshold, inpaint=False, segment=True, patch=patch)
        mask = mask != 0
    
    # Dilate the mask to ensure surrounding regions of cosmic rays are also masked
    if dilation_iterations != 0:  
        print('[INFO] [deepCR] Dilation of deepCR mask...')
        mask = binary_dilation(mask, iterations=dilation_iterations)

    # hdu_mask = fits.PrimaryHDU(np.array(mask*1, dtype=np.int32), hdu.header)

    # Copy the original image and set the masked regions (cosmic rays) to NaN
    data[mask] = np.nan
    
    print('[INFO] [deepCR] Interpolated deepCR mask...')
    # Interpolate over the masked regions to fill them with suitable values
    kernel = Gaussian2DKernel(x_stddev=1)
    data_i = interpolate_replace_nans(data, kernel)

    mask = hdu_mask.data == 0
    data_i[mask] = np.nan

    data_i = np.array(data_i, dtype=np.float32)
    hdu_i = fits.PrimaryHDU(data_i, hdu.header)
    
    return(hdu_i)