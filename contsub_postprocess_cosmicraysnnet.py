from deepCR import deepCR
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from scipy import ndimage
import numpy as np
import glob
import os 

def load_mask(inputdir):
    filename = glob.glob('%s/*mask.fits' %inputdir)[0]
    hdu = fits.open(filename)[0]
    return(hdu.data==0) 

def interpolate_masked(image_masked, x_stddev=1):
    '''
    Interpolate NaN pixels in the masked image using a convolution.

    Parameters:
    image_masked: ndarray
        Masked 2D image data.
    x_stddev: float
        Standard deviation for Gaussian kernel.

    Returns:
    image_interpolated: ndarray
        Interpolated image.
    '''

    # Create a Gaussian kernel for the convolution
    kernel = Gaussian2DKernel(x_stddev=x_stddev)

    # Interpolate the NaN pixels in the image using a convolution
    image_interpolated = interpolate_replace_nans(image_masked, kernel)
    
    return(image_interpolated)

def cosmicray_finder_nnet(input_filename, output_filename, dilation_iterations=5, threshold=0.25, patch=1024,
                          model_path='/Users/abarnes/opt/anaconda3/lib/python3.9/site-packages/learned_models/mask/ACS-WFC-F606W.pth'):
    """
    Remove cosmic rays from a FITS image using the deepCR model and save the cleaned image to a specified output filename.
    
    Parameters:
    - input_filename (str): Path to the input FITS file containing the image with cosmic rays.
    - output_filename (str): Path to the output FITS file to save the cleaned image.
    - model_path (str): Path to the deepCR model used for cosmic ray removal. Default is 'ACS-WFC-F606W.pth'.
    - dilation_iterations (int): Number of dilation iterations to expand the mask. Default is 5.
    - threshold (float): Threshold value for the deepCR model to detect cosmic rays. Default is 0.25.
    
    Returns:
    None: The function writes the cleaned image to a specified FITS file.
    """
    
    # Load the FITS file and extract image data and header
    hdu = fits.open(input_filename)[0]
    image = np.array(hdu.data, dtype=float32)

    # Look elsewhere for model 
    if not os.path.isfile(model_path):
        model_path = '/lustre/opsw/work/abarnes/applications/anaconda3/lib/python3.11/site-packages/learned_models/mask/ACS-WFC-F606W.pth'

    # Create an instance of deepCR with specified model configuration
    mdl = deepCR(mask=model_path, device="CPU")

    print('[INFO] [deepCR] Running deepCR...')
    # Apply the model to the input image to detect cosmic rays and generate a mask
    if patch==None: 
        mask = mdl.clean(image, threshold=threshold, inpaint=False)
    else: 
        print('[INFO] [deepCR] Running with patch=%i' %patch)
        mask = mdl.clean(image, threshold=threshold, inpaint=False, segment=True, patch=patch)
    
    print('[INFO] [deepCR] Dilation of deepCR mask...')
    # Dilate the mask to ensure surrounding regions of cosmic rays are also masked
    mask = ndimage.binary_dilation(mask, iterations=dilation_iterations)
    hdu_mask = fits.PrimaryHDU(np.array(mask*1, dtype=np.int32), hdu.header)

    # Copy the original image and set the masked regions (cosmic rays) to NaN
    image_ = image.copy()
    image_[mask] = np.nan
    
    print('[INFO] [deepCR] Interpolated deepCR mask...')
    # Interpolate over the masked regions to fill them with suitable values
    interpolated_image = interpolate_masked(image_)
    interpolated_image = interpolate_masked(interpolated_image)

    mask = load_mask(os.path.dirname(output_filename))
    interpolated_image[mask] = np.nan
    interpolated_image = np.array(interpolated_image, dtype=np.float32)

    print('[INFO] [deepCR] Saving the deepCR mask...')    
    # Write the cleaned image to the output FITS file
    hdu_mask.writeto(output_filename.replace('.fits', '_masknnet.fits'), overwrite=True)

    hdu_interpolated_image = fits.PrimaryHDU(interpolated_image, hdu.header)
    hdu_interpolated_image.writeto(output_filename, overwrite=True)
    
    return()