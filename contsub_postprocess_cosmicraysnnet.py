from deepCR import deepCR
from astropy.io import fits
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from scipy import ndimage
import numpy as np
import glob

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

def cosmicray_finder_nnet(input_filename, output_filename, dilation_iterations=5, threshold=0.25,
                          model_path='ACS-WFC-F606W-2-32'):
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
    image = hdu.data

    # Create an instance of deepCR with specified model configuration
    mdl = deepCR(mask=model_path, device="CPU")

    # Apply the model to the input image to detect cosmic rays and generate a mask
    mask = mdl.clean(image, threshold=threshold, inpaint=False)
    
    # Dilate the mask to ensure surrounding regions of cosmic rays are also masked
    mask = ndimage.binary_dilation(mask, iterations=dilation_iterations)
    hdu_mask = fits.PrimaryHDU(mask*1, hdu.header)

    # Copy the original image and set the masked regions (cosmic rays) to NaN
    image_ = image.copy()
    image_[mask] = np.nan
    
    # Interpolate over the masked regions to fill them with suitable values
    interpolated_image = interpolate_masked(image_)
    interpolated_image = interpolate_masked(interpolated_image)

    mask = load_mask(os.path.dirname(output_filename))
    interpolated_image[mask] = np.nan

    # Create a new HDU with the interpolated image and the original header
    hdu_interpolated_image = fits.PrimaryHDU(interpolated_image, hdu.header)
    
    # Write the cleaned image to the output FITS file
    hdu_mask.writeto(output_filename.replace('.fits', '_masknnet.fits'), overwrite=True)
    hdu_interpolated_image.writeto(output_filename, overwrite=True)
    
    return()