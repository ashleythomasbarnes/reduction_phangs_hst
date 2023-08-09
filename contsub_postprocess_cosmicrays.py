import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import SqrtStretch
from astropy.visualization.mpl_normalize import ImageNormalize
from photutils.aperture import CircularAperture
from photutils.profiles import RadialProfile
from photutils.detection import DAOStarFinder
from photutils.segmentation import SourceFinder
from photutils.detection import find_peaks
from astropy.stats import sigma_clipped_stats
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.modeling import models, fitting
from astropy.nddata import Cutout2D
import numpy as np
from tqdm.auto import tqdm
import gc
import os 
from astropy.io import fits

def cosmicray_finder(hdu):

    image = hdu.data.copy()
    
    print("[INFO] Calculating mean, median, and standard deviation of the image")
    # Calculate the mean, median and standard deviation of the image
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    image = image-median

    print("[INFO] Setting threshold for star detection")
    # Set the threshold for star detection
    threshold = (5.0 * std)

    print("[INFO] Initializing DAOStarFinder and finding sources")
    # Use DAOStarFinder to find stars in the image
    daofind = DAOStarFinder(fwhm=3.0, threshold=threshold) 
    sources = daofind(image)
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))
    
    # Prepare a mask for valid sources
    mask = []
    
    print("[INFO] Processing radial profiles of sources")
    # Check the radial profile of each source and add it to the mask if it doesn't contain negative values
    for i in tqdm(range(len(positions)), desc='Processing radial profiles of sources'):
        rp = RadialProfile(image, positions[i], np.arange(15))
        if (rp.profile<0).any() != True:
            mask.append(i)
        
    print("[INFO] Removing masked sources")
    # Remove masked sources
    sources.remove_rows(mask)

    # Get the position of the remaining sources
    positions = np.transpose((sources['xcentroid'], sources['ycentroid']))

    print("[INFO] Found positions of cosmic rays")
    return positions


def gaussian_fitted_data(cutout, source_positions, size=20):
    '''
    Fit 2D Gaussian model to cutout data.

    Parameters:
    cutout: ndarray
        2D image cutout.
    source_positions: list
        Positions of the sources.
    size: int
        Size of the cutout.

    Returns:
    p: Gaussian2D model
        Fitted 2D Gaussian model.
    '''

    # Create a 2D Gaussian model
    y, x = np.mgrid[:cutout.shape[0], :cutout.shape[1]]
    p_init = models.Gaussian2D(amplitude=cutout.data.max(), x_mean=cutout.data.shape[1]/2, y_mean=cutout.data.shape[0]/2)

    # Fit the data using astropy.modeling
    fit_p = fitting.LevMarLSQFitter()
    p = fit_p(p_init, x, y, cutout.data)
        
    return p


def mask_ellipse(image, x_centre, y_centre, x_fwhm, y_fwhm):
    '''
    Apply an elliptical mask to the image.

    Parameters:
    image: ndarray
        2D image data.
    x_centre, y_centre: float
        The x and y coordinates of the center of the ellipse.
    x_fwhm, y_fwhm: float
        Full-width at half-maximum along x and y axes.

    Returns:
    image_masked: ndarray
        Image after applying the mask.
    '''

    # Create the same size array as image
    y, x = np.indices(image.shape)

    # Generate an elliptical mask based on input parameters
    mask = ((x - x_centre) ** 2 / x_fwhm ** 2 + (y - y_centre) ** 2 / y_fwhm ** 2) < 1

    # Apply the mask to the image
    image_masked = np.where(mask, np.nan, image)

    return image_masked


def create_circular_mask(h, w, center=None, radius=None):

    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if radius is None: # use the smallest distance between the center and image walls
        radius = min(center[0], center[1], w-center[0], h-center[1])

    Y, X = np.ogrid[:h, :w]
    dist_from_center = np.sqrt((X - center[0])**2 + (Y-center[1])**2)

    mask = dist_from_center <= radius
    return (mask)


def interpolate_masked(image_masked, x_stddev=0.75):
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


def subtract_intp_2d_gaussian_fitted_data(image_fit, image_sub, source_positions, size_cut=50, size_fit=20):
    '''
    Subtract interpolated 2D Gaussian fitted data from image.

    Parameters:
    image_fit, image_sub: ndarray
        2D image data.
    source_positions: list
        Positions of the sources.
    size_cut: int
        Size for cutting the image.
    size_fit: int
        Size for fitting the Gaussian.

    Returns:
    image_int: ndarray
        Image after subtraction and interpolation.
    '''

    image_int = image_sub.copy()
    
    # Loop over each source position
    for source_position in tqdm(source_positions, desc='Processing sources'):
        # Make a cutout around the source position
        cutout = Cutout2D(image_fit.data, source_position, size=size_cut)
        cutout.data[np.isnan(cutout.data)] = 0
        
        p = gaussian_fitted_data(cutout, source_position, size=size_fit)

        x_mean, ymean = cutout.to_original_position((p.x_mean.value, p.y_mean.value))
        
        image_masked = mask_ellipse(image_int, x_mean, ymean, p.x_fwhm, p.y_fwhm)
        
    image_int = interpolate_masked(image_masked)
        
    return(image_int)


def subtract_intp_cut_data(hdu, source_positions, size_mask=3):
    '''
    Subtract interpolated cutout data from image.

    Parameters:
    image_sub: ndarray
        2D image data.
    source_positions: list
        Positions of the sources.
    size_mask: int
        Size for masking

    Returns:
    image_int: ndarray
        Image after subtraction and interpolation.
    '''

    image_int = hdu.data.copy()
    image_masked = hdu.data.copy()
    h, w = image_masked.shape[:2]

    # for source_position in tqdm(source_positions, desc='Processing interpolation of masked sources'):
        # x1 = int(source_position[1] - size_mask)
        # x2 = int(source_position[1] + size_mask)
        # y1 = int(source_position[0] - size_mask)
        # y2 = int(source_position[0] + size_mask)
        # image_masked[x1:x2, y1:y2] = np.nan
    
    for source_position in tqdm(source_positions, desc='Processing interpolation of masked sources'):
        
        circular_mask = create_circular_mask(h, w, source_position, size_mask)
        image_masked[circular_mask] = np.nan
        
    image_int = interpolate_masked(image_masked)
    hdu_int = fits.PrimaryHDU(image_int, hdu.header)

    # Deleting image_int and image_masked and collected unreferenced objects to free memory
    del image_int, image_masked
    _ = gc.collect()
    
    return(hdu_int)

def load_mask(inputdir):
    filename = glob.glob('%s/*mask.fits' %inputdir)[0]
    hdu = fits.open(filename)[0]
    return(hdu.data==0) 

def save_masked(hdu, output_filename):

    mask = load_mask(os.path.dirname(output_filename))
    hdu.data[mask] = np.nan

    hdu.writeto(output_filename, overwrite=True)
    print(f"[INFO] Making and saved as {output_filename}")
    return()