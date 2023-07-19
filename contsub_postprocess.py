import os
from astropy.io import fits
from astropy.table import Table
import astropy.wcs as wcs
import astropy.units as u
from astropy.convolution import Gaussian2DKernel, convolve, convolve_fft, interpolate_replace_nans
from astropy.stats import mad_std
import numpy as np
from radio_beam import Beam
from scipy.ndimage.morphology import binary_dilation
from reproject import reproject_interp
from matplotlib import pyplot as plt

from astropy.stats import SigmaClip
from photutils.segmentation import detect_threshold, detect_sources
from photutils.utils import circular_footprint
from photutils.background import Background2D
import numpy as np
from astropy.io import fits

def process_halpha_units(input_filename, output_filename):
    """
    Processes the Halpha FITS file by changing the units to erg/s/cm^2.

    Args:
        input_filename (str): Path to the input Halpha FITS file.
        output_filename (str): Path to save the processed FITS file.

    Returns:
        astropy.io.fits.ImageHDU: Processed Halpha FITS image HDU.
    """
    print(f"[INFO] Processing {input_filename}")

    # Open the input FITS file and get the primary HDU
    hdu_hst = fits.open(input_filename)[0]

    # Get the necessary header keywords for scaling and conversion
    photflam = hdu_hst.header['PHOTFLAM']
    photplam = hdu_hst.header['PHOTPLAM']
    photbw = hdu_hst.header['PHOTBW']

    print("[INFO] Scaling the data...")
    # Scale the data using photflam and photbw
    hdu_hst.data = hdu_hst.data * photflam * photbw

    # Update the BUNIT header keyword to indicate the new units
    hdu_hst.header['BUNIT'] = ('erg/s/cm2', '1e-20 erg/s/cm2')

    print("[INFO] Converting units...")
    # Convert the data units to erg/s/cm^2
    hdu_hst.data = hdu_hst.data * 1e20

    print(f"[INFO] Saving the processed FITS file as {output_filename}")
    # Save the processed FITS file
    hdu_hst.writeto(output_filename, overwrite=True)

    return hdu_hst


def process_halpha_background(hdu_hst, halpha_filename, sigma_clip_sigma=3.0, detect_threshold_nsigma=2.0, 
                        detect_sources_npixels=10, footprint_radius=10, box_size=(100, 100), filter_size=(51, 51), 
                        exclude_percentile=95, fill_value=0.0):
    '''
    Function to subtract the background from a FITS file and write the output and the background to new FITS files.
    
    Parameters:
    hdu_hst: HDUList
        Input HDUList from which the background is to be subtracted.
    halpha_filename: str
        The original filename of the file containing hdu_hst. This is used to generate the filenames for the output files.
    
    Optional Parameters:
    sigma_clip_sigma: float (default=3.0)
        The number of standard deviations to use for the sigma clipping threshold.
    detect_threshold_nsigma: float (default=2.0)
        The number of standard deviations above the background for a pixel to be considered a source.
    detect_sources_npixels: int (default=10)
        The number of connected pixels that an object needs to have to be considered a source.
    footprint_radius: int (default=10)
        The radius of the circular footprint to use for the source detection.
    box_size: tuple (default=(100, 100))
        The size of the box to use for the background estimation.
    filter_size: tuple (default=(51, 51))
        The size of the box to use for the background smoothing.
    exclude_percentile: float (default=95)
        Exclude sources from the background estimation that have values in the top percentile.
    fill_value: float (default=0.0)
        The value to use for pixels that are masked or outside of the image.

    Returns: None
    '''

    # Extracting the data from the hdu object
    data = hdu_hst.data

    # Defining the sigma clipping parameters
    sigma_clip = SigmaClip(sigma=sigma_clip_sigma, maxiters=10)

    # Calculating the detection threshold for the data
    threshold = detect_threshold(data, nsigma=detect_threshold_nsigma, sigma_clip=sigma_clip)

    # Detecting sources in the data above the threshold
    segment_img = detect_sources(data, threshold, npixels=detect_sources_npixels)

    # Creating a circular footprint for source detection
    footprint = circular_footprint(radius=footprint_radius)

    # Making a source mask using the footprint
    mask = segment_img.make_source_mask(footprint=footprint)

    # Creating a coverage mask to handle NaN values
    coverage_mask = ~np.isnan(data)

    # Creating a background estimator
    bkg = Background2D(data, box_size=box_size, filter_size=filter_size, mask=mask, 
                       exclude_percentile=exclude_percentile, fill_value=fill_value)

    # Getting the estimated background and subtracting it from the data
    data_bg = bkg.background
    data_bgsub = data - data_bg

    # Create a mask where True values represent NaN values in the original data
    mask = np.isnan(data)

    # Apply the mask to the background data, replacing detected positions (True) with NaNs
    data_bg[mask] = np.nan
    
    # Creating copies of the hdu object
    hdu_hst_bg = hdu_hst.copy()
    hdu_hst_bgsub = hdu_hst.copy()

    # Replacing the data in the copied hdu objects with the background and the background subtracted data
    hdu_hst_bg.data = data_bg
    hdu_hst_bgsub.data = data_bgsub

    # Defining the output filenames and writing the output to the files
    output_filename = halpha_filename.replace('_raw.fits', '_bg.fits')
    hdu_hst_bg.writeto(output_filename, overwrite=True)

    output_filename = halpha_filename.replace('_raw.fits', '_bgsub.fits')
    hdu_hst_bgsub.writeto(output_filename, overwrite=True)
    
    return(hdu_hst_bgsub)


def process_halpha_muse(input_muse_filename, output_muse_filename):
    """
    Process the Halpha MUSE file.

    Args:
        input_muse_filename (str): Filename of the input MUSE file.
        output_muse_filename (str): Filename of the output MUSE file.

    Returns:
        astropy.io.fits.ImageHDU: HDU containing the Halpha flux data.
    """
    print(f"Processing Halpha MUSE file: {input_muse_filename}")  # Print a message indicating the file being processed

    # Open the MUSE file and extract the HDU with Halpha flux data
    hdu_muse_ha = fits.open(input_muse_filename)['HA6562_FLUX']

    # Create a new HDU with the Halpha data and header
    hdu_muse_ha = fits.PrimaryHDU(hdu_muse_ha.data, hdu_muse_ha.header)

    # Write the new HDU to the output file
    hdu_muse_ha.writeto(output_muse_filename, overwrite=True)

    # Return the Halpha flux HDU
    return hdu_muse_ha


#[old code]
# def process_halpha_muse(input_halpha_filename, input_nebulae_mask_filename, input_nebulae_catalog_filename,
#                         input_nebulae_mask_hst_filename, hdu_hst, galaxy):
#     """
#     Processes the Halpha MUSE FITS file by reprojecting the nebulae mask.

#     Args:
#         input_halpha_filename (str): Path to the input Halpha MUSE FITS file.
#         input_nebulae_mask_filename (str): Path to the input nebulae mask FITS file.
#         input_nebulae_catalog_filename (str): Path to the input nebulae catalog FITS file.
#         output_nebulae_mask_filename (str): Path to save the processed nebulae mask FITS file.

#     Returns:
#         None
#     """
#     print(f"Processing Halpha MUSE file: {input_halpha_filename}")
#     hdu_muse = fits.open(input_halpha_filename)[1]

#     print(f"Processing nebulae mask: {input_nebulae_mask_filename}")
#     hdu_nebmask_muse = fits.open(input_nebulae_mask_filename)[0]

#     print(f"Reading nebulae catalog: {input_nebulae_catalog_filename}")
#     table_nebcat = Table.read(input_nebulae_catalog_filename)
#     table_nebcat = table_nebcat[table_nebcat['gal_name'] == galaxy.upper()]

#     if not os.path.isfile(input_nebulae_mask_hst_filename):
#         print("Reprojecting nebulae mask...")
#         data_nebmask_hst, _ = reproject.reproject_interp(hdu_nebmask_muse, hdu_hst.header, order=0, parallel=True)
#         hdu_nebmask_hst = fits.PrimaryHDU(data_nebmask_hst, hdu_hst.header)
#         hdu_nebmask_hst.writeto(output_nebulae_mask_filename, overwrite=True)
#         print(f"Processing complete. Saving the processed nebulae mask as {output_nebulae_mask_filename}")
#     else:
#         print("Opening existing nebulae mask...")
#         hdu_nebmask_hst = fits.open(input_nebulae_mask_hst_filename)[0]

#     return(hdu_muse, hdu_nebmask_muse, table_nebcat, hdu_nebmask_hst)


def regrid(hdu_input, hdu_template, output_filename=None, conserve_flux=True):
    """
    Reprojects an input FITS image to match the WCS of a template FITS image.

    Args:
        hdu_input (astropy.io.fits.ImageHDU): Input FITS image HDU.
        hdu_template (astropy.io.fits.ImageHDU): Template FITS image HDU.
        output_filename (str, optional): Path to save the reprojected image as a new FITS file.
                                         Defaults to None.
        conserve_flux (bool, optional): Flag to conserve flux during reprojection. 
                                       Defaults to True.

    Returns:
        astropy.io.fits.ImageHDU: Reprojected FITS image HDU.
    """
    print("[INFO] Reprojecting the input image to match the template WCS...")

    # Extract the WCS information from the input and template headers
    wcs_input = wcs.WCS(hdu_input.header)
    wcs_template = wcs.WCS(hdu_template.header)

    # Calculate the pixel scale for input and template images
    pixscale_input = wcs.utils.proj_plane_pixel_area(wcs_input.celestial)
    pixscale_template = wcs.utils.proj_plane_pixel_area(wcs_template.celestial)

    # Reproject the input image to match the template WCS
    print("[INFO] Performing image reprojection...")
    data_output = reproject_interp(hdu_input, hdu_template.header, order=0, parallel=True)[0]
    hdu_output = fits.PrimaryHDU(data_output, hdu_template.header)
    print("[INFO] Image reprojection complete.")

    if conserve_flux:
        # Scale the output data to conserve flux (only if template pixel size is smaller than input)
        print(f"[INFO] Scaling the output data to conserve flux with factor {(pixscale_template / pixscale_input):.2f}")
        hdu_output.data = hdu_output.data * (pixscale_template / pixscale_input)
        print("[INFO] Flux scaling complete.")

    if output_filename is not None:
        # Save the reprojected image to a new FITS file
        print(f"[INFO] Saving the reprojected image to: {output_filename}")
        hdu_output.writeto(output_filename, overwrite=True)
        print("[INFO] Image saved successfully.")

    print("[INFO] Reprojection process completed.")
    return hdu_output


def smooth_image_with_beam(input_hdu, initial_resolution, desired_resolution, output_filename=None):
    """
    Smooths an input image with a beam kernel to match a desired resolution.

    Args:
        input_hdu (astropy.io.fits.ImageHDU): Input FITS image HDU.
        initial_resolution (astropy.units.Quantity): Initial resolution of the input image.
        desired_resolution (astropy.units.Quantity): Desired resolution for smoothing.
        output_file (str, optional): Path to save the smoothed image as a new FITS file.
                                    Defaults to None.

    Returns:
        astropy.io.fits.ImageHDU: Smoothed FITS image HDU.
    """
    print("[INFO] Smoothing the image with a beam kernel...")
    
    # Create a WCS object from the input HDU header
    wcs_ = wcs.WCS(input_hdu.header)

    # Calculate the pixel scale in degrees
    pixscale = wcs.utils.proj_plane_pixel_area(wcs_.celestial) ** 0.5 * u.deg
    print(f"[INFO] Pixel scale: {pixscale.to('arcsec'):.2f} arcsec")

    # Define the initial and desired beams
    initial_beam = Beam(initial_resolution)
    desired_beam = Beam(desired_resolution)

    print(f"[INFO] Initial Resolution: {initial_resolution.to('arcsec'):.2f} arcsec")
    print(f"[INFO] Desired Resolution: {desired_resolution.to('arcsec'):.2f} arcsec")
    
    # Create the convolution kernel
    convolution_kernel = desired_beam.deconvolve(initial_beam).as_kernel(pixscale)
    print("[INFO] Convolution kernel created.")

    # Convolve the image with the kernel to smooth it
    print("[INFO] Performing image convolution...")
    smoothed_data = convolve_fft(input_hdu.data, convolution_kernel, preserve_nan=True)
    print("[INFO] Image convolution complete.")

    output_hdu = fits.PrimaryHDU(smoothed_data, input_hdu.header)

    if output_filename is not None:
        # Save the smoothed image to a new FITS file
        print(f"[INFO] Saving the smoothed image to: {output_filename}")
        output_hdu.writeto(output_filename, overwrite=True)
        print("[INFO] Image saved successfully.")

    print("[INFO] Smoothing process completed.")
    return output_hdu


def scale_data_2D(data):
    '''
    Function to scale 2D data to the range [0, 1], ignoring NaN values.
    
    Parameters:
    data : ndarray
        The 2D array to be scaled.
    
    Returns:
    scaled_data : ndarray
        The scaled 2D array.
    '''
    # Define a mask for NaN values
    mask = np.isnan(data)
    
    # Get the minimum and maximum values from the data, ignoring NaNs
    data_min = np.nanmin(data)
    data_max = np.nanmax(data)
    
    # Create a copy of the data to avoid modifying the input array
    scaled_data = data.copy()
    
    # Scale the data using the min and max, applying the mask to ignore NaNs
    scaled_data[~mask] = (data[~mask] - data_min) / (data_max - data_min)
    
    # Preserve the NaNs in the scaled data
    scaled_data[mask] = np.nan

    return scaled_data


def save_diff_ratio_smoothed_image(hdu_muse_regrid, hdu_hst, hdu_hst_smoothed, output_ratio_filename, output_diff_filename, output_ratio_anchored_filename, output_diff_anchored_filename):
    """
    Calculates the ratio of a MUSE regridded image to a smoothed HST image,
    saves the ratio image and the anchored HST image as FITS files.

    Args:
        hdu_muse_regrid (astropy.io.fits.ImageHDU): MUSE regridded image HDU.
        hdu_hst (astropy.io.fits.ImageHDU): HST image HDU.
        hdu_hst_smoothed (astropy.io.fits.ImageHDU): Smoothed HST image HDU.
        output_ratio_filename (str): Path to save the ratio image as a FITS file.
        output_anchored_filename (str): Path to save the anchored HST image as a FITS file.

    Returns:
        Tuple: Contains the ratio smoothed image HDU and the anchored HST image HDU.
    """

    # Calculate the ratio of the MUSE regridded image to the smoothed HST image
    # ratio_smooth = hdu_muse_regrid.data / hdu_hst_smoothed.data
    
    # Scale data for ratio between 0 and 1
    data_hst_scaled = scale_data_2D(hdu_hst_smoothed.data) 
    data_muse_scaled = scale_data_2D(hdu_muse_regrid.data)
    
    # Ratio will not work without - needs to be fixed
    scale_factor = np.abs(np.nanpercentile(hdu_hst.data, 0.1))

    ratio_smooth = hdu_muse_regrid.data/(hdu_hst_smoothed.data+scale_factor)
    # ratio_smooth = data_muse_scaled/data_hst_scaled

    diff_smooth = hdu_muse_regrid.data - hdu_hst_smoothed.data

    # Create a copy of the smoothed HST image HDU and update its data with the ratio
    hdu_ratio_smooth = hdu_hst_smoothed.copy()
    hdu_diff_smooth = hdu_hst_smoothed.copy()

    hdu_ratio_smooth.data = ratio_smooth
    hdu_diff_smooth.data = diff_smooth

    # Save the ratio smoothed image to the output file
    hdu_ratio_smooth.writeto(output_ratio_filename, overwrite=True)
    hdu_diff_smooth.writeto(output_diff_filename, overwrite=True)

    print(f"[INFO] Ratio smoothed image saved as {output_ratio_filename}")
    print(f"[INFO] Diff smoothed image saved as {output_diff_filename}")

    # Create a copy of the HST image HDU and update its data with the anchored values using the ratio
    hdu_hst_ratio_anchored = hdu_hst.copy()
    hdu_hst_diff_anchored = hdu_hst.copy()

    hdu_hst_ratio_anchored.data = (hdu_hst.data+scale_factor) * ratio_smooth
    hdu_hst_diff_anchored.data = hdu_hst.data + diff_smooth

    # Save the anchored HST image to the output file
    hdu_hst_ratio_anchored.writeto(output_ratio_anchored_filename, overwrite=True)
    hdu_hst_diff_anchored.writeto(output_diff_anchored_filename, overwrite=True)

    print(f"[INFO] Anchored HST image saved as {output_ratio_anchored_filename}")
    print(f"[INFO] Anchored HST image saved as {output_diff_anchored_filename}")

    return hdu_ratio_smooth, hdu_diff_smooth, hdu_hst_ratio_anchored, hdu_hst_diff_anchored


def process_intnegs_anchored_image(hdu_hst_anchored, output_filename):
    """
    Process the anchored HST image by replacing negative values with NaNs and
    performing interpolation using Gaussian kernels.

    Args:
        hdu_hst_anchored (astropy.io.fits.ImageHDU): Anchored HST image HDU.
        output_filename (str): Filename for the processed anchored HST image.

    Returns:
        astropy.io.fits.ImageHDU: Processed anchored HST image HDU.
    """
    # Create a copy of the anchored HST image HDU
    hdu_hst_anchored_intnegs = hdu_hst_anchored.copy()

    # Replace negative values with NaNs
    mask = hdu_hst_anchored_intnegs.data <= 0
    mask = binary_dilation(mask, iterations=2)
    hdu_hst_anchored_intnegs.data[mask] = np.nan

    # Perform interpolation using Gaussian kernels
    kernel = Gaussian2DKernel(x_stddev=1)
    hdu_hst_anchored_intnegs.data = interpolate_replace_nans(hdu_hst_anchored_intnegs.data, kernel)
    hdu_hst_anchored_intnegs.data = interpolate_replace_nans(hdu_hst_anchored_intnegs.data, kernel)

    # Save the processed anchored HST image to the output file
    hdu_hst_anchored_intnegs.writeto(output_filename, overwrite=True)
    print(f"[INFO] Anchored HST image with negative values processed and saved as {output_filename}")

    return hdu_hst_anchored_intnegs


def create_2d_hist_and_fit(hdu_input1, hdu_input2, hdu_input3, output_filename, showplot=False):
    '''
    This function creates a 2D histogram from two input HDUs (Header Data Unit). 
    It then fits a line to the data and applies the fitted line to modify the 
    second HDU. The modified HDU is then written to a new FITS file.

    Parameters:
    hdu_input1: astropy.io.fits.hdu.image.PrimaryHDU
        The first input HDU - hdu_muse_regrid
    hdu_input2: astropy.io.fits.hdu.image.PrimaryHDU
        The second input HDU - hdu_hst_smoothed
    hdu_input3: astropy.io.fits.hdu.image.PrimaryHDU
        The third input HDU which will be modified by the function.
    output_filename: str
        The name of the output FITS file to which the modified HDU will be written.

    Returns:
    hdu_fit_anchored: astropy.io.fits.hdu.image.PrimaryHDU
        The second input HDU modified by the line of best fit.
    '''

    # Filter out NaN values from the data
    valid_indices = np.isfinite(hdu_input1.data) & np.isfinite(hdu_input2.data)
    x_data = hdu_input1.data[valid_indices]
    y_data = hdu_input2.data[valid_indices]

    # Calculate a line of best fit for the data
    slope, intercept = np.polyfit(x_data, y_data, 1)
    x_fit = np.linspace(np.min(x_data), np.max(x_data), 100)
    y_fit = slope * x_fit + intercept

    if showplot: 
        # Create a 2D histogram plot using the filtered data
        fig, ax = plt.subplots(figsize=(8, 8))
        hist = ax.hist2d(x_data, y_data, bins=50, cmap='turbo', norm=LogNorm())
        cbar = plt.colorbar(hist[3], ax=ax, label='Counts')
        ax.set_xlabel('Input1 Data')
        ax.set_ylabel('Input2 Data')
        ax.set_title('2D Histogram Plot')
        ax.plot(x_fit, y_fit, color='red', linewidth=2, label=f'y = {slope:.2f}x + {intercept:.2f}')
        ax.legend()

    # Apply the calculated line of best fit to the second input data and save it as a new FITS file
    hdu_fit_anchored = hdu_input3.copy()
    hdu_fit_anchored.data = (hdu_fit_anchored.data - intercept) / slope
    hdu_fit_anchored.writeto(output_filename, overwrite=True)

    # Return the modified HDU object
    return(hdu_fit_anchored)


def process_anchored_fit_image(hdu_hst_anchored, output_filename, sigma=5):
    """
    Process the anchored HST image by replacing negative values with NaNs and
    performing interpolation using Gaussian kernels.

    Args:
        hdu_hst_anchored (astropy.io.fits.ImageHDU): Anchored HST image HDU.
        output_filename (str): Filename for the processed anchored HST image.

    Returns:
        astropy.io.fits.ImageHDU: Processed anchored HST image HDU.
    """
    # Create a copy of the anchored HST image HDU
    hdu_hst_anchored_intnegs = hdu_hst_anchored.copy()

    # Replace negative values with NaNs
    mask = hdu_hst_anchored_intnegs.data <= 0
    std = mad_std(hdu_hst_anchored_intnegs.data[mask], ignore_nan=True)
    mask = hdu_hst_anchored_intnegs.data <= -std*sigma
    mask = binary_dilation(mask, iterations=2)
    hdu_hst_anchored_intnegs.data[mask] = np.nan

    # Perform interpolation using Gaussian kernels
    kernel = Gaussian2DKernel(x_stddev=1)
    hdu_hst_anchored_intnegs.data = interpolate_replace_nans(hdu_hst_anchored_intnegs.data, kernel)
    hdu_hst_anchored_intnegs.data = interpolate_replace_nans(hdu_hst_anchored_intnegs.data, kernel)

    # Save the processed anchored HST image to the output file
    hdu_hst_anchored_intnegs.writeto(output_filename, overwrite=True)
    print(f"[INFO] Anchored HST image with negative values processed and saved as {output_filename}")

    return hdu_hst_anchored_intnegs



def add_noise_to_image(hdu_input, output_filename):
    """
    Adds random noise to the input FITS image.

    Args:
        hdu_input (astropy.io.fits.ImageHDU): Input FITS image HDU.
        output_filename (str): Path to save the output FITS file with added noise.

    Returns:
        astropy.io.fits.ImageHDU: FITS image with added noise.
    """
    # Calculate the RMS of the input image data
    rms_1 = mad_std(hdu_input.data, ignore_nan=True)
    rms_2 = mad_std(hdu_input.data[hdu_input.data < 5 * rms_1], ignore_nan=True)

    # Generate random noise with the same shape as the input image
    noise = np.random.normal(0, rms_2, hdu_input.data.shape)

    # Create a copy of the input HDU and add the noise to the data
    hdu_wnoise = hdu_input.copy()
    hdu_wnoise.data = hdu_input.data + noise

    # Save the image with added noise to the output file
    hdu_wnoise.writeto(output_filename, overwrite=True)
    print(f"[INFO] Image with added noise saved as {output_filename}")

    return hdu_wnoise
