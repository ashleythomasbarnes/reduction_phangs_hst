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


def process_halpha_units(input_filename, output_filename):
    """
    Processes the Halpha FITS file by changing the units to erg/s/cm^2.

    Args:
        input_filename (str): Path to the input Halpha FITS file.
        output_filename (str): Path to save the processed FITS file.

    Returns:
        None
    """
    print(f"Processing {input_filename}")
    hdu_hst = fits.open(input_filename)[0]

    photflam = hdu_hst.header['PHOTFLAM']
    photplam = hdu_hst.header['PHOTPLAM']
    photbw = hdu_hst.header['PHOTBW']

    print("Scaling the data...")
    hdu_hst.data = hdu_hst.data * photflam * photbw

    hdu_hst.header['BUNIT'] = ('erg/s/cm2', '1e-20 erg/s/cm2')

    print("Converting units...")
    hdu_hst.data = hdu_hst.data * 1e20

    print(f"Saving the processed FITS file as {output_filename}")
    hdu_hst.writeto(output_filename, overwrite=True)
    
    return(hdu_hst)


def process_halpha_muse(input_halpha_filename, input_nebulae_mask_filename, input_nebulae_catalog_filename,
                        input_nebulae_mask_hst_filename, hdu_hst, galaxy):
    """
    Processes the Halpha MUSE FITS file by reprojecting the nebulae mask.

    Args:
        input_halpha_filename (str): Path to the input Halpha MUSE FITS file.
        input_nebulae_mask_filename (str): Path to the input nebulae mask FITS file.
        input_nebulae_catalog_filename (str): Path to the input nebulae catalog FITS file.
        output_nebulae_mask_filename (str): Path to save the processed nebulae mask FITS file.

    Returns:
        None
    """
    print(f"Processing Halpha MUSE file: {input_halpha_filename}")
    hdu_muse = fits.open(input_halpha_filename)[1]

    print(f"Processing nebulae mask: {input_nebulae_mask_filename}")
    hdu_nebmask_muse = fits.open(input_nebulae_mask_filename)[0]

    print(f"Reading nebulae catalog: {input_nebulae_catalog_filename}")
    table_nebcat = Table.read(input_nebulae_catalog_filename)
    table_nebcat = table_nebcat[table_nebcat['gal_name'] == galaxy.upper()]

    if not os.path.isfile(input_nebulae_mask_hst_filename):
        print("Reprojecting nebulae mask...")
        data_nebmask_hst, _ = reproject.reproject_interp(hdu_nebmask_muse, hdu_hst.header, order=0, parallel=True)
        hdu_nebmask_hst = fits.PrimaryHDU(data_nebmask_hst, hdu_hst.header)
        hdu_nebmask_hst.writeto(output_nebulae_mask_filename, overwrite=True)
        print(f"Processing complete. Saving the processed nebulae mask as {output_nebulae_mask_filename}")
    else:
        print("Opening existing nebulae mask...")
        hdu_nebmask_hst = fits.open(input_nebulae_mask_hst_filename)[0]

    return(hdu_muse, hdu_nebmask_muse, table_nebcat, hdu_nebmask_hst)


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


def save_ratio_smoothed_image(hdu_muse_regrid, hdu_hst_smoothed, output_ratio_filename, output_anchored_filename):
    """
    Calculates the ratio of a MUSE regridded image to a smoothed HST image,
    saves the ratio image and the anchored HST image as FITS files.

    Args:
        hdu_muse_regrid (astropy.io.fits.ImageHDU): MUSE regridded image HDU.
        hdu_hst_smoothed (astropy.io.fits.ImageHDU): Smoothed HST image HDU.
        output_filename (str): Path to save the ratio image as a FITS file.
    """
    ratio_smooth = hdu_muse_regrid.data / hdu_hst_smoothed.data

    hdu_ratio_smooth = hdu_hst_smoothed.copy()
    hdu_ratio_smooth.data = ratio_smooth

    hdu_ratio_smooth.writeto(output_ratio_filename, overwrite=True)
    print(f"[INFO] Ratio smoothed image saved as {output_ratio_filename}")

    hdu_hst_anchored = hdu_hst_smoothed.copy()
    hdu_hst_anchored.data = hdu_hst_smoothed.data * ratio_smooth
    hdu_hst_anchored.writeto(output_anchored_filename, overwrite=True)
    print(f"[INFO] Anchored HST image saved as {output_anchored_filename}")
    
    return(hdu_ratio_smooth, hdu_hst_anchored)


def process_anchored_image(hdu_hst_anchored, output_filename):
    """
    Process the anchored HST image by replacing negative values with NaNs and
    performing interpolation using Gaussian kernels.

    Args:
        hdu_hst_anchored (astropy.io.fits.ImageHDU): Anchored HST image HDU.
        galaxy (str): Name of the galaxy.

    Returns:
        astropy.io.fits.ImageHDU: Processed anchored HST image HDU.
    """
    hdu_hst_anchored_intnegs = hdu_hst_anchored.copy()
    mask = hdu_hst_anchored_intnegs.data <= 0
    mask = binary_dilation(mask, iterations=2)
    hdu_hst_anchored_intnegs.data[mask] = np.nan

    kernel = Gaussian2DKernel(x_stddev=1)
    hdu_hst_anchored_intnegs.data = interpolate_replace_nans(hdu_hst_anchored_intnegs.data, kernel)
    hdu_hst_anchored_intnegs.data = interpolate_replace_nans(hdu_hst_anchored_intnegs.data, kernel)
    hdu_hst_anchored_intnegs.writeto(output_filename, overwrite=True)
    print(f"[INFO] Anchored HST image with negative values processed and saved as {output_filename}")

    return (hdu_hst_anchored_intnegs)


def add_noise_to_image(hdu_input, output_filename):
    """
    Adds random noise to the input FITS image.

    Args:
        hdu_input (astropy.io.fits.ImageHDU): Input FITS image HDU.
        output_filename (str): Path to save the output FITS file with added noise.

    Returns:
        astropy.io.fits.ImageHDU: FITS image with added noise.
    """
    rms_1 = mad_std(hdu_input.data, ignore_nan=True)
    rms_2 = mad_std(hdu_input.data[hdu_input.data < 5 * rms_1], ignore_nan=True)
    noise = np.random.normal(0, rms_2, hdu_input.data.shape)

    hdu_wnoise = hdu_input.copy()
    hdu_wnoise.data = hdu_input.data + noise
    hdu_wnoise.writeto(output_filename, overwrite=True)
    print(f"[INFO] Image with added noise saved as {output_filename}")

    return hdu_wnoise