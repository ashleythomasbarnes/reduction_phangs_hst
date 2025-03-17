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


def get_smooth(hdu, initial_resolution, desired_resolution):
    
    # Create a WCS object from the input HDU header
    wcs_ = wcs.WCS(hdu.header)

    # Calculate the pixel scale in degrees
    pixscale = wcs.utils.proj_plane_pixel_area(wcs_.celestial) ** 0.5 * u.deg
    print(f"[INFO] Pixel scale: {pixscale.to('arcsec'):.2f} arcsec")

    # Define the initial and desired beams
    initial_beam = Beam(initial_resolution)
    desired_beam = Beam(desired_resolution)

    print(f"[INFO] Initial Resolution: {initial_resolution.to('arcsec'):.2f} arcsec")
    print(f"[INFO] Desired Resolution: {desired_resolution.to('arcsec'):.2f} arcsec")
    
    # Create the convolution kernel
    convolution_beam = (desired_resolution.to('arcsec')**2 - initial_resolution.to('arcsec')**2)**0.5
    convolution_kernel = desired_beam.deconvolve(initial_beam).as_kernel(pixscale)
    print(f"[INFO] Convolution kernel: {convolution_beam.to('arcsec'):.2f} arcsec")

    # Convolve the image with the kernel to smooth it
    print("[INFO] Performing image convolution...")
    smoothed_data = convolve_fft(hdu.data, convolution_kernel, preserve_nan=True, allow_huge=True)
    print("[INFO] Image convolution complete.")

    output_hdu = fits.PrimaryHDU(np.array(smoothed_data, dtype=np.float32), hdu.header)

    print("[INFO] Smoothing process completed.")
    return(output_hdu)


def get_regrid(hdu_input, hdu_template, output_filename=None, conserve_flux=True, order='bilinear'):

    print("[INFO] Reprojecting the input image to match the template WCS...")

    # Extract the WCS information from the input and template headers
    header_input = hdu_input.header
    header_template = hdu_template.header
    wcs_input = wcs.WCS(header_input)
    wcs_template = wcs.WCS(header_template)

    # Calculate the pixel scale for input and template images
    pixscale_input = wcs.utils.proj_plane_pixel_area(wcs_input.celestial)
    pixscale_template = wcs.utils.proj_plane_pixel_area(wcs_template.celestial)

    # Reproject the input image to match the template WCS
    print("[INFO] Performing image reprojection...")
    data_output = reproject_interp(hdu_input, header_template, order=order)[0]
    print("[INFO] Image reprojection complete.")

    keys = ['NAXIS', 'NAXIS1', 'NAXIS2', 'CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']
    for key in keys: 
        if key in header_template.keys():
            header_input[key] = header_template[key]
    hdu_output = fits.PrimaryHDU(data_output, header_input)

    if conserve_flux:
        # Scale the output data to conserve flux 
        print(f"[INFO] Scaling the output data to conserve flux with factor {(pixscale_template / pixscale_input):.2f}")
        hdu_output.data = hdu_output.data * (pixscale_template / pixscale_input)
        hdu_output.data = np.array(hdu_output.data, dtype=np.float32)
        print("[INFO] Flux scaling complete.")

    print("[INFO] Reprojection process completed.")
    return(hdu_output)