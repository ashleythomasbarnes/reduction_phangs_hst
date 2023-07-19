# Continuum subtraction and scaling for HST Ha data using MUSE data

## Code 1: run_pipeline_contsub.ipynb

The continuum_subtraction.py code performs continuum subtraction on HST data and post-processes the resulting images. It takes user inputs such as the galaxy name and file paths for the H-alpha, continuum, and MUSE data. The code generates several output files, including continuum-subtracted H-alpha and scaled continuum images. It also processes the H-alpha units and MUSE data, performs regridding, smoothing, and generates ratio images. Additionally, the code applies various transformations to the anchored HST image, such as intensity negations and adding white noise. The output files are saved in the specified output directory.

## Code 2: run_pipeline_contsub_qa.ipynb

The qa_analysis.py code performs quality analysis (QA) on HST and MUSE data. It takes user inputs such as the galaxy name, file paths for H-alpha, continuum, MUSE, nebulae mask, and nebulae catalog. The code reads FITS files, creates a dictionary of FITS data, and processes the MUSE catalog. It calculates the flux values for nebulae in the HST and MUSE data and generates QA analysis plots. The output of the code includes QA analysis plots that provide insights into the quality of the HST and MUSE data.

These codes are part of the reduction_phangs_hst package and utilize the astropy and warnings packages for various data processing and analysis tasks.

## Code 1: Continuum Subtraction

This script is designed to process images of galaxies in the H-alpha emission line using data from both the Hubble Space Telescope (HST) and MUSE (Multi Unit Spectroscopic Explorer) on the Very Large Telescope (VLT).

1. **User-Defined Inputs**: The first thing you'll see in the script are some user-defined inputs. These are variables that you need to update based on the specific data you're using.

    - `galaxy`: The name of the galaxy that you're working with. You would replace `'ngc628'` with the name of your chosen galaxy.
    - `halpha_inputfilename`: The path to the file containing the H-alpha data.
    - `cont1_inputfilename`, `cont2_inputfilename`: The paths to the files containing the continuum data.
    - `input_muse_filename`: The path to the MUSE data file.

2. **Library Imports**: This part of the script imports the necessary Python libraries and modules that will be used later in the script.

3. **Output Directory Creation**: The script extracts some substrings from the input file names to use in creating a directory where output files will be stored.

4. **Continuum Subtraction**: The next section of the script is dedicated to subtracting the continuum from the H-alpha emission line. The resulting files are saved into the output directory.

The Tasks here are as follows: 

The above script performs a continuum subtraction on an emission line image. Continuum subtraction is a technique used in astrophysics to isolate the light from a specific astronomical object in a field of view. In an image, the light from the object of interest (an emission line) is typically mixed with the "background" light (the continuum). The purpose of this script is to isolate the light from the emission line by subtracting the continuum.

The script starts by importing necessary modules, such as numpy and astropy.io for handling FITS files, and defines three functions. FITS (Flexible Image Transport System) is a digital file format useful for storing, transmitting, and manipulating scientific data.

The first function, `perform_continuum_subtraction_onecont`, subtracts a single continuum image from an emission line image. The function retrieves the photflam (photometric flux density) and photplam (pivot wavelength) for both the H-alpha emission line and the continuum images from their headers. It then calculates the AB magnitude zero point for both images, scales the continuum data, subtracts the scaled continuum from the H-alpha image, and writes the output to a file.

The second function, `perform_continuum_subtraction_twocont`, performs a similar task but for two continuum images. It retrieves the photflam and photplam for both continuum images and calculates their zero points. Then, it calculates weights for each continuum image based on their distances from the H-alpha pivot wavelength. The function then scales both continuum images using the calculated weights and combines them. It subtracts this combined continuum from the H-alpha image and writes the output to a file.

The third function, `continuum_subtraction`, is a wrapper function for the two previous functions. It opens the H-alpha and continuum FITS files, replaces any zero values with NaNs to prevent any computational issues (like division by zero), and calls the appropriate continuum subtraction function based on whether one or two continuum images are provided. It also prints some informational messages about the files being processed.

This script would be used in the context of reducing observational data, where the observer has obtained images of the sky in various filters, including some centered on specific emission lines (like H-alpha) and some meant to capture the broader spectral continuum. By performing this continuum subtraction, the observer can isolate the light specifically emitted at the H-alpha wavelength, which can provide information about the physical conditions in the observed objects. This kind of analysis is common in many areas of astrophysics, including studies of nebulae, galaxies, and the interstellar medium.

**The next part of the work flows is as follows:** 

5. **Post-Processing**: The next few blocks are involved in post-processing the H-alpha data, including unit conversion, background subtraction, and MUSE data processing. At each step, the processed data is saved into the output directory.

6. **Image Smoothing**: The HST image is then smoothed with the desired resolution. This involves determining the initial and desired resolutions, and applying a smoothing function to the HST data.

7. **Difference and Ratio Images**: This part of the script creates difference and ratio images from the HST and MUSE data, then saves the output. 

8. **Intensity Anomalies**: Next, the script deals with 'intensity negations' in the anchored HST images, which may arise due to negative pixel values.

9. **2D Histogram and Fit**: After that, the script creates a 2D histogram and fit for the provided HDU objects, and saves the modified second HDU object.

10. **Cosmic Ray Removal**: The last part of the script involves finding and removing cosmic rays from the images using a module dedicated to this task.

The Tasks here are as follows: 

`process_halpha_units(input_filename, output_filename)` This function opens a FITS (Flexible Image Transport System) file, which is a format commonly used in astronomy for storing image data, and changes the unit of measurement for the image data to erg/s/cm^2. It uses keywords from the file's header to apply the necessary scaling and conversion to the data. The file is then saved with the updated data and units. The function returns an `ImageHDU` object, which represents the primary data structure of the FITS file.

`process_halpha_background(...)` This function subtracts the background from a FITS file. It begins by extracting the data from the provided HDU (Header Data Unit) object, which is part of a FITS file. Then it sets up a sigma clipping object and calculates the detection threshold for the data. Sources (regions of the image where the data value exceeds the threshold) are then detected, and a mask is created that covers these sources. A mask is also created to cover regions of the data where the value is `NaN` (Not a Number). A background estimator is then created and used to calculate an estimate of the background, which is subtracted from the original data. The resulting background-subtracted data is then written to a new FITS file.

`process_halpha_muse(input_muse_filename, output_muse_filename)` This function processes a MUSE (Multi Unit Spectroscopic Explorer) file, which is another type of data file used in astronomy. The function opens the file and extracts the Halpha flux data. This data is then written to a new FITS file.

`regrid(hdu_input, hdu_template, output_filename=None, conserve_flux=True)` This function reprojects an input image to match the World Coordinate System (WCS) of a template image. WCS is a system used in astronomy that defines the physical coordinates for each pixel in an image. This reprojection can involve a change in the size, shape, and/or orientation of the input image.  The function optionally scales the output data to conserve flux (total light energy), a process that's necessary when the reprojection involves a change in pixel size. The function returns the reprojected image as an `ImageHDU` object, and optionally saves it to a new FITS file.

`smooth_image_with_beam(input_hdu, initial_resolution, desired_resolution, output_filename=None)`

`save_diff_ratio_smoothed_image(...)` function calculates the ratio of a MUSE regridded image to a smoothed HST image. It saves both the ratio and the difference between images as new FITS files. It also creates anchored images, which is the original image multiplied by the ratio or added by the difference, and saves them as well.

`process_intnegs_anchored_image(...)` function processes an anchored image by replacing negative values with NaN and then performs interpolation on these NaN values using a Gaussian kernel.

`create_2d_hist_and_fit(...)` creates a 2D histogram from two input HDUs, fits a line to the data, and applies the fitted line to modify the second HDU. It then writes this modified HDU to a new FITS file.

The last function `process_anchored_fit_image(...)` is very similar to `process_intnegs_anchored_image(...)`. It processes the anchored image by replacing the negative values with NaN and performs interpolation on these NaN values using a Gaussian kernel. However, this version uses the standard deviation calculated from the Median Absolute Deviation (MAD) to generate a mask to identify negative outliers rather than a simple less than or equal to zero condition.



To run this script on your computer, follow these steps:

1. **Setup Python Environment**: Make sure you have Python installed. If not, install Python (version 3.7 or later is recommended). This script uses several external libraries (`os`, `astropy`, `warnings`, and the specific module `reduction_phangs_hst`), so you may need to install them with pip or conda. For example, `pip install astropy`.

2. **Get the Data**: You will need to have the correct data files for the galaxy you want to analyze. The data should be in `.fits` format and should include H-alpha data, two continuum data files, and a MUSE data file. Update the file paths in the user input section to point to these files.

3. **Get the `reduction_phangs_hst` Module**: This script uses functions from the `reduction_phangs_hst` module, which is not a standard Python library. This might be a module provided by your research team or mentor. If you don't have this module, you'll need to obtain it and ensure it's in your Python path.

4. **Run the Script**: Once everything is set up, you can run the script in a Python environment. This could be through an integrated development environment (IDE) like PyCharm or VSCode, a notebook interface like Jupyter, or directly in a command line interface.

5. **Check the Results**: The script will output several processed `.fits` files into the output directory, which you can analyze as needed.

Remember that this script is specific to the data and tasks involved in this particular galaxy image analysis, and may need to be modified to suit different data or analysis goals.

## Prerequisites

- Python 3
- `reduction_phangs_hst` package
- `astropy` package
- `warnings` package

## Usage

1. Update the user inputs section in the code with the desired values for the following variables:
   - `galaxy`: Define the galaxy name.
   - `halpha_inputfilename`: Path to the H-alpha input file.
   - `cont2_inputfilename`: Path to the second continuum input file.
   - `cont1_inputfilename`: Path to the first continuum input file.
   - `input_muse_filename`: Path to the MUSE input file.

2. Run the code.

## Output

The code generates the following output files in the specified output directory:
- Continuum-subtracted H-alpha image: `{output_dir}/{galaxy}_halpha_raw.fits`
- Scaled continuum image: `{output_dir}/{galaxy}_cont_raw.fits`
- Processed H-alpha unit image: `{output_dir}/{galaxy}_halpha.fits`
- Processed H-alpha MUSE image: `{output_dir}/{galaxy}_musehalpha.fits`
- Regridded MUSE image: `{output_dir}/{galaxy}_musehalpha_regrid.fits`
- Smoothed HST image: `{halpha_filename.replace('_raw.fits', '_smoothed.fits')}`
- Ratio image: `{output_dir}/{galaxy}_hstmuseratio.fits`
- Anchored HST image: `{output_dir}/{galaxy}_anchored.fits`
- Anchored HST image with intensity negations: `{output_dir}/{galaxy}_anchored_intnegs.fits`
- Anchored HST image with added white noise: `{output_dir}/{galaxy}_anchored_wnoise.fits`
- Anchored HST image with intensity negations and added white noise: `{output_dir}/{galaxy}_anchored_intnegs_wnoise.fits`

## Code 2: QA Analysis

This code performs quality analysis (QA) on HST and MUSE data.

## Prerequisites

- Python 3
- `reduction_phangs_hst` package
- `astropy` package
- `warnings` package

## Usage

1. Update the user inputs section in the code with the desired values for the following variables:
   - `galaxy`: Define the galaxy name.
   - `halpha_inputfilename`: Path to the H-alpha input file.
   - `cont2_inputfilename`: Path to the second continuum input file.
   - `cont1_inputfilename`: Path to the first continuum input file.
   - `input_muse_filename`: Path to the MUSE input file.
   - `input_nebulae_mask_filename`: Path to the nebulae mask file.
   - `input_nebulae_catalog_filename`: Path to the nebulae catalog file.

2. Run the code.

## Output

The code generates the following output in the specified output directory:
- QA analysis plots.
