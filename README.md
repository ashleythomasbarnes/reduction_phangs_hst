# Continuum subtraction and scaling for HST Ha data using MUSE data

## Code 1: run_pipeline_contsub.ipynb

The continuum_subtraction.py code performs continuum subtraction on HST data and post-processes the resulting images. It takes user inputs such as the galaxy name and file paths for the H-alpha, continuum, and MUSE data. The code generates several output files, including continuum-subtracted H-alpha and scaled continuum images. It also processes the H-alpha units and MUSE data, performs regridding, smoothing, and generates ratio images. Additionally, the code applies various transformations to the anchored HST image, such as intensity negations and adding white noise. The output files are saved in the specified output directory.

## Code 2: run_pipeline_contsub_qa.ipynb

The qa_analysis.py code performs quality analysis (QA) on HST and MUSE data. It takes user inputs such as the galaxy name, file paths for H-alpha, continuum, MUSE, nebulae mask, and nebulae catalog. The code reads FITS files, creates a dictionary of FITS data, and processes the MUSE catalog. It calculates the flux values for nebulae in the HST and MUSE data and generates QA analysis plots. The output of the code includes QA analysis plots that provide insights into the quality of the HST and MUSE data.

These codes are part of the reduction_phangs_hst package and utilize the astropy and warnings packages for various data processing and analysis tasks.

## Code 1: Continuum Subtraction

This code performs continuum subtraction on HST data and post-processes the resulting images.

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

