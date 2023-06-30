# run_pipeline_contsub.ipynb
## Continuum Subtraction

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

# run_pipeline_contsub_qa.ipynb
## QA Analysis

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

