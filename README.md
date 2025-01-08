# Phangs HST continuum subtraction pipeline

## Overview

The reduction_phangs_hst repository provides tools and pipelines for reducing and analyzing Hubble Space Telescope (HST) observations as part of the PHANGS-HST project. It includes workflows for continuum subtraction, extinction correction, and emission line analysis of HST narrowband observations, with optional integration of MUSE data for calibration.

See Chandar et al. 2025 paper for more details.  

[![zenodo](image.png)](https://zenodo.org/records/14610187)

## Features
*	Continuum Subtraction
    *	Advanced subtraction of stellar continuum from narrowband HST images (e.g., F658N for H-alpha).
    *	Optional integration with MUSE data for background and [NII] contamination correction.
    *	Simplified continuum subtraction without MUSE data.
*	Extinction Correction:
    *	Applies forground corrections based on input extinction values.
*	Modular Pipeline:
    *	Customizable for vario us scientific goals, including emission line analysis and star formation studies.

## Directory Structure and Requirements

The pipeline expects a specific directory structure for input and output files. Ensure your data is organized accordingly before running the pipeline.

### Expected Input Directory Structure

For example, the current setup is the following: 

```rootdir = '/Users/abarnes/Dropbox/work/Smallprojects/galaxies/data_hstha/'```

### Example galaxy: ic5332

```
# HST data 
rootdir/ic5332/hst/
    - ic5332_uvis_f658n_exp_drc_sci.fits    # Narrowband science image
    - ic5332_uvis_f555w_exp_drc_sci.fits    # Broadband filter image (F555W)
    - ic5332_uvis_f814w_exp_drc_sci.fits    # Broadband filter image (F814W)
    - ic5332_uvis_f658n_err_drc_wht.fits    # Narrowband error image
    - ic5332_uvis_f555w_err_drc_wht.fits    # Broadband error image (F555W)
    - ic5332_uvis_f814w_err_drc_wht.fits    # Broadband error image (F814W)

# MUSE data 
# Note that these MUSE images are only needed if running for with MUSE 
# See below for different run cases... 
rootdir/ic5332/muse/
    - IC5332-0.87asec_UVIS_F658N.fits       # Synthetic F658N image from MUSE
    - IC5332-0.87asec_UVIS_F555W.fits       # Synthetic F555W image from MUSE
    - IC5332-0.87asec_UVIS_F814W.fits       # Synthetic F814W image from MUSE
    - IC5332_starmask.fits                  # Star mask
    - IC5332_nebmask.fits                   # Nebula mask
    - IC5332-0.87asec_MAPS.fits             # MUSE MAPS file

# HST Filter transmission curves
rootdir/../data_misc/hst_filters/
    - HST_ACS_WFC.F550M.dat
    - HST_ACS_WFC.F555W.dat
    - HST_ACS_WFC.F658N.dat
    - HST_ACS_WFC.F814W.dat
    - HST_WFC3_UVIS1.F555W.dat
    - HST_WFC3_UVIS1.F657N.dat
    - HST_WFC3_UVIS1.F658N.dat
    - HST_WFC3_UVIS1.F814W.dat
    - (and other relevant HST filter throughput files)
```

### Output Files

Processed outputs (e.g., continuum-subtracted images, calibrated H-alpha maps) will be saved in directories within the same rootdir under (`rootdir/ic5332/hst_contsub/`). The pipeline will create subdirectories for specific results if they do not already exist.

## Installation and Dependencies

### Dependencies

This project requires Python 3.x and several scientific libraries for data handling, analysis, and visualization. Below is the complete list of required Python packages:

Core Dependencies: 
*	Astropy: FITS file handling, WCS manipulation, modeling, and statistical tools
    *	astropy.io.fits
    *	astropy.wcs
    *	astropy.units
    *	astropy.stats
    *	astropy.convolution
    *	astropy.modeling
    *	astropy.table
*	Numpy: Numerical computations (numpy)
*	Scipy: Advanced data processing and optimization (scipy.ndimage, scipy.optimize)
*	Matplotlib: Visualization (matplotlib.pyplot, matplotlib.ticker)
*	Photutils: Background subtraction and source extraction (photutils.background)
*	Reproject: Reprojecting astronomical images (reproject)
*	Radio_beam: Beam manipulation for radio astronomy data (radio_beam.Beam)
*	Synphot: Synthetic photometry tools (synphot.SpectralElement, synphot.units)
*	Extinction: Extinction curve calculation (extinction)
*	Aplpy: Astronomical data visualization (aplpy)
*	DeepCR: Cosmic ray cleaning for astronomical data (deepCR)

Other Utilities
*	Colorcet: Colormap utilities (colorcet)
*	Glob: File search with patterns (glob)
*	Warnings: Suppressing warnings (warnings)
*	Datetime: Handling date and time (datetime)
*	OS: Operating system utilities (os)

Additional Features

These packages are used for advanced functionalities, such as Gaussian convolution, background matching, and handling of data at different resolutions:
*	Astropy.convolution: Convolution and replacing NaNs
*	Astropy.stats.mad_std: Median absolute deviation for robust statistics
*	Scipy.ndimage: Binary dilation and closing
*	Extinction: Calculating dust extinction

For specific modules (e.g., deepCR), additional installation steps may be required. Refer to the module documentation for detailed instructions.

Let me know if you need further assistance!

## Usage

###  Running the Continuum Subtraction Pipeline with MUSE

1.	Open the Jupyter Notebook:

    ```jupyter notebook hstha_pipeline/run_pipeline/run_pipeline_contsub_withmuse.ipynb```


2.	Follow the notebook’s step-by-step instructions to:
    *	Load the HST narrowband and broadband images.
    *	Apply continuum subtraction using MUSE synthetic photometry.
    *	Perform background and [NII] corrections.
    *	Generate corrected H-alpha maps.

### Running the Continuum Subtraction Pipeline without MUSE

For cases where MUSE data is unavailable, a simplified continuum subtraction pipeline is available:

1.	Open the notebook:

```jupyter notebook hstha_pipeline/run_pipeline/run_pipeline_contsub_withoutmuse.ipynb```


2.	Follow the instructions to:
    *	Load HST narrowband and broadband images.
    *	Perform continuum subtraction without applying MUSE-based corrections.

### Custom usage 

For custom use cases, you could edit or create versions of the following pipeline commands: 

#### With MUSE
``` 
    def initial_pipeline(self): 

        self.load_data_hst()
        self.load_data_muse()
        self.load_data_muse_stars()
        self.load_data_muse_neb()
        self.load_data_muse_halpha()
        self.load_data_hst_inv()
        self.get_resolution()
        self.get_bandpass_info()
        self.get_sampletable_info()
        self.get_extinction()
```
```
    def continuum_subtraction_pipeline_withMUSE(self):
        
        time_start = time.time()

        self.get_clean_paths()
        self.load_and_preprocess_errors()
        self.preprocess_hst_data()

        self.apply_covmask()
        self.convert_units()
        self.smooth_and_regrid()
        self.background_correction_from_muse()
        self.correct_extinction()
        self.subtract_continuum_all()

        self.convert_to_physical_units()
        self.NII_correction_from_muse()
        self.generate_map_plots()
        self.clean_headers()
        self.convert_units_to_arcsec2()
        # self.convert_to_cdelt()
        self.save_fit_tables()

        time_end = time.time()
        time_elapsed = time_end - time_start
        print(f"Run time: {time_elapsed/60} mins")
```

#### Without MUSE
```
    def initial_pipeline_noMUSE(self): 

        self.load_data_hst()
        self.load_data_hst_inv()
        self.get_bandpass_info()
        self.get_sampletable_info()
        self.get_extinction()
```
```
    def continuum_subtraction_pipeline_noMUSE(self):

        time_start = time.time()

        self.get_clean_paths()
        self.load_and_preprocess_errors()
        self.preprocess_hst_data()
        self.apply_covmask()
        self.convert_units()
        self.correct_extinction()
        self.subtract_continuum_witherr_hst()
        self.convert_to_physical_units()
        self.generate_map_plots()
        self.clean_headers()
        self.convert_units_to_arcsec2()
        # self.convert_to_cdelt()

        time_end = time.time()
        time_elapsed = time_end - time_start
        print(f"Run time: {time_elapsed/60} mins")
```



## Contributing

Contributions are welcome! If you find a bug, have a question, or want to propose new features:
1.	Open an issue on GitHub.
2.	Fork the repository and create a pull request with your changes.

## License

MIT

## Contact

For any inquiries, please reach out via GitHub or email at ashley.barnes@eso.org.