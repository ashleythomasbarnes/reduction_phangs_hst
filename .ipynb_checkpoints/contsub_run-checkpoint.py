#!/usr/bin/env python
# coding: utf-8

# Imports
import os
import astropy.units as u
from astropy.io import fits
from reduction_phangs_hst import contsub, contsub_misc, contsub_postprocess, contsub_postprocess_cosmicrays, contsub_postprocess_cosmicraysnnet
import warnings
import glob
warnings.filterwarnings('ignore')

def run_pipeline(start_again=False,
                run_contsub = False,
                run_contsub_wmuse = False,
                run_contsub_nomuse = False,
                run_cosmics = False,
                run_cosmicsnnet = False,
                galaxy = 'ngc4321',
                halpha_filter = 'f657n',
                cont1_filter = 'f555w',
                cont2_filter = 'f814w',
                inputdir_hst = '../hst/',
                inputdir_muse = '../muse/',
                outputdir = '../hst_contsub/', 
                threshold=0.9, 
                dilation_iterations=0):


    ###### ----------- Following will run automatically 
    ###### ----------- Little ot no user input needed
    halpha_inputfilename = glob.glob('%s%s_*%s*.fits' %(inputdir_hst, galaxy, halpha_filter))[0]
    cont1_inputfilename = glob.glob('%s%s_*%s*.fits' %(inputdir_hst, galaxy, cont1_filter))[0]
    cont2_inputfilename = glob.glob('%s%s_*%s*.fits' %(inputdir_hst, galaxy, cont2_filter))[0]

    # Create output directory for continuum subtraction
    output_dir = '%s/%s_%s_%s/' % (outputdir, halpha_filter, cont1_filter, cont2_filter) 

    # Set the output filename for the continuum-subtracted emission line image
    # Set the output filename for the scaled continuum image
    halpha_filename = '%s/%s_halpha_raw.fits' % (output_dir, galaxy)  
    cont_filename = '%s/%s_cont_raw.fits'  % (output_dir, galaxy)  

    # Define MUSE input files for flux postprocessing
    # If not here will overwrite commands to run...
    try: 
        input_muse_filename = glob.glob('%s*_MAPS.fits' %inputdir_muse)[0]
    except:
        print('[INFO] No MUSE file found..')
        print('[INFO] Will *not* run flux anchoring to MUSE data.')
        run_contsub_wmuse = False
        run_contsub_nomuse = True
        run_cosmics = False
        run_cosmicsnnet = False 

    # Remove the output directory and create the output directory if it doesn't exist
    if start_again:
        contsub_misc.remove_directory(output_dir, safety=False)  
    contsub_misc.create_directory(output_dir)  

    if run_contsub: 

        # Run continuum subtraction
        # Perform continuum subtraction on the halpha_inputfilename using the cont1_inputfilename and cont2_inputfilename
        contsub.continuum_subtraction(halpha_inputfilename, cont1_inputfilename, halpha_filename, cont_filename, cont2_inputfilename)

        # Define the output filename by replacing '_raw.fits' with '.fits'
        output_filename = halpha_filename.replace('_raw.fits', '.fits')

        # Process the Halpha units, remove nans and save the result to the output file
        hdu_hst = contsub_postprocess.process_halpha_units(halpha_filename, output_filename)
        hdu_hst = contsub_postprocess.remove_nan_padding(hdu_hst)
        hdu_hst.writeto(output_filename, overwrite=True)

        # Get mask of mosaic map
        contsub_postprocess.get_mask(output_filename)

        # Call the 'process_halpha_background' function, passing in the hdu_hst and halpha_filename. 
        # This function processes the H-alpha image to remove the background. 
        # The result is a new HDU (Header Data Unit), with the background subtracted from the original H-alpha image.
        hdu_hst_bgsub = contsub_postprocess.process_halpha_background(hdu_hst, halpha_filename)
        hdu_hst_bgsub = contsub_postprocess.mask_stars(hdu_hst_bgsub, halpha_filename.replace('_raw.fits', '_bgsub.fits'))

    if run_contsub_wmuse: 

        # Process the Halpha MUSE file and save the result to the output file
        output_muse_filename = '%s/%s_musehalpha.fits' % (output_dir, galaxy)
        hdu_muse = contsub_postprocess.process_halpha_muse(input_muse_filename, output_muse_filename)

        # Perform regridding of the MUSE data using the HST data and save the result to the output file
        output_muse_filename = '%s/%s_musehalpha_regrid.fits' % (output_dir, galaxy)
        hdu_muse_regrid = contsub_postprocess.regrid(hdu_muse, hdu_hst_bgsub, output_muse_filename)

        # Set the initial resolution for HST observations
        initial_resolution = 0.1 * u.arcsec
        desired_resolution = float(contsub_misc.extract_substring(input_muse_filename, pattern=r'\d+\.\d+')) * u.arcsec

        # Smooth the HST image with the desired beam resolution and save the result to the output file
        output_filename = halpha_filename.replace('_raw.fits', '_bgsub_smoothed.fits')
        hdu_hst_bgsub_smoothed = contsub_postprocess.smooth_image_with_beam(hdu_hst_bgsub, initial_resolution, desired_resolution, output_filename)

        # Perform regridding of the HST data using the MUSE data and save the result to the output file
        output_filename = halpha_filename.replace('_raw.fits', '_bgsub_smoothed_regrid.fits')
        hdu_hst_bgsub_smoothed_regrid = contsub_postprocess.regrid(hdu_hst_bgsub_smoothed, hdu_muse, output_filename)

        # These are just for pretty images
        # Save the difference and ratio images, and smoothed HST image, and obtain the resulting HDUs
        # Generate the output filename for the difference and ratio image 
        # Generate the output filename for the anchored HST images
        output_ratio_filename = halpha_filename.replace('_raw.fits', '_bgsub_ratio.fits')
        output_diff_filename = halpha_filename.replace('_raw.fits', '_bgsub_diff.fits')
        output_ratio_anchored_filename = halpha_filename.replace('_raw.fits', '_bgsub_ratio_anchored.fits')
        output_diff_anchored_filename = halpha_filename.replace('_raw.fits', '_bgsub_diff_anchored.fits')
        output_ratio_anchored_intnegs_filename = halpha_filename.replace('_raw.fits', '_bgsub_ratio_anchored_intnegs.fits')
        output = contsub_postprocess.save_diff_ratio_smoothed_image(hdu_muse_regrid, hdu_hst_bgsub, hdu_hst_bgsub_smoothed, 
                                                                   output_ratio_filename, output_diff_filename, 
                                                                   output_ratio_anchored_filename, output_diff_anchored_filename)
        hdu_ratio_smooth, hdu_diff_smooth, hdu_hst_bgsub_ratio_anchored, hdu_hst_bgsub_diff_anchored = output
        hdu_hst_bgsub_ratio_anchored_intnegs = contsub_postprocess.process_intnegs_anchored_image(hdu_hst_bgsub_ratio_anchored, output_ratio_anchored_filename)

        # Create a 2D histogram and fit for the provided HDU objects, and save the modified second HDU object
        # Define the output filename
        output_fit_anchored_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored.fits')
        hdu_hst_bgsub_fit_anchored = contsub_postprocess.create_2d_hist_and_fit(hdu_muse, hdu_hst_bgsub_smoothed_regrid, hdu_hst_bgsub, output_fit_anchored_filename)

    #     # Process the anchored fit image using the defined function, and save the result
    #     # Define the filename for the processed anchored fit image
    #     hdu_hst_bgsub_fit_anchored_intnegs_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs.fits')
    #     hdu_hst_bgsub_fit_anchored_intnegs = contsub_postprocess.process_anchored_fit_image(hdu_hst_bgsub_fit_anchored, hdu_hst_bgsub_fit_anchored_intnegs_filename)

    #     # Smooth the HST image with the desired beam resolution and save the result to the output file
    #     output_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_smoothed.fits')
    #     hdu_hst_bgsub_fit_anchored_intnegs_smoothed = contsub_postprocess.smooth_image_with_beam(hdu_hst_bgsub_fit_anchored_intnegs, initial_resolution, desired_resolution, output_filename)

    #     # These are for the cosmic ray subtraction
    #     # Save the difference and ratio images, and smoothed HST image, and obtain the resulting HDUs
    #     # Generate the output filename for the difference and ratio image 
    #     # Generate the output filename for the anchored HST images
    #     output_ratio_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_ratio.fits')
    #     output_diff_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_diff.fits')
    #     output_ratio_anchored_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_ratio_anchored.fits')
    #     output_diff_anchored_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_diff_anchored.fits')
    #     output = contsub_postprocess.save_diff_ratio_smoothed_image(hdu_muse_regrid, 
    #                                                                 hdu_hst_bgsub_fit_anchored_intnegs, 
    #                                                                 hdu_hst_bgsub_fit_anchored_intnegs_smoothed, 
    #                                                                 output_ratio_filename, 
    #                                                                 output_diff_filename, 
    #                                                                 output_ratio_anchored_filename,
    #                                                                 output_diff_anchored_filename)
    #     _, _, _, hdu_hst_bgsub_diff_anchored = output

    # if run_cosmics:

    #     # Use the cosmicray_finder function from the contsub_postprocess_cosmicrays module to locate 
    #     # the positions of cosmic rays in the hdu_hst_bgsub_diff_anchored image
    #     cosmicray_positions = contsub_postprocess_cosmicrays.cosmicray_finder(hdu_hst_bgsub_diff_anchored)

    #     # Subtract the cosmic rays from the hdu_hst_bgsub_fit_anchored_intnegs image using the subtract_intp_cut_data 
    #     # function from the contsub_postprocess_cosmicrays module
    #     hdu_hst_bgsub_fit_anchored_nocosmic = contsub_postprocess_cosmicrays.subtract_intp_cut_data(hdu_hst_bgsub_fit_anchored_intnegs, cosmicray_positions)

    #     # Define the filename for the processed anchored fit image
    #     hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_nocosmic.fits')
    #     contsub_postprocess_cosmicrays.save_masked(hdu_hst_bgsub_fit_anchored_nocosmic, hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_filename)

    #     # The following commented out code is for processing an image hdu_hst_bgsub_ratio_anchored_intnegs in the same manner
    #     # It first subtracts the cosmic rays, then defines a filename for the processed image, 
    #     # and finally writes the processed image to the defined filename
    #     hdu_hst_bgsub_ratio_anchored_nocosmic = contsub_postprocess_cosmicrays.subtract_intp_cut_data(hdu_hst_bgsub_ratio_anchored_intnegs, cosmicray_positions)
    #     hdu_hst_bgsub_ratio_anchored_intnegs_nocosmic_filename = halpha_filename.replace('_raw.fits', '_ratio_anchored_intnegs_nocosmic.fits')
    #     contsub_postprocess_cosmicrays.save_masked(hdu_hst_bgsub_ratio_anchored_nocosmic, hdu_hst_bgsub_ratio_anchored_intnegs_nocosmic_filename)

    # if run_cosmicsnnet:

    #     hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_nocosmic.fits')
    #     hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_nnet_filename = halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs_nocosmic_nnet.fits')

    #     contsub_postprocess_cosmicraysnnet.cosmicray_finder_nnet(hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_filename, hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_nnet_filename, 
    #                                                          threshold=threshold, dilation_iterations=dilation_iterations)

    # if run_contsub_nomuse:

    #     input_filename = halpha_filename.replace('_raw.fits', '_bgsub.fits')
    #     output_filename = halpha_filename.replace('_raw.fits', '_bgsub_intnegs.fits')
    #     hdu_hst_bgsub_ratio_anchored_intnegs = contsub_postprocess.process_anchored_fit_image(fits.open(input_filename)[0], output_filename)

    #     hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_filename = halpha_filename.replace('_raw.fits', '_bgsub_intnegs.fits')
    #     hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_nnet_filename = halpha_filename.replace('_raw.fits', '_bgsub_intnegs_nnet.fits')

    #     contsub_postprocess_cosmicraysnnet.cosmicray_finder_nnet(hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_filename, hdu_hst_bgsub_fit_anchored_intnegs_nocosmic_nnet_filename)