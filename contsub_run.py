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

    halpha_inputfilename = glob.glob('%s%s_*%s*.fits' %(inputdir_hst, galaxy, halpha_filter))[0]
    cont1_inputfilename = glob.glob('%s%s_*%s*.fits' %(inputdir_hst, galaxy, cont1_filter))[0]
    cont2_inputfilename = glob.glob('%s%s_*%s*.fits' %(inputdir_hst, galaxy, cont2_filter))[0]

    output_dir = '%s/%s_%s_%s/' % (outputdir, halpha_filter, cont1_filter, cont2_filter) 
    halpha_filename = '%s/%s_halpha_raw.fits' % (output_dir, galaxy)  
    cont_filename = '%s/%s_cont_raw.fits'  % (output_dir, galaxy)  

    # Check for MUSE data... 
    try: 
        input_muse_filename = glob.glob('%s*_MAPS.fits' %inputdir_muse)[0]
        input_muse_starmask_filename = glob.glob('%s*_starmask.fits' %inputdir_muse)[0]
    except:
        print('[INFO] No MUSE file found..')
        print('[INFO] Will *not* run flux anchoring to MUSE data.')
        run_contsub_wmuse = False
        run_contsub_nomuse = True
        run_cosmics = False
        run_cosmicsnnet = False 

    # Check/Make directory
    if start_again:
        contsub_misc.remove_directory(output_dir, safety=False)  
    contsub_misc.create_directory(output_dir)  

    ######
    if run_contsub: 

        # Run continuum subtraction
        contsub.continuum_subtraction(halpha_inputfilename, cont1_inputfilename, halpha_filename, cont_filename, cont2_inputfilename)
        
        output_filename = halpha_filename.replace('_raw.fits', '.fits')
        hdu_hst = contsub_postprocess.process_halpha_units(halpha_filename, output_filename)
        hdu_hst = contsub_postprocess.remove_nan_padding(hdu_hst)
        hdu_hst.writeto(output_filename, overwrite=True)

        # Get mask of mosaic map
        contsub_postprocess.get_mask(output_filename)

        # Background subtraction 
        hdu_hst_bgsub = contsub_postprocess.process_halpha_background(hdu_hst, halpha_filename)
        hdu_hst_bgsub.writeto(halpha_filename.replace('_raw.fits', '_bgsub.fits'), overwrite=True)

        # Mask stars
        hdu_muse_starmask = fits.open(input_muse_starmask_filename)[0]
        hdu_muse_starmask_regrid = contsub_postprocess.regrid(hdu_muse_starmask, hdu_hst_bgsub, conserve_flux=False, order='nearest-neighbor')
        hdu_hst_bgsub = contsub_postprocess.mask_stars(hdu_hst_bgsub, hdu_muse_starmask_regrid, halpha_filename.replace('_raw.fits', '_bgsub_starmasked.fits'))


    if run_contsub_wmuse: 

        # Load MUSE
        hdu_muse = contsub_postprocess.process_halpha_muse(input_muse_filename, '%s/%s_musehalpha.fits' % (output_dir, galaxy))
        hdu_muse_regrid = contsub_postprocess.regrid(hdu_muse, hdu_hst_bgsub, '%s/%s_musehalpha_regrid.fits' % (output_dir, galaxy))

        # Set the initial resolution for HST observations
        initial_resolution = 0.1 * u.arcsec
        desired_resolution = float(contsub_misc.extract_substring(input_muse_filename, pattern=r'\d+\.\d+')) * u.arcsec
        hdu_hst_bgsub_smoothed = contsub_postprocess.smooth_image_with_beam(hdu_hst_bgsub, initial_resolution, desired_resolution, halpha_filename.replace('_raw.fits', '_bgsub_smoothed.fits'))
        hdu_hst_bgsub_smoothed_regrid = contsub_postprocess.regrid(hdu_hst_bgsub_smoothed, hdu_muse, halpha_filename.replace('_raw.fits', '_bgsub_smoothed_regrid.fits'))

        # Anchoring via difference and ratio - these are just for pretty images...
        output = contsub_postprocess.save_diff_ratio_smoothed_image(hdu_muse_regrid, 
                                                                    hdu_hst_bgsub, 
                                                                    hdu_hst_bgsub_smoothed, 
                                                                    halpha_filename.replace('_raw.fits', '_bgsub_ratio.fits'), 
                                                                    halpha_filename.replace('_raw.fits', '_bgsub_diff.fits'), 
                                                                    halpha_filename.replace('_raw.fits', '_bgsub_ratio_anchored.fits'), 
                                                                    halpha_filename.replace('_raw.fits', '_bgsub_diff_anchored.fits'))
        _, _, hdu_hst_bgsub_ratio_anchored, _ = output
        hdu_hst_bgsub_ratio_anchored_intnegs = contsub_postprocess.process_intnegs_anchored_image(hdu_hst_bgsub_ratio_anchored, halpha_filename.replace('_raw.fits', '_bgsub_ratio_anchored_intnegs.fits'))

        # Anchoring via fitting - interpolate > 5sigma negatives
        hdu_hst_bgsub_fit_anchored = contsub_postprocess.create_2d_hist_and_fit(hdu_muse, 
                                                                                hdu_hst_bgsub_smoothed_regrid, 
                                                                                hdu_hst_bgsub, 
                                                                                hdu_muse_starmask, 
                                                                                halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored.fits'))
        hdu_hst_bgsub_fit_anchored_intnegs = contsub_postprocess.process_anchored_fit_image(hdu_hst_bgsub_fit_anchored, halpha_filename.replace('_raw.fits', '_bgsub_fit_anchored_intnegs.fits'))

    # if run_cosmics:

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