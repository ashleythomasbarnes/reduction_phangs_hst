#!/usr/bin/env python
# coding: utf-8

# ----------------------------
###### User defined inputs 
###### Update as instructed

# Remove everything and start again? 
start_again = True

# What to run?
run_contsub = True
run_contsub_wmuse = True
run_cosmics = True
run_cosmicsnnet = True 

# Define the galaxy
galaxy = 'ngc1566'

# Define the filters
halpha_filter = 'f658n'
cont1_filter = 'f555w'
cont2_filter = 'f814w'

# Define the directories
inputdir_hst = '../hst/'
inputdir_muse = '../muse/'
outputdir = '../hst_contsub/'

###### End of user inputs
# ----------------------------

from reduction_phangs_hst import contsub_run
contsub_run.run_pipeline(start_again=start_again,
                            run_contsub = run_contsub,
                            run_contsub_wmuse = run_contsub_wmuse,
                            run_cosmics = run_cosmics,
                            run_cosmicsnnet = run_cosmicsnnet,
                            galaxy = galaxy,
                            halpha_filter = halpha_filter,
                            cont1_filter = cont1_filter,
                            cont2_filter = cont2_filter,
                            inputdir_hst = inputdir_hst,
                            inputdir_muse = inputdir_muse,
                            outputdir = outputdir)