from tools_imports import *

class PyHSTHAContSub:
    def __init__(self, galaxy, rootdir, filters, instruments, 
                 galaxy_muse=None, rootdir_bp=None, sampletable_path=None, 
                 cr_threshold=0.9, cr_dilation_iterations=0):

        self.galaxy = galaxy
        self.rootdir = rootdir

        if galaxy not in rootdir: 
            self.rootdir = os.path.join(rootdir, galaxy)
            self.rootdir = self.rootdir+'/'

        self.filter_narrow = filters[0]
        self.filter_broad1 = filters[1]
        self.filter_broad2 = filters[2]

        self.instrument_narrow = instruments[0]
        self.instruments_broad1 = instruments[1]
        self.instruments_broad2 = instruments[2]

        if galaxy_muse is None:
            self.galaxy_muse = galaxy
        else: 
            self.galaxy_muse = galaxy_muse
        
        if rootdir_bp is None:
            self.rootdir_bp = '/Users/abarnes/Dropbox/work/Smallprojects/galaxies/data_misc/hst_filters/'
        else:
            self.rootdir_bp = rootdir_bp

        if sampletable_path is None:
            self.sampletable_path = '/Users/abarnes/Dropbox/work/Smallprojects/galaxies/data_misc/sample_table/phangs_sample_table_v1p6.fits'
        else: 
            self.sampletable_path = sampletable_path

        self.hdu_data = {}
        self.resolutions = {}
        self.bandpass = {}
        self.extcorr_factor = {}
        self.fit_bgcorr = {}
        self.fit_NII = {}

        self.cr_threshold = cr_threshold
        self.cr_dilation_iterations = cr_dilation_iterations


    def load_data_hst(self):
        # Load HDU files and store them in self.hdu_data
        self.hdu_data['hst_narrow'] = get_hdu(self.rootdir, f"hst/{self.galaxy.lower()}*_*{self.instrument_narrow.lower()}_*{self.filter_narrow.lower()}*_exp_drc_sci.fits")
        self.hdu_data['hst_broad1'] = get_hdu(self.rootdir, f"hst/{self.galaxy.lower()}*_*{self.instruments_broad1.lower()}_*{self.filter_broad1.lower()}*_exp_drc_sci.fits")
        self.hdu_data['hst_broad2'] = get_hdu(self.rootdir, f"hst/{self.galaxy.lower()}*_*{self.instruments_broad2.lower()}_*{self.filter_broad2.lower()}*_exp_drc_sci.fits")

    def load_data_hst_inv(self):
        # Load HDU files and store them in self.hdu_data
        self.hdu_data['hst_narrow_inv'] = get_hdu(self.rootdir, f"hst/{self.galaxy.lower()}*_*{self.instrument_narrow.lower()}_*{self.filter_narrow.lower()}*_err_drc_wht.fits")
        self.hdu_data['hst_broad1_inv'] = get_hdu(self.rootdir, f"hst/{self.galaxy.lower()}*_*{self.instruments_broad1.lower()}_*{self.filter_broad1.lower()}*_err_drc_wht.fits")
        self.hdu_data['hst_broad2_inv'] = get_hdu(self.rootdir, f"hst/{self.galaxy.lower()}*_*{self.instruments_broad2.lower()}_*{self.filter_broad2.lower()}*_err_drc_wht.fits")
        
    def load_data_muse(self):
        # Load HDU files and store them in self.hdu_data
        self.hdu_data['muse_narrow'] = get_hdu(self.rootdir, f"muse/{self.galaxy_muse.upper()}*_{self.instrument_narrow.upper()}_*{self.filter_narrow.upper()}.fits")
        self.hdu_data['muse_broad1'] = get_hdu(self.rootdir, f"muse/{self.galaxy_muse.upper()}*_{self.instruments_broad1.upper()}_*{self.filter_broad1.upper()}.fits")
        self.hdu_data['muse_broad2'] = get_hdu(self.rootdir, f"muse/{self.galaxy_muse.upper()}*_{self.instruments_broad2.upper()}_*{self.filter_broad2.upper()}.fits")
  
    def load_data_muse_halpha(self):
        self.hdu_data['muse_DAP'] = get_hdu(self.rootdir, f'muse/{self.galaxy_muse.upper()}*_MAPS.fits', 'all')
        self.hdu_data['muse_halpha'] = self.hdu_data['muse_DAP']['HA6562_FLUX']
        self.hdu_data['muse_halpha'] = fits.PrimaryHDU(self.hdu_data['muse_halpha'].data, self.hdu_data['muse_halpha'].header)
        self.hdu_data['muse_halpha'].header.set('EQUINOX',  2000.0, '[yr] Equinox of equatorial coordinates') 

    def load_data_muse_stars(self):
        self.hdu_data['muse_stars'] = get_hdu(self.rootdir, f"muse/{self.galaxy_muse.upper()}_starmask.fits")

    def load_data_muse_neb(self):
        self.hdu_data['muse_neb'] = get_hdu(self.rootdir, f"muse/{self.galaxy_muse.upper()}_nebmask.fits")

    def get_resolution(self):
        # Load and calculate resolutions
        self.resolutions['hst'] = 0.07 * u.arcsec
        _, file_muse_f65Xn = get_hdu(self.rootdir, f"muse/{self.galaxy_muse.upper()}*_{self.instrument_narrow.upper()}_*{self.filter_narrow.upper()}.fits", return_filename=True, print_filename=False)
        self.resolutions['muse'] = float(file_muse_f65Xn.split('asec')[0].split('-')[-1]) * u.arcsec

    def get_bandpass_info(self):
        # Load bandpass information
        self.bandpass_all = get_bandpassinfo(self.rootdir_bp)

        self.bandpass['narrow'] = self.bandpass_all[f'{self.instrument_narrow.upper()}_{self.filter_narrow.upper()}']
        self.bandpass['broad1'] = self.bandpass_all[f'{self.instruments_broad1.upper()}_{self.filter_broad1.upper()}']
        self.bandpass['broad2'] = self.bandpass_all[f'{self.instruments_broad2.upper()}_{self.filter_broad2.upper()}']

    def get_sampletable_info(self):
        # Load sample table information
        self.sampletable = QTable.read(self.sampletable_path, format='fits')
        self.sampletable_galaxy = self.sampletable[self.sampletable['name'] == self.galaxy_muse.lower()]

    def get_extinction(self):

        # E(B-V) from the sample table
        self.ebv = self.sampletable_galaxy['mwext_sf11'].value[0] 
        self.R = 3.1 # standard for Galaxy
        self.Av = self.ebv * self.R

        self.extcorr_factor['narrow'] = extinction.remove(extinction.ccm89(np.array([self.bandpass['narrow']['pivot']]), self.Av, self.R), 1)
        self.extcorr_factor['broad1'] = extinction.remove(extinction.ccm89(np.array([self.bandpass['broad1']['pivot']]), self.Av, self.R), 1)
        self.extcorr_factor['broad2'] = extinction.remove(extinction.ccm89(np.array([self.bandpass['broad2']['pivot']]), self.Av, self.R), 1)

    def get_clean_paths(self):
        # Make paths
        clean_paths(self.rootdir)

  # Step 1: Load and preprocess errors
    def load_and_preprocess_errors(self):

        self.hdu_data['hst_narrow_err'] = conv_inverse_variance_to_error(self.hdu_data['hst_narrow_inv'])
        self.hdu_data['hst_broad1_err'] = conv_inverse_variance_to_error(self.hdu_data['hst_broad1_inv'])
        self.hdu_data['hst_broad2_err'] = conv_inverse_variance_to_error(self.hdu_data['hst_broad2_inv'])

   # Step 2: Convert units and handle NaNs
    def preprocess_hst_data(self):

        for key in ['hst_narrow', 'hst_broad1', 'hst_broad2']:
            self.hdu_data[key] = get_nanzeros(self.hdu_data[key])
        self.hdu_data['hst_narrow'] = remove_nan_padding(self.hdu_data['hst_narrow'])

        for key in ['hst_broad1', 'hst_broad2']:
            self.hdu_data[key] = get_regrid(self.hdu_data[key], self.hdu_data['hst_narrow'])
        
        if 'hst_narrow_err' in self.hdu_data:
            for key in ['hst_narrow_err', 'hst_broad1_err', 'hst_broad2_err']:
                self.hdu_data[key] = get_regrid(self.hdu_data[key], self.hdu_data['hst_narrow'])

    # Step 3: Apply masks
    def apply_covmask(self):

        self.hdu_data['hst_narrow'], self.hdu_data['hst_broad1'], self.hdu_data['hst_broad2'] = get_covmask(
            self.hdu_data['hst_narrow'], self.hdu_data['hst_broad1'], self.hdu_data['hst_broad2'])

    # Step 4: Convert units
    def convert_units(self):

        for key in ['hst_narrow', 'hst_broad1', 'hst_broad2']:
            self.hdu_data[key] = get_electrons_2_ergcm2sA(self.hdu_data[key])
        
        if 'hst_narrow_err' in self.hdu_data:
            for key in ['hst_narrow_err', 'hst_broad1_err', 'hst_broad2_err']:
                photflam = self.hdu_data[key.replace('_err', '')].header['PHOTFLAM']
                self.hdu_data[key] = get_electrons_2_ergcm2sA(self.hdu_data[key], photflam)

        if 'muse_narrow' in self.hdu_data:
            for muse_key, band_key in [('muse_narrow', 'narrow'), ('muse_broad1', 'broad1'), ('muse_broad2', 'broad2')]:
                pivot = self.bandpass[band_key]['pivot']
                self.hdu_data[muse_key] = get_Jy_2_ergcm2sA(self.hdu_data[muse_key], pivot)

   # Step 5: Smooth and regrid to MUSE resolution
    def smooth_and_regrid(self):

        for key in ['hst_narrow', 'hst_broad1', 'hst_broad2']:
            sm_key = f'{key}_sm'
            smre_key = f'{sm_key}re'
            self.hdu_data[sm_key] = get_smooth(self.hdu_data[key], self.resolutions['hst'], self.resolutions['muse'])
            self.hdu_data[smre_key] = get_regrid(self.hdu_data[sm_key], self.hdu_data[key.replace('hst', 'muse')])

    # Step 6: Background correct flux to MUSE
    def background_correction_from_muse(self):

        for key, band_key in [('hst_narrow', 'narrow'), ('hst_broad1', 'broad1'), ('hst_broad2', 'broad2')]:
            
            muse_key = f'muse_{band_key}'
            smre_key = f'{key}_smre'
            an_key = f'{key}_bgcorr'
            smrean_key = f'{key}_smrean'

            output = get_anchoring_offset(self.hdu_data[muse_key], self.hdu_data[smre_key], self.hdu_data[key], self.hdu_data['muse_stars'], band_key, self.rootdir)
            self.hdu_data[an_key], self.hdu_data[smrean_key], self.fit_bgcorr[band_key] = output

    def correct_extinction(self):
        """
        Apply extinction correction to HDU data for narrowband, broadband, and error HDUs.
        """
        for key, band in zip(
            ['hst_narrow', 'hst_broad1', 'hst_broad2', 
            'hst_narrow_bgcorr', 'hst_broad1_bgcorr', 'hst_broad2_bgcorr',
            'muse_narrow', 'muse_broad1', 'muse_broad2', 
            'hst_narrow_err', 'hst_broad1_err', 'hst_broad2_err'], 
            ['narrow', 'broad1', 'broad2'] * 4  # Repeat for each type
        ):
            if key in self.hdu_data:
                self.hdu_data[key].data *= self.extcorr_factor[band]

    # Step 8: Continuum subtraction
    def subtract_continuum_noerr_muse(self):
        
        output_muse = get_contsub(
            self.hdu_data['muse_narrow'], self.hdu_data['muse_broad1'], self.hdu_data['muse_broad2'],
            self.bandpass['narrow']['pivot'], self.bandpass['broad1']['pivot'], self.bandpass['broad2']['pivot'])
        (self.hdu_data['muse_contsub'], self.hdu_data['muse_cont']) = output_muse

    def subtract_continuum_noerr_hst(self):

        output_hst = get_contsub(
            self.hdu_data['hst_narrow'], self.hdu_data['hst_broad1'], self.hdu_data['hst_broad2'],
            self.bandpass['narrow']['pivot'], self.bandpass['broad1']['pivot'], self.bandpass['broad2']['pivot'])
        (self.hdu_data['hst_contsub'], self.hdu_data['hst_cont']) = output_hst

    def subtract_continuum_witherr_hst(self):

        output_hst = get_contsub_err(
            self.hdu_data['hst_narrow'], self.hdu_data['hst_broad1'], self.hdu_data['hst_broad2'],
            self.hdu_data['hst_narrow_err'], self.hdu_data['hst_broad1_err'], self.hdu_data['hst_broad2_err'],
            self.bandpass['narrow']['pivot'], self.bandpass['broad1']['pivot'], self.bandpass['broad2']['pivot'])
        (self.hdu_data['hst_contsub'], self.hdu_data['hst_cont']), (self.hdu_data['hst_contsub_err'], self.hdu_data['hst_cont_err']) = output_hst

    def subtract_continuum_witherr_hst_bgcorr(self):

        output_hst = get_contsub_err(
            self.hdu_data['hst_narrow_bgcorr'], self.hdu_data['hst_broad1_bgcorr'], self.hdu_data['hst_broad2_bgcorr'],
            self.hdu_data['hst_narrow_err'], self.hdu_data['hst_broad1_err'], self.hdu_data['hst_broad2_err'],
            self.bandpass['narrow']['pivot'], self.bandpass['broad1']['pivot'], self.bandpass['broad2']['pivot'])
        (self.hdu_data['hst_contsub_bgcorr'], self.hdu_data['hst_cont_bgcorr']), (self.hdu_data['hst_contsub_bgcorr_err'], self.hdu_data['hst_cont_bgcorr_err']) = output_hst

    def convert_to_physical_units(self):

        # Get the photometric bandwidth
        photbw = self.bandpass['narrow']['rectwidth']

        # Define the keys to process
        keys_in = ['hst_contsub', 'hst_contsub_err', 'hst_cont', 'hst_cont_err',
            'hst_contsub_bgcorr', 'hst_contsub_bgcorr_err', 'hst_cont_bgcorr', 'hst_cont_bgcorr_err',
            'muse_contsub', 'muse_cont', 'hst_narrow_bgcorr']
        keys_out = ['hst_contsub', 'hst_contsub_err', 'hst_cont', 'hst_cont_err',
            'hst_contsub_bgcorr', 'hst_contsub_bgcorr_err', 'hst_cont_bgcorr', 'hst_cont_bgcorr_err',
            'muse_contsub', 'muse_cont', 'hst_narrow_bgcorr_f']

        # Perform the conversion
        for key_in, key_out in zip(keys_in, keys_out):
            if key_in in self.hdu_data:
                self.hdu_data[key_out] = get_ergcm2sA_2_ergcm2s(self.hdu_data[key_in], photbw)

    def correct_halpha_errors(self):
        """
        Correct errors for anchored H-alpha flux using scaling factors.
        """
        # Extract scaling factor from fit_halpha_an
        # Correct errors for H-alpha and anchored H-alpha
        slope_halpha = float(self.fit_NII['fit_halpha']['slope_bins'].value[0])
        self.hdu_data['hst_contsub_halpha_err'] = self.hdu_data['hst_contsub_bgcorr_err'].copy() 
        self.hdu_data['hst_contsub_halpha_err'].data = self.hdu_data['hst_contsub_halpha_err'].data / slope_halpha

    def NII_correction_from_muse(self, do_errors=True):
        """
        Anchor H-alpha flux to MUSE data and compute scaling factors.
        """
        # Anchor H-alpha flux
        output = get_anchoring_slope(self.hdu_data['muse_halpha'], self.hdu_data['muse_contsub'], self.hdu_data['hst_contsub_bgcorr'], self.hdu_data['muse_stars'], 'halpha', self.rootdir)
        self.hdu_data['hst_contsub_halpha'], self.hdu_data['muse_contsub_halpha'], self.fit_NII['fit_halpha'] = output

        if do_errors:
            self.correct_halpha_errors()    

    def generate_map_plots(self):
        """
        Generate and save map plots for the scaled H-alpha data.
        """

        if 'hst_contsub_halpha' in self.hdu_data:
            make_plots_map(self.hdu_data['hst_contsub_halpha'], self.galaxy, 'hst_halpha', self.rootdir)
        else: 
            make_plots_map(self.hdu_data['hst_contsub'], self.galaxy, 'hst_contsub', self.rootdir)


    def clean_headers(self):
        """
        Clean the headers of the scaled H-alpha HDUs.
        """
        for key in ['hst_contsub_halpha', 'hst_contsub_halpha_err', 
                    'hst_contsub', 'hst_contsub_err', 
                    'muse_halpha', 'muse_contsub_halpha']:
            
            if key in self.hdu_data:
                self.hdu_data[key] = clean_header(self.hdu_data[key])


    def convert_units_to_arcsec2(self):
        """
        Update units of H-alpha HDUs from per pixel to per arcsecond².
        """
        # Convert scaled H-alpha and error HDUs to arcseconds²
        for key in ['hst_contsub_halpha', 'hst_contsub_halpha_err', 
                    'hst_contsub', 'hst_contsub_err', 
                    'muse_halpha', 'muse_contsub_halpha']:
            if key in self.hdu_data:
                key_as = f'{key}_as'
                self.hdu_data[key_as] = convert_perpix_to_perarcsec(self.hdu_data[key])


    def convert_to_cdelt(self):
        """
        Convert PCi_j to CDELTi 
        """
        for key_ in ['hst_contsub_halpha', 'hst_contsub_halpha_err', 
                    'hst_contsub', 'hst_contsub_err', 
                    'muse_halpha', 'muse_contsub_halpha']:

            for key in [key_, f'{key_}_as']:

                header = self.hdu_data[key].header
                PC1_1 = header['PC1_1']
                PC2_2 = header['PC2_2']
                CDELT1 = header['CDELT1']
                CDELT2 = header['CDELT2']

                if 'PC2_1' not in header:
                    header['PC2_1'] = 0
                    header['PC1_2'] = 0

                PC2_1 = header['PC2_1']
                PC1_2 = header['PC1_2']

                if CDELT1 == 1:
                    CDELT1 = PC1_1
                    CDELT2 = PC2_2
                    
                    PC1_2 = PC1_2/PC1_1
                    PC2_1 = PC2_1/PC2_2
                    PC1_1 = 1
                    PC2_2 = 1

                header['PC1_1'] = PC1_1
                header['PC1_2'] = PC1_2
                header['PC2_1'] = PC2_1
                header['PC2_2'] = PC2_2
                header['CDELT1'] = CDELT1
                header['CDELT2'] = CDELT2

                self.hdu_data[key].header = header

    def save_hdu_files(self, compress=False, save_all=False, save_nomuse=False, save_ha_only=False):
        """
        Save HDU files to disk.

        Parameters:
        compress (bool): Whether to compress the FITS files when saving.
        all (bool): Whether to save all files or just the key files.
        """

        if not (save_ha_only | save_all | save_nomuse):
            print('Not saving anything...')
            return

        elif save_all: 
            # Define all files to save if `all` is True
            save_map = {
                f'{self.galaxy}_hst_ha.fits': 'hst_contsub_halpha',
                f'{self.galaxy}_hst_ha_err.fits': 'hst_contsub_halpha_err',
                f'{self.galaxy}_hst_ha_as2.fits': 'hst_contsub_halpha_as',
                f'{self.galaxy}_hst_ha_err_as2.fits': 'hst_contsub_halpha_err_as',
                f'{self.galaxy}_muse_ha_as2.fits': 'muse_halpha_as',
                f'{self.galaxy}_muse_contsub_ha_as2.fits': 'muse_contsub_halpha_as',
                f'{self.galaxy}_hst_{self.filter_broad1}_smre.fits': 'hst_broad1_smre',
                f'{self.galaxy}_hst_{self.filter_narrow}_smre.fits': 'hst_narrow_smre',
                f'{self.galaxy}_hst_{self.filter_broad2}_smre.fits': 'hst_broad2_smre',
                f'{self.galaxy}_hst_{self.filter_broad1}.fits': 'hst_broad1',
                f'{self.galaxy}_hst_{self.filter_narrow}.fits': 'hst_narrow',
                f'{self.galaxy}_hst_{self.filter_broad2}.fits': 'hst_broad2',
                f'{self.galaxy}_muse_{self.filter_broad1}.fits': 'muse_broad1',
                f'{self.galaxy}_muse_{self.filter_narrow}.fits': 'muse_narrow',
                f'{self.galaxy}_muse_{self.filter_broad2}.fits': 'muse_broad2',
                f'{self.galaxy}_hst_{self.filter_broad1}_bgcorr.fits': 'hst_broad1_bgcorr',
                f'{self.galaxy}_hst_{self.filter_narrow}_bgcorr.fits': 'hst_narrow_bgcorr',
                f'{self.galaxy}_hst_{self.filter_narrow}_bgcorr_f.fits': 'hst_narrow_bgcorr_f',
                f'{self.galaxy}_hst_{self.filter_broad2}_bgcorr.fits': 'hst_broad2_bgcorr',
                f'{self.galaxy}_muse_{self.filter_narrow}_contsub.fits': 'muse_contsub',
                f'{self.galaxy}_hst_{self.filter_narrow}_bgcorr_contsub.fits': 'hst_contsub_bgcorr',
                f'{self.galaxy}_hst_{self.filter_narrow}_contsub.fits': 'hst_contsub',
                f'{self.galaxy}_muse_{self.filter_narrow}_cont.fits': 'muse_cont',
                f'{self.galaxy}_hst_{self.filter_narrow}_bgcorr_cont.fits': 'hst_cont_bgcorr',
                f'{self.galaxy}_hst_{self.filter_narrow}_cont.fits': 'hst_cont',
                f'{self.galaxy}_muse_ha.fits': 'muse_halpha',
                f'{self.galaxy}_muse_contsub_ha.fits': 'muse_contsub_halpha',
            }

        elif save_ha_only: 
            # Define files to save if `all` is False (default)
            save_map = {
                f'{self.galaxy}_hst_ha.fits': 'hst_contsub_halpha',
                f'{self.galaxy}_hst_ha_err.fits': 'hst_contsub_halpha_err',
                f'{self.galaxy}_hst_ha_as2.fits': 'hst_contsub_halpha_as',
                f'{self.galaxy}_hst_ha_err_as2.fits': 'hst_contsub_halpha_err_as',
                f'{self.galaxy}_muse_ha_as2.fits': 'muse_halpha_as',
                f'{self.galaxy}_muse_contsub_ha_as2.fits': 'muse_contsub_halpha_as',
            }

        elif save_nomuse:
            save_map = {
                f'{self.galaxy}_hst_contsub.fits': 'hst_contsub',
                f'{self.galaxy}_hst_contsub_err.fits': 'hst_contsub_err',
            }

        # Iterate and save files
        for filename, hdu_key in save_map.items():
            if hdu_key in self.hdu_data:
                # Ensure `filepath` is passed as the filename argument
                write_hdu(self.hdu_data[hdu_key], self.rootdir, filename, compress=compress)

    def save_fit_tables(self):
        """
        Save FITS tables for background correction offsets and scaling slopes.
        """
        # Save background correction offset fits
        save_fittables_offsets(
            self.fit_bgcorr['broad2'], 
            self.fit_bgcorr['narrow'], 
            self.fit_bgcorr['broad1'], 
            self.rootdir
        )
        
        # Save slope fits for H-alpha corrections
        save_fittables_slope(
            self.fit_NII['fit_halpha'],
            self.rootdir
        )

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

    def initial_pipeline_noMUSE(self): 

        self.load_data_hst()
        self.load_data_hst_inv()
        self.get_bandpass_info()
        self.get_sampletable_info()
        self.get_extinction()

    def subtract_continuum_all(self):

        self.subtract_continuum_noerr_muse()
        self.subtract_continuum_noerr_hst()
        self.subtract_continuum_witherr_hst()
        self.subtract_continuum_witherr_hst_bgcorr()

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