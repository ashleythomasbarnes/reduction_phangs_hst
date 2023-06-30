import os
import aplpy
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import concurrent.futures
from tqdm import tqdm
from astropy.io import fits
from reproject import reproject_interp


def process_muse_catalouge(input_nebulae_catalog_filename, galaxy):
    """
    Processes the MUSE nebulae catalog for a specific galaxy.

    Args:
        input_nebulae_catalog_filename (str): Path to the input nebulae catalog FITS file.
        galaxy (str): Name of the galaxy.

    Returns:
        astropy.table.Table: Table containing the processed nebulae catalog.
    """
    print(f"Reading nebulae catalog: {input_nebulae_catalog_filename}")
    table_nebcat = Table.read(input_nebulae_catalog_filename)  # Read the nebulae catalog from the specified file
    table_nebcat = table_nebcat[table_nebcat['gal_name'] == galaxy.upper()]  # Filter the catalog for the given galaxy name

    return table_nebcat


def read_fits_files(input_dir, input_nebulae_mask_filename, galaxy):
    """
    Reads all .fits files in the specified directory and returns them as a dictionary.

    Args:
        input_dir (str): Path to the directory containing .fits files.
        input_nebulae_mask_hst_filename (str): Path to the input nebulae mask FITS file.
        galaxy (str): Name of the galaxy.

    Returns:
        dict: Dictionary containing the FITS data, where the keys are the file names.
    """
    fits_files = [file for file in os.listdir(input_dir) if file.endswith('.fits')]  # Get all .fits files in the directory
    fits_dict = {}

    for fits_file in fits_files:
        fits_path = os.path.join(input_dir, fits_file)  # Construct the full path to the .fits file
        hdu = fits.open(fits_path)[0]  # Open the .fits file and access the primary HDU
        fits_key = fits_file.replace('.fits', '').replace('%s_' % galaxy, '')
        fits_dict[fits_key] = hdu  # Add the FITS data to the dictionary with the file name as the key

    hdu = fits.open(input_nebulae_mask_filename)[0] 
    fits_dict['nebmask_muse'] = hdu  # Add the FITS data to the dictionary with the file name as the key
    fits_dict['nebmask_muse'].writeto('%s/%s_nebmask_muse.fits' % (input_dir, galaxy), overwrite=True)

    input_nebulae_mask_hst_filename = '%s/%s_nebmask_hst.fits' % (input_dir, galaxy)
    if not os.path.isfile(input_nebulae_mask_hst_filename):
        print("Reprojecting nebulae mask...")
        data_nebmask_hst, _ = reproject_interp(fits_dict['nebmask_muse'], fits_dict['halpha'].header, order=0, parallel=True)
        hdu_nebmask_hst = fits.PrimaryHDU(data_nebmask_hst, fits_dict['halpha'].header)
        hdu_nebmask_hst.writeto(input_nebulae_mask_hst_filename, overwrite=True)
        fits_dict['nebmask_hst'] = hdu_nebmask_hst  # Add the FITS data to the dictionary with the key 'hdu_nebmask_hst'
        print(f"Processing complete. Saving the processed nebulae mask as {input_nebulae_mask_hst_filename}")
    else:
        print("Opening existing nebulae mask...")
        hdu_nebmask_hst = fits.open(input_nebulae_mask_hst_filename)[0]
        fits_dict['nebmask_hst'] = hdu_nebmask_hst  # Add the FITS data to the dictionary with the key 'hdu_nebmask_hst'
        
    return fits_dict


def process_nebulae_flux(hdu_nebmask_hst, hdu_nebmask_muse, hdu_ha_hst, hdu_ha_hst_anchored, hdu_ha_muse):
    """
    Calculates the flux for each nebula using concurrent processing.

    Args:
        hdu_nebmask_hst (astropy.io.fits.PrimaryHDU): HDU for the nebulae mask in HST coordinates.
        hdu_nebmask_muse (astropy.io.fits.PrimaryHDU): HDU for the nebulae mask in MUSE coordinates.
        hdu_ha_hst (astropy.io.fits.PrimaryHDU): HDU for the Halpha data in HST coordinates.
        hdu_ha_hst_anchored (astropy.io.fits.PrimaryHDU): HDU for the anchored Halpha data in HST coordinates.
        hdu_ha_muse (astropy.io.fits.PrimaryHDU): HDU for the Halpha data in MUSE coordinates.

    Returns:
        tuple: Three NumPy arrays containing the flux values for HST, anchored HST, and MUSE, respectively.
    """
    ids_neb = np.unique(hdu_nebmask_muse.data)
    ids_neb = ids_neb[1:]

    flux_hst = np.ones(len(ids_neb))
    flux_hst_anchored = np.ones(len(ids_neb))
    flux_muse = np.ones(len(ids_neb))

    def process_neb(i):
        mask_hst = hdu_nebmask_hst.data == ids_neb[i]
        mask_muse = hdu_nebmask_muse.data == ids_neb[i]

        flux_hst = np.nansum(hdu_ha_hst.data[mask_hst])
        flux_hst_anchored = np.nansum(hdu_ha_hst_anchored.data[mask_hst])
        flux_muse = np.nansum(hdu_ha_muse.data[mask_muse])

        return flux_hst, flux_hst_anchored, flux_muse

    # Initialize a ThreadPoolExecutor
    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_neb, i) for i in range(len(ids_neb))]

        # Iterate over the completed futures
        for i, future in tqdm(enumerate(concurrent.futures.as_completed(futures)), total=len(futures), desc="Processing Nebs", unit="neb"):
            flux_hst[i], flux_hst_anchored[i], flux_muse[i] = future.result()

    return flux_hst, flux_hst_anchored, flux_muse


def plot_fits_data(fits_dict, galaxy):
    """
    Plots the FITS data using APLpy.

    Args:
        fits_dict (dict): Dictionary containing the FITS data.
        galaxy (str): Name of the galaxy.

    Returns:
        None
    """
    fig = plt.figure(figsize=(12, 5))

    ax1 = aplpy.FITSFigure(fits_dict['musehalpha_regrid'], subplot=(1, 4, 1), figure=fig)
    ax2 = aplpy.FITSFigure(fits_dict['halpha'], subplot=(1, 4, 2), figure=fig)
    ax3 = aplpy.FITSFigure(fits_dict['halpha_smoothed'], subplot=(1, 4, 3), figure=fig)
    ax4 = aplpy.FITSFigure(fits_dict['halpha_hstmuseratio'], subplot=(1, 4, 4), figure=fig)

    vmin, vmax = np.nanpercentile(fits_dict['musehalpha_regrid'].data, (2, 99))
    ax1.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    ax2.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    ax3.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    ax4.show_colorscale(vmin=np.nanmin(fits_dict['halpha_hstmuseratio'].data),
                        vmax=np.nanmax(fits_dict['halpha_hstmuseratio'].data),
                        cmap='turbo')

    for ax in [ax1, ax2, ax3, ax4]:
        ax.recenter(24.1721149, 15.7806457, 0.0333138)
        ax.axis_labels.hide()
        ax.tick_labels.hide()

    ax4.add_colorbar()
    fig.tight_layout()
    fig.savefig('./qa/%s_halpha_hstmuse_fluxcomp.pdf' % galaxy, dpi=300, bbox_inches='tight')
    plt.close('all')

    return(None)


def plot_flux_comparison(flux_muse, flux_hst, flux_hst_anchored, galaxy):
    """
    Plots a flux comparison between MUSE and HST flux values in two panels.

    Args:
        flux_muse (numpy.ndarray): Array of MUSE flux values.
        flux_hst (numpy.ndarray): Array of HST flux values.
        flux_hst_anchored (numpy.ndarray): Array of anchored HST flux values.
        galaxy (str): Name of the galaxy.

    Returns:
        None
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 5), sharex=True, sharey=True)

    # Plotting the first panel
    ax1.scatter(flux_muse, flux_hst, ec='none', fc='grey', s=25)
    ax1.scatter(flux_muse, flux_hst, ec='none', fc='white', s=12)
    ax1.scatter(flux_muse, flux_hst, alpha=0.1, s=25, ec='none', fc='C0')

    ax1.plot([1e-20, 1e20], np.array([1e-20, 1e20]), c='grey', ls='--', label='y=x')
    ax1.plot([1e-20, 1e20], np.array([1e-20, 1e20]) * 3, c='grey', ls=':', label='y=3x')
    ax1.plot([1e-20, 1e20], np.array([1e-20, 1e20]) / 3, c='grey', ls=':', label='y=x/3')

    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.set_xlim(flux_muse[flux_muse > 0].min(), flux_muse[flux_muse > 0].max())
    ax1.set_ylim(flux_hst[flux_hst > 0].min(), flux_hst[flux_hst > 0].max())

    ax1.grid(alpha=0.3, linestyle=':')
    ax1.legend(loc='best')

    ax1.set_xlabel('MUSE flux, (erg/s/cm2)')
    ax1.set_ylabel('HST flux, (erg/s/cm2)')

    # Plotting the second panel
    ax2.scatter(flux_muse, flux_hst_anchored, ec='none', fc='grey', s=25)
    ax2.scatter(flux_muse, flux_hst_anchored, ec='none', fc='white', s=12)
    ax2.scatter(flux_muse, flux_hst_anchored, alpha=0.1, s=25, ec='none', fc='C0')

    ax2.plot([1e-20, 1e20], np.array([1e-20, 1e20]), c='grey', ls='--', label='y=x')
    ax2.plot([1e-20, 1e20], np.array([1e-20, 1e20]) * 3, c='grey', ls=':', label='y=3x')
    ax2.plot([1e-20, 1e20], np.array([1e-20, 1e20]) / 3, c='grey', ls=':', label='y=x/3')

    ax2.set_xscale('log')
    ax2.set_yscale('log')

    ax2.set_xlim(flux_muse[flux_muse > 0].min(), flux_muse[flux_muse > 0].max())
    ax2.set_ylim(flux_hst_anchored[flux_hst_anchored > 0].min(), flux_hst_anchored[flux_hst_anchored > 0].max())

    ax2.grid(alpha=0.3, linestyle=':')
    ax2.legend(loc='best')

    ax2.set_xlabel('MUSE flux, (erg/s/cm2)')
    ax2.set_ylabel('HST anchored flux, (erg/s/cm2)')

    fig.savefig(f'./qa/{galaxy}_flux_musehstcomp.pdf', dpi=300, bbox_inches='tight')
    plt.close('all')

    return(None)