import os
import aplpy
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import concurrent.futures
from tqdm.auto import tqdm
from astropy.io import fits
from reproject import reproject_interp
from scipy.optimize import curve_fit

plt.style.use('paper')


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
        data_nebmask_hst, _ = reproject_interp(fits_dict['nebmask_muse'], fits_dict['halpha'].header, order='nearest-neighbor')
        hdu_nebmask_hst = fits.PrimaryHDU(data_nebmask_hst, fits_dict['halpha'].header)
        # hdu_nebmask_hst = fits.PrimaryHDU(np.array(data_nebmask_hst, dtype=np.int32), fits_dict['halpha'].header)
        hdu_nebmask_hst.writeto(input_nebulae_mask_hst_filename, overwrite=True)
        fits_dict['nebmask_hst'] = hdu_nebmask_hst  # Add the FITS data to the dictionary with the key 'hdu_nebmask_hst'
        print(f"Processing complete. Saving the processed nebulae mask as {input_nebulae_mask_hst_filename}")
    else:
        print("Opening existing nebulae mask...")
        hdu_nebmask_hst = fits.open(input_nebulae_mask_hst_filename)[0]
        fits_dict['nebmask_hst'] = hdu_nebmask_hst  # Add the FITS data to the dictionary with the key 'hdu_nebmask_hst'
        
    return fits_dict


def process_nebulae_flux(hdu_nebmask, hdu, table_nebcat):
    """
    Calculates the flux for each nebula using concurrent processing.

    Args:
        hdu_nebmask_hst (astropy.io.fits.PrimaryHDU): HDU for the nebulae mask
        hdu_nebmask (astropy.io.fits.PrimaryHDU): HDU 
    Returns:
    
    """
    # ids_neb = np.unique(hdu_nebmask.data)
    # ids_neb = ids_neb[~np.isnan(ids_neb)]
    # ids_neb = ids_neb[ids_neb!=-1]
    
    ids_neb = table_nebcat['region_ID']
    
    flux = np.ones(len(ids_neb))
    
    for i in tqdm(range(len(ids_neb))):
        mask = hdu_nebmask.data == ids_neb[i]
        flux[i] = np.nansum(hdu.data[mask])

    return flux


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

    vmin, vmax = np.nanpercentile(fits_dict['halpha_smoothed'].data, (2, 99))
    ax1.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    ax2.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    ax3.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    vmin, vmax = np.nanpercentile(fits_dict['halpha_hstmuseratio'].data, (2, 99))
    ax4.show_colorscale(vmin=vmin, vmax=vmax, cmap='turbo')

    for ax in [ax1, ax2, ax3, ax4]:
        ax.recenter(24.1721149, 15.7806457, 0.0333138)
        ax.axis_labels.hide()
        ax.tick_labels.hide()

    ax1.add_colorbar(location='top', axis_label_text='MUSE Intensity')
    ax2.add_colorbar(location='top', axis_label_text='HST Intensity')
    ax3.add_colorbar(location='top', axis_label_text='HST smoothed Int.')
    ax4.add_colorbar(location='top', axis_label_text='Ratio (MUSE/HST)')

    fig.tight_layout()
    fig.savefig('./qa/%s_halpha_hstmuse_fluxcomp_matchedcolor.pdf' % galaxy, dpi=300, bbox_inches='tight')
    fig.savefig('./qa/%s_halpha_hstmuse_fluxcomp_matchedcolor.png' % galaxy, dpi=300, bbox_inches='tight')
    plt.close('all')

    fig = plt.figure(figsize=(12, 5))

    ax1 = aplpy.FITSFigure(fits_dict['musehalpha_regrid'], subplot=(1, 4, 1), figure=fig)
    ax2 = aplpy.FITSFigure(fits_dict['halpha'], subplot=(1, 4, 2), figure=fig)
    ax3 = aplpy.FITSFigure(fits_dict['halpha_smoothed'], subplot=(1, 4, 3), figure=fig)
    ax4 = aplpy.FITSFigure(fits_dict['halpha_hstmuseratio'], subplot=(1, 4, 4), figure=fig)

    vmin, vmax = np.nanpercentile(fits_dict['musehalpha_regrid'].data, (2, 99))
    ax1.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    vmin, vmax = np.nanpercentile(fits_dict['halpha_smoothed'].data, (2, 99))
    ax2.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    ax3.show_colorscale(vmin=vmin, vmax=vmax, cmap='magma', stretch='log')
    vmin, vmax = np.nanpercentile(fits_dict['halpha_hstmuseratio'].data, (2, 99))
    ax4.show_colorscale(vmin=vmin, vmax=vmax, cmap='turbo')

    for ax in [ax1, ax2, ax3, ax4]:
        ax.recenter(24.1721149, 15.7806457, 0.0333138)
        ax.axis_labels.hide()
        ax.tick_labels.hide()

    ax1.add_colorbar(location='top', axis_label_text='MUSE Intensity')
    ax2.add_colorbar(location='top', axis_label_text='HST Intensity')
    ax3.add_colorbar(location='top', axis_label_text='HST smoothed Int.')
    ax4.add_colorbar(location='top', axis_label_text='Ratio (MUSE/HST)')

    fig.tight_layout()
    fig.savefig('./qa/%s_halpha_hstmuse_fluxcomp.pdf' % galaxy, dpi=300, bbox_inches='tight')
    fig.savefig('./qa/%s_halpha_hstmuse_fluxcomp.png' % galaxy, dpi=300, bbox_inches='tight')
    plt.close('all')


def plot_flux_comparison(flux_muse, flux_hst, output_filename, showplot=True):
    """
    Plots a flux comparison between MUSE and HST flux values in two panels.

    Args:
        flux_muse (numpy.ndarray): Array of MUSE flux values.
        flux_hst (numpy.ndarray): Array of HST flux values.
        galaxy (str): Name of the galaxy.

    Returns:
        None
    """
    fig, ax1 = plt.subplots(1, 1, figsize=(5, 5))

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
    
    fig.savefig(output_filename, dpi=300, bbox_inches='tight')
    
    if not showplot: 
        plt.close('all')


def make_histogram(hdu, bins=100, plot=True, percentilecut=True):
    """
    Generate histogram data from a FITS HDU image.
    
    Parameters:
    - hdu (HDU): The input HDU.
    - bins (int): Number of bins for the histogram.
    - plot (bool): Whether to plot the histogram.
    
    Returns:
    - hist (tuple): Histogram values and bin edges.
    """
    data = hdu.data.flatten()
    if percentilecut:
        data = data[data>np.nanpercentile(data, 0.01)]
        data = data[data<np.nanpercentile(data, 99.99)]
    hist = np.histogram(data[~np.isnan(data)], bins=bins)
    
    if plot:
        plt.hist(data[~np.isnan(data)], bins=bins)
        plt.xlabel("Pixel Value")
        plt.ylabel("Frequency")
        plt.title("Histogram of data")
        plt.yscale('log')
        plt.show()
        
    return hist

def get_statistics(hist):
    """
    Get statistics (mean, median, std, etc.) from histogram data.
    
    Parameters:
    - hist (tuple): Histogram values and bin edges.
    
    Returns:
    - stats (dict): Dictionary of statistics.
    """
    values = hist[0]
    bin_edges = hist[1]
    
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    mean = np.average(bin_centers, weights=values)
    std = np.sqrt(np.average((bin_centers - mean)**2, weights=values))
    median = bin_centers[np.searchsorted(np.cumsum(values), np.sum(values)/2)]
    
    stats = {"mean": mean, "median": median, "std": std}
    
    return stats

def gaussian(x, A, mu, sigma):
    """Gaussian function."""
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def fit_gaussian_to_histogram(hist, stats):
    
    """
    Fit a Gaussian function to histogram data.
    
    Parameters:
    - hist (tuple): Histogram values and bin edges.
    - stats (dict): Statistics of the histogram data (mean, std).
    
    Returns:
    - popt (tuple): Optimized parameters of the Gaussian fit.
    """
    
    values = hist[0]
    bin_edges = hist[1]
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_centers = np.linspace(bin_centers.min(), bin_centers.max(), len(values))
    
    p0 = [np.max(values), stats['mean'], stats['std']]  # Initial guess
    popt, _ = curve_fit(gaussian, bin_centers, values, p0=p0)
    
    # Plot the histogram and its fit
    # fig = plt.figure(figsize=(15,5))
    fig, axes = plt.subplots(1,2,gridspec_kw={'width_ratios':[2,1]}, figsize=(15,5), sharey=True)
    
    for ax in axes:
        
        ax.hist(bin_centers, bins=bin_edges, weights=values, alpha=0.4, label='Data: mean=%0.2f med=%0.2f std=%0.2f' %(stats['mean'], stats['median'], stats['std']))
        ax.plot(bin_centers, values, alpha=1, c='C0', lw=1, ds='steps-mid')
        ax.plot(bin_centers, gaussian(bin_centers, *popt), 'k--', label='Fit: mean=%0.2f std=%0.2f' %(popt[1], popt[2]))
        ax.plot([popt[1],popt[1]], [0.9, values.max()*1.5], 'k:', alpha=0.4, label='Fit mean')

        ax.set_xlabel("Flux [erg/s/cm2/pixel]")
        ax.set_yscale('log')
        ax.set_ylim([0.9, values.max()*1.5])
        ax.grid()
         
    axes[0].set_ylabel("Frequency")
    axes[0].legend()
    axes[0].set_xlim([bin_centers.min(), bin_centers.max()])
    # axes[1].set_xlim([-100, 100])
    axes[1].set_xlim([popt[1]-popt[2]*5,popt[1]+popt[2]*5])
    
    plt.tight_layout()
    
    return popt, fig

def run_histogram(hdu, outputfile=None, bins=1000):
    hist = make_histogram(hdu, bins=bins, plot=False)
    stats = get_statistics(hist)
    fit_params, fig = fit_gaussian_to_histogram(hist, stats)
    
    if outputfile is not None: 
        fig.savefig(outputfile)


def get_scaling(hdu):
    hdu_scaled = hdu.copy()
    data = hdu_scaled.data
    data = (data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    data = np.log10(data)
    # data = (data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    hdu_scaled.data = data
    return(hdu_scaled)

def plot_map(fits_dict, outputfile, showplots=True, norm=False):
    """
    Plots the FITS data using APLpy.

    Args:
        fits_dict (dict): Dictionary containing the FITS data.
        galaxy (str): Name of the galaxy.

    Returns:
        None
    """
    fig = plt.figure(figsize=(15, 15))

    ax = ['']*9
    
    keys = ['musehalpha_regrid', 
            'halpha_bgsub_smoothed', 
            'halpha', 
            'halpha_bgsub', 
            'halpha_bgsub_fit_anchored', 
            'halpha_bgsub_fit_anchored_intnegs', 
            'halpha_bgsub_fit_anchored_intnegs_nocosmic', 
            'halpha_bgsub_fit_anchored_intnegs_nocosmic_nnet', 
            'halpha_bgsub_fit_anchored_intnegs_ratio']
    
    for i in tqdm(range(len(keys))):
    
        if norm:
            hdu = get_scaling(fits_dict[keys[i]])
            vmin, vmax = np.nanpercentile(hdu.data, (2, 99.9))
        else: 
            hdu = fits_dict[keys[i]]
            vmin, vmax = np.nanpercentile(fits_dict['halpha_bgsub_smoothed'].data, (2, 99))
        
        ax[i] = aplpy.FITSFigure(hdu, subplot=(3, 3, i+1), figure=fig)
        ax[i].show_colorscale(vmin=vmin, vmax=vmax, cmap='turbo', stretch='linear')
        
        ax[i].axis_labels.hide()
        ax[i].tick_labels.hide()

        ax[i].add_colorbar(location='top', axis_label_text='%s' %keys[i])

    fig.tight_layout()
    
    if norm: 
        outputfile = outputfile.replace('.pdf','_norm.pdf')
        fig.savefig(outputfile, dpi=300, bbox_inches='tight')
    else:
        fig.savefig(outputfile, dpi=300, bbox_inches='tight')
    
    if not showplots:
        plt.close('all')
        
    return()