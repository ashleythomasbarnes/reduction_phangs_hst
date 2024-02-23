from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import warnings 
import colorcet
import matplotlib as mpl
import aplpy
from scipy.optimize import curve_fit

warnings.filterwarnings('ignore')

from tools_contsub_misc import * 
from tools_contsub_anchoring import * 
from tools_contsub_smoothregrid import * 

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.weight"] = "bold"
plt.rcParams["axes.labelweight"] = "bold"
plt.rcParams["xtick.direction"] = "in"
plt.rcParams["ytick.direction"] = "in" 
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True
plt.rcParams['xtick.minor.visible'] = True
plt.rcParams['ytick.minor.visible'] = True


def make_plots_diff(hdu_hst, hdu_muse, hdu_muse_stars, fit, filter='',  rootdir='./', appdir='hst_contsub/', make_plots=True, save_hdu=False):

    fit_slope = np.float32(fit['slope_bins'][0])
    # fit_slope = 1

    ratio = fits.PrimaryHDU(hdu_hst.data.copy()/hdu_muse.data.copy(), hdu_muse.header)
    diff = fits.PrimaryHDU((hdu_hst.data.copy()/fit_slope)-hdu_muse.data.copy(), hdu_muse.header)
    diff_ratio = fits.PrimaryHDU(((hdu_hst.data.copy()/fit_slope)-hdu_muse.data.copy())/hdu_muse.data.copy(), hdu_muse.header)

    ratio.data[~np.isfinite(ratio.data)] = np.nan 
    diff.data[~np.isfinite(diff.data)] = np.nan 
    diff_ratio.data[~np.isfinite(diff.data)] = np.nan 

    mask_stars = hdu_muse_stars.data!=0
    ratio.data[mask_stars] = np.nan
    diff.data[mask_stars] = np.nan
    diff_ratio.data[mask_stars] = np.nan
    
    if save_hdu: 
        ratio.writeto('ratio_%s.fits' %filter, overwrite=True)

    if make_plots: 

        fig = plt.figure(figsize=(15, 10))
        ax1 = fig.add_subplot(2, 3, 1)
        ax2 = fig.add_subplot(2, 3, 2)
        ax3 = fig.add_subplot(2, 3, 3)
        ax4 = fig.add_subplot(2, 3, 4)
        ax5 = fig.add_subplot(2, 3, 5)
        ax6 = fig.add_subplot(2, 3, 6)

        # Maps
        cmap = plt.cm.get_cmap('turbo', 11)
        img1 = ax1.imshow(ratio.data, vmin=0.5, vmax=1.5, cmap=cmap, origin='lower')
        pmin, pmax = np.nanpercentile(diff.data, [1,99])
        img2 = ax2.imshow(diff.data, vmin=pmin, vmax=pmax, cmap=cmap, origin='lower')
        pmin, pmax = np.nanpercentile(diff_ratio.data, [1,99])
        img3 = ax3.imshow(diff_ratio.data, vmin=pmin, vmax=pmax, cmap=cmap, origin='lower')

        cax1 = plt.colorbar(img1, ax=ax1, pad=0.01)
        cax2 = plt.colorbar(img2, ax=ax2, pad=0.01)
        cax3 = plt.colorbar(img3, ax=ax3, pad=0.01)

        cax1.set_label('%s HST/MUSE' %filter)
        cax2.set_label('%s HST*f - MUSE' %filter)
        cax3.set_label('%s (HST*f - MUSE)/MUSE' %filter)

        for ax in [ax1, ax2, ax3]:
            ax.set_yticks([])
            ax.set_xticks([])

        # Scatter
        ydata1 = ratio.data.flatten()
        ydata2 = diff.data.flatten()
        ydata3 = diff_ratio.data.flatten()

        for ax, ydata in zip([ax4,ax5,ax6], [ydata1, ydata2, ydata3]):

            xdata = hdu_muse.data.flatten().copy()
            ydata_ = ydata.copy()

            xmask = np.isfinite(xdata)
            xdata = xdata[xmask]
            ydata_ = ydata_[xmask]
            ymask = np.isfinite(ydata_)
            xdata = xdata[ymask]
            ydata_ = ydata_[ymask]

            mask1 = (xdata>np.nanpercentile(xdata,[1]))&(xdata<np.nanpercentile(xdata,[99]))
            mask2 = (ydata_>np.nanpercentile(ydata_,[1]))&(ydata_<np.nanpercentile(ydata_,[99]))
            mask = mask1&mask2
            xdata = xdata[mask]
            ydata_ = ydata_[mask]

            ax.scatter(xdata, ydata_, c='k', alpha=0.02, s=1, rasterized=True)

            bin_values = get_bins(xdata, ydata_, 20, equal_spaced=True)
            ax.scatter(bin_values[0], bin_values[1], ec='C0', fc='none', s=50, rasterized=True)
            ax.plot(bin_values[0], bin_values[1], c='C0', ls='--', rasterized=True)

            ax.set_xlim(np.nanpercentile(xdata,[0.5,99.5]))
            ax.set_ylim(np.nanpercentile(ydata_,[0.5,99.5]))

            ax.set_xlabel('MUSE')

        ax4.set_ylabel('%s HST/MUSE' %filter)
        ax5.set_ylabel('%s HST - MUSE' %filter)
        ax6.set_ylabel('%s (HST - MUSE)/MUSE' %filter)

        plt.tight_layout(w_pad=0.15, h_pad=0.15)
        fig.savefig(rootdir+appdir+'/figs/diffs_%s.png' %filter, bbox_inches='tight')
        plt.close('all')
        

def make_plots_muse_comp(hdu_muse_ha_contsub, hdu_muse_ha, matched=True, rootdir='./', appdir='hst_contsub/'):

    data1 = hdu_muse_ha_contsub.data
    data2 = hdu_muse_ha.data

    data1[data1==0] = np.nan
    data2[data2==0] = np.nan

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    if matched: 
        vmin = np.nanpercentile(data1, 0.1)
        vmax = np.nanpercentile(data1, 98)

        ax1.imshow(data1, vmin=vmin, vmax=vmax, origin='lower', cmap='inferno')
        img = ax2.imshow(data2, vmin=vmin, vmax=vmax, origin='lower', cmap='inferno')

    else:

        data1 = np.sqrt(data1)
        data2 = np.sqrt(data1)

        vmin1 = np.nanpercentile(data1, 0.1)
        vmax1 = np.nanpercentile(data1, 98)

        vmin2 = np.nanpercentile(data2, 0.1)
        vmax2 = np.nanpercentile(data2, 98)

        ax1.imshow(data1, vmin=vmin1, vmax=vmax1, origin='lower', cmap='inferno')
        ax2.imshow(data2, vmin=vmin2, vmax=vmax2, origin='lower', cmap='inferno')

    ax1.set_yticks([])
    ax1.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xticks([])

    bbox = dict(boxstyle='round', fc="w", ec="k")
    ax1.text(0.5, 0.95, 'MUSE Ha contsub', transform=ax1.transAxes, va='top', ha='center', bbox=bbox)
    ax2.text(0.5, 0.95, 'MUSE Ha', transform=ax2.transAxes, va='top', ha='center', bbox=bbox)

    plt.tight_layout()
    if matched: 
        fig.savefig(rootdir+appdir+'figs/comp_muse_ha_matchedcolor.png', bbox_inches='tight', dpi=150)
    else: 
        fig.savefig(rootdir+appdir+'figs/comp_muse_ha.png', bbox_inches='tight', dpi=150)
    plt.close('all')


def make_plots_hst_comp(hdu_hst_an_halpha, hdu_hst_halpha, matched=True, rootdir='./', appdir='hst_contsub/'):

    data1 = hdu_hst_an_halpha.data
    data2 = hdu_hst_halpha.data

    data1[data1==0] = np.nan
    data2[data2==0] = np.nan

    fig = plt.figure(figsize=(10, 5))
    ax1 = fig.add_subplot(1, 2, 1)
    ax2 = fig.add_subplot(1, 2, 2)

    if matched: 
        vmin = np.nanpercentile(data1, 0.1)
        vmax = np.nanpercentile(data1, 98)

        ax1.imshow(data1, vmin=vmin, vmax=vmax, origin='lower', cmap='inferno')
        img = ax2.imshow(data2, vmin=vmin, vmax=vmax, origin='lower', cmap='inferno')

    else:

        vmin1 = np.nanpercentile(data1, 0.1)
        vmax1 = np.nanpercentile(data1, 98)

        vmin2 = np.nanpercentile(data2, 0.1)
        vmax2 = np.nanpercentile(data2, 98)

        ax1.imshow(data1, vmin=vmin1, vmax=vmax1, origin='lower', cmap='inferno')
        ax2.imshow(data2, vmin=vmin2, vmax=vmax2, origin='lower', cmap='inferno')

    ax1.set_yticks([])
    ax1.set_xticks([])
    ax2.set_yticks([])
    ax2.set_xticks([])

    bbox = dict(boxstyle='round', fc="w", ec="k")
    ax1.text(0.5, 0.95, 'HST Ha contsub anchored', transform=ax1.transAxes, va='top', ha='center', bbox=bbox)
    ax2.text(0.5, 0.95, 'MUSE Ha contsub w/o anchored', transform=ax2.transAxes, va='top', ha='center', bbox=bbox)

    plt.tight_layout()
    if matched: 
        fig.savefig(rootdir+appdir+'figs/comp_hst_ha_matchedcolor.png', bbox_inches='tight', dpi=150)
    else: 
        fig.savefig(rootdir+appdir+'figs/comp_hst_ha.png', bbox_inches='tight', dpi=150)
    plt.close('all')


def make_plots_fluxsubmasks(hdu1, hdu2, hdu_neb, filter='', rootdir='./', appdir='hst_contsub/'):

    data_neb = hdu_neb.data.copy()
    ids = np.unique(data_neb)
    ids.sort()
    ids = ids[1:]

    flux_1 = np.ones(len(ids))
    flux_2 = np.ones(len(ids))

    for i in range(len(ids)): 
        mask = data_neb == ids[i]
        flux_1[i] = np.nansum(hdu1.data[mask])
        flux_2[i] = np.nansum(hdu2.data[mask])

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.scatter(flux_2[flux_2!=0], flux_2[flux_2!=0]/flux_1[flux_2!=0], fc='none', ec='C0')

    ax.grid(True, ls=':', color='k', alpha=0.2)
    ax.set_ylim(0, 2)
    ax.set_xscale('log')

    plt.tight_layout()
    fig.savefig(rootdir+appdir+'/figs/flux_nebmasks_%s.png' %filter, bbox_inches='tight', dpi=300)
    plt.close('all')

    median = np.nanmedian(flux_2[flux_2!=0]/flux_1[flux_2!=0])
    print('[INFO] Median flux ratio: %0.3f' %(median))


def make_plots_map(hdu, galaxy, filter, rootdir='./', appdir='hst_contsub/', smooth=True, smooth_factor=3, resolution = 0.07*u.arcsec):

    cmap1 = plt.cm.binary(np.linspace(0., 1, 16))
    cmap2 = colorcet.cm.fire(np.linspace(0, 1, 256))
    cmaplist = np.vstack((cmap1, cmap2))
    cmap = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', cmaplist)
    cmap.set_under(cmap(0))
    cmap.set_over(cmap(-1))
    cmap.set_bad(color=cmap(0))


    hdu = remove_nan_padding(hdu)
    hdu = get_smooth(hdu, resolution, resolution*smooth_factor)
    data = hdu.data

    data[data==0] = np.nan
    data = (data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    data[~np.isfinite(data)] = np.nan
    data = np.log10(data)
    data[~np.isfinite(data)] = np.nan
    data = (data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    data = data**0.001
    data = (data-np.nanmin(data))/(np.nanmax(data)-np.nanmin(data))
    hdu.data = data

    fig = plt.figure(figsize=(10, 10))
    ax = aplpy.FITSFigure(hdu, figure=fig)

    vmin = np.nanpercentile(data, 1)
    vmax = np.nanpercentile(data, 99.99)
    ax.show_colorscale(vmin=vmin, vmax=vmax, cmap=cmap, interpolation='none')

    bbox = dict(boxstyle='round', fc="w", ec="k")
    ax.add_label(0.05, 0.95, galaxy,  ha='left', va='top', size=20, bbox = bbox, relative=True)
    ax.ticks.set_color('black')

    ax.tick_labels.set_xformat('hh:mm:ss')
    ax.tick_labels.set_yformat('hh:mm:ss')

    ax_plot = fig.get_axes()[0]
    ax_plot.grid(True, ls=':', color='k', alpha=0.3)

    fig.savefig(rootdir+appdir+'/figs/%s_map_%s.png' %(galaxy, filter), bbox_inches='tight')


def make_histogram(data, bins=1000, plot=True, percentilecut=True):
    """
    Generate histogram data from a FITS HDU image.
    
    Parameters:
    - hdu (HDU): The input HDU.
    - bins (int): Number of bins for the histogram.
    - plot (bool): Whether to plot the histogram.
    
    Returns:
    - hist (tuple): Histogram values and bin edges.
    """
    data = data.flatten()
    if percentilecut:
        data = data[data>np.nanpercentile(data, 0.01)]
        data = data[data<np.nanpercentile(data, 99.99)]
    hist = np.histogram(data[~np.isnan(data)], bins=bins)
    
    if plot:
        plt.hist(data[~np.isnan(data)], bins=bins)
        plt.xlabel("Pixel Value")
        plt.ylabel("Frequency")
        plt.title("Histogram of data")
        # plt.yscale('log')
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
        ax.grid(True, ls=':', color='k', alpha=0.2)
         
    axes[0].set_ylabel("Frequency")
    axes[0].legend()
    axes[0].set_xlim([bin_centers.min(), bin_centers.max()])
    # axes[1].set_xlim([-100, 100])
    axes[1].set_xlim([popt[1]-popt[2]*4,popt[1]+popt[2]*4])
    
    plt.tight_layout()
    
    return popt, fig

def make_plots_histogram(hdu, hdu_muse_neb_re, hdu_muse_stars_re, filter='',  rootdir='./', appdir='hst_contsub/'):

    mask1 = hdu_muse_neb_re.data.copy() == -1
    mask2 = hdu_muse_stars_re.data.copy() != 1
    mask = ~mask1 & mask2
    mask = mask1 & mask2

    data = hdu.data.copy()
    data_ = data.copy()
    data_[~mask] = np.nan

    hist = make_histogram(data_, bins=1000, plot=False)
    stats = get_statistics(hist)
    fit_params, fig = fit_gaussian_to_histogram(hist, stats)
    
    axes = fig.get_axes()
    bbox = dict(boxstyle='round', fc="w", ec="k")
    axes[0].text(0.98, 0.75, filter, transform=axes[0].transAxes, va='top', ha='right', bbox=bbox)

    plt.tight_layout(w_pad=0.15, h_pad=0.15)
    fig.savefig(rootdir+appdir+'/figs/histo_%s.png' %filter, bbox_inches='tight')
    plt.close('all')