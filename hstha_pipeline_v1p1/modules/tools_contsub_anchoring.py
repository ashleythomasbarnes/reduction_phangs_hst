import astropy.wcs as wcs
import numpy as np
from matplotlib import pyplot as plt
from astropy.modeling import models, fitting
from astropy.table import Table, vstack 
from matplotlib.ticker import (MultipleLocator, LogLocator)
import warnings 
from scipy.optimize import curve_fit
warnings.filterwarnings('ignore')


def save_fittables_offsets(fit_f555w, fit_f65Xn, fit_f814w, rootdir='./', appdir='hst_contsub/'):
    
    table_fits = vstack([fit_f555w, fit_f65Xn, fit_f814w])
    table_fits.write(rootdir+appdir+'/table_fits_filteranchor.fits', overwrite=True)
    table_fits.write(rootdir+appdir+'/table_fits_filteranchor.csv', overwrite=True)


def save_fittables_slope(fit_halpha, fit_an_halpha, rootdir='./', appdir='hst_contsub/'):
    
    table_fits = vstack([fit_halpha, fit_an_halpha])
    table_fits.write(rootdir+appdir+'/table_fits_haanchor.fits', overwrite=True)
    table_fits.write(rootdir+appdir+'/table_fits_haanchor.csv', overwrite=True)
    
    
def get_bins(data1, data2, num_bins, min_val=None, max_val=None, equal_spaced=True):

    # Calculate min and max if not provided
    if min_val is None:
        min_val = np.nanmin(data1)
    if max_val is None:
        max_val = np.nanmax(data1)

    mask = (data1>min_val) & (data1<max_val)
    data1 = data1.copy()[mask]
    data2 = data2.copy()[mask]

    # Generate bin edges
    # for equal_spaced bins = space spacing
    if equal_spaced:
        bin_edges = np.linspace(min_val, max_val, num_bins + 1)
    # for equal_number bins = same number of points
    else: 
        percentiles = np.linspace(0, 100, num_bins+1)
        bin_edges = np.percentile(data1, percentiles)

    # Calculate mean values within each bin
    binned_values1 = np.zeros(num_bins)
    binned_values2 = np.zeros(num_bins)

    for i in range(num_bins):

        bin_mask = ((data1>=bin_edges[i]) & (data1<=bin_edges[i+1]))
        # bin_mask = (data1>=bin_edges[i])
        bin_ids = np.where(bin_mask)

        values_in_bin1 = data1[bin_ids]
        values_in_bin2 = data2[bin_ids]

        binned_values1[i] = np.nanmedian(values_in_bin1)
        binned_values2[i] = np.nanmedian(values_in_bin2)

    return binned_values1, binned_values2

def get_anchoring_offset(hdu1, hdu2, hdu3, hdu_stars, filter='', rootdir='./', appdir='hst_contsub/', make_plots=True):

    ### 
    hdu1 = hdu1.copy()
    hdu2 = hdu2.copy()
    hdu3 = hdu3.copy()
    data1 = hdu1.data.copy()
    data2 = hdu2.data.copy()
    data3 = hdu3.data.copy()

    # Mask zeros 
    mask_zero1 = data1==0
    mask_zero2 = data2==0
    data1[(mask_zero1|mask_zero2)] = np.nan
    data2[(mask_zero1|mask_zero2)] = np.nan

    # Mask with starmask 
    mask_stars = hdu_stars.data!=0
    data1[mask_stars] = np.nan
    data2[mask_stars] = np.nan

    data1 = data1.flatten()
    data2 = data2.flatten()

    valid_indices = np.isfinite(data1) & np.isfinite(data2)
    data1 = data1[valid_indices]
    data2 = data2[valid_indices]

    # Mask to only lowest value points 
    x_per = np.percentile(data1, [0.01, 99])
    y_per = np.percentile(data2, [0.01, 99])

    x_mask = (data1>x_per[0])&(data1<x_per[1])
    y_mask = (data2>y_per[0])&(data2<y_per[1])

    data1 = data1[x_mask&y_mask]
    data2 = data2[x_mask&y_mask]

    # Get bins with equal number of points in each bin 
    min_val, max_val = np.percentile(data1, [0, 65]) 
    bin_values = get_bins(data1, data2, 25, equal_spaced=False, min_val=min_val, max_val=max_val)

    # Fit binned data
    model_poly = models.Polynomial1D(degree=1)
    fitter_poly = fitting.LinearLSQFitter() 
    best_fit_poly_bins = fitter_poly(model_poly, bin_values[0], bin_values[1])
    intercept_bins, slope_bins = best_fit_poly_bins.parameters
    
    def func_fixed(x, b):
        # For fixed slope - y = x + b
        return x + b

    best_fitfixed_bins, _ = curve_fit(func_fixed, bin_values[0], bin_values[1])
    interceptfixed_bins = best_fitfixed_bins[0]

    x_fit = np.linspace(-1e3, 1e3, 10000)
    y_fit_bins = slope_bins * x_fit + intercept_bins
    ###

    # Extract the WCS information from the input and template headers
    wcs1 = wcs.WCS(hdu1.header)
    wcs3 = wcs.WCS(hdu3.header)
    pixscale1 = wcs.utils.proj_plane_pixel_area(wcs1.celestial)
    pixscale3 = wcs.utils.proj_plane_pixel_area(wcs3.celestial)

    pixscale_ratio = (pixscale3 / pixscale1)
    hdu3.data = hdu3.data - (intercept_bins*pixscale_ratio) # HST full resolution 
    hdu2.data = hdu2.data - (intercept_bins) # HST smoothed

    fit = [filter, slope_bins, intercept_bins, intercept_bins*pixscale_ratio, interceptfixed_bins]
    table_fit = Table(np.array(fit), names=['filter', 'slope_bins', 'intercept_lowres', 'intercept_highres', 'interceptfixed_bins'])
    #### 

    if make_plots: 

        fig = plt.figure(figsize=(12, 5))
        ax1 = fig.add_subplot(1, 4, 1)
        ax2 = fig.add_subplot(1, 4, 2)
        ax3 = fig.add_subplot(1, 4, 3)
        ax4 = fig.add_subplot(1, 4, 4)

        #data

        for ax in [ax1,ax2,ax3, ax4]:
            ax.scatter(data1, data2, c='k', alpha=0.01, s=1, rasterized=True)
            ax.scatter(bin_values[0], bin_values[1], fc='none', ec='C0', alpha=1, s=30, zorder=5)
            ax.plot(bin_values[0], bin_values[1], c='C0', alpha=1, zorder=5)

            # fits 
            offset = (0 - (intercept_bins)) / slope_bins
            ax.plot(x_fit, y_fit_bins, color='C0', linewidth=2, linestyle=':', label=f'y = {slope_bins:.4f}x + {intercept_bins:.4g}')
            ax.plot(x_fit, x_fit, 'k', linewidth=2, linestyle=':', label=f'y = x')
            ax.plot([offset, offset], [-100,0], color='C0', linewidth=2, linestyle=':', label=f'Offset = {offset:.4g}')

            ax.set_xlabel('Flux density (MUSE) [erg/s/cm-2/A/pix]')
            ax.set_ylabel('Flux density (HST smoothed, regrid) [erg/s/cm-2/A/pix]')
            ax.grid(True, ls=':', color='k', alpha=0.2)

        ax1.legend(title=filter, loc='upper left', fontsize=8)
        ax1.set_xlim(np.nanpercentile(data1, [0,75]))
        ax1.set_ylim(np.nanpercentile(data2, [0,75]))

        ax2.set_xlim(np.nanpercentile(data1, [0,30]))
        ax2.set_ylim(np.nanpercentile(data2, [0,30]))

        offsety = intercept_bins + (offset*slope_bins)
        ax3.set_xlim([offset-10, offset+10])
        ax3.set_ylim(offsety-10, offsety+10)

        ax4.set_xlim(np.nanpercentile(data1[data1>0], [0.01,99.99]))
        ax4.set_ylim(np.nanpercentile(data2[data2>0], [0.01,99.99]))

        ax4.set_xscale('log')
        ax4.set_yscale('log')

        plt.tight_layout()
        fig.savefig(rootdir+appdir+'/figs/fit_%s.png' %filter, bbox_inches='tight')
        plt.close('all')

    return(hdu3, hdu2, table_fit)

def get_anchoring_slope(hdu1, hdu2, hdu3, hdu_neb, filter='', rootdir='./', appdir='hst_contsub/', make_plots=True):

    ### 
    hdu1 = hdu1.copy()
    hdu2 = hdu2.copy()
    hdu3 = hdu3.copy()
    data1 = hdu1.data.copy()
    data2 = hdu2.data.copy()
    data3 = hdu3.data.copy()

    # Mask zeros 
    mask_zero1 = data1==0
    mask_zero2 = data2==0
    data1[(mask_zero1&mask_zero2)] = np.nan
    data2[(mask_zero1&mask_zero2)] = np.nan

    # Mask with nebmask 
    mask_neb = hdu_neb.data==-1
    data1[mask_neb] = np.nan
    data2[mask_neb] = np.nan

    valid_indices = np.isfinite(data1) & np.isfinite(data2)
    data1 = data1[valid_indices]
    data2 = data2[valid_indices]

    # Mask to only lowest value points 
    x_per = np.percentile(data1, [0.1, 99.9])
    y_per = np.percentile(data2, [0.1, 99.9])

    x_mask = (data1>x_per[0])&(data1<x_per[1])
    y_mask = (data2>y_per[0])&(data2<y_per[1])

    data1 = data1[x_mask&y_mask]
    data2 = data2[x_mask&y_mask]

    # Get bins with equal number of points in each bin 
    min_val, max_val = np.percentile(data1, [10, 90]) 
    # bin_values = get_bins(data1, data2, 20, equal_spaced=True, min_val=min_val, max_val=max_val)
    bin_values = get_bins(data1, data2, 20, equal_spaced=False, min_val=min_val, max_val=max_val)

    # Fit binned data
    model_poly = models.Polynomial1D(degree=1)
    fitter_poly = fitting.LinearLSQFitter() 
    best_fit_poly_bins = fitter_poly(model_poly, bin_values[0], bin_values[1])
    intercept_bins, slope_bins = best_fit_poly_bins.parameters

    def func_fixed(x, a):
        # For fixed offeset - y = xa + 0
        return a*x

    best_fitfixed_bins, _ = curve_fit(func_fixed, bin_values[0], bin_values[1])
    slopefixed_bins = best_fitfixed_bins[0]

    x_fit = np.linspace(-1e3, np.nanmax(data2), 10000)
    y_fit_bins = slope_bins * x_fit + intercept_bins
    ###

    # Extract the WCS information from the input and template headers
    wcs1 = wcs.WCS(hdu1.header)
    wcs3 = wcs.WCS(hdu3.header)
    pixscale1 = wcs.utils.proj_plane_pixel_area(wcs1.celestial)
    pixscale3 = wcs.utils.proj_plane_pixel_area(wcs3.celestial)

    pixscale_ratio = (pixscale3 / pixscale1)
    # hdu3.data = (hdu3.data - (intercept_bins*pixscale_ratio)) / slope_bins # HST full resolution 
    # hdu2.data = (hdu2.data - (intercept_bins)) / slope_bins # HST smoothed
    hdu3.data = hdu3.data / slope_bins # HST full resolution 
    hdu2.data = hdu2.data / slope_bins # HST smoothed

    fit = [filter, slope_bins, intercept_bins, intercept_bins*pixscale_ratio, slopefixed_bins]
    table_fit = Table(np.array(fit), names=['filter', 'slope_bins', 'intercept_lowres', 'intercept_highres', 'slopefixed_bins'])
    #### 

    if make_plots: 

        fig = plt.figure(figsize=(10, 5))
        ax1 = fig.add_subplot(1, 2, 1)
        ax2 = fig.add_subplot(1, 2, 2)

        for ax in [ax1,ax2]:

            #data
            ax.scatter(data1, data2, c='k', alpha=0.01, s=1, rasterized=True)

            #bins 
            ax.scatter(bin_values[0], bin_values[1], fc='none', ec='C0', alpha=1, s=30, zorder=5)
            ax.plot(bin_values[0], bin_values[1], c='C0', alpha=1, zorder=5)

            # fits 
            offset = (0 - (intercept_bins)) / slope_bins
            ax.plot(x_fit, y_fit_bins, color='C0', linewidth=2, linestyle=':', label=f'y = {slope_bins:.4f}x + {intercept_bins:.4g}')
            ax.plot(x_fit, x_fit, 'k', linewidth=2, linestyle=':', label=f'y = x')
            ax.plot([offset, offset], [-100,0], color='C0', linewidth=2, linestyle=':', label=f'Offset = {offset:.4g}')

            for f in [0.1,0.5,2,10]:
                ax.plot(x_fit, x_fit*f, 'k', linewidth=2, linestyle='--', alpha=0.1)

            ax.set_xlabel('Flux (MUSE Ha) [erg/s/cm-2/pix]')
            ax.set_ylabel('Flux (MUSE contsub) [erg/s/cm-2/pix]')
            ax.legend(title=filter, loc='upper left')
            ax.grid(True, ls=':', color='k', alpha=0.2)

        ax1.set_xlim(np.nanpercentile(data1, [0,99]))
        ax1.set_ylim(np.nanpercentile(data2, [0,99]))

        ax2.set_xlim(np.nanpercentile(data1[data1>0], [0.01,99.99]))
        ax2.set_ylim(np.nanpercentile(data2[data2>0], [0.01,99.99]))

        ax2.set_xscale('log')
        ax2.set_yscale('log')

        plt.tight_layout()
        fig.savefig(rootdir+appdir+'/figs/fit_%s.png' %filter, bbox_inches='tight')
        plt.close('all')

    return(hdu3, hdu2, table_fit)

# def get_anchoring_slope(hdu1, hdu2, hdu3, hdu_neb, filter='', rootdir='./', appdir='hst_contsub/', make_plots=True):

#     ### 
#     hdu1 = hdu1.copy()
#     hdu2 = hdu2.copy()
#     hdu3 = hdu3.copy()
#     data1 = hdu1.data.copy()
#     data2 = hdu2.data.copy()
#     data3 = hdu3.data.copy()

#     # Mask zeros 
#     mask_zero1 = data1==0
#     mask_zero2 = data2==0
#     data1[(mask_zero1&mask_zero2)] = np.nan
#     data2[(mask_zero1&mask_zero2)] = np.nan

#     # Mask with nebmask 
#     mask_neb = hdu_neb.data==-1
#     data1[mask_neb] = np.nan
#     data2[mask_neb] = np.nan

#     valid_indices = np.isfinite(data1) & np.isfinite(data2)
#     data1 = data1[valid_indices]
#     data2 = data2[valid_indices]

#     # Mask to only lowest value points 
#     # x_per = np.percentile(data1, [0, 100])
#     # y_per = np.percentile(data2, [0, 100])
#     x_per = np.percentile(data1, [25, 75])
#     y_per = np.percentile(data2, [25, 75])

#     x_mask = (data1>x_per[0])&(data1<x_per[1])
#     y_mask = (data2>y_per[0])&(data2<y_per[1])

#     data1 = data1[x_mask&y_mask]
#     data2 = data2[x_mask&y_mask]

#     # Get bins with equal number of points in each bin 
#     bin_values = get_bins(data1, data2, 20, equal_spaced=True)

#     # Fit binned data
#     model_poly = models.Polynomial1D(degree=1)
#     fitter_poly = fitting.LinearLSQFitter() 
#     best_fit_poly_bins = fitter_poly(model_poly, bin_values[0], bin_values[1])
#     intercept_bins, slope_bins = best_fit_poly_bins.parameters

#     x_fit = np.linspace(-1e3, np.nanmax(data2), 10000)
#     y_fit_bins = slope_bins * x_fit + intercept_bins
#     ###

#     # Extract the WCS information from the input and template headers
#     wcs1 = wcs.WCS(hdu1.header)
#     wcs3 = wcs.WCS(hdu3.header)
#     pixscale1 = wcs.utils.proj_plane_pixel_area(wcs1.celestial)
#     pixscale3 = wcs.utils.proj_plane_pixel_area(wcs3.celestial)

#     pixscale_ratio = (pixscale3 / pixscale1)
#     offset1 = (0 - (intercept_bins*pixscale_ratio)) / slope_bins
#     offset2 = (0 - (intercept_bins)) / slope_bins

#     hdu3.data = (hdu3.data - (intercept_bins*pixscale_ratio)) / slope_bins # HST full resolution 
#     hdu2.data = (hdu2.data - (intercept_bins*pixscale_ratio)) / slope_bins # HST smoothed

#     fit = [filter, slope_bins, intercept_bins, intercept_bins*pixscale_ratio, offset1, offset2]
#     table_fit = Table(np.array(fit), names=['filter', 'slope_bins', 'intercept_lowres', 'intercept_highres', 'offset_lowres', 'offset_highres'])
#     #### 

#     if make_plots: 

#         fig = plt.figure(figsize=(10, 5))
#         ax1 = fig.add_subplot(1, 2, 1)
#         ax2 = fig.add_subplot(1, 2, 2)

#         for ax in [ax1,ax2]:

#             #data
#             ax.scatter(data1, data2, c='k', alpha=0.01, s=1, rasterized=True)

#             #bins 
#             ax.scatter(bin_values[0], bin_values[1], fc='none', ec='C0', alpha=1, s=30, zorder=5)
#             ax.plot(bin_values[0], bin_values[1], c='C0', alpha=1, zorder=5)

#             # fits 
#             offset = (0 - (intercept_bins)) / slope_bins
#             ax.plot(x_fit, y_fit_bins, color='C0', linewidth=2, linestyle=':', label=f'y = {slope_bins:.4f}x + {intercept_bins:.4g}')
#             ax.plot(x_fit, x_fit, 'k', linewidth=2, linestyle=':', label=f'y = x')
#             ax.plot([offset, offset], [-100,0], color='C0', linewidth=2, linestyle=':', label=f'Offset = {offset:.4g}')

#             for f in [0.5,2]:
#                 ax.plot(x_fit, x_fit*f, 'k', linewidth=2, linestyle='--', alpha=0.1)

#             ax.set_xlabel('Flux (MUSE Ha) [erg/s/cm-2/pix]')
#             ax.set_ylabel('Flux (MUSE contsub) [erg/s/cm-2/pix]')
#             ax.legend(title=filter, loc='upper left')
#             ax.grid(True, ls=':', color='k', alpha=0.2)

#         ax2.set_xscale('log')
#         ax2.set_yscale('log')

#         plt.tight_layout()
#         fig.savefig(rootdir+appdir+'/figs/fit_%s.png' %filter, bbox_inches='tight')
#         plt.close('all')

#     return(hdu3, hdu2, table_fit)

# def get_anchoring_slope(hdu1, hdu2, hdu3, hdu_neb, filter='', rootdir='./', appdir='hst_contsub/', make_plots=True):

#     ### 
#     hdu1 = hdu1.copy()
#     hdu2 = hdu2.copy()
#     hdu3 = hdu3.copy()
#     data1 = hdu1.data.copy()
#     data2 = hdu2.data.copy()
#     data3 = hdu3.data.copy()

#     # Mask zeros 
#     mask_zero1 = data1==0
#     mask_zero2 = data2==0
#     data1[(mask_zero1&mask_zero2)] = np.nan
#     data2[(mask_zero1&mask_zero2)] = np.nan

#     # Mask with nebmask 
#     mask_neb = hdu_neb.data==-1
#     data1[mask_neb] = np.nan
#     data2[mask_neb] = np.nan

#     data_neb = hdu_neb.data
#     mids = list(np.unique(data_neb))
#     mids.remove(-1)
#     data1_nebs = []
#     data2_nebs = []
#     for i, mid in enumerate(mids): 
#         data1_nebs += [np.nansum(data1[data_neb==mid])]
#         data2_nebs += [np.nansum(data2[data_neb==mid])]

#     data1_nebs = np.array(data1_nebs)
#     data2_nebs = np.array(data2_nebs)

#     # Get bins with equal number of points in each bin 
#     min_val, max_val = np.nanpercentile(data1_nebs, [0.1, 99.9])
#     bin_values = get_bins(data1_nebs, data2_nebs, 20, equal_spaced=False, min_val=min_val, max_val=max_val)

#     # Fit binned data
#     model_poly = models.Polynomial1D(degree=1)
#     fitter_poly = fitting.LinearLSQFitter() 
#     best_fit_poly_bins = fitter_poly(model_poly, bin_values[0], bin_values[1])
#     intercept_bins, slope_bins = best_fit_poly_bins.parameters

#     x_fit = np.linspace(-1e3, np.nanmax(data1_nebs)*1.2, 10000)
#     y_fit_bins = slope_bins * x_fit + intercept_bins
#     ##

#     # Extract the WCS information from the input and template headers
#     wcs1 = wcs.WCS(hdu1.header)
#     wcs3 = wcs.WCS(hdu3.header)
#     pixscale1 = wcs.utils.proj_plane_pixel_area(wcs1.celestial)
#     pixscale3 = wcs.utils.proj_plane_pixel_area(wcs3.celestial)

#     pixscale_ratio = (pixscale3 / pixscale1)
#     offset1 = (0 - (intercept_bins*pixscale_ratio)) / slope_bins
#     offset2 = (0 - (intercept_bins)) / slope_bins
#     hdu3.data = hdu3.data / slope_bins # HST full resolution 
#     hdu2.data = hdu2.data / slope_bins # HST smoothed

#     fit = [filter, slope_bins, intercept_bins, intercept_bins*pixscale_ratio, offset1, offset2]
#     table_fit = Table(np.array(fit), names=['filter', 'slope_bins', 'intercept_lowres', 'intercept_highres', 'offset_lowres', 'offset_highres'])
#     #### 

#     if make_plots: 

#         fig = plt.figure(figsize=(10, 5))
#         ax1 = fig.add_subplot(1, 2, 1)
#         ax2 = fig.add_subplot(1, 2, 2)

#         for ax in [ax1,ax2]:

#             #data
#             ax.scatter(data1_nebs, data2_nebs, c='k', alpha=0.75, s=2, rasterized=True)

#             #bins 
#             ax.scatter(bin_values[0], bin_values[1], fc='none', ec='C0', alpha=1, s=30, zorder=5)
#             ax.plot(bin_values[0], bin_values[1], c='C0', alpha=1, zorder=5)

#             # fits 
#             offset = (0 - (intercept_bins)) / slope_bins
#             ax.plot(x_fit, y_fit_bins, color='C0', linewidth=2, linestyle='--', label=f'y = {slope_bins:.4f}x + {intercept_bins:.4g}', alpha=0.5)
#             ax.plot(x_fit, x_fit, 'k', linewidth=2, linestyle='--', label=f'y = x', alpha=0.5)

#             ax.set_xlabel('Flux (MUSE Ha) [erg/s/cm-2/pix]')
#             ax.set_ylabel('Flux (MUSE contsub) [erg/s/cm-2/pix]')
#             ax.legend(title=filter, loc='upper left')
#             ax.grid(True, ls=':', color='k', alpha=0.2)

#         ax2.set_xscale('log')
#         ax2.set_yscale('log')

#         plt.tight_layout()
#         fig.savefig(rootdir+appdir+'/figs/fit_%s.png' %filter, bbox_inches='tight')
#         plt.close('all')

#     return(hdu3, hdu2, table_fit)

# def get_anchoring_slope(hdu1, hdu2, hdu3, hdu_neb, filter='', rootdir='./', appdir='hst_contsub/', make_plots=True):

#     ### 
#     hdu1 = hdu1.copy()
#     hdu2 = hdu2.copy()
#     hdu3 = hdu3.copy()
#     data1 = hdu1.data.copy()
#     data2 = hdu2.data.copy()
#     data3 = hdu3.data.copy()

#     # Mask zeros 
#     mask_zero1 = data1==0
#     mask_zero2 = data2==0
#     data1[(mask_zero1&mask_zero2)] = np.nan
#     data2[(mask_zero1&mask_zero2)] = np.nan

#     # Mask with nebmask 
#     mask_neb = hdu_neb.data==-1
#     data1[mask_neb] = np.nan
#     data2[mask_neb] = np.nan

#     valid_indices = np.isfinite(data1) & np.isfinite(data2)
#     data1 = data1[valid_indices]
#     data2 = data2[valid_indices]

#     # Mask to only lowest value points 
#     x_per = np.percentile(data1, [0, 100])
#     y_per = np.percentile(data2, [0, 100])

#     x_mask = (data1>x_per[0])&(data1<x_per[1])
#     y_mask = (data2>y_per[0])&(data2<y_per[1])

#     data1 = data1[x_mask&y_mask]
#     data2 = data2[x_mask&y_mask]

#     # Get bins with equal number of points in each bin 
#     bin_values = get_bins(data1, data2, 20, equal_spaced=True)

#     # Fit binned data
#     model_poly = models.Polynomial1D(degree=1)
#     fitter_poly = fitting.LinearLSQFitter() 
#     best_fit_poly_bins = fitter_poly(model_poly, bin_values[0], bin_values[1])
#     intercept_bins, slope_bins = best_fit_poly_bins.parameters

#     x_fit = np.linspace(-1e3, np.nanmax(data2), 10000)
#     y_fit_bins = slope_bins * x_fit + intercept_bins
#     ###

#     # Extract the WCS information from the input and template headers
#     wcs1 = wcs.WCS(hdu1.header)
#     wcs3 = wcs.WCS(hdu3.header)
#     pixscale1 = wcs.utils.proj_plane_pixel_area(wcs1.celestial)
#     pixscale3 = wcs.utils.proj_plane_pixel_area(wcs3.celestial)

#     pixscale_ratio = (pixscale3 / pixscale1)
#     offset1 = (0 - (intercept_bins*pixscale_ratio)) / slope_bins
#     offset2 = (0 - (intercept_bins)) / slope_bins
#     hdu3.data = hdu3.data / slope_bins # HST full resolution 
#     hdu2.data = hdu2.data / slope_bins # HST smoothed

#     fit = [filter, slope_bins, intercept_bins, intercept_bins*pixscale_ratio, offset1, offset2]
#     table_fit = Table(np.array(fit), names=['filter', 'slope_bins', 'intercept_lowres', 'intercept_highres', 'offset_lowres', 'offset_highres'])
#     #### 

#     if make_plots: 

#         fig = plt.figure(figsize=(10, 5))
#         ax1 = fig.add_subplot(1, 2, 1)
#         ax2 = fig.add_subplot(1, 2, 2)

#         for ax in [ax1,ax2]:

#             #data
#             ax.scatter(data1, data2, c='k', alpha=0.01, s=1, rasterized=True)

#             #bins 
#             ax.scatter(bin_values[0], bin_values[1], fc='none', ec='C0', alpha=1, s=30, zorder=5)
#             ax.plot(bin_values[0], bin_values[1], c='C0', alpha=1, zorder=5)

#             # fits 
#             offset = (0 - (intercept_bins)) / slope_bins
#             ax.plot(x_fit, y_fit_bins, color='C0', linewidth=2, linestyle=':', label=f'y = {slope_bins:.4f}x + {intercept_bins:.4g}')
#             ax.plot(x_fit, x_fit, 'k', linewidth=2, linestyle=':', label=f'y = x')
#             ax.plot([offset, offset], [-100,0], color='C0', linewidth=2, linestyle=':', label=f'Offset = {offset:.4g}')

#             ax.set_xlabel('Flux (MUSE Ha) [erg/s/cm-2/pix]')
#             ax.set_ylabel('Flux (MUSE contsub) [erg/s/cm-2/pix]')
#             ax.legend(title=filter, loc='upper left')
#             ax.grid(True, ls=':', color='k', alpha=0.2)

#         ax2.set_xscale('log')
#         ax2.set_yscale('log')

#         plt.tight_layout()
#         fig.savefig(rootdir+appdir+'/figs/fit_%s.png' %filter, bbox_inches='tight')
#         plt.close('all')

#     return(hdu3, hdu2, table_fit)