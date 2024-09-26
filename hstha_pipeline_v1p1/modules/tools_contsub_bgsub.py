from photutils.background import Background2D
from scipy.ndimage import binary_dilation
from astropy.stats import mad_std
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

import warnings 
warnings.filterwarnings('ignore')

def get_bgmask(data, plot=False):

    rms = mad_std(data, ignore_nan=True)
    rms = mad_std(data[data<rms*2], ignore_nan=True)
    mask_high = data > rms*5
    mask_low = data > rms
    mask = binary_dilation(mask_high, mask=mask_low, iterations=-1)

    if plot: 
        # make plot showing the mask and data 
        fig, ax = plt.subplots(1, 2, figsize=(10, 5))
        ax[0].imshow(data, origin='lower', cmap='turbo', vmin=0, vmax=rms*5)
        ax[1].imshow(mask, origin='lower', cmap='Greys', vmin=0, vmax=1)

    return mask, rms

def get_bgsub(hdu, box_size=(50, 50), filter_size=(25, 25), plot=False, want_bg=False):

    data = hdu.data
    hdr = hdu.header

    mask, rms = get_bgmask(data, plot=False)

    bkg = Background2D(data, box_size=box_size, filter_size=filter_size, mask=mask)

    data_bg = bkg.background
    data_bgsub = data - data_bg
    coverage_mask = np.isnan(data)
    data_bg[coverage_mask] = np.nan

    if plot:
        # make plot showing the mask, data, data_bg and data_bgsub
        fig, ax = plt.subplots(1, 4, figsize=(20, 5))
        ax = ax.flatten()
        ax[0].imshow(data, origin='lower', cmap='turbo', vmin=0, vmax=rms*3)
        ax[1].imshow(mask, origin='lower', cmap='Greys', vmin=0, vmax=1)
        ax[2].imshow(data_bg, origin='lower', cmap='turbo', vmin=np.nanmin(data_bg), vmax=np.nanmax(data_bg))
        ax[3].imshow(data_bgsub, origin='lower', cmap='turbo', vmin=0, vmax=rms*3)

    hdu_bg = fits.PrimaryHDU(data_bg, header=hdr)
    hdu_bgsub = fits.PrimaryHDU(data_bgsub, header=hdr)

    if want_bg:
        return hdu_bg, hdu_bgsub
    else:
        return hdu_bgsub
