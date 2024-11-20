import numpy as np
from astropy.io import fits

def conv_inverse_variance_to_error(hdu):

    data = hdu.data.copy()
    hdr = hdu.header.copy()

    err = np.sqrt(1/data) 
    hdu_err = fits.PrimaryHDU(err, header=hdr)

    return hdu_err