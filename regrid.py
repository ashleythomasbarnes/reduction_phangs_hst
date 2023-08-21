from astropy.io import fits
from reproject import reproject_interp

def adjust_central_reference_pixel(input_filename, output_filename, shift_x, shift_y):
    """
    Adjusts the central reference pixel of a FITS image by shifting its coordinates.

    Args:
        input_filename (str): Path to the input FITS file.
        output_filename (str): Path to save the output FITS file.
        shift_x (float): Amount to shift in the x-direction.
        shift_y (float): Amount to shift in the y-direction.
    """
    print(f"[INFO] Adjusting the central reference pixel of {input_filename}")
    
    # Open the input FITS file
    hdu = fits.open(input_filename)[0]
    hdu_tmp = fits.open(input_filename)[0]

    # Update the CRPIX1 and CRPIX2 values
    hdu_tmp.header['CRPIX1'] += shift_x
    hdu_tmp.header['CRPIX2'] += shift_y

    # Reproject the image with the updated header
    print("[INFO] Reprojecting the image...")
    array, _ = reproject_interp(hdu_tmp, hdu.header, parallel=True)
    hdu_shifted = fits.PrimaryHDU(np.array(array, dtype=np.float32), hdu.header)

    # Save the shifted FITS file
    print(f"[INFO] Saving the shifted FITS file as {output_filename}")
    hdu_shifted.writeto(output_filename, overwrite=True)

    return None