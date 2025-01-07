import numpy as np
from astropy.io import fits
from astropy.io import ascii
from scipy.integrate import simps
from glob import glob
from matplotlib import pyplot as plt
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy import constants as c
import matplotlib.colors as mcolors
import os
from astropy.table import Table, vstack


lightspeed=c.c.to(u.Angstrom/u.s)

from convolve_filter import readfilter,filter_flux,write_image_from_cube

#googledrive='/Users/00102515/Google Drive/Shared drives/PHANGS/Archive/'
#MUSEcubedir=googledrive+'MUSE/Dr2.2/copt/datacubes/'
PHANGSdrive='/Volumes/Mythfit/PHANGS/'
MUSEcubedir=PHANGSdrive+'MUSE/datacubes/'
filterdir='./'
filters=['F555W','F657N','F658N','F814W']
out_path='./'

cubefiles=glob(MUSEcubedir+'*.fits')
galaxies=[]
galfiles=[]
for file in cubefiles:
    temp=(file.split('/'))[-1].split('.fits')
    galfiles.append(temp[0])
    galaxies.append(((temp[0].split('-'))[0]))      
galaxies=np.array(galaxies)
galfiles=np.array(galfiles)
ind=(np.argsort(galaxies))
galaxies=galaxies[ind]
galfiles=galfiles[ind]

i=3
cubefile=MUSEcubedir+galfiles[i]+'.fits'
galaxy=galaxies[i]
cubedat,cubehd=fits.getdata(cubefile, header=True)
cubeunits=1e-20*u.erg/(u.s*u.cm**2*u.Angstrom)
wave=np.arange(cubehd['NAXIS3'])*cubehd['CD3_3']+cubehd['CRVAL3']
for HSTfilter in filters:
   filterImage=np.zeros(np.shape(cubedat)[1:3])*np.NAN
   filtwave,filt=readfilter('HST_WFC3_UVIS1.'+HSTfilter+'.dat',filterdir=filterdir)
   for j in range((np.shape(cubedat))[-1]):
      spec=cubedat[:,:,j]
      if np.count_nonzero(np.isfinite(spec))>0:
        filterImage[:,j],lam0=filter_flux(spec,wave*u.Angstrom,filt,filtwave)
   outfile=out_path+galaxy+'/'+galfiles[i]+'_'+HSTfilter+'.fits'
   write_image_from_cube(filterImage*1000.0,cubehd,target=galaxy,BUNIT='Jansky',filtername=HSTfilter,outfile=outfile)
