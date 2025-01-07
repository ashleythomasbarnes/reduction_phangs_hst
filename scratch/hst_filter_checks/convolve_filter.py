import numpy as np
import fnmatch
from scipy.integrate import simps
from astropy import units as u
from astropy import constants as c
from astropy.io import fits
from astropy.wcs import WCS

lightspeed=c.c.to(u.Angstrom/u.s)

def readfilter(filtername,filterdir='/Users/00102515/Work/filters/'):
    try:
       waveF,filt=np.loadtxt(filterdir+filtername,unpack=True)
    except:
       print(filterdir+filtername+" not a SVO filter file.")
       return -1
    return waveF*u.Angstrom, filt
            
@u.quantity_input(waveS=u.Angstrom)
def filter_flux(spectrum, waveS,filt,waveF,Flam=False,calibspec='AB',effwave=True):
    newfilt=np.interp(waveS,waveF,filt,left=0.0,right=0.0)
    lam_eff=simps(filt)/simps(filt/waveF**2)*waveF.unit
    if calibspec=='AB':
        temp=np.nan_to_num(np.multiply(spectrum.transpose(),newfilt*waveS))
        flux=simps(temp,waveS)/simps(newfilt/waveS,waveS)/lightspeed.value
        funit=temp.unit/u.Hz
#    elif calibspec=='flat':
 # ToDo
    if Flam:
       flux=flux*(lightspeed.value)/(lam_eff.value)**2
       funit=funit*(u.Hz/lam_eff.unit)
    if effwave:
        return flux.transpose()*funit,lam_eff
    else:
        return flux.transpose()*funit

def write_image_from_cube(filterImage,cubehd,target='', BUNIT='',filtername='', outfile='Image.fits'):
    newhead=cubehd
    newhead['NAXIS']=2
    for key in list(newhead):
      if fnmatch.fnmatch(key,'*3*'):
        del newhead[key]
    newImage=fits.PrimaryHDU(data=filterImage,header=newhead)
    newImage.header['WCSAXES']=2
    if len(BUNIT)==0:
      try:
        newImage.header['BUNIT']=filterImage.unit.to_string('console')
      except:
        newImage.header['BUNIT']=''
    else:
      newImage.header['BUNIT']=BUNIT
    if len(target) > 0:
      newImage.header['TARGNAME']=target
    if len(filtername) > 0:
      newImage.header['FILTER']=filtername
    newImage.writeto(outfile,overwrite=True)
 
