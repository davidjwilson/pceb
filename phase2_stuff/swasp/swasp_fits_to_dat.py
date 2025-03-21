import pyfits as fits #incase I want to use it on work desktop
#from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os

#turns SWASP fits files into dat files. 


star='UZ_Sex'
#for fits_file in files:
#	if fits_file[-4:] == 'mxlo':		
hdulist = fits.open('1SWASPJ102834.89-000029.1.fits')
scidata = hdulist[1].data
t=scidata['TMID']
for fluxnum in ['FLUX1', 'FLUX2', 'FLUX3']:
  plt.errorbar(t, scidata[fluxnum], yerr=scidata[fluxnum+'_ERR'],  label=fluxnum, marker='o', ls ='none', capsize=0)
plt.legend()    
    
    
"""
f = scidata['FLUX'][0]
w = [scidata['WAVELENGTH'][0]+i*scidata['DELTAW'][0] for i in xrange(len(f))]
#e = scidata['SIGMA'][0]
#print w
w, f, = np.array(w), np.array(f)
w, f = w[(f>0)], f[(f>0)]

hdulist.close()
plt.figure(star+'_iue')
plt.plot(w, f, drawstyle='steps-mid')
#plt.plot(w,e)


#fl=open('iue/dat_files/'+star+'_'+str(fits_file[0:2])+'_iue.dat', 'w')
#for j in xrange(len(w)):
#  fl.write((str(w[j])+'   '+str(f[j])+'\n'))
"""
plt.show()


	
