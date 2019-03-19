#import pyfits as fits #incase I want to use it on work desktop
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
import os

#turns IUE fits files into dat files. 

stars = os.listdir('iue')
for star in stars:
	if star != 'dat_files':
		files = os.listdir('iue/'+star)
		for fits_file in files:
			if fits_file[-4:] == 'mxlo':		
				hdulist = fits.open('iue/'+star+'/'+fits_file)
				scidata = hdulist[1].data
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
				
				
				fl=open('iue/dat_files/'+star+'_'+str(fits_file[0:2])+'_iue.dat', 'w')
				for j in xrange(len(w)):
				  fl.write((str(w[j])+'   '+str(f[j])+'\n'))
				
				plt.show()
			
		
	
