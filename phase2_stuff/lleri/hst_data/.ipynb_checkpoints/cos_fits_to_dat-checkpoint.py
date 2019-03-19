import numpy as np
import matplotlib.pyplot as plt
import os
#import pyfits as fits
import astropy.io.fits as fits
from astropy.convolution import convolve, Box1DKernel
import glob

def boxcar(flux,factor):
    smoothed_flux=convolve(flux,Box1DKernel(factor))
    return smoothed_flux

def ensure_dir(d):
    if not os.path.exists(d):
        os.makedirs(d)

def filewriter(wavelength, flux, error, dq, save_path=os.getcwd()+'/', filename=''): #writes spectra to file in format (wavelength   flux   fluxerror)
    ensure_dir(save_path)
    fl=open((save_path+filename),'w') #makes a file to put the values of the spectra in
    for w, f, e, q in zip(wavelength, flux, error, dq):
        fl.write('%f %g %g %i\n'%(w,f,e,q))

def plot_spectrum(wavelength, flux, error, dq, rootname):
    flux = boxcar(flux, 5)
    plt.figure(rootname)
    plt.plot(wavelength, flux, '0.5')
    wavelength_dq, flux_dq = wavelength[dq==0], flux[dq==0]
    plt.plot(wavelength_dq, flux_dq, label = 'spectrum')
    plt.plot(wavelength, error, label='error')  
    plt.legend()
    plt.show()

stars = ['ll_eri']

for k in range(len(stars)):
    star=stars[k]
    print(star)
    data_loc=''#'../c25_hst/'+star+'/hst/'
    for spectrum in glob.glob(data_loc+'*x1dsum.fit*'):

        hdul=fits.open(spectrum)

        #details of observation from headers
        date=hdul[1].header['DATE-OBS']
        cenwave = hdul[0].header['CENWAVE']
        arm = hdul[0].header['DETECTOR']

        #data
        data=hdul[1].data      
        w1, f1, e1, dq1 = data['wavelength'][0], data['flux'][0], data['error'][0], data['dq'][0]
        w2, f2, e2, dq2 = data['wavelength'][1], data['flux'][1], data['error'][1], data['dq'][1]
        w = np.concatenate((w2, w1))
        f = np.concatenate((f2, f1))
        e = np.concatenate((e2, e1))
        dq = np.concatenate((dq2, dq1))
        
        #tidy up
        hdul.close()

        #write to file
        filename=(star+'_COS_'+arm+'_'+str(cenwave)+'_'+date+'.dat')
        save_path = data_loc+'spectra/'  
        filewriter(w, f, e, dq, save_path = save_path, filename = filename  )

        #make a plot
        plot_spectrum(w, f, e, dq, filename)
