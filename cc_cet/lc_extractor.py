#!/usr/bin/python
import numpy as np
import matplotlib.pyplot as plt
import glob
import astropy.io.fits as fits
import os
from astropy.table import Table
from astropy.io import ascii

"""
Script to extract a lightcurve from HST/COS corrtag files
Requires astropy, matplotlib, numpy. 
Saves each FP_POS exposure separately (for FUV exposures)
For each exposure the counts from the A and B  segments (FUV) or A, B and C traces (NUV), if present, are combined.
Airglow from Lymman alpha and ~1300A Oi is removed.
Error is photon noise only. 
Optional: Plots lightcurve. 

Usage: call the function lc_maker()

lc_extractor.lc_maker(file_path=os.getcwd()+'/', save_path=os.getcwd()+'/lightcurves/', bin_time=1., plot=True, save_file='long')

Arguments: 
	- file_path = string, where your corrtag files are. Default is the curret directory.
	- save_path = sring, where you want the output to be saves. 
	Default is a new "lightcurves" directory in the current directory
	- bin_time = float, time in s to bin the lightcurve to. Default is 1.0s.
	- qual_check = boolean, masks out flagged pixels. Default is True.
	- plot = boolean, makes a plot of the combined lightcurve. Default is True.
    - save_file = file name to save the lightcurves as. Default is 'long':
        - 'long' = targname_band_expstart_rootname_bintime_lc.dat'
        - 'short' = rootname_bintime_lc.dat
	
Outputs: 
	- Lightcurve of each exposure saved as [exposure rootname]_[bintime]s_lc.dat.
	Lightcurves saved as time(s since MJD=0) counts(s-1) error(s-1). 

"""

def region_mask(x, y, slope, intercept, height):
	mask = (y > slope*x+intercept-height/2.) & (y < slope*x+intercept+height/2.)
	return mask 

def filewriter(time, counts, error, save_path, filename): 
    # writes lightcurves to dat files
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    savdat = Table([time, counts, error], names=['TIME', 'COUNTS', 'ERROR'])
    ascii.write(savdat, save_path+filename, format='ecsv', overwrite=True)

def extract_counts(header, data, bins, qual_check, seg):
    slope = header['SP_SLP_'+seg]
    sp_intercept = header['SP_LOC_'+seg]
    sp_height = float(header['SP_HGT_'+seg])

    #background regions
    bk1_intercept = header['B_BKG1_'+seg]
    bk1_height = float(header['B_HGT1_'+seg])
    bk2_intercept = header['B_BKG1_'+seg]
    bk2_height = float(header['B_HGT1_'+seg])

    #data 
    x = data['XCORR']
    y = data['YCORR']
    time = data['TIME']
    w = data['WAVELENGTH']
    dq = data['DQ']

    #mask out flagged pixels
    if qual_check == True:
        x, y, time, w = x[dq==0], y[dq==0], time[dq==0], w[dq==0]

    #mask out airglow from lyman alpha and oi
    wave_mask = (w < 1214.)|(w > 1217.)&(w < 1301.)|(w > 1307.)
    x, y, time = x[wave_mask], y[wave_mask], time[wave_mask]

    #extract lightcurve from spectrum
    sp_mask = region_mask(x, y, slope, sp_intercept, sp_height)
    sp_lc = np.histogram(time[sp_mask], bins)
    t,sp_counts = sp_lc[1][:-1], sp_lc[0]

    #background
    bk1_mask = region_mask(x, y, slope, bk1_intercept, bk1_height)
    bk1_lc = np.histogram(time[bk1_mask], bins)
    bk2_mask = region_mask(x, y, slope, bk2_intercept, bk2_height)
    bk2_lc = np.histogram(time[bk2_mask], bins)
    bk_counts = (bk1_lc[0]+bk2_lc[0])*(sp_height/(bk1_height+bk2_height)) #sum background counts and normalise to spectrum area

    #background subtraction 
    counts_bksub = sp_counts - bk_counts 

    return t, counts_bksub


def lc_maker(file_path=os.getcwd()+'/', save_path=os.getcwd()+'/lightcurves/', bin_time=1., qual_check=True, plot=True, save_file='long'):

    #find the corrtag files, and end the script if there aren't any
    tag_files = glob.glob(file_path+'*corrtag*')
    if len(tag_files) == 0:
        print ('There are no corrtag files in file_path :(.')
        os._exit(1)

    #find all rootnames
    rootnames = np.array([], dtype=str)
    bands = np.array([], dtype=str)
    exp_times = np.array([], dtype=float)
    stars = np.array([], dtype=str)
    start_times = np.array([], dtype=float)
    for tag in tag_files:
        with fits.open(tag) as hdul:
            rootname = hdul[0].header['ROOTNAME']
            if rootname not in rootnames:
                rootnames = np.append(rootnames,rootname)
                bands = np.append(bands, hdul[0].header['DETECTOR'])
                exp_times = np.append(exp_times, hdul[1].header['EXPTIME'])
                stars = np.append(stars, hdul[0].header['TARGNAME'])
                start_times = np.append(start_times, hdul[1].header['EXPSTART'])
  
    for rootname, band, exp_time, star, start_time in zip(rootnames, bands, exp_times, stars, start_times):
        nbins = int(exp_time/bin_time)
        counts = np.zeros(nbins)
        
        if band == 'FUV':

            #checks if both segments are available 
            segs = ['a', 'b']
            if (file_path+rootname+'_corrtag_a.fits') not in tag_files:
                segs = ['b']
            if (file_path+rootname+'_corrtag_b.fits') not in tag_files:
                segs = ['a']

            for seg in segs:
                tag_file = rootname+'_corrtag_'+seg+'.fits'
                seg = seg.upper() #header keywords are uppercase
                
                hdul = fits.open(file_path+tag_file)
                header = hdul[1].header
                data = hdul[1].data
                hdul.close()
                            
                t, counts_bksub = extract_counts(header, data , nbins, qual_check, seg)                
                counts += counts_bksub
            
        
        elif band == 'NUV':

            tag_file = rootname+'_corrtag.fits'
            segs = ['A', 'B', 'C']
            hdul = fits.open(file_path+tag_file)
            header = hdul[1].header
            data = hdul[1].data
            hdul.close()
            
            for seg in segs:
                t, counts_bksub = extract_counts(header, data, nbins, qual_check, seg)
                counts += counts_bksub


        #calculate photon noise
        error = counts**0.5 

        #convert time to absolute time
        t_adg = t + (start_time*86400.)

        #convert to counts s-1 
        counts_sec = counts/bin_time
        error_sec = error/bin_time
        
        if save_file == 'long':            
            save_file = star+'_'+band+'_'+str(start_time)+'_'+rootname+'_'+str(bin_time)+'s_lc.ecsv'
        elif save_file =='short':
            save_file = rootname+'_'+str(bin_time)+'s_lc.ecsv'
        
        filewriter(t_adg, counts_sec, error_sec, save_path, save_file)
        if plot == True:
            plot_lc(star, rootname, t_adg, counts_sec, error_sec, bin_time)

        
        
        
def plot_lc(star, rootname, time, counts, error, bin_time):
	plt.figure(star+'_'+rootname+'_'+str(bin_time)+'s')
	plt.subplots_adjust(top=0.99, right =0.99)
	plt.errorbar(time-time[0], counts, yerr = error, ls='none', marker='o')
	plt.xlabel('Time (s)', size=20)
	plt.ylabel('Counts (s$^{-1}$)', size=20)
	plt.show()
	
#test
#lc_maker(file_path = '0845_lc_test/', save_path = '0845_lc_test/test/', bin_time=5.)