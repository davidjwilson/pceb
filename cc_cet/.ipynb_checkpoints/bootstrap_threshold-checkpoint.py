'''
Script written to estimate a 0.1% false alarm probability significance threshold
via the bootsrapt method discussed in VanderPlas 2018, Bell et al. 2019

Written by Joseph Guidry

Last updated: Feb 10, 2020
'''

import sys
import numpy as np
from astropy.stats import bootstrap
from astropy.timeseries import LombScargle as ls
import astropy.coordinates as coord
from astropy.time import Time
import warnings
warnings.filterwarnings('ignore')
import os 
from astropy.utils import NumpyRNGContext
import multiprocessing as mp
from joblib import Parallel, delayed
import astropy.units as u


#############################################################
##
##  Progress Bar Code. Code from Stack Overflow:
##  https://stackoverflow.com/questions/3002085/python-to-print-out-status-bar-and-percentage
##
#############################################################

## Provide the interation counter (count=int)
## and the action being performed (action=string)
def progress_bar(count,total,action):
    sys.stdout.write('\r')
    sys.stdout.write(action)
    sys.stdout.write("[%-20s] %d%%  %d/%d" % ('='*int((count*20/total)),\
                                            count*100/total,\
                                            count,total))
    sys.stdout.flush()
    return

# Function to Calculate Lomb-Scargle Periodogram (LSP) with Amplitude Units
def calc_lsp(time,flux,freqs):
    lsp = ls(time,flux).power(freqs,normalization='psd')
    norm_lsp = np.sqrt(abs(4.0*(lsp/len(time))))
    return norm_lsp

# Define a separate max function if parallelizing
def calc_lsp_max(time,flux,freqs,idx,bnum,act):
    lsp = ls(time,flux).power(freqs,normalization='psd')
    norm_lsp = np.sqrt(abs(4.0*(lsp/len(time))))
    progress_bar(idx+1,bnum,act)
    return np.max(norm_lsp)

# Define punction to perform bootstrapping and calculate the False Alarm Probability significance threshold
def calc_threshold(time,flux,bootnum=10000,percentile=99.9,parallel=False):
    # Define frequency array to run the LSP on
    # Oversample the baseline by a factor of 2
    fnyq = (1/(20*u.s.to(u.day)))
    fres = 1. / (2 * (np.max(time) - np.min(time)))
    farr = np.linspace(1e-6,fnyq,int(fnyq/fres))

    # Bootstrap the time-series 10,000 times with replacement
    # and specify the random seed for replicability
    with NumpyRNGContext(1):
        bootstrapped_lcs = bootstrap(flux,bootnum=bootnum)

    # Initialize progress bar
    action = 'Performing Bootstrapping...' # Progress bar message
    progress_bar(0,bootnum,action)

    # If the operation is to be parallelized:
    if parallel:
    # Calculate the maximum value for the periodogram of each bootstrapped light curve
        max_vals = Parallel(n_jobs=4)(delayed(calc_lsp_max)(time,bootstrapped_lcs[i],farr,i,bootnum,action) for i in range(bootnum))

    else:
        # Define a list to append the maximum values into
        max_vals = []
        # Calculate the maximum value for the periodogram of each bootstrapped light curve
        for i in range(bootnum):
            max_vals.append(np.max(calc_lsp(time,bootstrapped_lcs[i],farr)))
            progress_bar(i+1,bootnum,action)

    # Estimate the false alarm probability
    fap = np.percentile(max_vals,percentile)
    # Compute the mean amplitude of the original periodogram
    og_mean_amp = np.mean(calc_lsp(time,flux,farr))

    print('\n')
    print('--------------------------------------------------------------------------')
    print("The 0.1 % False Alarm Probability threshold: {:.4f} %".format(fap*1e2))
    print("This is equal to {:.4f} times the original peridogram's amplitude".format(fap/og_mean_amp))

    # Make python talk to you to let you know the script is finished
    # os.system("say 'Your bootstrapping routine is finished running.'")

    # Return the 0.1% FAP in both percent and units of the original periodogram's mean amplitude
    return fap*1e2, fap/og_mean_amp

######### Example snippet code to run on a light curves in the Google Drive tarred folder ##########
# import pandas as pd
from astropy.table import Table
lcs = ['lightcurves/cc__cos_ldlc01_10s_fppos_norm.ecsv', 'lightcurves/cc__cos_ldlc51_10s_fppos_norm.ecsv']
for lc in lcs:
    print(lc)
    data = Table.read(lc)
    obj_name = 'CC_Cet'
    time, flux = data['JD'], data['FLUX']-1
    # import matplotlib.pyplot as plt

    # plt.scatter(time,flux)
    # plt.show()
    # lc_fname = 'WDJ171251p78-191550p28_lcdata.csv'
    # obj_name = 'WD' + ' ' + lc_fname[2:7] + lc_fname[12:17]
    # lc = pd.read_csv('../bootstrap_ztf_lcs/' + lc_fname)
    # time = lc.bjd_sec.values
    # flux = lc.flux.values

    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
    print("%%%%%%%%     Commencing bootstrapping for", obj_name,"    %%%%%%%%")
    # Estimate the false alarm probability lever in units of percent and <A>
    fap_percent, fap_amp = calc_threshold(time,flux,bootnum=25,percentile=99.9,parallel=True)


