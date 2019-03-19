import matplotlib.pyplot as plt
import numpy as np
import os
from astropy import constants as const
from astropy import units as u
from astropy.analytic_functions import blackbody_lambda
import phot

#makes the "flaring" spectrum for the PCEBs, ie wd + 9000K blackbody

#uz sex first as there is a deadline

models =os.listdir('models_adj')
d_U=np.array([-8.,-8,-2.8,-8,-8,-8,-2.8,-8,-8,-8,-8], dtype=float) #change in u mag if flaring


tab=np.genfromtxt('pceb_bop.csv', names=True, delimiter=',', dtype=None)
#tab = tab[tab['Name'] not in ['V727 Car', 'WD 0137-349']]
for k in xrange(len(tab)):
        t = tab[k]
        if t['Name'] not in ['WD 0137-349']:
	  plt.figure(t['Name'], figsize=(8,6))
	  plt.subplots_adjust(top=0.99, bottom=0.1, left=0.1, right=0.99)
	  plt.xscale('log')
	  plt.yscale('log')
	  for m in models:
	    if m[0:5] == str(t['wd_teff']):
	      mod=m
	  mw, mf = np.loadtxt('models_adj/'+mod, unpack=True)

	  plt.plot(mw, mf)


	# blackbody
	  #w_bb=np.linspace(mw[0], mw[-1], 1000)
	  f_bb=blackbody_lambda(mw, 9000)
	#  plt.plot(w_bb, f_bb*1.E-14, 'r')
	  #B_abs = 11.048 #from Bergeron's tables for logg = 8, t =17000
	  #d = 38. #distance
	  #B_ap = B_abs +(5*np.log10(d/10.))
	  #V_mdwarf= 12.8 +(5*np.log10(d/10.))
	  #print V_mdwarf
	  #U_mag = V_mdwarf + 1.661  + 1.222
	  #print U_mag
	  #print B_ap
	  #mod_B = np.median(mf[(mw >4450.-20) & (mw < 3543+20.)
	  #mod_B = 
	  #U_mag=t['B'] +1.22 
	  f1=phot.return_fluxes(t['U_calc']+d_U[k], 0, 'u', 'sdss')[0]
	  #print f1
	  bb_slice  = (mw >3543.-20) & (mw < 3543.+20.)
	  f_bb = f_bb/(np.median(f_bb[bb_slice])/f1)
	  plt.plot(mw, f_bb, 'r')
	  
	  f_comb = mf + f_bb
	  plt.plot(mw, f_comb, 'g')
	  
	  name = t['Name'].replace(' ','_')
	  fl=open('combined_models/'+name+'_dU'+str(d_U[0])+'.dat', 'w')
	  for j in xrange(len(mw)):
	    fl.write((str(mw[j])+'   '+str(f_comb[j])+'\n'))
	  #plt.plot(mw, f_comb, 'c')
	
plt.show()