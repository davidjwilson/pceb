import matplotlib.pyplot as plt
import numpy as np
import phot as phot
import os

#plots photometry of mdwarfs


band_names = np.array(['fuv', 'nuv', 'u', 'g', 'r', 'i', 'z', 'B','V', 'j', 'h', 'k'], dtype=str)
survey = np.array(['galex', 'galex', 'sdss', 'sdss', 'sdss', 'sdss', 'sdss', 'Johnson', 'Johnson', '2mass', '2mass', '2mass'], dtype=str)
bands = np.array([1516., 2267.,3543., 4770., 6231., 7625., 9134., 4450., 5510.,12200., 16300., 21900.], dtype=float)

distance=[110,38,238,37.3,102,50,36,191,133,214,136,243]

sp_type=['M4', 'M3', 'M3', 'M4', 'L8', 'M3', 'M4', 'M3', 'M3', 'M5', 'M3', 'M3']

models =os.listdir('models')
mdwarf_models = os.listdir('sdss_mdwarf')

plt.figure('pcebs', figsize=(16,12))
plt.subplots_adjust(top=0.99, bottom=0.01, left=0.01, right=0.99)
tab=np.genfromtxt('pceb_bop.csv', names=True, delimiter=',', dtype=None)
for k in xrange(len(tab)):
        t = tab[k]
	plt.subplot(4,3, k+1)
	
	#iue
	iue_files = os.listdir('iue/dat_files')
	name = t['Name'].replace(' ','_')
	if name+'_sw_iue.dat' in iue_files:
	  w, f = np.loadtxt('iue/dat_files/'+name+'_sw_iue.dat', unpack=True)
	  plt.plot(w,f,'g')
	if name+'_lw_iue.dat' in iue_files:
	  w, f = np.loadtxt('iue/dat_files/'+name+'_lw_iue.dat', unpack=True)
	  plt.plot(w,f,'g')
	
	#photometry
	w=[]
	f=[]
#	if np.isnan(t[
	for i in xrange(len(band_names)):
		v=t[band_names[i]]
		if np.isnan(v)==False:
		#	w.append(bands[i])
			f1=phot.return_fluxes(t[band_names[i]], 0, band_names[i], survey[i])
			f.append(f1[0])
			w.append(f1[2])
	#print w, f
	plt.scatter(w,f, color='r')
	yup=max(f)+0.1*max(f)
	ydn=min(f)-0.1*min(f)
	plt.ylim(ydn, yup)
	plt.xscale('log')
	plt.yscale('log')
	for m in models:
	  if m[0:5] == str(t['wd_teff']):
	    mod=m
	mw, mf = np.loadtxt('models/'+mod, unpack=True)
	mf *= (10./distance[k])**2.
	if np.isnan(t['fuv'])==False:
	  gal=[]
	  for i in xrange(2):
	    f_slice = (mw > bands[i]-10.) & (mw < bands[i]+10.)
	    gal.append(np.median(mf[f_slice])/f[i])
	  scale=np.mean(gal)
	if name+'_sw_iue.dat' in iue_files: #UZ Sex galax data is poor, scale to iue instead
	  m_slice= (mw > bands[0]-10.) & (mw < bands[0]+10.)
	  iue_w, iue_f = np.loadtxt('iue/dat_files/'+name+'_sw_iue.dat', unpack=True)
	  iue_slice = (iue_w > bands[0]-10.) & (iue_w < bands[0]+10.)
	  scale = np.median(mf[m_slice])/np.median(iue_f[iue_slice])
	mf /= scale
	
	   

	
	fl=open('models_adj/'+mod[:-4]+'_adj.dat', 'w')
	for j in xrange(len(mw)):
	  if mw[j] != mw[j-1]:
            fl.write((str(mw[j])+'   '+str(mf[j])+'\n'))
	"""
	for m in mdwarf_models:
	  if m[5:7] == sp_type[k]:
	    mod2=m

	mdw, mdf = np.loadtxt('sdss_mdwarf/'+mod2, unpack=True)
	mdf *= (0.01/distance[k])**2.
	"""
	for i in xrange(len(band_names)):
	  plt.annotate(band_names[i], (bands[i],min(f))) 
	plt.annotate(t['Name'], (w[-3],f[1]))
	plt.xlim(1000., max(w)+0.1*max(w))
	plt.plot(mw, mf)
	#plt.plot(mdw, mdf, 'r')
	if t['Name'] == 'LM com': #test, plots the sdss spectrum
	  dsw, dsf = np.loadtxt('mwdd_opt_spectra/spect.V_LMCom.sdss.txt', unpack=True, skiprows=2, delimiter=',')
	  plt.plot(dsw, dsf*1.7E11, 'r')
	if t['Name'] == 'WD 1339+606': #test, plots the sdss spectrum
	  dsw, dsf = np.loadtxt('mwdd_opt_spectra/spect.WD1339+606.sdss.txt', unpack=True, skiprows=2, delimiter=',')
	  plt.plot(dsw, dsf*1.5E11, 'r')
	
	
	
	
plt.show()
			
		
 
