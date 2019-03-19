import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gs
#plots pft and crts lighcurves of ux sex


gs1=gs.GridSpec(1,3)
#gs2=gs.GridSpec(3,2)

star='UZ Sex'

plt.figure(star, figsize=(10,6))
plt.subplots_adjust(top=0.99, bottom=0.1, left=0.1, right=0.99, wspace=0.05)
col='#d62728'
col2='#1f77bf'


#plt.subplot(gs1[0])
plt.minorticks_on()
plt.ylabel('Normalised Flux', size=20)
t, f, e = np.loadtxt('ptf/UZSex_ptf.dat', unpack=True)
t= t + (2454904.6748572-2400000.5)
plt.errorbar(t, f/np.median(f), yerr = e/np.median(f), color=col, label='PTF', marker='o', mec=col, ls ='none', capsize=0)
plt.legend()


#plt.subplot(gs1[0,-2:])
#plt.yticks(visible=False)
plt.minorticks_on()
plt.xlabel('Date (MJD)', size=20)
t, f, e = np.loadtxt('crts/uzsex_crts.csv', unpack=True, usecols=(5, 1,2), delimiter=',', skiprows=1)
plt.errorbar(t, f/np.median(f), yerr = e/np.median(f), color=col2, label='CRTS', marker='o', mec=col2, ls ='none', capsize=0)
#plt.legend(numpoints=0)
#plt.annotate('CRTS', (0.9,0.9), xycoords='figure fraction', size=20) 
plt.legend()
plt.xlim(53400,56490)
plt.ylim(0.971, 1.029)

plt.show()