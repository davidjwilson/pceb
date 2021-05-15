import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
from astropy.convolution import convolve, Box1DKernel
import glob

cpath = '/media/david/5tb_storage1/pceb_data/eg_uma/'
x = glob.glob('{}*x1dsum.fits'.format(cpath))[0]

smooth=5
data = fits.getdata(x, 1)
wc = np.array([], dtype=float)
fc = np.array([], dtype=float)
ec = np.array([], dtype=float)
for dt in data[::-1]:
    wi, fi, ei, dq = dt['WAVELENGTH'], dt['FLUX'], dt['ERROR'], dt['DQ']
   # mask = (fi>0) & (dq == 0) & (wi < 1213) | (wi > 1217) & (fi>0) & (dq == 0) 
    mask = (dq==0)
    wi, fi, ei = wi[mask], fi[mask], ei[mask]
    wc = np.concatenate((wc, wi))
    fc = np.concatenate((fc, fi))
    ec = np.concatenate((ec, ei))    

fc = convolve(fc,Box1DKernel(smooth))
ec = convolve(ec,Box1DKernel(smooth))/(smooth**0.5)

# plt.plot(wc, fc)

mask = (wc < 1210) | (wc > 1220) 
wc1, fc1, ec1 = wc[mask], fc[mask], ec[mask]

plt.figure()
plt.plot(wc1, fc1)

linestab = np.genfromtxt('linelist.csv', dtype=None, names=True, delimiter=',',encoding=None)

lines = linestab[linestab['ISM'] =='n']
ism = linestab[linestab['ISM'] =='y']
# print(ism)
[plt.axvline(line, ls='--', c='C1') for line in lines['wavelength']]
[plt.axvline(line, ls='--', c='C2') for line in ism['wavelength']]


plt.show()