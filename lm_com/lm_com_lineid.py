#simple script to plot an lm com spectrum and find the lines

import numpy as np
import matplotlib.pyplot as plt
import glob
from astropy.convolution import convolve, Box1DKernel
from astropy.table import Table

def boxcar(flux,factor):
    smoothed_flux=convolve(flux,Box1DKernel(factor))
    return smoothed_flux

ism = Table.read('../ism_lines.csv')

o_lines=[1152.150,1302.170,1304.860,1306.030]
si_lines=[1190.416,1193.292,1194.500,1197.394,1246.740,1248.426,
1250.091,1250.436,1251.164,1260.422,1264.738,
1305.592,1309.276,1346.884,
1348.543,1350.072,1350.516,1350.656,1352.635,1353.721,
1526.707,1533.431,3854.758,3857.112,3863.690,4129.219,
4132.059,5042.430,5057.394,6348.864,6373.132,
1140.546,1141.579,1142.285,1144.309,1144.959,1154.998,
1155.959,1156.782,1158.101,1160.252,1161.579,1206.500,
1206.555,1294.545,1296.726,1298.892,1301.149,1303.323,
1341.458,1342.389,1365.253,1417.237,1393.775,1402.770]
c_lines= [1334.530,1335.660,1335.708, 1174.930,1175.260,1175.590,1175.710,1175.987,1176.370]


spectra = glob.glob('*.dat')


plt.figure('lm_com')

w, f, e, dq = np.loadtxt(spectra[0], unpack=True)
w, f = w[dq==0], f[dq==0]
f = boxcar(f,5)
ly_mask = (w < 1214.5)|(w > 1217.)
plt.plot(w[ly_mask],f[ly_mask])
#ism
[plt.axvline(line,ls='--', c='k') for line in ism['rest_lambda']]
#o
[plt.axvline(line,ls='--', c='g') for line in o_lines]
#si
[plt.axvline(line,ls='--', c='r') for line in si_lines]
#c
[plt.axvline(line,ls='--', c='b') for line in c_lines]

plt.xlim(w[0],w[-1])



plt.show()