import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as fits
import glob
from astropy.table import Table
from astropy.io import ascii
from astropy.convolution import convolve, Box1DKernel
import astropy.units as u

x1 = []

def on_key(event):
    global x1
    if event.key == 'w':
        x1.append(event.xdata)
        print('%.3f' %event.xdata)
        plt.close()
    if event.key == 'e':
        x1.append(-1)
        print('%.3f' %event.xdata)
        plt.close()


"""
Plots each line in EG Uma in turn to identify the line position. If no lines, marks it with -1

"""
data = Table.read('wd_removed_eguma.ecsv')
w, f, e = data['WAVELENGTH'], data['FLUX'], data['ERROR']
#f = convolve(f,Box1DKernel(5))
#e = convolve(e,Box1DKernel(5))/(5**0.5)

linelist = Table.read('/home/david/work/muscles/FUV_linelist.csv')
mask = (linelist['Wavelength'] > w[0]) & (linelist['Wavelength'] < w[-1]) & (linelist['Likelihood to measure'] == 'High')
linelist = linelist[mask]

for row in linelist:
    line= row['Wavelength']
    mask = (w > line-2) & (w <line+2)
    fig =plt.figure(row['Ion'])
    plt.step(w[mask], f[mask], where='mid')
    plt.axvline(line, c='C1', ls='--')
    dv= -140*u.km/u.s
    line_shift = dv.to(u.AA, equivalencies=u.doppler_optical(line*u.AA))
    plt.axvline(line_shift.value, c='C3', ls='--')
    [plt.axvline(lines, c='C2', ls='--', zorder=-10, alpha=0.5) for lines in linelist['Wavelength']]
    plt.xlim(line-2, line+2)
    cid = fig.canvas.mpl_connect('key_press_event',on_key)
    plt.show()
    
savedat = Table([linelist['Wavelength'], x1], names=['FILENAME', 'Xs'])
ascii.write(savedat, 'lines_to_measure.ecsv', format='ecsv', overwrite=True)
    