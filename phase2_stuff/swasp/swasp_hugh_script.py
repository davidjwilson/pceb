from astropy.io import fits

# Hugh's swasp script

def readfile(fname):
	if str(fname.split('/')[-1]).find('WASP')!=-1:
			a=fits.open(fname)
			if mag==0:
				fluxmed=(np.median(a[1].data.TAMFLUX2[:][a[1].data.TAMFLUX2[:]>0]))
				#Getting flux and repositioning at 1
				flux=(a[1].data.TAMFLUX2[:][a[1].data.TAMFLUX2[:]>0])/fluxmed
				flux_errs=(a[1].data.TAMFLUX2_ERR[:][a[1].data.TAMFLUX2[:]>0])/(fluxmed)*flux
				time=(a[1].data.TMID[a[1].data.TAMFLUX2[:]>0]).astype(float)
				Lc=np.column_stack(((time/86400.0)+3005.5, flux, flux_errs))
				#Sorting by TMID
			else:
				#Leaving in magnitudes...
				mags=WaspMag(a[1].data.TAMFLUX2[:][a[1].data.TAMFLUX2[:]>0])
				magerr=WaspMagErr(a[1].data.TAMFLUX2_ERR[:][a[1].data.TAMFLUX2[:]>0], a[1].data.TAMFLUX2[:][a[1].data.TAMFLUX2[:]>0])
				time=(a[1].data.TMID[a[1].data.TAMFLUX2[:]>0]).astype(float)
				Lc=np.column_stack(((time/86400.)+3, mags, magerr))
			Lc=Lc[np.argsort(Lc[:, 0])]
			return Lc
	
lc = readfile('1SWASPJ102834.89-000029.1.fits')
print lc
