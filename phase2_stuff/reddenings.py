"""A module containing Reddening transformations
The standard reddening law is the Rv=3.1 Fitzpatrick (1999)

	cardelli -- Rv=3.1 Cardelli et al. (1989) reddening law
	deredden -- Rv=x Fitzpatrick (1999) reddening law
	deredden_fix -- Rv=3.1 Fitzpatrick (1999) reddening law

This module, automatically loads VPHAS+ and SDSS passbands into
astSED objects
"""


from numpy import *
from scipy import interpolate
from astLib import astSED
import urllib2
import pyfits
from astLib import astWCS
from astLib import astCoords
from subprocess import *
import string
class passbands:
        """This class loads all the passbands of several telescope/surveys:
	The original files can be found here:
	http://svo2.cab.inta-csic.es/theory/fps3/index.php?mode=browse
	
	They are stored on disk here:
	/storage/astro2/phsmav/data2/filters/"""
	
	passbands_home = '/storage/astro2/phsmav/data2/filters/'
	# loading GALEX passbands #
	fuv_band_galex = astSED.Passband(passbands_home+'galex/GALEX_GALEX.FUV.dat')
	nuv_band_galex = astSED.Passband(passbands_home+'galex/GALEX_GALEX.NUV.dat')
	# loading iphas passbands #
	r_band_iphas = astSED.Passband(passbands_home+'iphas/r.dat')
	i_band_iphas = astSED.Passband(passbands_home+'iphas/i.dat')
	h_band_iphas = astSED.Passband(passbands_home+'iphas/ha_wfc.dat')	
	# loading vphas passbands #
	u_band_vphas = astSED.Passband(passbands_home+'vphas/u_modified_atm_tel.txt')
	g_band_vphas = astSED.Passband(passbands_home+'vphas/g_modified_atm_tel.txt')
	r_band_vphas = astSED.Passband(passbands_home+'vphas/r_modified_atm_tel.txt')
	i_band_vphas = astSED.Passband(passbands_home+'vphas/i_modified_atm_tel.txt')
	h_band_vphas = astSED.Passband(passbands_home+'vphas/h_modified_atm_tel.txt')
	# loading atlas passbands #
	u_band_atlas = astSED.Passband(passbands_home+'atlas/u_modified_atm_tel.txt')
	g_band_atlas = astSED.Passband(passbands_home+'atlas/g_modified_atm_tel.txt')
	r_band_atlas = astSED.Passband(passbands_home+'atlas/r_modified_atm_tel.txt')
	i_band_atlas = astSED.Passband(passbands_home+'atlas/i_modified_atm_tel.txt')
	z_band_atlas = astSED.Passband(passbands_home+'atlas/z_SDSS.txt')      
	# loading SDSS passbands #
	u_band_sdss = astSED.Passband(passbands_home+'sdss/u.dat')
	g_band_sdss = astSED.Passband(passbands_home+'sdss/g.dat')
	r_band_sdss = astSED.Passband(passbands_home+'sdss/r.dat')
	i_band_sdss = astSED.Passband(passbands_home+'sdss/i.dat')
	z_band_sdss = astSED.Passband(passbands_home+'sdss/z.dat')
	# loading APASS filters #
	B_band_apass = astSED.Passband(passbands_home+'apass/APASS_B.dat')        
	V_band_apass = astSED.Passband(passbands_home+'apass/APASS_V.dat')
	g_band_apass = astSED.Passband(passbands_home+'apass/APASS-sdss_g.dat')    
	r_band_apass = astSED.Passband(passbands_home+'apass/APASS-sdss_r.dat')
	i_band_apass = astSED.Passband(passbands_home+'apass/APASS-sdss_i.dat')	
	# loading Johnson passbands #
	B_band_johnson = astSED.Passband(passbands_home+'johnson/johnson_B.dat')	
	V_band_johnson = astSED.Passband(passbands_home+'johnson/johnson_V.dat')
	# loading Bessell passbands #
	B_band_bessell = astSED.Passband(passbands_home+'bessell/bessell_B.dat') 
	V_band_bessell = astSED.Passband(passbands_home+'bessell/bessell_V.dat')	  
	# loading CTIO Johnson-Kron-Cousins #
	B_band_ctio = astSED.Passband(passbands_home+'ctio/johnson_B_4202-1050.dat')	     
	V_band_ctio = astSED.Passband(passbands_home+'ctio/johnson_V_5475-1000.dat')
	R_band_ctio = astSED.Passband(passbands_home+'ctio/cousins_R_6400-1450.dat')	     
	I_band_ctio = astSED.Passband(passbands_home+'ctio/cousins_I_8118-1415.dat')	
	# loading APASS filters #
	u_band_OAJ = astSED.Passband(passbands_home+'OAJ/JPLUS_01.dat')	   
	g_band_OAJ = astSED.Passband(passbands_home+'OAJ/JPLUS_06.dat')
	r_band_OAJ = astSED.Passband(passbands_home+'OAJ/JPLUS_08.dat')    
	h_band_OAJ = astSED.Passband(passbands_home+'OAJ/JPLUS_09.dat')
	i_band_OAJ = astSED.Passband(passbands_home+'OAJ/JPLUS_10.dat')   
	z_band_OAJ = astSED.Passband(passbands_home+'OAJ/JPLUS_12.dat')
	# loading Gaia filters #
	Gp_band_gaia = astSED.Passband(passbands_home+'gaia/Gp.dat')       
	Bp_band_gaia = astSED.Passband(passbands_home+'gaia/Bp.dat')
	Rp_band_gaia = astSED.Passband(passbands_home+'gaia/Rp.dat') 
	# loading Hipparcos/Tycho-2 filters #
	Hp_band_tycho = astSED.Passband(passbands_home+'tycho/Hipparcos_Hipparcos.Hp_MvB.dat') 
	Bt_band_tycho = astSED.Passband(passbands_home+'tycho/TYCHO_TYCHO.B_MvB.dat')
	Vt_band_tycho = astSED.Passband(passbands_home+'tycho/TYCHO_TYCHO.V_MvB.dat')	  
	# loading 2MASS passbands #
	j_band_2mass = astSED.Passband(passbands_home+'2mass/J.dat')
	h_band_2mass = astSED.Passband(passbands_home+'2mass/H.dat')
	k_band_2mass = astSED.Passband(passbands_home+'2mass/Ks.dat')
	# loading UKIDSS filters #
	z_band_ukidss = astSED.Passband(passbands_home+'ukidss/UKIRT_UKIDSS.Z.dat')	     
	y_band_ukidss = astSED.Passband(passbands_home+'ukidss/UKIRT_UKIDSS.Y.dat')
	j_band_ukidss = astSED.Passband(passbands_home+'ukidss/UKIRT_UKIDSS.J.dat')    
	h_band_ukidss = astSED.Passband(passbands_home+'ukidss/UKIRT_UKIDSS.H.dat')
	k_band_ukidss = astSED.Passband(passbands_home+'ukidss/UKIRT_UKIDSS.K.dat')	
	# loading VISTA passbands #
	z_band_vista = astSED.Passband(passbands_home+'vista/Paranal_VISTA.Z.dat')
	y_band_vista = astSED.Passband(passbands_home+'vista/Paranal_VISTA.Y.dat')
	j_band_vista = astSED.Passband(passbands_home+'vista/Paranal_VISTA.J.dat')
	h_band_vista = astSED.Passband(passbands_home+'vista/Paranal_VISTA.H.dat')
	k_band_vista = astSED.Passband(passbands_home+'vista/Paranal_VISTA.Ks.dat')   
	# loading WISE passbands #
	w1_band_wise = astSED.Passband(passbands_home+'wise/W1.dat')
	w2_band_wise = astSED.Passband(passbands_home+'wise/W2.dat')
	w3_band_wise = astSED.Passband(passbands_home+'wise/W3.dat')
	w4_band_wise = astSED.Passband(passbands_home+'wise/W4.dat')
	# loading Spitzer passbands #
	i1_band_spitzer = astSED.Passband(passbands_home+'spitzer/Spitzer_IRAC.I1.dat')
	i2_band_spitzer = astSED.Passband(passbands_home+'spitzer/Spitzer_IRAC.I2.dat')
	i3_band_spitzer = astSED.Passband(passbands_home+'spitzer/Spitzer_IRAC.I3.dat')
	i4_band_spitzer = astSED.Passband(passbands_home+'spitzer/Spitzer_IRAC.I4.dat')
	m1_band_spitzer = astSED.Passband(passbands_home+'spitzer/Spitzer_MIPS.24mu.dat')	
	# loading Pan-STARRS passbands #
	g_band_ps1 = astSED.Passband(passbands_home+'panstarrs/PS1.g.dat')
	r_band_ps1 = astSED.Passband(passbands_home+'panstarrs/PS1.r.dat')
	i_band_ps1 = astSED.Passband(passbands_home+'panstarrs/PS1.i.dat')
	z_band_ps1 = astSED.Passband(passbands_home+'panstarrs/PS1.z.dat')
	y_band_ps1 = astSED.Passband(passbands_home+'panstarrs/PS1.y.dat')      





def cardelli(spectrum, color_exc):
	"""This function dereddens OPTICAL/NIR (1.1 micron^-1 < x < 3.3 micron^-1) 
	spectra, by using the relationship given in Cardelli, Clayton and Mathis, 
	ApJ 345:245, 1989.
	Output: Corrected spectrum, Extinction
	"""
	Rv= 3.1
	Av= Rv * color_exc
	x = (10000./spectrum[:,0])
	y= (x-1.82)
	a0,a1,a2,a3,a4,a5,a6,a7 = 1,0.17699,0.50447,0.02427,0.72085,0.01979,0.77530,0.32999
	b1,b2,b3,b4,b5,b6,b7 = 1.41338,2.28305,1.07233,5.38434,0.62251,5.30260,2.09002
	a_x = a0 + a1*y - a2*(y**2) - a3*(y**3) + a4*(y**4) + a5*(y**5) - a6*(y**6) + a7*(y**7)
	b_x = b1*y + b2*(y**2) + b3*(y**3) - b4*(y**4) - b5*(y**5) + b6*(y**6) - b7*(y**7)
	A_lambda = (a_x + b_x/Rv)*Av
	deredden_flux = spectrum[:,1]*power(10., 0.4*A_lambda)
	A = asarray((spectrum[:,0], deredden_flux))
	S= swapaxes(A,0,1)
	return S, str(Av)

def deredden(spectrum, ebv, tot_to_sel=3.1):
       """This function dereddens OPTICAL/NIR (1.1 micron^-1 < x < 3.3 micron^-1) 
       spectra, by using the relationship given in Fitzpatrick, 1999, PASP, 111, 63.
       R dependent curve"""
       R= tot_to_sel
       c2 = -0.824 + 4.717 / R
       c1 = 2.03 - 3.007*c2
       x0 = 4.596
       gamma = 0.99
       c3 = 3.23
       c4 = 0.41
       
       Rv= 3.1
       Av= R * ebv
       x = (10000./spectrum[:,0])      
       anchors=[0.000,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846]
       curve=[0.000, 0, 0, 0, 0, 0 , 0, 0, 0]
       curve[1] = 0.265*R/Rv
       curve[2] = 0.829*R/Rv
       curve[3] = -0.426+1.0044*R
       curve[4] = -0.050 + 1.0016*R
       curve[5] = 0.701 + 1.0016*R
       curve[6] = 1.208 + 1.0032*R -0.00033*R**2
       curve[7] = c1 + c2*anchors[7] + c3*pow(anchors[7],2)/(pow(anchors[7]**2 - x0**2, 2) + pow(gamma*anchors[7],2)) + R
       curve[8] = c1 + c2*anchors[8] + c3*pow(anchors[8],2)/(pow(anchors[8]**2 - x0**2, 2) + pow(gamma*anchors[8],2)) + R 

       spl = interpolate.splrep(anchors,curve, k=1)
       A_lambda= interpolate.splev(x,spl)*ebv
       deredden_flux = spectrum[:,1]*power(10., 0.4*A_lambda)
       A = array((spectrum[:,0], deredden_flux))
       S= swapaxes(A,0,1)
       return S, str(Av)

def deredden_fix(spectrum, ebv):
       """This function dereddens OPTICAL/NIR (1.1 micron^-1 < x < 3.3 micron^-1) 
       spectra, by using the relationship given in Fitzpatrick, 1999, PASP, 111, 63.
       """
       Rv= 3.1
       Av= Rv * ebv
       x = (10000./spectrum[:,0])      
       anchors=[0.000,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846]
       curve=[0.000,0.265,0.829,2.688,3.055,3.806,4.315,6.265,6.591]
       spl = interpolate.splrep(anchors,curve,k=1)
       A_lambda= interpolate.splev(x,spl)*ebv
       deredden_flux = spectrum[:,1]*power(10., 0.4*A_lambda)
       A = asarray((spectrum[:,0], deredden_flux))
       S= swapaxes(A,0,1)
       return S, str(Av)

def extinction_law(wavelengths, law="Fitzpatrick", tot_to_sel=3.1, ebv=1, out="A_lambda"):
	"""Plots the extinction law for an input reddening law, as:
	A_lambda = Rv - E(lambda-V)/E(B-V).
	
	Input:
	Wavelength: array like wavelength range
	
	if law=Fitzpatrick:
		uses the Fitzpatrick 1999 law
	elif law=Cardelli:
		uses the Cardelli et al 1989 law
	tot_to_sel=3.1: float indicating the total to selective ratio
		it is fixed to 3.1 for law=Cardelli
	
	Outputs:
	if out="A_lambda":
		it returns an array A_lambda/E(B-V) vs 10E+5/lambda (micron^-1)
	elif out="E_lambda":
		it returns an array E(lambda-V)/E(B-V)  vs 10E+5/lambda (micron^-1)
				
	"""
	X = wavelengths.copy()
	
	if law=="Fitzpatrick":
		R= tot_to_sel
		c2 = -0.824 + 4.717 / R
		c1 = 2.03 - 3.007*c2
		x0 = 4.596
		gamma = 0.99
		c3 = 3.23
		c4 = 0.41
	
		Rv= 3.1
		Av= R * ebv
		x = (10000./X)      
		anchors=[0.000,0.377,0.820,1.667,1.828,2.141,2.433,3.704,3.846]
		curve=[0.000, 0, 0, 0, 0, 0 , 0, 0, 0]
		curve[1] = 0.265*R/Rv
		curve[2] = 0.829*R/Rv
		curve[3] = -0.426+1.0044*R
		curve[4] = -0.050 + 1.0016*R
		curve[5] = 0.701 + 1.0016*R
		curve[6] = 1.208 + 1.0032*R -0.00033*R**2
		curve[7] = c1 + c2*anchors[7] + c3*pow(anchors[7],2)/(pow(anchors[7]**2 - x0**2, 2) + pow(gamma*anchors[7],2)) + R
		curve[8] = c1 + c2*anchors[8] + c3*pow(anchors[8],2)/(pow(anchors[8]**2 - x0**2, 2) + pow(gamma*anchors[8],2)) + R 

		spl = interpolate.splrep(anchors,curve, k=1)
		if ebv!=0:
			A_lambda = (interpolate.splev(x,spl)*ebv)
		elif ebv==0:
			A_lambda = (interpolate.splev(x,spl)*ebv)		
		if out=="A_lambda":
			extinction_curve = swapaxes(vstack((x,A_lambda)), 0, 1)
			return extinction_curve
		elif out=="E_lambda":
			E_lambda = A_lambda - R
			extinction_curve = swapaxes(vstack((x,E_lambda)), 0, 1)
			return extinction_curve	
	elif law=="Cardelli":
		Rv= 3.1
		Av= Rv * ebv
		x = (10000./X)
		y= (x-1.82)
		a0,a1,a2,a3,a4,a5,a6,a7 = 1,0.17699,0.50447,0.02427,0.72085,0.01979,0.77530,0.32999
		b1,b2,b3,b4,b5,b6,b7 = 1.41338,2.28305,1.07233,5.38434,0.62251,5.30260,2.09002
		a_x = a0 + a1*y - a2*(y**2) - a3*(y**3) + a4*(y**4) + a5*(y**5) - a6*(y**6) + a7*(y**7)
		b_x = b1*y + b2*(y**2) + b3*(y**3) - b4*(y**4) - b5*(y**5) + b6*(y**6) - b7*(y**7)
		A_lambda = (a_x + b_x/Rv)*Av
		if out=="A_lambda":
			extinction_curve = swapaxes(vstack((x,A_lambda)), 0, 1)
			return extinction_curve
		elif out=="E_lambda":
			E_lambda = A_lambda/ebv - R
			extinction_curve = swapaxes(vstack((x,E_lambda)), 0, 1)
			return extinction_curve
					
		
def compute_extinction(passband, survey, law="Fitzpatrick", tot_to_sel=3.1, ebv=1.):
	""" Computes extinction in the passbands of several surveys.
	Allowed passbands are: 
	iphas: r,i,ha
	vphas: u,g,r,i,ha
	sdss: u,g,r,i,z
	"""
	if survey!="iphas" and survey!="vphas" and survey!="sdss" and survey!="2mass"  and survey!="wise":
		print 'not yet implemented'
	else:
		passbands()
		X = arange(0, 120000,0.1)
		A_lambda = extinction_law(X, law=law, tot_to_sel=tot_to_sel, ebv=ebv, out="A_lambda")
		A_lambda = astSED.SED(10000/A_lambda[:,0], A_lambda[:,1])	
		if survey=="vphas":
			if passband=="u":
				passband=passbands.u_band_vphas
				A_mag = A_lambda.calcFlux(passband)*ebv
				return A_mag
			elif passband=="g":
				passband=passbands.g_band_vphas		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="r":
				passband=passbands.r_band_vphas		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="i":
				passband=passbands.i_band_vphas		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="ha":
				passband=passbands.h_band_vphas	  		 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
		elif survey=="sdss":
			if passband=="u":
				passband=passbands.u_band_sdss
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag
			elif passband=="g":
				passband=passbands.g_band_sdss		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="r":
				passband=passbands.r_band_sdss		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="i":
				passband=passbands.i_band_sdss		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="z":
				passband=passbands.z_band_sdss		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag
		if survey=="iphas":	
			if passband=="r":
				passband=passbands.r_band_iphas		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="i":
				passband=passbands.i_band_iphas		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="ha":
				passband=passbands.h_band_iphas	  		 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag
		if survey=="2mass":	
			if passband=="j":
				passband=passbands.j_band_2mass		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="h":
				passband=passbands.h_band_2mass		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="k":
				passband=passbands.k_band_2mass	  		 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag
		elif survey=="wise":
			if passband=="w1":
				passband=passbands.w1_band_wise
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag
			elif passband=="w2":
				passband=passbands.w2_band_wise		  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="w3":
				passband=passbands.w3_band_wise	  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag		
			elif passband=="w4":
				passband=passbands.w4_band_wise	  	 
				A_mag = A_lambda.calcFlux(passband)*ebv   	  
				return A_mag	







def sfd_web(x,y,system='gal',radius=2.0):
    """On-line interrogator of Schlegel, Finkbeiner & Davis (1998) dust maps.
    Inputs are coordinates in the specified system.
    system: 'gal', 'equ j2000'
    radius: 2, to 10 deg"""
    
    if system == 'equ j2000':
       system = 'equ+j2000'
    
    url = 'http://irsa.ipac.caltech.edu/cgi-bin/DUST/nph-dust?locstr=%.2f+%.2f+%s&regSize=%.2f' % (x,y,system,radius)
    #response = urllib2.urlopen(url)
    html = urllib2.urlopen(url).readlines()
    #h = array(html.replace(' ','').split('\n'))
    
    EBV_SF =   float(html[22].replace(' ','')[:-6]) #float(h[argwhere(h=='<refPixelValueSandF>')[0][0]+1][:-5])
    dEBV_SF =  float(h[argwhere(h=='<stdSandF>')[0][0]+1][:-5]) 	  #float(h[argwhere(h=='<stdSandF>')[0][0]+1][:-5])
    EBV_SFD =  float(html[25].replace(' ','')[:-6])   #float(h[argwhere(h=='<refPixelValueSFD>')[0][0]+1][:-5])
    dEBV_SFD = float(h[argwhere(h=='<stdSFD>')[0][0]+1][:-5])		  #float(h[argwhere(h=='<stdSFD>')[0][0]+1][:-5])
    
    return EBV_SF, EBV_SFD

def sfd(x,y,system='gal'):
    """OInterrogator of Schlegel, Finkbeiner & Davis (1998) dust maps.
    Inputs are coordinates in the specified system.
    system: 'gal', 'equ j2000'"""
    if system == 'equ j2000':
        x,y = astCoords.convertCoords('J2000', 'GALACTIC', x,y, 2000.)
    ipath="/storage/astro2/phsmav/data2/ism/dust_dir/sfd98/"
    command=str("/storage/astro2/phsmav/data2/ism/dust_dir/sfd98/dust_getval %.2f %.2f ipath=%s" % (x,y,ipath))
    p = Popen(command, shell=True, stdin=PIPE, stdout=PIPE).communicate()
    Ebv=string.split(p[0])[5]
    return float(Ebv)

def sfd_old(x,y,system='gal'):
   """Measure E(B-V) from Schlegel, Finkbeiner & Davis (1998) dust maps
   default system: gal
   other: equ j2000
   """
   
   Nwcs = astWCS.WCS('/storage/astro2/phsmav/data2/ism/dust_dir/maps/SFD_dust_4096_ngp.fits')
   F = pyfits.open('/storage/astro2/phsmav/data2/ism/dust_dir/maps/SFD_dust_4096_ngp.fits',memmap=True)
   N = F[0].data
   F.close()
   
   Swcs = astWCS.WCS('/storage/astro2/phsmav/data2/ism/dust_dir/maps/SFD_dust_4096_sgp.fits')
   F = pyfits.open('/storage/astro2/phsmav/data2/ism/dust_dir/maps/SFD_dust_4096_sgp.fits',memmap=True)
   S = F[0].data
   F.close()	
   
   if hasattr(x, "__len__"):
      m,m_int, sigma = zeros((len(x))),zeros((len(x))),zeros((len(x)))
      if system=='equ j2000':
         gal = array([astCoords.convertCoords('J2000', 'GALACTIC', x[i],y[i], 2000.) for i in range(0,len(x))])
         x, y = gal[:,0], gal[:,1]
      for i in range(0,len(x)):
      
   	 n = Nwcs.wcs2pix(x[i],y[i])
   	 s = Swcs.wcs2pix(x[i],y[i])
   	 nC = array([[int(n[0]),int(n[1])],[int(n[0])+1,int(n[1])],[int(n[0])+1,int(n[1])+1],[int(n[0]),int(n[1])+1]])
         sC = array([[int(s[0]),int(s[1])],[int(s[0])+1,int(s[1])],[int(s[0])+1,int(s[1])+1],[int(s[0]),int(s[1])+1]])
         nC = nC[logical_and(nC[:,0]>0, nC[:,1]>0)]
         sC = sC[logical_and(sC[:,0]>0, sC[:,1]>0)]
	 nV = array([N[nC[0][0],nC[0][1]],N[nC[1][0],nC[1][1]],N[nC[2][0],nC[2][1]],N[nC[3][0],nC[3][1]]])
   	 sV = array([S[sC[0][0],sC[0][1]],S[sC[1][0],sC[1][1]],S[sC[2][0],sC[2][1]],S[sC[3][0],sC[3][1]]])
   	 V = hstack([nV,sV])
	 nt = (n[0]-int(n[0]))
	 ng = (n[1]-int(n[1]))
	 nv = (1-nt)*(1-ng)*nV[0]+nt*(1-ng)*nV[1]+nt*ng*nV[2]+(1-nt)*ng*nV[3]
	 st = (s[0]-int(s[0]))
	 sg = (s[1]-int(s[1]))
	 sv = (1-st)*(1-sg)*sV[0]+st*(1-ng)*sV[1]+st*sg*sV[2]+(1-st)*sg*sV[3]
	 m_i = array([nv,sv])
	 m_int[i] = mean(m_i[where(m_i>0)])
	 m[i] = mean(V[where(V>0)])
   	 sigma[i] = std(V[where(V>0)])     
   else:
   	 x, y = astCoords.convertCoords('J2000', 'GALACTIC', x,y, 2000.)
	 n = Nwcs.wcs2pix(x,y)
   	 s = Swcs.wcs2pix(x,y)
   	 nC = array([[int(n[0]),int(n[1])],[int(n[0])+1,int(n[1])],[int(n[0])+1,int(n[1])+1],[int(n[0]),int(n[1])+1]])
         sC = array([[int(s[0]),int(s[1])],[int(s[0])+1,int(s[1])],[int(s[0])+1,int(s[1])+1],[int(s[0]),int(s[1])+1]])
         nC = nC[logical_and(nC[:,0]>0, nC[:,1]>0)]
	 sC = sC[logical_and(sC[:,0]>0, sC[:,1]>0)]
	 nV = array([N[nC[0][0],nC[0][1]],N[nC[1][0],nC[1][1]],N[nC[2][0],nC[2][1]],N[nC[3][0],nC[3][1]]])
   	 sV = array([S[sC[0][0],sC[0][1]],S[sC[1][0],sC[1][1]],S[sC[2][0],sC[2][1]],S[sC[3][0],sC[3][1]]])
   	 V = hstack([nV,sV])
	 nt = (n[0]-int(n[0]))
	 ng = (n[1]-int(n[1]))
	 nv = (1-nt)*(1-ng)*nV[0]+nt*(1-ng)*nV[1]+nt*ng*nV[2]+(1-nt)*ng*nV[3]
	 st = (s[0]-int(s[0]))
	 sg = (s[1]-int(s[1]))
	 sv = (1-st)*(1-sg)*sV[0]+st*(1-ng)*sV[1]+st*sg*sV[2]+(1-st)*sg*sV[3]
	 m_i = array([nv,sv])
	 m_int = mean(m_i[where(m_i>0)])
	 m= mean(V[where(V>0)])
   	 sigma = std(V[where(V>0)])
   
   return m, m_int, sigma
    
