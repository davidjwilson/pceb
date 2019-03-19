"""Module for useful photometric transformations

	add_DM	-- compute apparent magnitudes for a given distance and reddening
	
	
"""

import reddenings
from numpy import *
from astLib import astSED
from pylab import *

def add_DM(M, passband, survey, dist, ebv, law="Fitzpatrick", tot_to_sel=3.1):
	"""Shift to cluster distance modulus (m-M)
	Inputs:
	M: intrinsic magnitude
	passband, surveys: allowed passbands and surveys are defined in reddenings.compute_extinction() 
	ebv: E(B-V)
	dist: distance in pc
	law: either Fitzpatrick (1999) or Cardelli et al (1989) laws
	R_v: 3.1"""

	Am = reddenings.compute_extinction(passband, survey, law=law, tot_to_sel=tot_to_sel, ebv=ebv)
	m = M + 5*log10(dist) - 5 + Am
	return m

def vega2AB(mag,err,passband):
	Vega = astSED.VEGA
	b1 = astSED.Passband('/storage/astro2/phsmav/data2/filters/poss/poss1_b.dat')
	b2 = astSED.Passband('/storage/astro2/phsmav/data2/filters/poss/poss2_b.dat')
	r1 = astSED.Passband('/storage/astro2/phsmav/data2/filters/poss/poss1_r.dat')
	r2 = astSED.Passband('/storage/astro2/phsmav/data2/filters/poss/poss2_r.dat')
	if passband=='b1':
		f0 = Vega.calcFlux(b1)
		flux,eflux = f0*10**(-0.4*mag),-0.4*f0*err*10**(-0.4*mag)
		ab = []
		for j in range(0, len(flux)):
			ab.append(astSED.flux2Mag(flux[j],eflux[j],b1))
		return array(ab)
	elif passband=='r1':
		f0 = Vega.calcFlux(r1)
		flux,eflux = f0*10**(-0.4*mag),-0.4*f0*err*10**(-0.4*mag)
		ab = []
		for j in range(0, len(flux)):
			ab.append(astSED.flux2Mag(flux[j],eflux[j],b1))
		return array(ab)
	elif passband=='b2':
		f0 = Vega.calcFlux(b2)
		flux,eflux = f0*10**(-0.4*mag),-0.4*f0*err*10**(-0.4*mag)
		ab = []
		for j in range(0, len(flux)):
			ab.append(astSED.flux2Mag(flux[j],eflux[j],b1))
		return array(ab)
	elif passband=='r2':
		f0 = Vega.calcFlux(r2)
		flux,eflux = f0*10**(-0.4*mag),-0.4*f0*err*10**(-0.4*mag)
		ab = []
		for j in range(0, len(flux)):
			ab.append(astSED.flux2Mag(flux[j],eflux[j],b1))
		return array(ab)


def average_errors(color, err_color, mag, err_mag, m1,m2,dm):
	r = arange(m1,m2+dm,dm)
	i = 0
	C,err_C,M,err_M = [], [], [], []
	while i<len(r)-1:
		c = color[logical_and(mag>=r[i],mag<=r[i+1])]
		c = mean(c[where(isnan(c)==False)])
		#m = color[logical_and(mag>=r[i],mag<=r[i+1])]
		#m = mean(m[where(isnan(m)==False)])
		m = r[i]+dm/2.		
		delta_c = err_color[logical_and(mag>=r[i],mag<=r[i+1])]
		delta_c = mean(delta_c[where(isnan(delta_c)==False)])
		delta_m = err_mag[logical_and(mag>=r[i],mag<=r[i+1])]
		delta_m = mean(delta_m[where(isnan(delta_m)==False)])
		C.append(c),err_C.append(delta_c),M.append(m),err_M.append(delta_m) 
		i = i +1
	C,err_C,M,err_M = array(C),array(err_C),array(M),array(err_M)
	return C, err_C, M, err_M



def reduced_proper_motion(mag, e_mag, PM, e_PM):
	""" Reduced proper motions in the input Band
	PM: is the total propermotion in arcsec/yr"""
	Hm = mag + 5*log10(PM) + 5
	e_Hm = sqrt(e_mag**2+(5.*e_PM/(log(10)*PM))**2)
	return Hm, e_Hm

def return_fluxes(mag, mag_err, passband, survey="vphas",units='erg/cm^2/s/A'):
	"""Transforming magnitudes into fluxes.
	Avaliable surveys, passbands:
	iphas: r, i, ha
	atlas, atlas-AB: u,g,r,i,z
	gaia: gp, bp, rp
	tycho, hipparcos: hp, bt, vt
	vphas: u, g, r, i, ha
	sdss: u,g,r,i,z
	panstarrs: g,r,i,y,z
	galex: fuv, nuv
	apass: b,v,g,r,i
	2mass: j,h,k
	ukidds: y,z,j,h,k
	wise: w1,w2,w3,w4
	spitzer: i1,i2,i3,i4
	Johnson: B,V
	Bessell: B,V
	ctio: B,V,R,I
	
	Output: tuple(flux,flux_error,wavelength)
	"""
	
	vega = astSED.VegaSED()
	if survey=="iphas":	
		if passband=="r":
			v_r = vega.calcFlux(reddenings.passbands.r_band_iphas)
			f,e_f = v_r*10**(-mag/2.5),abs(v_r*10**(-(mag+mag_err)/2.5)-v_r*10**(-mag/2.5))
			lambda_f = reddenings.passbands.r_band_iphas.effectiveWavelength()
			return f,e_f,lambda_f		
		elif passband=="i":
			v_i = vega.calcFlux(reddenings.passbands.i_band_iphas)
			f,e_f = v_i*10**(-mag/2.5),abs(v_i*10**(-(mag+mag_err)/2.5)-v_i*10**(-mag/2.5))
			lambda_f = reddenings.passbands.i_band_iphas.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="ha":
			v_ha = vega.calcFlux(reddenings.passbands.h_band_iphas)
			f,e_f = v_ha*10**(-mag/2.5),abs(v_ha*10**(-(mag+mag_err)/2.5)-v_ha*10**(-mag/2.5))
			lambda_f = reddenings.passbands.h_band_iphas.effectiveWavelength()
			return f,e_f,lambda_f	
	elif survey=="atlas":	
		if passband=="u":
			v_u = vega.calcFlux(reddenings.passbands.u_band_atlas)
			f,e_f = v_u*10**(-mag/2.5),abs(v_u*10**(-(mag+mag_err)/2.5)-v_u*10**(-mag/2.5))
			lambda_f = reddenings.passbands.u_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="g":
			v_g = vega.calcFlux(reddenings.passbands.g_band_atlas)
			f,e_f = v_g*10**(-mag/2.5),abs(v_g*10**(-(mag+mag_err)/2.5)-v_g*10**(-mag/2.5))
			lambda_f = reddenings.passbands.g_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="r":
			v_r = vega.calcFlux(reddenings.passbands.r_band_atlas)
			f,e_f = v_r*10**(-mag/2.5),abs(v_r*10**(-(mag+mag_err)/2.5)-v_r*10**(-mag/2.5))
			lambda_f = reddenings.passbands.r_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f		
		elif passband=="i":
			v_i = vega.calcFlux(reddenings.passbands.i_band_atlas)
			f,e_f = v_i*10**(-mag/2.5),abs(v_i*10**(-(mag+mag_err)/2.5)-v_i*10**(-mag/2.5))
			lambda_f = reddenings.passbands.i_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="z":
			v_z = vega.calcFlux(reddenings.passbands.z_band_atlas)
			f,e_f = v_z*10**(-mag/2.5),abs(v_z*10**(-(mag+mag_err)/2.5)-v_z*10**(-mag/2.5))
			lambda_f = reddenings.passbands.z_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="gaia":	
		if passband=="gp":
			v_gp = vega.calcFlux(reddenings.passbands.Gp_band_gaia)
			f,e_f = v_gp*10**(-mag/2.5),abs(v_gp*10**(-(mag+mag_err)/2.5)-v_gp*10**(-mag/2.5))
			lambda_f = reddenings.passbands.Gp_band_gaia.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="bp":
			v_bp = vega.calcFlux(reddenings.passbands.Bp_band_gaia)
			f,e_f = v_bp*10**(-mag/2.5),abs(v_bp*10**(-(mag+mag_err)/2.5)-v_bp*10**(-mag/2.5))
			lambda_f = reddenings.passbands.Bp_band_gaia.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="rp":
			v_rp = vega.calcFlux(reddenings.passbands.Rp_band_gaia)
			f,e_f = v_rp*10**(-mag/2.5),abs(v_rp*10**(-(mag+mag_err)/2.5)-v_rp*10**(-mag/2.5))
			lambda_f = reddenings.passbands.Rp_band_gaia.effectiveWavelength()
			return f,e_f,lambda_f	
	
	
	elif survey=="tycho" or survey=="hipparcos":	
		if passband=="hp":
			v_hp = vega.calcFlux(reddenings.passbands.Hp_band_tycho)
			f,e_f = v_hp*10**(-mag/2.5),abs(v_hp*10**(-(mag+mag_err)/2.5)-v_hp*10**(-mag/2.5))
			lambda_f = reddenings.passbands.Hp_band_tycho.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="bt":
			v_bt = vega.calcFlux(reddenings.passbands.Bt_band_tycho)
			f,e_f = v_bt*10**(-mag/2.5),abs(v_bt*10**(-(mag+mag_err)/2.5)-v_bt*10**(-mag/2.5))
			lambda_f = reddenings.passbands.Bt_band_tycho.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="vt":
			v_vt = vega.calcFlux(reddenings.passbands.Vt_band_tycho)
			f,e_f = v_vt*10**(-mag/2.5),abs(v_vt*10**(-(mag+mag_err)/2.5)-v_vt*10**(-mag/2.5))
			lambda_f = reddenings.passbands.Vt_band_tycho.effectiveWavelength()
			return f,e_f,lambda_f	
	elif survey=="atlas-AB":	
		if passband=="u":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.u_band_atlas)
			lambda_f = reddenings.passbands.u_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="g":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.g_band_atlas)
			lambda_f = reddenings.passbands.g_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="r":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.r_band_atlas)
			lambda_f = reddenings.passbands.r_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="i":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.i_band_atlas)
			lambda_f = reddenings.passbands.i_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="z":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.z_band_atlas)
			lambda_f = reddenings.passbands.z_band_atlas.effectiveWavelength()
			return f,e_f,lambda_f			
	elif survey=="vphas":	
		if passband=="u":
			v_u = vega.calcFlux(reddenings.passbands.u_band_vphas)
			f,e_f = v_u*10**(-mag/2.5),abs(v_u*10**(-(mag+mag_err)/2.5)-v_u*10**(-mag/2.5))
			lambda_f = reddenings.passbands.u_band_vphas.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="g":
			v_g = vega.calcFlux(reddenings.passbands.g_band_vphas)
			f,e_f = v_g*10**(-mag/2.5),abs(v_g*10**(-(mag+mag_err)/2.5)-v_g*10**(-mag/2.5))
			lambda_f = reddenings.passbands.g_band_vphas.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="r":
			v_r = vega.calcFlux(reddenings.passbands.r_band_vphas)
			f,e_f = v_r*10**(-mag/2.5),abs(v_r*10**(-(mag+mag_err)/2.5)-v_r*10**(-mag/2.5))
			lambda_f = reddenings.passbands.r_band_vphas.effectiveWavelength()
			return f,e_f,lambda_f		
		elif passband=="i":
			v_i = vega.calcFlux(reddenings.passbands.i_band_vphas)
			f,e_f = v_i*10**(-mag/2.5),abs(v_i*10**(-(mag+mag_err)/2.5)-v_i*10**(-mag/2.5))
			lambda_f = reddenings.passbands.i_band_vphas.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="ha":
			v_ha = vega.calcFlux(reddenings.passbands.h_band_vphas)
			f,e_f = v_ha*10**(-mag/2.5),abs(v_ha*10**(-(mag+mag_err)/2.5)-v_ha*10**(-mag/2.5))
			lambda_f = reddenings.passbands.h_band_vphas.effectiveWavelength()
			return f,e_f,lambda_f	
	elif survey=="sdss":
		if passband=="u":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.u_band_sdss)
			lambda_f = reddenings.passbands.u_band_sdss.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="g":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.g_band_sdss)
			lambda_f = reddenings.passbands.g_band_sdss.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="r":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.r_band_sdss)
			lambda_f = reddenings.passbands.r_band_sdss.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="i":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.i_band_sdss)
			lambda_f = reddenings.passbands.i_band_sdss.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="z":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.z_band_sdss)
			lambda_f = reddenings.passbands.z_band_sdss.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="panstarrs":
		if passband=="g":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.g_band_ps1)
			lambda_f = reddenings.passbands.g_band_ps1.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="r":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.r_band_ps1)
			lambda_f = reddenings.passbands.r_band_ps1.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="i":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.i_band_ps1)
			lambda_f = reddenings.passbands.i_band_ps1.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="y":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.y_band_ps1)
			lambda_f = reddenings.passbands.y_band_ps1.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="z":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.z_band_ps1)
			lambda_f = reddenings.passbands.z_band_ps1.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="galex":
		if passband=="fuv":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.fuv_band_galex)
			lambda_f = reddenings.passbands.fuv_band_galex.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="nuv":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.nuv_band_galex)
			lambda_f = reddenings.passbands.nuv_band_galex.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="apass":
		if passband=="b":
			v_b = vega.calcFlux(reddenings.passbands.B_band_apass)
			f,e_f = v_b*10**(-mag/2.5),abs(v_b*10**(-(mag+mag_err)/2.5)-v_b*10**(-mag/2.5))
			lambda_f = reddenings.passbands.B_band_apass.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="v":
			v_v = vega.calcFlux(reddenings.passbands.V_band_apass)
			f,e_f = v_v*10**(-mag/2.5),abs(v_v*10**(-(mag+mag_err)/2.5)-v_v*10**(-mag/2.5))
			lambda_f = reddenings.passbands.V_band_apass.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="g":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.g_band_apass)
			lambda_f = reddenings.passbands.g_band_apass.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="r":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.r_band_apass)
			lambda_f = reddenings.passbands.r_band_apass.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="i":
			f,e_f = astSED.mag2Flux(mag,mag_err,reddenings.passbands.i_band_apass)
			lambda_f = reddenings.passbands.i_band_apass.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="2mass":
		if passband=="j":
			v_j = vega.calcFlux(reddenings.passbands.j_band_2mass)
			f,e_f = v_j*10**(-mag/2.5),abs(v_j*10**(-(mag+mag_err)/2.5)-v_j*10**(-mag/2.5))
			lambda_f = reddenings.passbands.j_band_2mass.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="h":
			v_h = vega.calcFlux(reddenings.passbands.h_band_2mass)
			f,e_f = v_h*10**(-mag/2.5),abs(v_h*10**(-(mag+mag_err)/2.5)-v_h*10**(-mag/2.5))
			lambda_f = reddenings.passbands.h_band_2mass.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="k":
			v_k = vega.calcFlux(reddenings.passbands.k_band_2mass)
			f,e_f = v_k*10**(-mag/2.5),abs(v_k*10**(-(mag+mag_err)/2.5)-v_k*10**(-mag/2.5))
			lambda_f = reddenings.passbands.k_band_2mass.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="ukidss":
		if passband=="y":
			v_y = vega.calcFlux(reddenings.passbands.y_band_ukidss)
			f,e_f = v_y*10**(-mag/2.5),abs(v_y*10**(-(mag+mag_err)/2.5)-v_y*10**(-mag/2.5))
			lambda_f = reddenings.passbands.y_band_ukidss.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="z":
			v_z = vega.calcFlux(reddenings.passbands.z_band_ukidss)
			f,e_f = v_z*10**(-mag/2.5),abs(v_z*10**(-(mag+mag_err)/2.5)-v_z*10**(-mag/2.5))
			lambda_f = reddenings.passbands.z_band_ukidss.effectiveWavelength()
		        return f,e_f,lambda_f
		elif passband=="j":
			v_j = vega.calcFlux(reddenings.passbands.j_band_ukidss)
			f,e_f = v_j*10**(-mag/2.5),abs(v_j*10**(-(mag+mag_err)/2.5)-v_j*10**(-mag/2.5))
			lambda_f = reddenings.passbands.j_band_ukidss.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="h":
			v_h = vega.calcFlux(reddenings.passbands.h_band_ukidss)
			f,e_f = v_h*10**(-mag/2.5),abs(v_h*10**(-(mag+mag_err)/2.5)-v_h*10**(-mag/2.5))
			lambda_f = reddenings.passbands.h_band_ukidss.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="k":
			v_k = vega.calcFlux(reddenings.passbands.k_band_ukidss)
			f,e_f = v_k*10**(-mag/2.5),abs(v_k*10**(-(mag+mag_err)/2.5)-v_k*10**(-mag/2.5))
			lambda_f = reddenings.passbands.k_band_ukidss.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="wise":
		if passband=="w1":
			v_w1 = vega.calcFlux(reddenings.passbands.w1_band_wise)
			f,e_f = v_w1*10**(-mag/2.5),abs(v_w1*10**(-(mag+mag_err)/2.5)-v_w1*10**(-mag/2.5))
			lambda_f = reddenings.passbands.w1_band_wise.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="w2":
			v_w2 = vega.calcFlux(reddenings.passbands.w2_band_wise)
			f,e_f = v_w2*10**(-mag/2.5),abs(v_w2*10**(-(mag+mag_err)/2.5)-v_w2*10**(-mag/2.5))
			lambda_f = reddenings.passbands.w2_band_wise.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="w3":
			v_w3 = vega.calcFlux(reddenings.passbands.w3_band_wise)
			f,e_f = v_w3*10**(-mag/2.5),abs(v_w3*10**(-(mag+mag_err)/2.5)-v_w3*10**(-mag/2.5))
			lambda_f = reddenings.passbands.w3_band_wise.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="w4":
			v_w4 = vega.calcFlux(reddenings.passbands.w4_band_wise)
			f,e_f = v_w4*10**(-mag/2.5),abs(v_w4*10**(-(mag+mag_err)/2.5)-v_w4*10**(-mag/2.5))
			lambda_f = reddenings.passbands.w4_band_wise.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="spitzer":
		if passband=="i1":
			v_i1 = vega.calcFlux(reddenings.passbands.i1_band_spitzer)
			f,e_f = v_i1*10**(-mag/2.5),abs(v_i1*10**(-(mag+mag_err)/2.5)-v_i1*10**(-mag/2.5))
			lambda_f = reddenings.passbands.i1_band_spitzer.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="i2":
			v_i2 = vega.calcFlux(reddenings.passbands.i2_band_spitzer)
			f,e_f = v_i2*10**(-mag/2.5),abs(v_i2*10**(-(mag+mag_err)/2.5)-v_i2*10**(-mag/2.5))
			lambda_f = reddenings.passbands.i2_band_spitzer.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="i3":
			v_i3 = vega.calcFlux(reddenings.passbands.i3_band_spitzer)
			f,e_f = v_i3*10**(-mag/2.5),abs(v_i3*10**(-(mag+mag_err)/2.5)-v_i3*10**(-mag/2.5))
			lambda_f = reddenings.passbands.i3_band_spitzer.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="i4":
			v_i4 = vega.calcFlux(reddenings.passbands.i4_band_spitzer)
			f,e_f = v_i4*10**(-mag/2.5),abs(v_i4*10**(-(mag+mag_err)/2.5)-v_i4*10**(-mag/2.5))
			lambda_f = reddenings.passbands.i4_band_spitzer.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="m1":
			v_m1 = vega.calcFlux(reddenings.passbands.m1_band_spitzer)
			f,e_f = v_m1*10**(-mag/2.5),abs(v_m1*10**(-(mag+mag_err)/2.5)-v_m1*10**(-mag/2.5))
			lambda_f = reddenings.passbands.m1_band_spitzer.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="Johnson":
		if passband=="V":
			v_v = vega.calcFlux(reddenings.passbands.V_band_johnson)
			f,e_f = v_v*10**(-mag/2.5),abs(v_v*10**(-(mag+mag_err)/2.5)-v_v*10**(-mag/2.5))
			lambda_f = reddenings.passbands.V_band_johnson.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="B":
			v_b = vega.calcFlux(reddenings.passbands.B_band_johnson)
			f,e_f = v_b*10**(-mag/2.5),abs(v_b*10**(-(mag+mag_err)/2.5)-v_b*10**(-mag/2.5))
			lambda_f = reddenings.passbands.B_band_johnson.effectiveWavelength()
			return f,e_f,lambda_f
	elif survey=="Bessell":
		if passband=="V":
			v_v = vega.calcFlux(reddenings.passbands.V_band_bessell)
			f,e_f = v_v*10**(-mag/2.5),abs(v_v*10**(-(mag+mag_err)/2.5)-v_v*10**(-mag/2.5))
			lambda_f = reddenings.passbands.V_band_bessell.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="B":
			v_b = vega.calcFlux(reddenings.passbands.B_band_bessell)
			f,e_f = v_b*10**(-mag/2.5),abs(v_b*10**(-(mag+mag_err)/2.5)-v_b*10**(-mag/2.5))
			lambda_f = reddenings.passbands.B_band_bessell.effectiveWavelength()
			return f,e_f,lambda_f			
	elif survey=="ctio":
		if passband=="V":
			v_v = vega.calcFlux(reddenings.passbands.V_band_ctio)
			f,e_f = v_v*10**(-mag/2.5),abs(v_v*10**(-(mag+mag_err)/2.5)-v_v*10**(-mag/2.5))
			lambda_f = reddenings.passbands.V_band_ctio.effectiveWavelength()
			return f,e_f,lambda_f	
		elif passband=="B":
			v_b = vega.calcFlux(reddenings.passbands.B_band_ctio)
			f,e_f = v_b*10**(-mag/2.5),abs(v_b*10**(-(mag+mag_err)/2.5)-v_b*10**(-mag/2.5))
			lambda_f = reddenings.passbands.B_band_ctio.effectiveWavelength()
			return f,e_f,lambda_f
		elif passband=="R":
			v_b = vega.calcFlux(reddenings.passbands.R_band_ctio)
			f,e_f = v_b*10**(-mag/2.5),abs(v_b*10**(-(mag+mag_err)/2.5)-v_b*10**(-mag/2.5))
			lambda_f = reddenings.passbands.R_band_ctio.effectiveWavelength()
			return f,e_f,lambda_f			
		elif passband=="I":
			v_b = vega.calcFlux(reddenings.passbands.I_band_ctio)
			f,e_f = v_b*10**(-mag/2.5),abs(v_b*10**(-(mag+mag_err)/2.5)-v_b*10**(-mag/2.5))
			lambda_f = reddenings.passbands.I_band_ctio.effectiveWavelength()
			return f,e_f,lambda_f			
			
			

def galex_mag2flux(mag):
	""" This gives the flux in cnts/sec
	Note: the magnitude must not include the zeropoint, i.e.
	mag = m -zp, where zp=20.08 for the NUV and zp=18.82 for the FUV """	
	flux = 10**(-mag*0.4)
	return flux

def aper_corr_galex(mag, band, aperture=3):
	""" Aperture correction for GALEX magnitudes
	from Morrissey et al. (2007)"""
	if band=="nuv":
		if aperture==1:
			m = mag - 2.09
			return m
		elif aperture==2:
			m = mag - 1.33
			return m
		elif aperture==3:
			m = mag - 0.59
			return m
		elif aperture==4:
			m = mag - 0.23
			return m
		elif aperture==5:
			m = mag - 0.13
			return m
		elif aperture==6:
			m = mag - 0.09
			return m
		elif aperture==7:
			m = mag - 0.07
			return m
	elif band=="fuv":
		if aperture==1:
			m = mag - 1.65
			return m
		elif aperture==2:
			m = mag - 0.96
			return m
		elif aperture==3:
			m = mag - 0.36
			return m
		elif aperture==4:
			m = mag - 0.15
			return m
		elif aperture==5:
			m = mag - 0.10
			return m
		elif aperture==6:
			m = mag - 0.09
			return m
		elif aperture==7:
			m = mag - 0.07
			return m

def zp_galex(mag, band):
	""" Zero points
	from Morrissey et al. (2007)"""
	if band=="nuv":
		m = 20.08 + mag
		return m
	elif band=="fuv":
		m = 18.82 + mag
		return m

def nonlinearity_galex(mag, band, zp=False):
	""" zp==False: Nonlinearity correction for GALEX magnitudes
	from Morrissey et al. (2007)
	zp==True: Nonlinearity correction for GALEX magnitudes
	from 
	"""	

	if zp==False:
		if band=="nuv":
			C0,C1,C2 = -0.314,1.365,-0.103 
			log_MR = -0.4*mag
			if type(log_MR) == ndarray:
				#log_MR[where(log_MR<0.)] = nan
				log_PR_plus = log_MR.copy()
				log_PR_plus[where(10**log_MR<100)] = (-C1 + sqrt(C1**2 - 4*C2*(C0-log_MR[where(10**log_MR<100)])))/(2*C2)
			else:
				if 10**log_MR<100:
					log_PR_plus = (-C1 + sqrt(C1**2 - 4*C2*(C0-log_MR)))/(2*C2)
				else:
					log_PR_plus = log_MR
			mag_corr = -2.5*log_PR_plus
			return mag_corr

		elif band=="fuv":
			C0,C1,C2 = -0.531,1.696,-0.225 
			log_MR = -0.4*mag
			if type(log_MR) == ndarray:
				#log_MR[where(log_MR<0.)] = nan
				log_PR_plus = log_MR.copy()
				log_PR_plus[where(10**log_MR<100)] = (-C1 + sqrt(C1**2 - 4*C2*(C0-log_MR[where(10**log_MR<100)])))/(2*C2)
			else:
				if 10**log_MR<100:
					log_PR_plus = (-C1 + sqrt(C1**2 - 4*C2*(C0-log_MR)))/(2*C2)
				else:
					log_PR_plus = log_MR
			mag_corr = -2.5*log_PR_plus
			return mag_corr
	if zp==True:
		if band=="nuv":
			C0,C1,C2 = 2.634,26.316,-245.329
			mag_corr = C0+sqrt(C1*mag+C2)
			return mag_corr

		elif band=="fuv":
			C0,C1,C2 = 5.371,20.00,-210.200 
			mag_corr = C0+sqrt(C1*mag+C2)
			return mag_corr



def load_intrinsic_magnitudes(tipo="MS", modelli="pickles98", classe="V", logg=None):
    """ This functions loads the absolute magnitudes for
    main sequence stars or white dwarfs.
    
    Allowed inputs:
    tipo: MS, DA, DB, BB
    modelli: pickles98, koester (MLalpha=0.8, mass-radius relation by Bergeron) 
    classe: V, IV, III, I, None
    logg: 5.25-9.5 (DA white dwarfs), 7-9 (DB white dwarfs) 
    
    main sequence models are stored here:
    /storage/astro2/phsmav/data/model_atmospheres/spectral_libraries/pickles/uvk/
    white dwarf models are stored here:
    /storage/astro2/phsmav/data2/model_atmospheres/DA/DA_SDSS_2014/
    /storage/astro2/phsmav/data2/model_atmospheres/DB/data/sdss/  
    """
    if tipo=="MS" and modelli=="pickles98":
        if classe=="V":
	    P = csv2rec('/storage/astro2/phsmav/data2/models/ms_colours/pickles/magnitudes_BC/pickles_magnitudes_0.00.csv')
            P = P[where(P['lumclass']==5)]
        elif classe=="IV":
	    P = csv2rec('/storage/astro2/phsmav/data2/models/ms_colours/pickles/magnitudes_BC/pickles_magnitudes_0.00.csv')
            P = P[where(P['lumclass']==4)]
        elif classe=="III":
	    P = csv2rec('/storage/astro2/phsmav/data2/models/ms_colours/pickles/magnitudes_BC/pickles_magnitudes_0.00.csv')
            P = P[where(P['lumclass']==3)]
        elif classe=="I":
	    P = csv2rec('/storage/astro2/phsmav/data2/models/ms_colours/pickles/magnitudes_BC/pickles_magnitudes_0.00.csv')
            P = P[where(P['lumclass']==1)]
        elif classe is None:
	    P = csv2rec('/storage/astro2/phsmav/data2/models/ms_colours/pickles/magnitudes_BC/pickles_magnitudes_0.00.csv')
    elif tipo=="DA":
        P = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_koester10_magnitudes.dat')
        if logg is not None:
	    P= P[where(P['logg']==logg)]
    elif tipo=="DB":
        P = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/db_koester10_magnitudes.dat')
        if logg is not None:
	    P= P[where(P['logg']==logg)]
    elif tipo=="BB":
        P = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/black-body_magnitudes.dat')
        if logg is not None:
	    P= P[where(P['logg']==logg)]
    return P







def galex_poss_intrinsic(color1,color2,classe,color='red'):
	""" Function plotting MS,DA sequences in the GALEX/APASS/SDSS color-color diagrams """
	#logg650 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_650_Av0.txt')
	logg700 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/DA/log_g_700_Av0.txt')
	logg750 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/DA/log_g_750_Av0.txt')
	logg800 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/DA/log_g_800_Av0.txt')
	logg850 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/DA/log_g_850_Av0.txt')
	logg900 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/DA/log_g_900_Av0.txt')
	#MS = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/MS/MS_Av0.txt',usecols=[1,2,3,4,5,6,7])
	#ZZlogg700 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/ZZCeti/log_g_7.00_Av0.txt')
	#ZZlogg750 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/ZZCeti/log_g_7.50_Av0.txt')
	#ZZlogg800 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/ZZCeti/log_g_8.00_Av0.txt')
	#ZZlogg850 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/ZZCeti/log_g_8.50_Av0.txt')
	#ZZlogg900 = loadtxt('/storage/astro2/phsmav/data2/ppmxl/tracks/ZZCeti/log_g_9.00_Av0.txt')	
	if color1=="fuv-b1" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,2],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,2],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,2],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,2],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,2],logg900[:,1], color=color,zorder=5)
	elif color1=="fuv-b2" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,3],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,3],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,3],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,3],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,3],logg900[:,1], color=color,zorder=5)
	elif color1=="fuv-r1" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,4],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,4],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,4],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,4],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,4],logg900[:,1], color=color,zorder=5)
	elif color1=="fuv-r2" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,5],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,5],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,5],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,5],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,5],logg900[:,1], color=color,zorder=5)
	if color1=="nuv-b1" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,6],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,6],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,6],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,6],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,6],logg900[:,1], color=color,zorder=5)
	elif color1=="nuv-b2" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,7],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,7],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,7],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,7],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,7],logg900[:,1], color=color,zorder=5)
	elif color1=="nuv-r1" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,8],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,8],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,8],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,8],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,8],logg900[:,1], color=color,zorder=5)
	elif color1=="nuv-r2" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,9],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,9],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,9],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,9],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,9],logg900[:,1], color=color,zorder=5)
	elif color1=="b1-r1" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,10],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,10],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,10],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,10],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,10],logg900[:,1], color=color,zorder=5)
	elif color1=="b2-r2" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,11],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,11],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,11],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,11],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,11],logg900[:,1], color=color,zorder=5)

	if color1=="nuv-r" and color2=="fuv-nuv" and classe=="ZZ":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(ZZlogg700[:,3],ZZlogg700[:,1], color=color,zorder=5)
		plot(ZZlogg750[:,3],ZZlogg750[:,1], color=color,zorder=5)
		plot(ZZlogg800[:,3],ZZlogg800[:,1], color=color,zorder=5)
		plot(ZZlogg850[:,3],ZZlogg850[:,1], color=color,zorder=5)
		plot(ZZlogg900[:,3],ZZlogg900[:,1], color=color,zorder=5)
	elif color1=="g-i" and color2=="nuv-g" and classe=="ZZ":
		#plot(logg650[:,4],logg650[:,2], color=color,zorder=5)
		plot(ZZlogg700[:,4],ZZlogg700[:,2], color=color,zorder=5)
		plot(ZZlogg750[:,4],ZZlogg750[:,2], color=color,zorder=5)
		plot(ZZlogg800[:,4],ZZlogg800[:,2], color=color,zorder=5)
		plot(ZZlogg850[:,4],ZZlogg850[:,2], color=color,zorder=5)
		plot(ZZlogg900[:,4],ZZlogg900[:,2], color=color,zorder=5)
	elif color1=="g-r" and color2=="nuv-g" and classe=="ZZ":
		#plot(logg650[:,5],logg650[:,2], color=color,zorder=5)
		plot(ZZlogg700[:,5],ZZlogg700[:,2], color=color,zorder=5)
		plot(ZZlogg750[:,5],ZZlogg750[:,2], color=color,zorder=5)
		plot(ZZlogg800[:,5],ZZlogg800[:,2], color=color,zorder=5)
		plot(ZZlogg850[:,5],ZZlogg850[:,2], color=color,zorder=5)
		plot(ZZlogg900[:,5],ZZlogg900[:,2], color=color,zorder=5)

	elif color1=="nuv-r" and color2=="fuv-nuv" and classe=="MS":
		plot(MS[:,3],MS[:,1], color=color,zorder=5, ls='dashed')
	elif color1=="g-r" and color2=="nuv-g" and classe=="MS":
		plot(MS[:,5],MS[:,2], color=color,zorder=5, ls='dashed')	
	elif color1=="g-i" and color2=="nuv-g" and classe=="MS":
		plot(MS[:,4],MS[:,2], color=color,zorder=5, ls='dashed')

def vista_intrinsic(color1,color2,classe='MS',reddening=[0,0],dist=10,color='red',ls='solid'):
        """ Vista intrinsic colours of Pickles 1998 stars. Only main sequence and giants for now"""
	P = csv2rec('/storage/astro2/phsmav/data2/models/ms_colours/pickles/magnitudes/pickles_magnitudes_vista_0.00.csv')	
	Av = P[where(P['lumclass']==5)]
	Aiii = P[where(P['lumclass']==3)]
	if color1=="j-h" and color2=="z-y" and classe=="MS":
		plot(Av['j']-Av['h']+reddening[0],Av['z']-Av['y']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="h-k" and color2=="j-h" and classe=="MS":
		plot(Av['h']-Av['k']+reddening[0],Av['j']-Av['h']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="j-h" and color2=="j" and classe=="MS":
		plot(Av['j']-Av['h']+reddening[0],Av['j']+reddening[1]+5*log(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="h-k" and color2=="h" and classe=="MS":
		plot(Av['h']-Av['k']+reddening[0],Av['h']+reddening[1]+5*log(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="j-k" and color2=="j" and classe=="MS":
		plot(Av['j']-Av['k']+reddening[0],Av['j']+reddening[1]+5*log(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	if color1=="j-h" and color2=="z-y" and classe=="GB":
		plot(Aiii['j']-Aiii['h']+reddening[0],Aiii['z']-Aiii['y']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="h-k" and color2=="j-h" and classe=="GB":
		plot(Aiii['h']-Aiii['k']+reddening[0],Aiii['j']-Aiii['h']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="j-h" and color2=="j" and classe=="GB":
		plot(Aiii['j']-Aiii['h']+reddening[0],Aiii['j']+reddening[1]+5*log(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="h-k" and color2=="h" and classe=="GB":
		plot(Aiii['h']-Aiii['k']+reddening[0],Aiii['h']+reddening[1]+5*log(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="j-k" and color2=="j" and classe=="GB":
		plot(Aiii['j']-Aiii['k']+reddening[0],Aiii['j']+reddening[1]+5*log(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)

def vphas_intrinsic(color1,color2,classe,reddening=[0,0],dist=10,color='red',ls='solid'):
	""" Function plotting MS,DA sequences in the VPHAS+ color-color diagrams 
	reddening = [E(color1),E(color2)]"""
	DA = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_koester10_magnitudes.dat')
	logg700 = DA[where(DA['logg']==7.00)]
	logg750 = DA[where(DA['logg']==7.50)]
	logg800 = DA[where(DA['logg']==8.00)]
	logg850 = DA[where(DA['logg']==8.50)] 
	logg900 = DA[where(DA['logg']==9.00)]
	DB = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/db_koester10_magnitudes.dat')

	
	DBlogg700 =  DB[where(DB['logg']==7.00)] 
	DBlogg750 =  DB[where(DB['logg']==7.50)]
	DBlogg800 =  DB[where(DB['logg']==8.00)]	  
	DBlogg850 =  DB[where(DB['logg']==8.50)]
	DBlogg900 =  DB[where(DB['logg']==9.00)]
	
	MS = load_intrinsic_magnitudes()
	
	if color1=="g-r" and color2=="u-g" and classe=="DA":
		plot(logg700['g_vphas']-logg700['r_vphas']+reddening[0],logg700['u_vphas']-logg700['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg750['g_vphas']-logg750['r_vphas']+reddening[0],logg750['u_vphas']-logg750['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg800['g_vphas']-logg800['r_vphas']+reddening[0],logg800['u_vphas']-logg800['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg850['g_vphas']-logg850['r_vphas']+reddening[0],logg850['u_vphas']-logg850['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg900['g_vphas']-logg900['r_vphas']+reddening[0],logg900['u_vphas']-logg900['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="g-r" and classe=="DA":
		plot(logg700['r_vphas']-logg700['i_vphas']+reddening[0],logg700['g_vphas']-logg700['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg750['r_vphas']-logg750['i_vphas']+reddening[0],logg750['g_vphas']-logg750['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg800['r_vphas']-logg800['i_vphas']+reddening[0],logg800['g_vphas']-logg800['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg850['r_vphas']-logg850['i_vphas']+reddening[0],logg850['g_vphas']-logg850['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg900['r_vphas']-logg900['i_vphas']+reddening[0],logg900['g_vphas']-logg900['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="r-Ha" and classe=="DA":
		plot(logg700['r_vphas']-logg700['i_vphas']+reddening[0],logg700['r_vphas']-logg700['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg750['r_vphas']-logg750['i_vphas']+reddening[0],logg750['r_vphas']-logg750['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg800['r_vphas']-logg800['i_vphas']+reddening[0],logg800['r_vphas']-logg800['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg850['r_vphas']-logg850['i_vphas']+reddening[0],logg850['r_vphas']-logg850['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg900['r_vphas']-logg900['i_vphas']+reddening[0],logg900['r_vphas']-logg900['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="g-r" and color2=="g" and classe=="DA":
		plot(logg700['g_vphas']-logg700['r_vphas']+reddening[0],logg700['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg750['g_vphas']-logg750['r_vphas']+reddening[0],logg750['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg800['g_vphas']-logg800['r_vphas']+reddening[0],logg800['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg850['g_vphas']-logg850['r_vphas']+reddening[0],logg850['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg900['g_vphas']-logg900['r_vphas']+reddening[0],logg900['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="r" and classe=="DA":
		plot(logg700['r_vphas']-logg700['i_vphas']+reddening[0],logg700['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg750['r_vphas']-logg750['i_vphas']+reddening[0],logg750['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg800['r_vphas']-logg800['i_vphas']+reddening[0],logg800['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg850['r_vphas']-logg850['i_vphas']+reddening[0],logg850['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg900['r_vphas']-logg900['i_vphas']+reddening[0],logg900['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	if color1=="g-r" and color2=="u-g" and classe=="DB":
		plot(DBlogg700['g_vphas']-DBlogg700['r_vphas']+reddening[0],DBlogg700['u_vphas']-DBlogg700['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg750['g_vphas']-DBlogg750['r_vphas']+reddening[0],DBlogg750['u_vphas']-DBlogg750['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg800['g_vphas']-DBlogg800['r_vphas']+reddening[0],DBlogg800['u_vphas']-DBlogg800['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg850['g_vphas']-DBlogg850['r_vphas']+reddening[0],DBlogg850['u_vphas']-DBlogg850['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg900['g_vphas']-DBlogg900['r_vphas']+reddening[0],DBlogg900['u_vphas']-DBlogg900['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="g-r" and classe=="DB":
		plot(DBlogg700['r_vphas']-DBlogg700['i_vphas']+reddening[0],DBlogg700['g_vphas']-DBlogg700['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg750['r_vphas']-DBlogg750['i_vphas']+reddening[0],DBlogg750['g_vphas']-DBlogg750['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg800['r_vphas']-DBlogg800['i_vphas']+reddening[0],DBlogg800['g_vphas']-DBlogg800['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg850['r_vphas']-DBlogg850['i_vphas']+reddening[0],DBlogg850['g_vphas']-DBlogg850['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg900['r_vphas']-DBlogg900['i_vphas']+reddening[0],DBlogg900['g_vphas']-DBlogg900['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="r-Ha" and classe=="DB":
		plot(DBlogg700['r_vphas']-DBlogg700['i_vphas']+reddening[0],DBlogg700['r_vphas']-DBlogg700['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg750['r_vphas']-DBlogg750['i_vphas']+reddening[0],DBlogg750['r_vphas']-DBlogg750['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg800['r_vphas']-DBlogg800['i_vphas']+reddening[0],DBlogg800['r_vphas']-DBlogg800['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg850['r_vphas']-DBlogg850['i_vphas']+reddening[0],DBlogg850['r_vphas']-DBlogg850['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg900['r_vphas']-DBlogg900['i_vphas']+reddening[0],DBlogg900['r_vphas']-DBlogg900['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	if color1=="g-r" and color2=="u-g" and classe=="MS":
		plot(MS['g_vphas']-MS['r_vphas']+reddening[0],MS['u_vphas']-MS['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="g-r" and classe=="MS":
		plot(MS['r_vphas']-MS['i_vphas']+reddening[0],MS['g_vphas']-MS['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="r-Ha" and classe=="MS":
		plot(MS['r_vphas']-MS['i_vphas']+reddening[0],MS['r_vphas']-MS['ha_vphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="g-r" and color2=="g" and classe=="MS":
		plot(MS['g_vphas']-MS['r_vphas']+reddening[0],MS['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="r" and classe=="MS":
		plot(MS['r_vphas']-MS['i_vphas']+reddening[0],MS['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)

def iphas_intrinsic(color1,color2,classe,reddening=[0,0],dist=10,color='red',ls='solid'):
	""" Function plotting MS,DA sequences in the VPHAS+ color-color diagrams 
	reddening = [E(color1),E(color2)]"""
	DA = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_koester10_magnitudes.dat')
	logg700 = DA[where(DA['logg']==7.00)]
	logg750 =  DA[where(DA['logg']==7.50)]
	logg800 =  DA[where(DA['logg']==8.00)]
	logg850 = DA[where(DA['logg']==8.50)] 
	logg900 =  DA[where(DA['logg']==9.00)]
	DB = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/db_koester10_magnitudes.dat')
	DBlogg700 = DB[where(DB['logg']==7.00)]
	DBlogg750 =  DB[where(DB['logg']==7.50)]
	DBlogg800 =  DB[where(DB['logg']==8.00)]
	DBlogg850 = DB[where(DB['logg']==8.50)] 
	DBlogg900 =  DB[where(DB['logg']==9.00)]	
	MS = loadtxt('/storage/astro2/phsmav/data2/iphas/tracks/iphas.txt',usecols=[1,2])
	xA0 = arange(0.029, 2.9, 0.05)
	yA0= -0.009+0.330*xA0-0.0455*xA0**2
	if color1=="r-i" and color2=="r-Ha" and classe=="DA":
		plot(logg700['r_iphas']-logg700['i_iphas']+reddening[0],logg700['r_iphas']-logg700['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg750['r_iphas']-logg750['i_iphas']+reddening[0],logg750['r_iphas']-logg750['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg800['r_iphas']-logg800['i_iphas']+reddening[0],logg800['r_iphas']-logg800['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg850['r_iphas']-logg850['i_iphas']+reddening[0],logg850['r_iphas']-logg850['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg900['r_iphas']-logg900['i_iphas']+reddening[0],logg900['r_iphas']-logg900['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="r" and classe=="DA":
		plot(logg700['r_iphas']-logg700['i_iphas']+reddening[0],logg700['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg750['r_iphas']-logg750['i_iphas']+reddening[0],logg750['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg800['r_iphas']-logg800['i_iphas']+reddening[0],logg800['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg850['r_iphas']-logg850['i_iphas']+reddening[0],logg850['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(logg900['r_iphas']-logg900['i_iphas']+reddening[0],logg900['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
	if color1=="r-i" and color2=="r-Ha" and classe=="DB":
		plot(DBlogg700['r_iphas']-DBlogg700['i_iphas']+reddening[0],DBlogg700['r_iphas']-DBlogg700['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg750['r_iphas']-DBlogg750['i_iphas']+reddening[0],DBlogg750['r_iphas']-DBlogg750['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg800['r_iphas']-DBlogg800['i_iphas']+reddening[0],DBlogg800['r_iphas']-DBlogg800['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg850['r_iphas']-DBlogg850['i_iphas']+reddening[0],DBlogg850['r_iphas']-DBlogg850['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg900['r_iphas']-DBlogg900['i_iphas']+reddening[0],DBlogg900['r_iphas']-DBlogg900['ha_iphas']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="r-i" and color2=="r" and classe=="DB":
		plot(DBlogg700['r_iphas']-DBlogg700['i_iphas']+reddening[0],DBlogg700['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg750['r_iphas']-DBlogg750['i_iphas']+reddening[0],DBlogg750['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg800['r_iphas']-DBlogg800['i_iphas']+reddening[0],DBlogg800['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg850['r_iphas']-DBlogg850['i_iphas']+reddening[0],DBlogg850['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)
		plot(DBlogg900['r_iphas']-DBlogg900['i_iphas']+reddening[0],DBlogg900['r_iphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3,ls=ls)	

	if color1=="r-i" and color2=="r-Ha" and classe=="MS":
		plot(MS[:,0]+reddening[0],MS[:,1]+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	if color1=="r-i" and color2=="r-Ha" and classe=="A0":
		plot(xA0,yA0, color=color,zorder=5,lw=0.3,ls='dashed')

def uvex_intrinsic(color1,color2,classe,reddening=[0,0],dist=10,color='red',ls='solid'):
	""" Function plotting MS,DA sequences in the VPHAS+ color-color diagrams 
	reddening = [E(color1),E(color2)]"""
	KDA = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_koester_uvex')
	KDB = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/db_koester_uvex')
	BDA = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_bergeron_uvex')
	if color1=="r-i" and color2=="r-Ha" and classe=="DA":
		plot(KDA['ri']+reddening[0],KDA['rha']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(BDA['ri']+reddening[0],BDA['rha']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="g-r" and color2=="u-g" and classe=="DA":
		plot(KDA['gr']+reddening[0],KDA['ug']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(BDA['gr']+reddening[0],BDA['ug']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="g-r" and color2=="He-r" and classe=="DA":
		plot(KDA['gr']+reddening[0],KDA['her']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
		plot(BDA['gr']+reddening[0],BDA['her']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)        
	
	elif color1=="r-i" and color2=="r-Ha" and classe=="DB":
		plot(KDB['ri']+reddening[0],KDB['rha']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="g-r" and color2=="u-g" and classe=="DA":
		plot(KDB['gr']+reddening[0],KDB['ug']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	elif color1=="g-r" and color2=="He-r" and classe=="DA":
		plot(KDB['gr']+reddening[0],KDB['her']+reddening[1], color=color,zorder=5,lw=0.3,ls=ls)
	


def atlas_intrinsic(color1,color2,classe,reddening=[0,0],dist=10,color='red'):
	""" Function plotting MS,DA sequences in the ATLAS color-color diagrams 
	reddening = [E(color1),E(color2)]"""
	DA = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_koester10_magnitudes.dat')
	logg700 = DA[where(DA['logg']==7.00)]
	logg750 =  DA[where(DA['logg']==7.50)]
	logg800 =  DA[where(DA['logg']==8.00)]
	logg850 = DA[where(DA['logg']==8.50)] 
	logg900 =  DA[where(DA['logg']==9.00)]
	
	if color1=="g-r" and color2=="u-g" and classe=="DA":
		plot(logg700['g_vphas']-logg700['r_vphas']+reddening[0],logg700['u_vphas']-logg700['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg750['g_vphas']-logg750['r_vphas']+reddening[0],logg750['u_vphas']-logg750['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg800['g_vphas']-logg800['r_vphas']+reddening[0],logg800['u_vphas']-logg800['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg850['g_vphas']-logg850['r_vphas']+reddening[0],logg850['u_vphas']-logg850['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg900['g_vphas']-logg900['r_vphas']+reddening[0],logg900['u_vphas']-logg900['g_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
	elif color1=="r-i" and color2=="g-r" and classe=="DA":
		plot(logg700['r_vphas']-logg700['i_vphas']+reddening[0],logg700['g_vphas']-logg700['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg750['r_vphas']-logg750['i_vphas']+reddening[0],logg750['g_vphas']-logg750['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg800['r_vphas']-logg800['i_vphas']+reddening[0],logg800['g_vphas']-logg800['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg850['r_vphas']-logg850['i_vphas']+reddening[0],logg850['g_vphas']-logg850['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg900['r_vphas']-logg900['i_vphas']+reddening[0],logg900['g_vphas']-logg900['r_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
	elif color1=="i-z" and color2=="r-i" and classe=="DA":
		plot(logg700['i_vphas']-logg700['z_atlas']+reddening[0],logg700['r_vphas']-logg700['i_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg750['i_vphas']-logg750['z_atlas']+reddening[0],logg750['r_vphas']-logg750['i_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg800['i_vphas']-logg800['z_atlas']+reddening[0],logg800['r_vphas']-logg800['i_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg850['i_vphas']-logg850['z_atlas']+reddening[0],logg850['r_vphas']-logg850['i_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
		plot(logg900['i_vphas']-logg900['z_atlas']+reddening[0],logg900['r_vphas']-logg900['i_vphas']+reddening[1], color=color,zorder=5,lw=0.3)
	elif color1=="g-r" and color2=="g" and classe=="DA":
		plot(logg700['g_vphas']-logg700['r_vphas']+reddening[0],logg700['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg750['g_vphas']-logg750['r_vphas']+reddening[0],logg750['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg800['g_vphas']-logg800['r_vphas']+reddening[0],logg800['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg850['g_vphas']-logg850['r_vphas']+reddening[0],logg850['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg900['g_vphas']-logg900['r_vphas']+reddening[0],logg900['g_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
	elif color1=="r-i" and color2=="r" and classe=="DA":
		plot(logg700['r_vphas']-logg700['i_vphas']+reddening[0],logg700['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg750['r_vphas']-logg750['i_vphas']+reddening[0],logg750['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg800['r_vphas']-logg800['i_vphas']+reddening[0],logg800['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg850['r_vphas']-logg850['i_vphas']+reddening[0],logg850['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
		plot(logg900['r_vphas']-logg900['i_vphas']+reddening[0],logg900['r_vphas']+reddening[1]+5*log10(dist)-5, color=color,zorder=5,lw=0.3)
	elif color1=="g-r" and color2=="Hg" and classe=="DA":
		plot(logg800['g_vphas']-logg800['r_vphas'],logg800['r_vphas']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_vphas']-logg800['r_vphas'],logg800['r_vphas']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_vphas']-logg800['r_vphas'],logg800['r_vphas']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	elif color1=="g-i" and color2=="Hg" and classe=="DA":
		plot(logg800['g_vphas']-logg800['i_vphas'],logg800['r_vphas']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		#plot(logg800['g_vphas']-logg800['i_vphas'],logg800['r_vphas']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		#plot(logg800['g_vphas']-logg800['i_vphas'],logg800['r_vphas']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	elif color1=="g-z" and color2=="Hg" and classe=="DA":
		plot(logg800['g_vphas']-logg800['z_atlas'],logg800['r_vphas']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_vphas']-logg800['z_atlas'],logg800['r_vphas']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_vphas']-logg800['z_atlas'],logg800['r_vphas']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	

def intrinsic_colours(color1,color2,classe,diagram='ccd',reddening=[0,0],dist=10,color='red'):
    """Colours are defined as lists, i.e.:
    ['r_apass','i_apass'], or ['g_apass'], 
    or reduced proper motion: mag + 5 +5*log10(Vt)-3.38
    Available bands: galex, sdss, apass, vphas, 
                     atlas, iphas, 2mass, ukidss, wise
    classe: DA, DB, MS,III
    diagram: ccd, cmd, rpm
    reddening: [ebv1,ebv2], or [ebv1,A_lambda]
    dist: pc    
    """
    if classe=="DA":
        M = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_koester10_magnitudes.dat')
	M = M[where(M['teff']<100000.)]
	logg700 = M[where(M['logg']==7.00)]
	logg750 = M[where(M['logg']==7.50)]
	logg800 = M[where(M['logg']==8.00)]
	logg850 = M[where(M['logg']==8.50)] 
	logg900 = M[where(M['logg']==9.00)]        
	if diagram=="ccd":
	    plot(logg700[color1[0]]-logg700[color1[1]]+reddening[0],logg700[color2[0]]-logg700[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)
 	    plot(logg750[color1[0]]-logg750[color1[1]]+reddening[0],logg750[color2[0]]-logg750[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)   
 	    plot(logg800[color1[0]]-logg800[color1[1]]+reddening[0],logg800[color2[0]]-logg800[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)   
	    plot(logg850[color1[0]]-logg850[color1[1]]+reddening[0],logg850[color2[0]]-logg850[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)    
	    plot(logg900[color1[0]]-logg900[color1[1]]+reddening[0],logg900[color2[0]]-logg900[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)    
        elif diagram=="cmd":
	    plot(logg700[color1[0]]-logg700[color1[1]]+reddening[0],logg700[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)
 	    plot(logg750[color1[0]]-logg750[color1[1]]+reddening[0],logg750[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)   
 	    plot(logg800[color1[0]]-logg800[color1[1]]+reddening[0],logg800[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)   
	    plot(logg850[color1[0]]-logg850[color1[1]]+reddening[0],logg850[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)    
	    plot(logg900[color1[0]]-logg900[color1[1]]+reddening[0],logg900[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)    
        elif diagram=="rpm":
 	    plot(logg800[color1[0]]-logg800[color1[1]],logg800[color2[0]]+5*log10(20)-3.38, color=color,zorder=2,lw=0.3)   
 	    plot(logg800[color1[0]]-logg800[color1[1]],logg800[color2[0]]+5*log10(40)-3.38, color=color,zorder=2,lw=0.3)   
 	    plot(logg800[color1[0]]-logg800[color1[1]],logg800[color2[0]]+5*log10(150)-3.38, color=color,zorder=2,lw=0.3)     
    elif classe=="DB":
        M = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/db_koester10_magnitudes.dat')
	logg700 = M[where(M['logg']==7.00)]
	logg750 = M[where(M['logg']==7.50)]
	logg800 = M[where(M['logg']==8.00)]
	logg850 = M[where(M['logg']==8.50)] 
        if diagram=="ccd":
	    plot(logg700[color1[0]]-logg700[color1[1]]+reddening[0],logg700[color2[0]]-logg700[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)
 	    plot(logg750[color1[0]]-logg750[color1[1]]+reddening[0],logg750[color2[0]]-logg750[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)   
 	    plot(logg800[color1[0]]-logg800[color1[1]]+reddening[0],logg800[color2[0]]-logg800[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)   
	    plot(logg850[color1[0]]-logg850[color1[1]]+reddening[0],logg850[color2[0]]-logg850[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)    
        elif diagram=="cmd":
	    plot(logg700[color1[0]]-logg700[color1[1]]+reddening[0],logg700[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)
 	    plot(logg750[color1[0]]-logg750[color1[1]]+reddening[0],logg750[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)   
 	    plot(logg800[color1[0]]-logg800[color1[1]]+reddening[0],logg800[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)   
	    plot(logg850[color1[0]]-logg850[color1[1]]+reddening[0],logg850[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)    
        elif diagram=="rpm":
 	    plot(logg800[color1[0]]-logg800[color1[1]],logg800[color2[0]]+5*log10(20)-3.38, color=color,zorder=2,lw=0.3)
 	    plot(logg800[color1[0]]-logg800[color1[1]],logg800[color2[0]]+5*log10(40)-3.38, color=color,zorder=2,lw=0.3)
 	    plot(logg800[color1[0]]-logg800[color1[1]],logg800[color2[0]]+5*log10(150)-3.38, color=color,zorder=2,lw=0.3)
    elif classe=="MS":
        M = load_intrinsic_magnitudes()
        if diagram=="ccd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]-M[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)
        elif diagram=="cmd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)
    elif classe=="IV":
        M = load_intrinsic_magnitudes(classe="IV")
        if diagram=="ccd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]-M[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)
        elif diagram=="cmd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)
    elif classe=="III":
        M = load_intrinsic_magnitudes(classe="III")
        if diagram=="ccd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]-M[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)
        elif diagram=="cmd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)

    elif classe=="I":
        M = load_intrinsic_magnitudes(classe="I")
        if diagram=="ccd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]-M[color2[1]]+reddening[1], color=color,zorder=2,lw=0.3)
        elif diagram=="cmd":
	    plot(M[color1[0]]-M[color1[1]]+reddening[0],M[color2[0]]+reddening[1]+5*log10(dist)-5, color=color,zorder=2,lw=0.3)



def galex_apass_intrinsic(color1,color2,classe,color='red'):
	""" Function plotting MS,DA sequences in the GALEX/APASS/SDSS color-color diagrams 
	classe: MS, DA
	"""
	
	#logg650 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_650_Av0.txt')
	DA = csv2rec('/storage/astro2/phsmav/data2/models/cooling_models/da_koester10_magnitudes.dat')
	logg700 = DA[where(DA['logg']==7.00)]
	logg750 = DA[where(DA['logg']==7.50)]
	logg800 = DA[where(DA['logg']==8.00)]
	logg850 = DA[where(DA['logg']==8.50)] 
	logg900 = DA[where(DA['logg']==9.00)]
	
	MS = load_intrinsic_magnitudes()
	#MS = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/MS/MS_Av0.txt',usecols=[1,2,3,4,5,6,7])
	if color1=="nuv-r" and color2=="fuv-nuv" and classe=="DA":
		plot(logg700['nuv']-logg700['r_apass'],logg700['fuv']-logg700['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg750['nuv']-logg750['r_apass'],logg750['fuv']-logg750['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg800['nuv']-logg800['r_apass'],logg800['fuv']-logg800['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg850['nuv']-logg850['r_apass'],logg850['fuv']-logg850['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg900['nuv']-logg900['r_apass'],logg900['fuv']-logg900['nuv'], color=color,zorder=5,lw=0.3)
	elif color1=="nuv-g" and color2=="fuv-nuv" and classe=="DA":
		plot(logg700['nuv']-logg700['g_apass'],logg700['fuv']-logg700['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg750['nuv']-logg750['g_apass'],logg750['fuv']-logg750['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg800['nuv']-logg800['g_apass'],logg800['fuv']-logg800['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg850['nuv']-logg850['g_apass'],logg850['fuv']-logg850['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg900['nuv']-logg900['g_apass'],logg900['fuv']-logg900['nuv'], color=color,zorder=5,lw=0.3)
	elif color1=="g-i" and color2=="nuv-g" and classe=="DA":
		plot(logg700['g_apass']-logg700['i_apass'],logg700['nuv']-logg700['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg750['g_apass']-logg750['i_apass'],logg750['nuv']-logg750['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['i_apass'],logg800['nuv']-logg800['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg850['g_apass']-logg850['i_apass'],logg850['nuv']-logg850['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg900['g_apass']-logg900['i_apass'],logg900['nuv']-logg900['g_apass'], color=color,zorder=5,lw=0.3)
	elif color1=="r-i" and color2=="nuv-r" and classe=="DA":
		plot(logg700['r_apass']-logg700['i_apass'],logg700['nuv']-logg700['r_apass'], color=color,zorder=5,lw=0.3)
		plot(logg750['r_apass']-logg750['i_apass'],logg750['nuv']-logg750['r_apass'], color=color,zorder=5,lw=0.3)
		plot(logg800['r_apass']-logg800['i_apass'],logg800['nuv']-logg800['r_apass'], color=color,zorder=5,lw=0.3)
		plot(logg850['r_apass']-logg850['i_apass'],logg850['nuv']-logg850['r_apass'], color=color,zorder=5,lw=0.3)
		plot(logg900['r_apass']-logg900['i_apass'],logg900['nuv']-logg900['r_apass'], color=color,zorder=5,lw=0.3)
	elif color1=="g-r" and color2=="nuv-g" and classe=="DA":
		plot(logg700['g_apass']-logg700['r_apass'],logg700['nuv']-logg700['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg750['g_apass']-logg750['r_apass'],logg750['nuv']-logg750['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['r_apass'],logg800['nuv']-logg800['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg850['g_apass']-logg850['r_apass'],logg850['nuv']-logg850['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg900['g_apass']-logg900['r_apass'],logg900['nuv']-logg900['g_apass'], color=color,zorder=5,lw=0.3)
	elif color1=="b-v" and color2=="nuv-v" and classe=="DA":
		plot(logg700['b_apass']-logg700['v_apass'],logg700['nuv']-logg700['v_apass'], color=color,zorder=5,lw=0.3)
		plot(logg750['b_apass']-logg750['v_apass'],logg750['nuv']-logg750['v_apass'], color=color,zorder=5,lw=0.3)
		plot(logg800['b_apass']-logg800['v_apass'],logg800['nuv']-logg800['v_apass'], color=color,zorder=5,lw=0.3)
		plot(logg850['b_apass']-logg850['v_apass'],logg850['nuv']-logg850['v_apass'], color=color,zorder=5,lw=0.3)
		plot(logg900['b_apass']-logg900['v_apass'],logg900['nuv']-logg900['v_apass'], color=color,zorder=5,lw=0.3)
	elif color1=="nuv-v" and color2=="fuv-nuv" and classe=="DA":
		plot(logg700['nuv']-logg700['v_apass'],logg700['fuv']-logg700['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg750['nuv']-logg750['v_apass'],logg750['fuv']-logg750['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg800['nuv']-logg800['v_apass'],logg800['fuv']-logg800['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg850['nuv']-logg850['v_apass'],logg850['fuv']-logg850['nuv'], color=color,zorder=5,lw=0.3)
		plot(logg900['nuv']-logg900['v_apass'],logg900['fuv']-logg900['nuv'], color=color,zorder=5,lw=0.3)	
	elif color1=="g-h" and color2=="nuv-g" and classe=="DA":
		plot(logg700['g_apass']-logg700['h_2mass'],logg700['nuv']-logg700['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg750['g_apass']-logg750['h_2mass'],logg750['nuv']-logg750['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['h_2mass'],logg800['nuv']-logg800['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg850['g_apass']-logg850['h_2mass'],logg850['nuv']-logg850['g_apass'], color=color,zorder=5,lw=0.3)
		plot(logg900['g_apass']-logg900['h_2mass'],logg900['nuv']-logg900['g_apass'], color=color,zorder=5,lw=0.3)	
	# reduced proper motions #
	elif color1=="nuv-g" and color2=="Hg" and classe=="DA":
		plot(logg800['nuv']-logg800['g_apass'],logg800['g_apass']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['nuv']-logg800['g_apass'],logg800['g_apass']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['nuv']-logg800['g_apass'],logg800['g_apass']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	elif color1=="nuv-r" and color2=="Hr" and classe=="DA":
		plot(logg800['nuv']-logg800['r_apass'],logg800['g_apass']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['nuv']-logg800['r_apass'],logg800['g_apass']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['nuv']-logg800['r_apass'],logg800['g_apass']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	elif color1=="g-r" and color2=="Hg" and classe=="DA":
		plot(logg800['g_apass']-logg800['r_apass'],logg800['g_apass']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['r_apass'],logg800['g_apass']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['r_apass'],logg800['g_apass']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	elif color1=="g-i" and color2=="Hg" and classe=="DA":
		plot(logg800['g_apass']-logg800['i_apass'],logg800['g_apass']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['i_apass'],logg800['g_apass']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['i_apass'],logg800['g_apass']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	elif color1=="g-h" and color2=="Hg" and classe=="DA":
		plot(logg800['g_apass']-logg800['h_2mass'],logg800['g_apass']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['h_2mass'],logg800['g_apass']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['g_apass']-logg800['h_2mass'],logg800['g_apass']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)
	elif color1=="b-v" and color2=="Hv" and classe=="DA":
		plot(logg800['b_apass']-logg800['v_apass'],logg800['v_apass']+ 5*log10(20)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['b_apass']-logg800['v_apass'],logg800['v_apass']+ 5*log10(40)-3.38, color=color,zorder=5,lw=0.3)
		plot(logg800['b_apass']-logg800['v_apass'],logg800['v_apass']+5*log10(150)-3.38, color=color,zorder=5,lw=0.3)



        
	elif color1=="nuv-r" and color2=="fuv-nuv" and classe=="ZZ":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(ZZlogg700[:,3],ZZlogg700[:,1], color=color,zorder=5)
		plot(ZZlogg750[:,3],ZZlogg750[:,1], color=color,zorder=5)
		plot(ZZlogg800[:,3],ZZlogg800[:,1], color=color,zorder=5)
		plot(ZZlogg850[:,3],ZZlogg850[:,1], color=color,zorder=5)
		plot(ZZlogg900[:,3],ZZlogg900[:,1], color=color,zorder=5)
	elif color1=="nuv-g" and color2=="fuv-nuv" and classe=="ZZ":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(ZZlogg700[:,2],ZZlogg700[:,1], color=color,zorder=5)
		plot(ZZlogg750[:,2],ZZlogg750[:,1], color=color,zorder=5)
		plot(ZZlogg800[:,2],ZZlogg800[:,1], color=color,zorder=5)
		plot(ZZlogg850[:,2],ZZlogg850[:,1], color=color,zorder=5)
		plot(ZZlogg900[:,2],ZZlogg900[:,1], color=color,zorder=5)
	elif color1=="g-i" and color2=="nuv-g" and classe=="ZZ":
		#plot(logg650[:,4],logg650[:,2], color=color,zorder=5)
		plot(ZZlogg700[:,4],ZZlogg700[:,2], color=color,zorder=5)
		plot(ZZlogg750[:,4],ZZlogg750[:,2], color=color,zorder=5)
		plot(ZZlogg800[:,4],ZZlogg800[:,2], color=color,zorder=5)
		plot(ZZlogg850[:,4],ZZlogg850[:,2], color=color,zorder=5)
		plot(ZZlogg900[:,4],ZZlogg900[:,2], color=color,zorder=5)
	elif color1=="r-i" and color2=="nuv-r" and classe=="ZZ":
		#plot(ZZlogg650[:,4],ZZlogg650[:,2], color=color,zorder=5)
		plot(ZZlogg700[:,6],ZZlogg700[:,3], color=color,zorder=5)
		plot(ZZlogg750[:,6],ZZlogg750[:,3], color=color,zorder=5)
		plot(ZZlogg800[:,6],ZZlogg800[:,3], color=color,zorder=5)
		plot(ZZlogg850[:,6],ZZlogg850[:,3], color=color,zorder=5)
		plot(ZZlogg900[:,6],ZZlogg900[:,3], color=color,zorder=5)
	elif color1=="g-r" and color2=="nuv-g" and classe=="ZZ":
		#plot(logg650[:,5],logg650[:,2], color=color,zorder=5)
		plot(ZZlogg700[:,5],ZZlogg700[:,2], color=color,zorder=5)
		plot(ZZlogg750[:,5],ZZlogg750[:,2], color=color,zorder=5)
		plot(ZZlogg800[:,5],ZZlogg800[:,2], color=color,zorder=5)
		plot(ZZlogg850[:,5],ZZlogg850[:,2], color=color,zorder=5)
		plot(ZZlogg900[:,5],ZZlogg900[:,2], color=color,zorder=5)


	# main sequence #
	elif color1=="nuv-r" and color2=="fuv-nuv" and classe=="MS":
		plot(MS['nuv']-MS['r_apass'],MS['fuv']-MS['nuv'], color=color,zorder=5, ls='dashed',lw=0.3)
	elif color1=="nuv-g" and color2=="fuv-nuv" and classe=="MS":
		plot(MS['nuv']-MS['r_apass'],MS['fuv']-MS['nuv'], color=color,zorder=5, ls='dashed',lw=0.3)
	elif color1=="g-r" and color2=="nuv-g" and classe=="MS":
		plot(MS['g_apass']-MS['r_apass'],MS['nuv']-MS['g_apass'], color=color,zorder=5, ls='dashed',lw=0.3)	
	elif color1=="g-i" and color2=="nuv-g" and classe=="MS":
		plot(MS['g_apass']-MS['i_apass'],MS['nuv']-MS['g_apass'], color=color,zorder=5, ls='dashed',lw=0.3)
	elif color1=="nuv-v" and color2=="fuv-nuv" and classe=="MS":
		plot(MS['nuv']-MS['v_apass'],MS['fuv']-MS['nuv'], color=color,zorder=5, ls='dashed',lw=0.3)
	elif color1=="b-v" and color2=="nuv-v" and classe=="MS":
		plot(MS['b_apass']-MS['v_apass'],MS['nuv']-MS['v_apass'], color=color,zorder=5, ls='dashed',lw=0.3)	
	elif color1=="g-h" and color2=="nuv-g" and classe=="MS":
		plot(MS['g_apass']-MS['h_2mass'],MS['nuv']-MS['g_apass'], color=color,zorder=5, ls='dashed',lw=0.3)	




def galex_sdss_intrinsic(color1,color2,classe,color='red'):
	""" Function plotting MS,DA sequences in the GALEX/APASS/SDSS color-color diagrams """
	#logg650 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_650_Av0.txt')
	logg700 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_700_Av0.txt')
	logg750 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_750_Av0.txt')
	logg800 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_800_Av0.txt')
	logg850 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_850_Av0.txt')
	logg900 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/DA/log_g_900_Av0.txt')
	MS = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/MS/MS_Av0.txt',usecols=[1,2,3,4,5,6,7])
	ZZlogg700 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/ZZCeti/log_g_7.00_Av0.txt')
	ZZlogg750 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/ZZCeti/log_g_7.50_Av0.txt')
	ZZlogg800 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/ZZCeti/log_g_8.00_Av0.txt')
	ZZlogg850 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/ZZCeti/log_g_8.50_Av0.txt')
	ZZlogg900 = loadtxt('/storage/astro2/phsmav/data2/galex/tracks/ZZCeti/log_g_9.00_Av0.txt')	
	if color1=="nuv-r" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,3],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,3],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,3],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,3],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,3],logg900[:,1], color=color,zorder=5)
	elif color1=="g-r" and color2=="u-g" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,5],logg700[:,7], color=color,zorder=5)
		plot(logg750[:,5],logg750[:,7], color=color,zorder=5)
		plot(logg800[:,5],logg800[:,7], color=color,zorder=5)
		plot(logg850[:,5],logg850[:,7], color=color,zorder=5)
		plot(logg900[:,5],logg900[:,7], color=color,zorder=5)
	elif color1=="nuv-g" and color2=="fuv-nuv" and classe=="WD":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(logg700[:,2],logg700[:,1], color=color,zorder=5)
		plot(logg750[:,2],logg750[:,1], color=color,zorder=5)
		plot(logg800[:,2],logg800[:,1], color=color,zorder=5)
		plot(logg850[:,2],logg850[:,1], color=color,zorder=5)
		plot(logg900[:,2],logg900[:,1], color=color,zorder=5)
	elif color1=="g-i" and color2=="nuv-g" and classe=="WD":
		#plot(logg650[:,4],logg650[:,2], color=color,zorder=5)
		plot(logg700[:,4],logg700[:,2], color=color,zorder=5)
		plot(logg750[:,4],logg750[:,2], color=color,zorder=5)
		plot(logg800[:,4],logg800[:,2], color=color,zorder=5)
		plot(logg850[:,4],logg850[:,2], color=color,zorder=5)
		plot(logg900[:,4],logg900[:,2], color=color,zorder=5)
	elif color1=="g-r" and color2=="nuv-g" and classe=="WD":
		#plot(logg650[:,5],logg650[:,2], color=color,zorder=5)
		plot(logg700[:,5],logg700[:,2], color=color,zorder=5)
		plot(logg750[:,5],logg750[:,2], color=color,zorder=5)
		plot(logg800[:,5],logg800[:,2], color=color,zorder=5)
		plot(logg850[:,5],logg850[:,2], color=color,zorder=5)
		plot(logg900[:,5],logg900[:,2], color=color,zorder=5)
	if color1=="nuv-r" and color2=="fuv-nuv" and classe=="ZZ":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(ZZlogg700[:,3],ZZlogg700[:,1], color=color,zorder=5)
		plot(ZZlogg750[:,3],ZZlogg750[:,1], color=color,zorder=5)
		plot(ZZlogg800[:,3],ZZlogg800[:,1], color=color,zorder=5)
		plot(ZZlogg850[:,3],ZZlogg850[:,1], color=color,zorder=5)
		plot(ZZlogg900[:,3],ZZlogg900[:,1], color=color,zorder=5)
	elif color1=="nuv-g" and color2=="fuv-nuv" and classe=="ZZ":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(ZZlogg700[:,2],ZZlogg700[:,1], color=color,zorder=5)
		plot(ZZlogg750[:,2],ZZlogg750[:,1], color=color,zorder=5)
		plot(ZZlogg800[:,2],ZZlogg800[:,1], color=color,zorder=5)
		plot(ZZlogg850[:,2],ZZlogg850[:,1], color=color,zorder=5)
		plot(ZZlogg900[:,2],ZZlogg900[:,1], color=color,zorder=5)
	elif color1=="g-r" and color2=="u-g" and classe=="ZZ":
		#plot(logg650[:,3],logg650[:,1], color=color,zorder=5)
		plot(ZZlogg700[:,5],ZZlogg700[:,7], color=color,zorder=5)
		plot(ZZlogg750[:,5],ZZlogg750[:,7], color=color,zorder=5)
		plot(ZZlogg800[:,5],ZZlogg800[:,7], color=color,zorder=5)
		plot(ZZlogg850[:,5],ZZlogg850[:,7], color=color,zorder=5)
		plot(ZZlogg900[:,5],ZZlogg900[:,7], color=color,zorder=5)
	elif color1=="g-i" and color2=="nuv-g" and classe=="ZZ":
		#plot(logg650[:,4],logg650[:,2], color=color,zorder=5)
		plot(ZZlogg700[:,4],ZZlogg700[:,2], color=color,zorder=5)
		plot(ZZlogg750[:,4],ZZlogg750[:,2], color=color,zorder=5)
		plot(ZZlogg800[:,4],ZZlogg800[:,2], color=color,zorder=5)
		plot(ZZlogg850[:,4],ZZlogg850[:,2], color=color,zorder=5)
		plot(ZZlogg900[:,4],ZZlogg900[:,2], color=color,zorder=5)
	elif color1=="g-r" and color2=="nuv-g" and classe=="ZZ":
		#plot(logg650[:,5],logg650[:,2], color=color,zorder=5)
		plot(ZZlogg700[:,5],ZZlogg700[:,2], color=color,zorder=5)
		plot(ZZlogg750[:,5],ZZlogg750[:,2], color=color,zorder=5)
		plot(ZZlogg800[:,5],ZZlogg800[:,2], color=color,zorder=5)
		plot(ZZlogg850[:,5],ZZlogg850[:,2], color=color,zorder=5)
		plot(ZZlogg900[:,5],ZZlogg900[:,2], color=color,zorder=5)

def ZZCeti_box(plane="ugr"):
	"""ZZ Ceti selection box in SDSS colou-colour-plane"""
	p1 = polyfit([-0.04,-0.15],[-0.7,0.8],deg=1)
	p2 = polyfit([-0.04,-0.27],[-0.7,0.8],deg=1)
	x = arange(-0.27,-0.05,0.01)
	y1 = polyval(p1,x)
	y2 = polyval(p2,x)
	plot(x[y1<=0.8],y1[y1<=0.8], ls='dashed',color='m')
	plot(x[y2<=0.8],y2[y2<=0.8], ls='dashed',color='m')
