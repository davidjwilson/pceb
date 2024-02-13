#!/bin/bash

#spectra

evselect table=m1_S002_ImagingEvts_events_gtifiltered.fit withspectrumset=yes spectrumset=EMOS1source_spectrum.fits energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN circle(25778.605,23870.203,600))'

evselect table=m2_S003_ImagingEvts_events_gtifiltered.fit withspectrumset=yes spectrumset=EMOS2source_spectrum.fits energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN circle(25778.605,23870.203,600))'

evselect table=pn_S001_ImagingEvts_events_gtifiltered.fit withspectrumset=yes  spectrumset=EPNsource_spectrum.fits energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN circle(25778.605,23870.203,600))'
#acceptchanrange=yes
#background

evselect table=m1_S002_ImagingEvts_events_gtifiltered.fit withspectrumset=yes spectrumset=EMOS1background_spectrum.fits energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN circle(23481.4,24343.3,1166.12))'

evselect table=m2_S003_ImagingEvts_events_gtifiltered.fit withspectrumset=yes spectrumset=EMOS2background_spectrum.fits energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=11999 expression='#XMMEA_EM && (PATTERN<=12) && ((X,Y) IN circle(23481.4,24343.3,1166.12))'

evselect table=pn_S001_ImagingEvts_events_gtifiltered.fit withspectrumset=yes spectrumset=EPNbackground_spectrum.fits energycolumn=PI spectralbinsize=5 withspecranges=yes specchannelmin=0 specchannelmax=20479 expression='(FLAG==0) && (PATTERN<=4) && ((X,Y) IN circle(23481.4,24343.3,1166.12))'

#backscale

backscale spectrumset=EMOS1source_spectrum.fits badpixlocation=m1_S002_ImagingEvts_events_gtifiltered.fit

backscale spectrumset=EMOS1background_spectrum.fits badpixlocation=m1_S002_ImagingEvts_events_gtifiltered.fit

backscale spectrumset=EMOS2source_spectrum.fits badpixlocation=m2_S003_ImagingEvts_events_gtifiltered.fit

backscale spectrumset=EMOS2background_spectrum.fits badpixlocation=m2_S003_ImagingEvts_events_gtifiltered.fit

backscale spectrumset=EPNsource_spectrum.fits badpixlocation=pn_S001_ImagingEvts_events_gtifiltered.fit

backscale spectrumset=EPNbackground_spectrum.fits badpixlocation=pn_S001_ImagingEvts_events_gtifiltered.fit

#RMF
rmfgen spectrumset=EMOS1source_spectrum.fits rmfset=EMOS1.rmf withenergybins=yes energymin=0.1 energymax=12.0 nenergybins=2400 

rmfgen spectrumset=EMOS2source_spectrum.fits rmfset=EMOS2.rmf withenergybins=yes energymin=0.1 energymax=12.0 nenergybins=2400 

rmfgen spectrumset=EPNsource_spectrum.fits rmfset=EPN.rmf withenergybins=yes energymin=0.1 energymax=12.0 nenergybins=2400 

#arf

arfgen spectrumset=EMOS1source_spectrum.fits arfset=EMOS1.arf withrmfset=yes rmfset=EMOS1.rmf badpixlocation=m1_S002_ImagingEvts_events_gtifiltered.fit detmaptype=psf

arfgen spectrumset=EMOS2source_spectrum.fits arfset=EMOS2.arf withrmfset=yes rmfset=EMOS2.rmf badpixlocation=m2_S003_ImagingEvts_events_gtifiltered.fit detmaptype=psf

arfgen spectrumset=EPNsource_spectrum.fits arfset=EPN.arf withrmfset=yes rmfset=EPN.rmf badpixlocation=pn_S001_ImagingEvts_events_gtifiltered.fit detmaptype=psf


#combine
epicspeccombine pha="EMOS1source_spectrum.fits EMOS2source_spectrum.fits EPNsource_spectrum.fits" bkg="EMOS1background_spectrum.fits EMOS2background_spectrum.fits EPNbackground_spectrum.fits" rmf="EMOS1.rmf EMOS2.rmf EPN.rmf" arf="EMOS1.arf EMOS2.arf EPN.arf" filepha="src_spectrum_grp.ds" filebkg="bkg_spectrum_grp.ds" filersp="response_grp.rmf" allowHEdiff=yes

#bin spectrum
specgroup spectrumset=src_spectrum_grp.ds groupedset=src_spectrum_bin25.ds mincounts=25 lastbin="setbad"

