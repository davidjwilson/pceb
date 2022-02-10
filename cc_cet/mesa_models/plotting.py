#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 23:44:43 2020

@author: odette
"""

from nugridpy import mesa as ms
import numpy as np
import matplotlib.pyplot as plt


path = '/home/eduardo/Odette/Other-Works/David_Wilson/PCEB/CC-Cet/MESA/star_plus_point_mass/'

data  = ms.history_data(path, slname = 'binary_history_1.0xGR.data')
data0 = ms.history_data(path, slname = 'binary_history_2.47xGR.data')


porb  = data.get('period_days')*24.
porb0 = data0.get('period_days')*24.

mdot  = data.get('lg_mtransfer_rate')
mdot0 = data0.get('lg_mtransfer_rate')

age  = data.get('age')
age0 = data0.get('age')


# plotting  Age and Porb
plt.plot(age0, porb0,  color = '#F08D1B', label=r'$\dot{J}_{\mathrm{tot}}$=2.47$\times\dot{J}_{\mathrm{GR}}$') 
plt.axhline(1.83479, color='#F08D1B', linestyle='dashed')
plt.plot(age, porb,  color='#1B70F0', label=r'$\dot{J}_{\mathrm{tot}}$=1.00$\times\dot{J}_{\mathrm{GR}}$') 
plt.axhline(1.89239, color='#1B70F0', linestyle='dashed')

# This is the tsd estimated by Schreiber and Ganscke (2003; Table 2) for CC-Cet -- log(tsd)=10.27
plt.axvline(10**(10.27), color='grey', linestyle='dashed')

plt.text(1.2e10, 1.6, '~1.83h', color='#F08D1B')
plt.text(1.2e10, 2.0, '~1.89h', color='#1B70F0')

plt.legend(loc='upper center')
plt.xlabel('age (yr)')
plt.ylabel(r'$P_{\mathrm{orb}}$ (h)')
plt.show()


# plotting Mdot versus Porb
plt.plot(porb0, mdot0, color = '#F08D1B', label=r'$\dot{J}_{\mathrm{tot}}$=2.47$\times\dot{J}_{\mathrm{GR}}$') 
plt.axvline(1.83479, color='#F08D1B', linestyle='dashed')
plt.plot(porb,mdot, color='#1B70F0', label=r'$\dot{J}_{\mathrm{tot}}$=1.00$\times\dot{J}_{\mathrm{GR}}$') 
plt.axvline(1.89239, color='#1B70F0', linestyle='dashed')

plt.legend(loc='upper right')
plt.xlabel('age (yr)')
plt.ylabel(r'$P_{\mathrm{orb}}$ (h)')
plt.show()
