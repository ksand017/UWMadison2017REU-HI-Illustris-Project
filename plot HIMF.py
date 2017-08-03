#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 13:11:06 2017

@author: sanderson

Files Needed:
    'Galaxies HI det table with lower limit 1e6.tex'
    (made from 'HI and Optic det table.py' script by changing upper and lower limits for HI detection)
    
Result:
    Plots the HIMF (Nuetral Hydrogen Mass Function) for all HI detections in specfied Illustris suite and snapshot
"""

#===================================================================================================================
from astropy.table import Table, Column
import numpy as np
from datetime import datetime
import illustris_python as il
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import plotly.plotly as py
import math
import random

startTime = datetime.now()


galaxy_HI_table = Table.read('Galaxies HI det table with lower limit 1e6.tex', format = 'latex')

galaxy_IDs = galaxy_HI_table['Galaxy IDs']
print(len(galaxy_IDs))
galaxy_HI_Masses = galaxy_HI_table['HI Mass 1e10[$M_\odot$]']
print(len(galaxy_HI_Masses))
galaxies_6 = []
galaxies_625 = []
galaxies_65 = []
galaxies_675 = []
galaxies_7 = []
galaxies_725 = []
galaxies_75 = []
galaxies_775 = []
galaxies_8 = []
galaxies_825 = []
galaxies_85 = []
galaxies_875 = []
galaxies_9 = []
galaxies_925 = []
galaxies_95 =[]
galaxies_975 = []
galaxies_10 = []
galaxies_1025 = []
galaxies_105 =[]
galaxies_1075 = []
galaxies_11 = []
histogram_data = []
masses = []
V_max = ((106.5)**3) 
for i in range(len(galaxy_IDs)):
    g = galaxy_IDs[i]
    h1mass = galaxy_HI_Masses[i]*(1e10)
    masses.append(h1mass)
    if ((h1mass >= (1e6)) & (h1mass < (10**6.25))):
        galaxies_6.append(g)
        histogram_data.append(6)
    elif ((h1mass >= (10**6.25)) & ((h1mass) < (10**6.5))):
        galaxies_625.append(g)
        histogram_data.append(6.25)
    elif ((h1mass >= (10**6.5)) & (h1mass < (10**6.75))):
        galaxies_65.append(g)
        histogram_data.append(6.5)
    elif ((h1mass >= (10**6.75)) & (h1mass < (1e7))):
        galaxies_675.append(g)
        histogram_data.append(6.75)
    elif ((h1mass >= (1e7)) & (h1mass < (10**7.25))):
        galaxies_7.append(g)
        histogram_data.append(7)
    elif ((h1mass >= (10**7.25)) & (h1mass < (10**7.5))):
        galaxies_725.append(g)
        histogram_data.append(7.25)
    elif ((h1mass >= (10**7.5)) & (h1mass < (10**7.75))):
        galaxies_75.append(g)
        histogram_data.append(7.5)
    elif ((h1mass >= (10**7.75)) & (h1mass < (1e8))):
        galaxies_775.append(g)
        histogram_data.append(7.75)
    elif ((h1mass >= (1e8)) & (h1mass < (10**8.25))):
        galaxies_8.append(g)
        histogram_data.append(8)
    elif ((h1mass >= (10**8.25)) & (h1mass < (10**8.5))):
        galaxies_825.append(g)
        histogram_data.append(8.25)
    elif ((h1mass >= (10**8.5)) & (h1mass < (10**8.75))):
        galaxies_85.append(g)
        histogram_data.append(8.5)
    elif ((h1mass >= (10**8.75)) & (h1mass < (1e9))):
        galaxies_875.append(g)
        histogram_data.append(8.75)
    elif ((h1mass >= (1e9)) & (h1mass < (10**9.25))):
        galaxies_9.append(g)
        histogram_data.append(9)
    elif ((h1mass >= (10**9.25)) & (h1mass < (10**9.5))):
        galaxies_925.append(g)
        histogram_data.append(9.25)
    elif ((h1mass >= (10**9.5)) & (h1mass < (10**9.75))):
        galaxies_95.append(g)
        histogram_data.append(9.5)
    elif ((h1mass >= (10**9.75)) & (h1mass < (1e10))):
        galaxies_975.append(g)
        histogram_data.append(9.75)
    elif ((h1mass >= (1e10)) & (h1mass < (10**10.25))):
        galaxies_10.append(g)
        histogram_data.append(10)
    elif ((h1mass >= (10**10.25)) & (h1mass < (10**10.5))):
        galaxies_1025.append(g)
        histogram_data.append(10.25)
    elif ((h1mass >= (10**10.5)) & (h1mass < (10**10.75))):
        galaxies_105.append(g)
        histogram_data.append(10.5)
    elif ((h1mass >= (10**10.75)) & (h1mass < (1e11))):
        galaxies_1075.append(g)
        histogram_data.append(10.75)

        
print("number of 10^6 HI galaxies =" + str(histogram_data.count(6)))
print("number of 10^6.25 HI galaxies =" + str(histogram_data.count(6.25)))
print("number of 10^6.5 HI galaxies =" + str(histogram_data.count(6.5)))
print("number of 10^6.75 HI galaxies =" + str(histogram_data.count(6.75)))
print("number of 10^7 HI galaxies =" + str(histogram_data.count(7)))
print("number of 10^7.25 HI galaxies =" + str(histogram_data.count(7.25)))
print("number of 10^7.5 HI galaxies =" + str(histogram_data.count(7.5)))
print("number of 10^7.75 HI galaxies =" + str(histogram_data.count(7.75)))
print("number of 10^8 HI galaxies =" + str(histogram_data.count(8)))
print("number of 10^8.25 HI galaxies =" + str(histogram_data.count(8.25)))
print("number of 10^8.5 HI galaxies =" + str(histogram_data.count(8.5)))
print("number of 10^8.75 HI galaxies =" + str(histogram_data.count(8.75)))
print("number of 10^9 HI galaxies =" + str(histogram_data.count(9)))
print("number of 10^9.25 HI galaxies =" + str(histogram_data.count(9.25)))
print("number of 10^9.5 HI galaxies =" + str(histogram_data.count(9.5)))
print("number of 10^9.75 HI galaxies =" + str(histogram_data.count(9.75)))
print("number of 10^10 HI galaxies =" + str(histogram_data.count(10)))
print("number of 10^10.25 HI galaxies =" + str(histogram_data.count(10.25)))
print("number of 10^10.5 HI galaxies =" + str(histogram_data.count(10.5)))
print("number of 10^10.75 HI galaxies =" + str(histogram_data.count(10.75)))

#Volume normalized HIMF plot
'''
mass_bins = [6,6.25,6.5,6.75,7,7.25,7.5,7.75,8,8.25,8.5,8.75,9,9.25,9.5,9.75,10,10.25,10.5,10.75]
frequency = [(len(galaxies_6)/V_max), (len(galaxies_625)/V_max), (len(galaxies_65)/V_max), (len(galaxies_675)/V_max), (len(galaxies_7)/V_max), (len(galaxies_725)/V_max), (len(galaxies_75)/V_max), (len(galaxies_775)/V_max), (len(galaxies_8)/V_max), (len(galaxies_825)/V_max), (len(galaxies_85)/V_max), (len(galaxies_875)/V_max), (len(galaxies_9)/V_max), (len(galaxies_925)/V_max), (len(galaxies_95)/V_max), (len(galaxies_975)/V_max), (len(galaxies_10)/V_max), (len(galaxies_1025)/V_max), (len(galaxies_105)/V_max), (len(galaxies_1075)/V_max)]
frequency_log = np.log10(frequency)
fig = plt.figure()
plt.plot(mass_bins, frequency_log, label = 'HIMF')
plt.title("HIMF")
plt.xlabel(" log($M_HI$/$M_\odot$)")
plt.ylabel("# galaxies")
'''

fig1 = plt.figure()
plt.hist(histogram_data, bins = 19)
plt.title("HIMF")
plt.xlabel(" log($M_HI$/$M_\odot$)")
plt.ylabel("# galaxies")

print(datetime.now()-startTime)
