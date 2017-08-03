#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 14:06:47 2017

@author: sanderson

Downloads Needed:
    Illustris_Python
    specified snapshot number Groupcat download(FOF & Subfind)
    
Files Needed:
    'Galaxies HI det table.tex' 
    'ids for groups with 4-5 optically detectable galaxies.tex'
    'ids for groups with 6-9 optically detectable galaxies.tex'
    'ids for groups with 10-19 optically detectable galaxies.tex'
    'ids for groups with 20-30 optically detectable galaxies.tex'
    (all created from the 'HI and Optic det table.py' script)
    
Result:
   makes the HOD(Halo Occupation Distribution) plot for all groups using luminosity as proxy for group mass 
"""
#===================================================================================================================
from astropy.table import Table, Column
import numpy as np
from datetime import datetime
import illustris_python as il
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
import random

startTime = datetime.now()


basePath = './Illustris-1_snap132' #directory to illustris_python and illustris files
fields1 = ['SubhaloPos', 'SubhaloGrNr', 'SubhaloSFR', 'SubhaloMass', 'SubhaloStellarPhotometrics'] #fields to be requested
subhalos = il.groupcat.loadSubhalos(basePath,132,fields=fields1) #Call illustris_python scrips to request specified fields
fields2 = ['GroupFirstSub', 'GroupNsubs', 'GroupPos', 'Group_M_Crit200', 'Group_R_Crit200', 'GroupMass'] #fields to be requested
halos = il.groupcat.loadHalos(basePath,132,fields=fields2) #Call illustris_python scrips to request specified fields



halo_id_table_45 = Table.read('ids for groups with 4-5 optically detectable galaxies_mag17.tex', format = 'latex')
halo_id_table_69 = Table.read('ids for groups with 6-9 optically detectable galaxies_mag17.tex', format = 'latex')
halo_id_table_1019 = Table.read('ids for groups with 10-19 optically detectable galaxies_mag17.tex', format = 'latex')
halo_id_table_2030 = Table.read('ids for groups with 20-30 optically detectable galaxies_mag17.tex', format = 'latex')

HI_det_table = Table.read('Galaxies HI det table.tex', format = 'latex')
galaxies_HI_det_table_1e10= Table.read('Galaxies HI det table HI mass fraction_1e10 lower limit.tex', format = 'latex')
galaxies_HI_det_table_1e95= Table.read('Galaxies HI det table HI mass fraction_1e9.5 lower limit.tex', format = 'latex')

print("minimum HI mass =" + str(min(galaxies_HI_det_table_1e10['HI Mass 1e10'])) + " (1e10)solar masses")

IDs_45 = halo_id_table_45['IDs']
IDs_69 = halo_id_table_69['IDs']
IDs_1019 = halo_id_table_1019['IDs']
IDs_2030 = halo_id_table_2030['IDs']

print(len(IDs_45))
print(len(IDs_69))
print(len(IDs_1019))
print(len(IDs_2030))

IDs = []
print(type(IDs))
for h in IDs_45:
    IDs.append(h)
for i in IDs_69:
    IDs.append(i)
for j in IDs_1019:
    IDs.append(j)   
for k in IDs_2030:
    IDs.append(k)

print(len(IDs))



HI_det = []
groupMasses = []
OpticCount = []
num_HI_det = []
HIcount = []
GroupMasses = []
group_magnitudes = []
#group_luminosities = []
optic_det_ids = []

#get list of optical detections and HI detections
for h in IDs:
    opticDet = []
    group_lum = []
    numOpticDet = 0
    numHIdet = 0
    groupMass = halos['GroupMass'][h]
    subhalo_idlist = np.arange(halos['GroupFirstSub'][h], halos['GroupFirstSub'][h] + halos['GroupNsubs'][h]) 
    #list of subhalos in that group
    for s in subhalo_idlist:
        mag = subhalos['SubhaloStellarPhotometrics'][s][5] #mag for galaxy in h
        if (mag < -18): #optical detection test
            numOpticDet+=1
            opticDet.append(mag)
            optic_det_ids.append(s)
            lum = (3.86e26)*(10**((4.77-mag)/2.5))
            group_lum.append(lum)
    group_luminosities = sum(group_lum)
#    group_lum = (3.86e26)*(10**((4.77-group_mag)/2.5))
    group_mass = (350/(3.86e26))*group_luminosities
    print(group_mass)
    OpticCount.append(numOpticDet)
    GroupMasses.append(group_mass)
    HI_det_for_h = np.where(galaxies_HI_det_table_1e10['GrNr'] == h)
    for gal in HI_det_for_h[0]:
        numHIdet+=1
        HI_det.append(galaxies_HI_det_table_1e10['Galaxy IDs'][gal])
    HIcount.append(numHIdet)

total_subs = len(optic_det_ids) + len(HI_det)
total_optics = sum(OpticCount)
total_H1 = sum(HIcount)
for i in optic_det_ids:
    if i in HI_det:
        total_subs-= 1
        total_optics-=1
for i in HI_det:
    if i in optic_det_ids:
        total_H1-=1
group_mass_logmsun = np.log10(GroupMasses)
galaxy_counts_optic = np.log10(OpticCount)
galaxy_counts_HI = np.log10(HIcount)

        

fig = plt.figure()
plt.plot(group_mass_logmsun, galaxy_counts_optic, '+', label = 'num Optical detected galaxies')
plt.plot(group_mass_logmsun, galaxy_counts_HI, 'd', fillstyle = 'none', label = 'num HI detected galaxies')

#plt.plot(groupMasses, opticDet, '+', label = 'num Optical detected galaxies')
#plt.plot(groupMasses, HIdet, 'd', fillstyle = 'none', label = 'num HI detected galaxies')
ax=fig.add_subplot(1,1,1)
plt.xlabel("Halo Mass log(M/[$M_\odot$])")
plt.ylabel("Number of Galaxies log(N)")


#plt.xlabel("Halo Mass[$1e10[$M_\odot$]/h$]")
#plt.ylabel("Number of Galaxies")
plt.title("Halo Occupation Number (using luminosity as proxy for group mass)")
plt.legend()
plt.show()
print("the total number of groups is: " + str(len(GroupMasses)))
print("the number of HI only detected galaxies is: " + str(total_H1))  
print("the number of optical only detected galaxies is: " + str(total_optics))

print(datetime.now()-startTime)
    

