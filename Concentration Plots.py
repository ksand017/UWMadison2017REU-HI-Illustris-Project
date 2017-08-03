#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 13:13:15 2017

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

Results:
    creates plot for the C31 concentration (see Hess & Wilcots 2013) as a function of bin membership
    creates plotof the fraction of HI detections as a function of bin membership
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
from operator import itemgetter, attrgetter

startTime = datetime.now()


basePath = './Illustris-1_snap132' #directory to illustris_python and illustris files
fields1 = ['SubhaloPos', 'SubhaloGrNr', 'SubhaloSFR', 'SubhaloMass', 'SubhaloStellarPhotometrics'] #fields to be requested
subhalos = il.groupcat.loadSubhalos(basePath,132,fields=fields1) #Call illustris_python scrips to request specified fields
fields2 = ['GroupFirstSub', 'GroupNsubs', 'GroupPos', 'Group_M_Crit200', 'Group_R_Crit200'] #fields to be requested
halos = il.groupcat.loadHalos(basePath,132,fields=fields2) #Call illustris_python scrips to request specified fields


group_table = Table.read('groupidsIll1_snap132.tex', format = 'latex')

halo_id_table_45 = Table.read('ids for groups with 4-5 optically detectable galaxies.tex', format = 'latex')
halo_id_table_69 = Table.read('ids for groups with 6-9 optically detectable galaxies.tex', format = 'latex')
halo_id_table_1019 = Table.read('ids for groups with 10-19 optically detectable galaxies.tex', format = 'latex')
halo_id_table_2030 = Table.read('ids for groups with 20-30 optically detectable galaxies.tex', format = 'latex')

#IDs for child galaxies
HI_det_table = Table.read('Galaxies HI det table.tex', format = 'latex')
galaxies_HI_det_table_1e10= Table.read('Galaxies HI det table HI mass fraction_1e10 lower limit.tex', format = 'latex')
galaxies_HI_det_table_1e95= Table.read('Galaxies HI det table HI mass fraction_1e9.5 lower limit.tex', format = 'latex')

#IDs for parent halos
IDs_45 = halo_id_table_45['IDs']
IDs_69 = halo_id_table_69['IDs']
IDs_1019 = halo_id_table_1019['IDs']
IDs_2030 = halo_id_table_2030['IDs']

print("number of optical groups =" + str(len(IDs_2030)))

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
 
length = 75000
half_length = length/2
neg_half_length = -1*length/2

    
#make sorted table for optical detections
def get_sorted_optic_table(IDs, binnumber):
    r_optic_xy = []
    optic_det = []
    for h in IDs:
        subhalo_idlist = np.arange(halos['GroupFirstSub'][h], halos['GroupFirstSub'][h] + halos['GroupNsubs'][h]) 
        #list of subhalos in that group
        grnr = h
        grpos_x = halos['GroupPos'][grnr][0]
        grpos_y = halos['GroupPos'][grnr][1]
        grpos_z = halos['GroupPos'][grnr][2]
        for s in subhalo_idlist:
            mag = subhalos['SubhaloStellarPhotometrics'][s][5] #mag for galaxy in h
            if (mag < -18): #optical detection test
                optic_det.append(s)
    for oid in optic_det:
        
        dx = halos['GroupPos'][grnr][0] - subhalos['SubhaloPos'][oid][0]
        dy = halos['GroupPos'][grnr][1] - subhalos['SubhaloPos'][oid][1]
        dz = halos['GroupPos'][grnr][2] - subhalos['SubhaloPos'][oid][2]
        
        if (dx <= neg_half_length):
            subhalos['SubhaloPos'][oid][0] = subhalos['SubhaloPos'][oid][0] - 75000
        elif (dx >= half_length):
            subhalos['SubhaloPos'][oid][0] = subhalos['SubhaloPos'][oid][0] + 75000
        
        if (dy <= neg_half_length):
            subhalos['SubhaloPos'][oid][1] = subhalos['SubhaloPos'][oid][1] - 75000
        elif (dy >= half_length):
            subhalos['SubhaloPos'][oid][1] = subhalos['SubhaloPos'][oid][1] + 75000
        
        if (dz <= neg_half_length):
            subhalos['SubhaloPos'][oid][2] = subhalos['SubhaloPos'][oid][2] - 75000
            #print("new z=" + str(subhalos['SubhaloPos'][subhalo][2]))
        elif (dz >= half_length):
            subhalos['SubhaloPos'][oid][2] = subhalos['SubhaloPos'][oid][2] + 75000
            #print("new z=" + str(subhalos['SubhaloPos'][subhalo][2]))
        x_optic = subhalos['SubhaloPos'][oid][0]-grpos_x
        y_optic = subhalos['SubhaloPos'][oid][1]-grpos_y
        z_optic = subhalos['SubhaloPos'][oid][2]-grpos_z
        r = math.sqrt(((x_optic**2)+(y_optic**2)))
        r_optic_xy.append(r)
    print("number of optical galaxies =" + str(len(optic_det)))    
    optic_table = Table([optic_det, r_optic_xy])
    

    optic_table.write("optic_table_unsorted_" + str(binnumber) + ".tex", format = 'latex', overwrite = True)
    sorted_optic_table = sorted(optic_table, key=itemgetter(1))
    return sorted_optic_table 
'''
sorted_op_table_45 = get_sorted_optic_table(IDs_45, 45)
sorted_op_table_69 = get_sorted_optic_table(IDs_69, 69)
sorted_op_table_1019 = get_sorted_optic_table(IDs_1019, 1019)
sorted_op_table_2030 = get_sorted_optic_table(IDs_2030, 2030)
'''
def get_sorted_HI_table(IDs, binnumber):
    HI_det = []
    r_HI_xy = []
    for h in IDs:
        HI_det_for_h = np.where(galaxies_HI_det_table_1e95['GrNr'] == h)
        for gal in HI_det_for_h[0]:
            HI_det.append(galaxies_HI_det_table_1e95['Galaxy IDs'][gal])
    for hid in HI_det:
        grnr1 = subhalos['SubhaloGrNr'][hid]
        grpos1_x = halos['GroupPos'][grnr1][0]
        grpos1_y = halos['GroupPos'][grnr1][1]
        grpos1_z = halos['GroupPos'][grnr1][2]
        dx = halos['GroupPos'][grnr1][0] - subhalos['SubhaloPos'][hid][0]
        dy = halos['GroupPos'][grnr1][1] - subhalos['SubhaloPos'][hid][1]
        dz = halos['GroupPos'][grnr1][2] - subhalos['SubhaloPos'][hid][2]
        #print("dz=" + str(dz))
        
        if (dx <= neg_half_length):
            subhalos['SubhaloPos'][hid][0] = subhalos['SubhaloPos'][hid][0] - 75000
        elif (dx >= half_length):
            subhalos['SubhaloPos'][hid][0] = subhalos['SubhaloPos'][hid][0] + 75000
        
        if (dy <= neg_half_length):
            subhalos['SubhaloPos'][hid][1] = subhalos['SubhaloPos'][hid][1] - 75000
        elif (dy >= half_length):
            subhalos['SubhaloPos'][hid][1] = subhalos['SubhaloPos'][hid][1] + 75000
        
        if (dz <= neg_half_length):
            subhalos['SubhaloPos'][hid][2] = subhalos['SubhaloPos'][hid][2] - 75000
            #print("new z=" + str(subhalos['SubhaloPos'][subhalo][2]))
        elif (dz >= half_length):
            subhalos['SubhaloPos'][hid][2] = subhalos['SubhaloPos'][hid][2] + 75000
            #print("new z=" + str(subhalos['SubhaloPos'][subhalo][2]))
        x_HI = subhalos['SubhaloPos'][hid][0]-grpos1_x
        y_HI = subhalos['SubhaloPos'][hid][1]-grpos1_y
        z_HI = subhalos['SubhaloPos'][hid][2]-grpos1_z
        r = math.sqrt(((x_HI**2)+(y_HI**2)))
        r_HI_xy.append(r)
    print("number of HI galaxies =" + str(len(HI_det)))
    HI_table = Table([HI_det, r_HI_xy])
    
    HI_table.write("HI_table_unsorted_" + str(binnumber) + ".tex", format = 'latex', overwrite = True)
    sorted_HI_table = sorted(HI_table, key=itemgetter(1))
    print(sorted_HI_table[:300])
    return sorted_HI_table

sorted_HI_table_45 = get_sorted_HI_table(IDs_45, 45)
#sorted_HI_table_69 = get_sorted_HI_table(IDs_69, 69)
#sorted_HI_table_1019 = get_sorted_HI_table(IDs_1019, 1019)
#sorted_HI_table_2030 = get_sorted_HI_table(IDs_2030, 2030)
'''
HI_radii_concentrations = []

op_l_45 = sorted_op_table_45
op_l_69 = sorted_op_table_69
op_l_1019 = sorted_op_table_1019
op_l_2030 = sorted_op_table_2030


HI_l_45 = sorted_HI_table_45
HI_l_69 = sorted_HI_table_69
HI_l_1019 = sorted_HI_table_1019
HI_l_2030 = sorted_HI_table_2030

#optical
def get_C_31(sorted_table):
    op_l = len(sorted_table)
    i_25 = int(round(0.25*op_l, 0))#placeholder for the radius in which 25% of galaxies fall in
    i_75 = int(round(0.75*op_l, 0))#placeholder for the radius in which 75% of galaxies fall in
    r_25 = sorted_table[i_25]['col1']
    r_75 = sorted_table[i_75]['col1']
    C_31 = r_75/r_25
    print(C_31)
    return C_31

op_radii_concentrations = [get_C_31(op_l_45), get_C_31(op_l_69), get_C_31(op_l_1019), get_C_31(op_l_2030)]

HI_radii_concentrations = [get_C_31(HI_l_45), get_C_31(HI_l_69), get_C_31(HI_l_1019), get_C_31(HI_l_2030)]

bin_numbers = [1, 2, 3, 4]
labels = ["4<=N<=5", "6<=N<=9", "10<=N<=19", "20<=N<=30"]

fig = plt.figure()
plt.plot(bin_numbers, op_radii_concentrations, '-', label = 'Optical detected galaxies')
plt.plot(bin_numbers, HI_radii_concentrations, '--', label = 'HI detected galaxies')
plt.xticks(bin_numbers, labels)
plt.title("Concentrations of galaxies")
plt.xlabel("Group Membership")
plt.ylabel("Group Concentration C31c")
plt.legend()
plt.show()


HI_fractions = [1071.0/1663.0, 719.0/1451.0, 670.0/1415.0, 265.0/625.0]
fig1 = plt.figure()
plt.plot(bin_numbers, HI_fractions, '-', label = 'Fraction of HI detections per bin')
plt.xticks(bin_numbers, labels)
plt.title('Fraction of HI detections per bin')
plt.xlabel("Group Membership")
plt.ylabel("Fraction of HI detections")
plt.show()
'''
print(datetime.now()-startTime)
    
    
    
#print(optic_table)
#optic_table.write("optic_table_unsorted_45.tex", format = 'latex')
#sorted_optic_table = sorted(optic_table, key=itemgetter(1))
#print(sorted_optic_table)
    
    
    
    
    
    
    
    
    