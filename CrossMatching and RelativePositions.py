#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 13:31:29 2017

@author: sanderson


Downloads Needed:
    Illustris_Python
    specified snapshot number Groupcat download(FOF & Subfind)
    specified snapshot number Snapshot dowanload
    
Files Needed:
    'Galaxies HI det table.tex' 
    'ids for groups with 4-5 optically detectable galaxies.tex'
    'ids for groups with 6-9 optically detectable galaxies.tex'
    'ids for groups with 10-19 optically detectable galaxies.tex'
    'ids for groups with 20-30 optically detectable galaxies.tex'
    (all created from the 'HI and Optic det table.py' script)
    
Result:
    Makes plots for the relative position of subhalos belonging to a member bin. 
    Plots can be random samples from bin or the entire bin.
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
fields2 = ['GroupFirstSub', 'GroupNsubs', 'GroupPos', 'Group_M_Crit200', 'Group_R_Crit200'] #fields to be requested
halos = il.groupcat.loadHalos(basePath,132,fields=fields2) #Call illustris_python scrips to request specified fields



halo_id_table_45 = Table.read('ids for groups with 4-5 optically detectable galaxies_mag17.tex', format = 'latex')
halo_id_table_69 = Table.read('ids for groups with 6-9 optically detectable galaxies_mag17.tex', format = 'latex')
halo_id_table_1019 = Table.read('ids for groups with 10-19 optically detectable galaxies_mag17.tex', format = 'latex')
halo_id_table_2030 = Table.read('ids for groups with 20-30 optically detectable galaxies_mag17.tex', format = 'latex')

HI_det_table = Table.read('Galaxies HI det table.tex', format = 'latex')
galaxies_HI_det_table_1e10= Table.read('Galaxies HI det table HI mass fraction_1e10 lower limit.tex', format = 'latex')
galaxies_HI_det_table_1e95= Table.read('Galaxies HI det table HI mass fraction_1e9.5 lower limit.tex', format = 'latex')
galaxies_HI_det_table_1e9= Table.read('Galaxies HI det table HI mass fraction_1e9 lower limit.tex', format = 'latex')

print("total number of HI detected galaxies in simulation =" + str(len(HI_det_table)))

IDs_45 = halo_id_table_45['IDs']
IDs_69 = halo_id_table_69['IDs']
IDs_1019 = halo_id_table_1019['IDs']
IDs_2030 = halo_id_table_2030['IDs']

#print(HI_det_table)
#print(halo_id_table_45)



length = 75000
half_length = length/2
neg_half_length = -1*length/2
random_IDs_45 = random.sample(IDs_45, 5)
random_IDs_69 = random.sample(IDs_69, 5)
random_IDs_1019 = random.sample(IDs_1019, 5)
random_IDs_2030 = random.sample(IDs_2030, 5)
optic_det = []
HI_det = []
x_optic = []
y_optic = []
z_optic = []
x_HI = []
y_HI = []
z_HI = []

def make_rando_plots(random_IDs, bin_number):
    
    for h in random_IDs:
        print('group id = ' + str(h))
        optic_det = []
        HI_det = []
        x_optic = []
        y_optic = []
        z_optic = []
        x_HI = []
        y_HI = []
        z_HI = []
        group_cm_x = halos['GroupPos'][h][0]
        group_cm_y = halos['GroupPos'][h][1]
        group_cm_z = halos['GroupPos'][h][2]
        subhalo_idlist = np.arange(halos['GroupFirstSub'][h], halos['GroupFirstSub'][h] + halos['GroupNsubs'][h]) 
        #list of subhalos in that group
        for s in subhalo_idlist:
            mag = subhalos['SubhaloStellarPhotometrics'][s][5] #mag for galaxy in h
            if (mag < -18): #optical detection test
                optic_det.append(s)
        HI_det_for_h = np.where(galaxies_HI_det_table_1e10['GrNr'] == h)
        for gal in HI_det_for_h[0]:
            HI_det.append(galaxies_HI_det_table_1e10['Galaxy IDs'][gal])
        for oid in optic_det:
            grnr = subhalos['SubhaloGrNr'][oid]
            grpos_x = halos['GroupPos'][grnr][0]
            grpos_y = halos['GroupPos'][grnr][1]
            grpos_z = halos['GroupPos'][grnr][2] 
            dx = halos['GroupPos'][grnr][0] - subhalos['SubhaloPos'][oid][0]
            dy = halos['GroupPos'][grnr][1] - subhalos['SubhaloPos'][oid][1]
            dz = halos['GroupPos'][grnr][2] - subhalos['SubhaloPos'][oid][2]
            #print("dz=" + str(dz))
                
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
            x_optic.append(subhalos['SubhaloPos'][oid][0]-grpos_x)
            y_optic.append(subhalos['SubhaloPos'][oid][1]-grpos_y)
            z_optic.append(subhalos['SubhaloPos'][oid][2]-grpos_z)
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
            x_HI.append(subhalos['SubhaloPos'][hid][0]-grpos1_x)
            y_HI.append(subhalos['SubhaloPos'][hid][1]-grpos1_y)
            z_HI.append(subhalos['SubhaloPos'][hid][2]-grpos1_z)
        
        total_subs = len(optic_det) + len(HI_det)
        grRcrit200 = halos['Group_R_Crit200'][h]
        for i in optic_det:
            if i in HI_det:
                total_subs-= 1
        num_HI_det = len(HI_det)
        fig = plt.figure()
        plt.plot(x_optic, y_optic, '.', label = 'Optical detected galaxies')
        plt.plot(x_HI, y_HI, 'o', fillstyle = 'none', label = 'HI detected galaxies')
        ax=fig.add_subplot(1,1,1)
        circ=plt.Circle((0,0), radius=grRcrit200, color='g', fill=False)
        ax.add_patch(circ)
        plt.title("relative positions for bin " + str(bin_number))
        plt.axis([-1000,1000,-1000,1000])
        plt.legend()
        plt.show()
        print("for bin " + str(bin_number) )
        print("the group number is: " + str(h))
        print("the total number of subhalos is: " + str(total_subs))
        print("the number of HI detected galaxies is: " + str(num_HI_det))
    
        
    
#make_rando_plots(random_IDs_45, 45)
#make_rando_plots(random_IDs_69, 69)
#make_rando_plots(random_IDs_1019, 1019) 
#make_rando_plots(random_IDs_2030, 2030)  
    
def plot_relative_positions(IDs, bin_number):  
    optic_det = []
    HI_det = []
    x_optic = []
    y_optic = []
    z_optic = []
    x_HI = []
    y_HI = []
    z_HI = []
    OnlyOpticCount = 0
    OnlyHIcount = 0
    bothCount = 0
    num_groups = len(IDs)
    for h in IDs:
        subhalo_idlist = np.arange(halos['GroupFirstSub'][h], halos['GroupFirstSub'][h] + halos['GroupNsubs'][h]) 
        #list of subhalos in that group
        for s in subhalo_idlist:
            mag = subhalos['SubhaloStellarPhotometrics'][s][5] #mag for galaxy in h
            if (mag < -18): #optical detection test
                optic_det.append(s)
        HI_det_for_h = np.where(galaxies_HI_det_table_1e95['GrNr'] == h)
        for gal in HI_det_for_h[0]:
            HI_det.append(galaxies_HI_det_table_1e95['Galaxy IDs'][gal])
            
    #print(HI_det)
    for oid in optic_det:
        grnr = subhalos['SubhaloGrNr'][oid]
        grpos_x = halos['GroupPos'][grnr][0]
        grpos_y = halos['GroupPos'][grnr][1]
        grpos_z = halos['GroupPos'][grnr][2]
        dx = halos['GroupPos'][grnr][0] - subhalos['SubhaloPos'][oid][0]
        dy = halos['GroupPos'][grnr][1] - subhalos['SubhaloPos'][oid][1]
        dz = halos['GroupPos'][grnr][2] - subhalos['SubhaloPos'][oid][2]
        #print("dz=" + str(dz))
        
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
        x_optic.append(subhalos['SubhaloPos'][oid][0]-grpos_x)
        y_optic.append(subhalos['SubhaloPos'][oid][1]-grpos_y)
        z_optic.append(subhalos['SubhaloPos'][oid][2]-grpos_z)
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
        x_HI.append(subhalos['SubhaloPos'][hid][0]-grpos1_x)
        y_HI.append(subhalos['SubhaloPos'][hid][1]-grpos1_y)
        z_HI.append(subhalos['SubhaloPos'][hid][2]-grpos1_z)
        
    total_subs = len(optic_det) + len(HI_det)
    print("num optical det =" +str(len(optic_det)))
    print("num HI det =" +str(len(HI_det)))
    for i in optic_det:
        if i in HI_det:
            bothCount+=1
            total_subs-= 1
        elif i not in HI_det:
            OnlyOpticCount+=1
    for i in HI_det:
        if i not in optic_det:
            OnlyHIcount+=1
    num_HI_det = len(HI_det)
            
    fig = plt.figure()
    plt.plot(x_optic, y_optic, '.', label = 'Optical detected galaxies')
    plt.plot(x_HI, y_HI, 'o', fillstyle = 'none', label = 'HI detected galaxies')
    ax=fig.add_subplot(1,1,1)
    circ=plt.Circle((0,0), radius=1000, color='g', fill=False)
    ax.add_patch(circ)
    plt.axis([-1100,1100,-1100,1100])
    plt.title("relative positions for bin " + str(bin_number))
    plt.legend()
    plt.show()
    print("for bin " + str(bin_number) )
    print("the total number of groups is: " + str(num_groups))
    print("the total number of subhalos is: " + str(total_subs))
    print("the number of HI only detected galaxies is: " + str(OnlyHIcount))
    print("the number of optical only detected galaxies is: " + str(OnlyOpticCount))
    print("the number of galaxies detected both optically and in HI is: " + str(bothCount))


    
plot_relative_positions(IDs_45, 45)
   
plot_relative_positions(IDs_69, 69)

plot_relative_positions(IDs_1019, 1019)

plot_relative_positions(IDs_2030, 2030)
      
print(datetime.now()-startTime)
