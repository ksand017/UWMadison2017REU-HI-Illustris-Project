#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 17 14:51:56 2017

@author: sanderson

seperates all galaxies produced from the MakeGalaxiesTable.py script in to tables based on 
their number of optically detected galaxies as well as their HI mass

Downloads Needed:
    Illustris_Python
    specified snapshot number Groupcat download(FOF & Subfind)
    specified snapshot number Snapshot dowanload
    
Files Needed:
    galaxies_idsIll1_snap132.tex created from the MakeGalaxiesTable.py script
    
Result:
    4 tables; optically detected galaxies binned by optical group membership
    1 table; all HI detections above the 1e6 and below 1e11 solar masses (parameters can be changed according to preference)
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


#group_table = Table.read('groupidsIll1_snap132.tex', format = 'latex')
galaxy_table = Table.read('galaxies_idsIll1_snap132.tex', format = 'latex')
'''
print(len(group_table))
groupm_2030 = np.where((group_table['Number of Subhalos'] >= 20) & (group_table['Number of Subhalos'] <= 30))
groupm_1019 = np.where((group_table['Number of Subhalos'] >= 10) & (group_table['Number of Subhalos'] <= 19 ))
groupm_69 = np.where((group_table['Number of Subhalos'] >= 6) & (group_table['Number of Subhalos'] <= 9 ))
groupm_45 = np.where((group_table['Number of Subhalos'] >= 4) & (group_table['Number of Subhalos'] <= 5 ))
print(len(groupm_45[0]))

print('number of subs in 2030 bin= ' + str(len(groupm_2030[0])))
print('number of subs in 1019 bin= ' + str(len(groupm_1019[0])))
print('number of subs in 69 bin= ' + str(len(groupm_69[0])))
print('number of subs in 45 bin= ' + str(len(groupm_45[0])))
'''
#group_ids = group_table['Group IDs']
galaxy_ids = galaxy_table['Sub IDs']

def make_galaxy_HI_det_table(galaxy_ids):
    subhalos_ids = []
    HI_mass = []
    Group_number = []
    for g in galaxy_ids:
        print(g)
        subhalo_gas_data = il.snapshot.loadSubhalo(basePath,132,str(g),'gas')
        ranger_gas = subhalo_gas_data['count']
        HI_mass_partial_sub = [] #list for HI masses of a single subhalo
        if (('Masses' in subhalo_gas_data) & ('NeutralHydrogenAbundance' in subhalo_gas_data)):
            
            for i in range(ranger_gas):
                gas_m = subhalo_gas_data['Masses'][i] #gas mass for a single cell within a single subhalo in 1e10 solar masses
                HI_abun = subhalo_gas_data['NeutralHydrogenAbundance'][i]#HI abundance for a single cell within a single subhalo
                HI_m = gas_m*HI_abun*(0.76) #HI mass for a single cell within a single subhalo 
                HI_mass_partial_sub.append(HI_m) #list of HI masses for all cells in a single subhalo
            HI_sum = sum(HI_mass_partial_sub) #total HI mass for a single sub
            HI_sum_detection_test = HI_sum*(1e10) #mass of HI in solar masses for a single sub
            if((HI_sum_detection_test >= (1e9)) & (HI_sum_detection_test <= (1e11))): #test to see if HI mass is detectable
                HI_mass.append(HI_sum)
                subhalos_ids.append(g)
                Group_number.append(subhalos['SubhaloGrNr'][g])
        print(datetime.now()-startTime)
    
    galaxies_HI_det_table = Table([subhalos_ids, HI_mass, Group_number], names = ('Galaxy IDs', 'HI Mass 1e10[$M_\odot$]', 'GrNr'))
    galaxies_HI_det_table.write('Galaxies HI det table HI mass fraction_1e9 lower limit.tex', format = 'latex')

make_galaxy_HI_det_table(galaxy_ids)

def grouphalo_optic_det_table(group_ids):
    subhalo_r_band_mag = []
    halo_ids_45 = []
    halo_ids_69 = []
    halo_ids_1019 = []
    halo_ids_2030 = []
    num_optical_members_45 = []
    num_optical_members_69 = []
    num_optical_members_1019 = []
    num_optical_members_2030 = []
    num_total_members_45 = []
    num_total_members_69 = []
    num_total_members_1019 = []
    num_total_members_2030 = []
    groups_luminosities_45 = []
    groups_luminosities_69 = []
    groups_luminosities_1019 = []
    groups_luminosities_2030 = []
    for h in group_ids:
        print("group id =" + str(h))
        mag = []
        group_lum = []
        subhalo_idlist = np.arange(halos['GroupFirstSub'][h], halos['GroupFirstSub'][h] + halos['GroupNsubs'][h])
        count = 0
        mark = len(subhalo_idlist)
        print("subhalo id list =" + str(subhalo_idlist))
        print("length of subhalo id list is " + str(mark))
        for s in subhalo_idlist:
            mag = subhalos['SubhaloStellarPhotometrics'][s][5]
            if mag > -17:
                continue
            else:
                count = count + 1
                lum = (3.86e26)*(10**((4.77-mag)/2.5))
                group_lum.append(lum)
                print("count is at " + str(count))
        group_luminosity = sum(group_lum)
        if ((count<=5) & (count>=4)):
            halo_ids_45.append(h)
            num_optical_members_45.append(count)
            num_total_members_45.append(len(subhalo_idlist))
            groups_luminosities_45.append(group_luminosity)
        elif ((count<=9) & (count>=6)):
            halo_ids_69.append(h)
            num_optical_members_69.append(count)
            num_total_members_69.append(len(subhalo_idlist))
            groups_luminosities_69.append(group_luminosity)
        elif ((count<=19) & (count>=10)):
            halo_ids_1019.append(h)
            num_optical_members_1019.append(count)
            num_total_members_1019.append(len(subhalo_idlist))
            groups_luminosities_1019.append(group_luminosity)
        elif ((count<=30) & (count>=20)):
            halo_ids_2030.append(h)
            num_optical_members_2030.append(count)
            num_total_members_2030.append(len(subhalo_idlist))
            groups_luminosities_2030.append(group_luminosity)
    '''
    for g in halo_ids:
         subhalo_idlist1 = np.arange(halos['GroupFirstSub'][g], halos['GroupFirstSub'][g] + halos['GroupNsubs'][g])
         for s1 in subhalo_idlist1:
             subhalos_ids.append(s1)
    '''
    num_45 = len(halo_ids_45)
    number_list_45 = range(0,num_45)
    halo_id_table_45 = Table([number_list_45, halo_ids_45, num_optical_members_45, num_total_members_45, groups_luminosities_45], names = ('no.','IDs', 'Num of optical galaxies', 'Num of total galaxies', 'Luminosity of Group'))           
    halo_id_table_45.write('ids for groups with 4-5 optically detectable galaxies_mag17.tex', format = 'latex')
    
    num_69 = len(halo_ids_69)
    number_list_69 = range(0,num_69)
    halo_id_table_69 = Table([number_list_69, halo_ids_69, num_optical_members_69, num_total_members_69, groups_luminosities_69], names = ('no.','IDs', 'Num of optical galaxies', 'Num of total galaxies', 'Luminosity of Group'))           
    halo_id_table_69.write('ids for groups with 6-9 optically detectable galaxies_mag17.tex', format = 'latex')
    
    num_1019 = len(halo_ids_1019)
    number_list_1019 = range(0,num_1019)
    halo_id_table_1019 = Table([number_list_1019, halo_ids_1019, num_optical_members_1019, num_total_members_1019, groups_luminosities_1019], names = ('no.','IDs', 'Num of optical galaxies', 'Num of total galaxies', 'Luminosity of Group'))           
    halo_id_table_1019.write('ids for groups with 10-19 optically detectable galaxies_mag17.tex', format = 'latex')
    
    num_2030 = len(halo_ids_2030)
    number_list_2030 = range(0,num_2030)
    halo_id_table_2030 = Table([number_list_2030, halo_ids_2030, num_optical_members_2030, num_total_members_2030, groups_luminosities_2030], names = ('no.','IDs', 'Num of optical galaxies', 'Num of total galaxies', 'Luminosity of Group'))           
    halo_id_table_2030.write('ids for groups with 20-30 optically detectable galaxies_mag17.tex', format = 'latex')
    
    #all_optic_det_table = Table([subhalos_ids, subhalo_r_band_mag], names = ('IDs', 'r-band magnitude'))
    #ptic_det_table_ids = Table([subhalos_ids], names = ('IDs'))
    #all_optic_det_table.write('all subhalo optical detections for' + str(bin_number) + '.tex', format='latex')                
    #halo_id_table.write('ids for groups with optically detectable galaxies' + str(bin_number) + '.tex', format = 'latex')
#grouphalo_optic_det_table(galaxy_ids)

print(datetime.now()-startTime)
