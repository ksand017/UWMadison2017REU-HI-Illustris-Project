# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
Dowloads Needed:
    Illustris_Python
    specified snapshot number Groupcat download(FOF & Subfind)
    
Files Needed:
    'Galaxies HI det table HI mass fraction_1e9 lower limit.tex' where the lower limit 1e9 can be changed depending on 
    preference 
    (created from the 'HI and Optic det table.py' script)
    
Result:
    A table of galaxies corresponding to a specified HI mass limit. Table includes IDs, HI masses, HI mass fractions, 
    and group numbers

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

galaxy_table = Table.read('Galaxies HI det table with lower limit 1e6.tex', format = 'latex')

galaxy_ids = galaxy_table['Galaxy IDs']
HI_mass = galaxy_table['HI Mass 1e10[$M_\odot$]']

print("maximum HI mass =" + str(max(HI_mass)) + " (1e10)solar masses")
print("minimum HI mass =" + str(min(HI_mass)) + " (1e10)solar masses")


subhalos_ids = []
Group_number = []
HI_mass_fraction = []
HI_Mass = []
ranger = len(galaxy_ids)
for i in range(ranger):
    Id = galaxy_ids[i]
    print(Id)
    HI_m = HI_mass[i]*1e10
    if((HI_m >= (10**(9.5))) & (HI_m <= (1e11))):
        HI_Mass.append(HI_mass[i])
        subhalos_ids.append(Id)
        total_mass = subhalos['SubhaloMass'][Id]
        grnr = subhalos['SubhaloGrNr'][Id]
        Group_number.append(grnr)
        HI_mass_fraction.append(HI_mass[i]/total_mass)

galaxies_HI_det_table = Table([subhalos_ids, HI_Mass, HI_mass_fraction, Group_number], names = ('Galaxy IDs', 'HI Mass 1e10', 'HI Mass Fraction', 'GrNr'))
galaxies_HI_det_table.write('Galaxies HI det table HI mass fraction_1e9.5 lower limit.tex', format = 'latex')
   

    



def make_galaxy_HI_det_table(galaxy_ids):
    subhalos_ids = []
    HI_mass = []
    Group_number = []
    HI_mass_fraction = []
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
            if((HI_sum_detection_test >= (1e10)) & (HI_sum_detection_test <= (1e11))): #test to see if HI mass is detectable
                HI_mass.append(HI_sum)
                subhalos_ids.append(g)
                subhalo_mass = subhalos['SubhaloMass'][g]
                HI_mass_fraction.append(HI_sum/subhalo_mass)
                Group_number.append(subhalos['SubhaloGrNr'][g])
        print(datetime.now()-startTime)
    
    galaxies_HI_det_table = Table([subhalos_ids, HI_mass, HI_mass_fraction, Group_number], names = ('Galaxy IDs', 'HI Mass 1e10', 'HI Mass Fraction', 'GrNr'))
    galaxies_HI_det_table.write('Galaxies HI det table HI mass fraction with lower limit 1e10.tex', format = 'latex')

#make_galaxy_HI_det_table(galaxy_ids)
    
    
print(datetime.now()-startTime)
