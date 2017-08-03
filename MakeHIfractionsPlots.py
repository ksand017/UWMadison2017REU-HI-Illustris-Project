#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 31 09:56:21 2017

@author: sanderson

Downloads Needed:
    Illustris_Python
    specified snapshot number Groupcat download(FOF & Subfind)
    
Files Needed:
    'Galaxies HI det table HI mass fraction.tex'
    (created from Make HI fraction Det Table.py)
    
Result:
    Plots HI mass fraction as a function of total mass for each galaxy
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

galaxies_HI_det_table= Table.read('Galaxies HI det table HI mass fraction.tex', format = 'latex')
galaxies_HI_det_table_1e6= Table.read('Galaxies HI det table HI mass fraction_1e6 lower limit.tex', format = 'latex')
galaxies_HI_det_table_1e9= Table.read('Galaxies HI det table HI mass fraction_1e9 lower limit.tex', format = 'latex')
galaxies_HI_det_table_1e10= Table.read('Galaxies HI det table HI mass fraction_1e10 lower limit.tex', format = 'latex')
galaxies_HI_det_table_1e95= Table.read('Galaxies HI det table HI mass fraction_1e9.5 lower limit.tex', format = 'latex')



def plot_HI_fractions(galaxy_table, low_limit):
    IDs = galaxy_table['Galaxy IDs']
    HI_mass = galaxy_table['HI Mass 1e10']
    print("For a HI mass lower limit of " + str(low_limit) + ":")
    print("maximum HI mass =" + str(max(HI_mass)) + " (1e10)solar masses")
    print("minimum HI mass =" + str(min(HI_mass)) + " (1e10)solar masses")
    
    HI_fractions = galaxy_table['HI Mass Fraction']
    total_mass = []
    for i in IDs:
        mass = subhalos['SubhaloMass'][i]
        total_mass.append(mass)
    print("number of HI detections =" + str(len(IDs)))
    HI_log_fractions = np.log10(HI_fractions)
    total_mass_log = np.log10(total_mass)
    plt.figure()
    plt.plot(total_mass_log, HI_fractions, '.')
    plt.title("HI Mass Fraction vs. total mass for Illustris-1 Simulation: HI lower limit =" + str(low_limit))
    plt.xlabel(" log($M_t$/$M_\odot$)")
    plt.ylabel(" M_HI$/$M_t$$M_\odot$")
    
    
    plt.show()
    
    fig1 = plt.figure()
    plt.hist(HI_log_fractions, bins = 19)
    plt.title("HI Mass Fractions")
    plt.xlabel(" log($M_HI$/$M_\odot$)")
    plt.ylabel("frequency")
    
    plt.show()
    
    print("maximum HI mass fraction =" + str(max(HI_fractions)))
    print("minimum HI mass fraction =" + str(min(HI_fractions)))
    
    
#plot_HI_fractions(galaxies_HI_det_table, 1e9)

#plot_HI_fractions(galaxies_HI_det_table_1e6, 1e6)

plot_HI_fractions(galaxies_HI_det_table_1e95, 10**9.5)
print(datetime.now()-startTime)