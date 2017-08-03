#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 12:29:15 2017

@author: sanderson


Downloads Needed:
    Illustris_Python
    specified snapshot number Groupcat download(FOF & Subfind)
    
Result:
    makes table of all subhalos in a simulation suite along with their positions, group number, and individual mass

"""


import illustris_python as il
from astropy.table import Table, Column
import numpy as np
from datetime import datetime
 
startTime = datetime.now()
basePath = './Illustris-1_snap132' #directory to illustris_python and illustris files
fields1 = ['SubhaloPos', 'SubhaloGrNr', 'SubhaloSFR', 'SubhaloMass', 'SubhaloStellarPhotometrics'] #fields to be requested
subhalos = il.groupcat.loadSubhalos(basePath,132,fields=fields1) #Call illustris_python scrips to request specified fields


print(subhalos.keys())
num_subs = 4389572 #check illustris website to find the number of subhalos needed for this simulation suite
subhalos_ID = range(0, num_subs) #make list of halo ids for specified illustris run
print(len(subhalos_ID))

print(datetime.now()-startTime)

all_subs = Table([subhalos_ID, subhalos['SubhaloPos'], subhalos['SubhaloGrNr'], subhalos['SubhaloMass']], names=('Sub IDs', 'Position(x,y,z)', 'Group Number', 'Mass'))
#make table of data for all groups in specified run

all_subs.write('galaxies_idsIll1_snap132.tex', format = 'latex')
#write table to file

print(datetime.now()-startTime)
#retrieve parent halo info and store in table named 'groupids.txt'print(all_groups['Number of Subhalos'][0])
