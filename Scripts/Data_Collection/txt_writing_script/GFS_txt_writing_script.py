#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:38:37 2019

@author: Rarrell

This is the script used to write the GFS_Variables.txt, GFS_TypeOfHeights.txt,
and GFS_Hieghts.txt
"""

#%% Import modules

import sys, os, warnings
import numpy as np
import pygrib

#%% Load a random GFS file and collect variable names, type of heights, and
#   heights.
gfs = pygrib.open('/Users/Rarrell/Desktop/Thesis_Research/Scripts/' +\
                  'Data_Collection/txt_writing_script/gfs_3_20190610_0000_006.grb2')

# Determine how many variables are in the grib file.
n = 0

for g in gfs:
    n = n + 1

# Initialize some variables
VariableName  = ['tmp'] * n
TypeOfHeights = ['tmp'] * n
Heights       = np.ones((n, )) * np.nan

# Fill the variables with values. 

# This must be called again otherwise the next loop is ignored.
gfs = pygrib.open('/Users/Rarrell/Desktop/Thesis_Research/Scripts/' +\
                  'Data_Collection/txt_writing_script/gfs_3_20190610_0000_006.grb2')

n = 0
for g in gfs:
    VariableName[n]  = g.name
    TypeOfHeights[n] = g.typeOfLevel
    Heights[n]         = int(g.level)
    n = n + 1

print(VariableName)
print(TypeOfHeights)
print(Heights)

gfs.close()

#%% Sort, and determine the unique values of the variables.

VariableName  = np.sort(VariableName)
TypeOfHeights = np.sort(TypeOfHeights)
Heights       = np.sort(Heights)

VarNameUniq      = np.unique(VariableName)
TypeOfHeightUniq = np.unique(TypeOfHeights)
HeightsUniq      = np.unique(Heights)

print(VarNameUniq)
print(TypeOfHeightUniq)
print(HeightsUniq)

#%% Write the variables names to a .txt file.
f = open('GFS_Variable_Names.txt', 'w')

for name in VarNameUniq:
    f.write(name + ',\n')

f.close()

#%% Write the type of heights to a .txt file.
f = open('GFS_TypeOfHeights.txt', 'w')

for typeOfHeight in TypeOfHeightUniq:
    f.write(typeOfHeight + ',\n')

f.close()

#%% Write the heights to a .txt file.
f = open('GFS_Height_Values.txt', 'w')

for height in HeightsUniq:
    f.write(str(height) + ',\n')

f.close()