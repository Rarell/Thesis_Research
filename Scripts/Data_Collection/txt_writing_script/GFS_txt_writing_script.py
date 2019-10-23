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

# This file is for operational use in codes. The new line add an empty column
# that makes the easier to read versions more difficult to use practically in
# codework.
f = open('GFS_VarNames.txt', 'w')

for name in VarNameUniq:
    if name == VarNameUniq[-1]:
        f.write(name)
    else:
        f.write(name + ',')

f.close()

#%% Write the type of heights to a .txt file.
f = open('GFS_Type_Of_Heights.txt', 'w')

for typeOfHeight in TypeOfHeightUniq:
    f.write(typeOfHeight + ',\n')

f.close()

# This file is for operational use in codes. The new line add an empty column
# that makes the easier to read versions more difficult to use practically in
# codework.
f = open('GFS_TypeOfHeights.txt', 'w')

for typeOfHeight in TypeOfHeightUniq:
    if typeOfHeight == TypeOfHeightUniq[-1]:
        f.write(typeOfHeight)
    else:
        f.write(typeOfHeight + ',')

f.close()

#%% Write the heights to a .txt file.
f = open('GFS_Height_Values.txt', 'w')

for height in HeightsUniq:
    f.write(str(height) + ',\n')

f.close()

# This file is for operational use in codes. The new line add an empty column
# that makes the easier to read versions more difficult to use practically in
# codework.
f = open('GFS_Heights.txt', 'w')

for height in HeightsUniq:
    if height == HeightsUniq[-1]:
        f.write(str(height))
    else:
        f.write(str(height) + ',')

f.close()

#%% Write a .txt file for the forecast time.
# Note, starting June 14, 2019 the GFS forecasts became a consistent
# 3 hour intervals from 000 to 384. Might modify this for those later.

# Initialize the forecast hours. It has 92 values (excluding 000).
FT = ['tmp'] * 92
n = 0

for i in range(3, 240+3, 3):
    if i < 10:
        FT[n] = '00' + str(i)
    elif (i < 100) & (i > 10):
        FT[n] = '0' + str(i)
    else:
        FT[n] = str(i)
    n = n + 1

for i in range(240+12, 384+12, 12):
    FT[n] = str(i)
    n = n + 1

print(FT)

f = open('GFS_Forecast_Times.txt', 'w')

for ForecastTime in FT:
    f.write(ForecastTime + '\n')
    
f.close()

# This file is for operational use in codes. The new line add an empty column
# that makes the easier to read versions more difficult to use practically in
# codework.
f = open('GFS_FH.txt', 'w')

for ForecastTime in FT:
    if ForecastTime == FT[-1]:
        f.write(ForecastTime)
    else:
        f.write(ForecastTime + ',')
    
f.close()