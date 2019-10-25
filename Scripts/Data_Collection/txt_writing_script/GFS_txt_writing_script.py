#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 13:38:37 2019

@author: Rarrell

This is the script used to write the GFS_Variables.txt, GFS_TypeOfHeights.txt,
  and GFS_Hieghts.txt files for the extract_GFS_variable.sh program.
"""

#%%
#####################################
### Import some libraries ###########
#####################################

import sys, os, warnings
import numpy as np
import pygrib

#%% 
#############################################
### Loada .grb2 and Obtain Information ######
#############################################

# Load a random GFS file and collect variable names, type of heights, and
#   heights
gfs = pygrib.open('./Scripts/Data_Collection/txt_writing_script/' +\
                  'gfs_3_20190610_0000_006.grb2')

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
gfs = pygrib.open('./Scripts/Data_Collection/txt_writing_script/' +\
                  'gfs_3_20190610_0000_006.grb2')

# Obtain all the variable names, type of heights, and heights in a GFS .grb2 file
n = 0
for g in gfs:
    VariableName[n]  = g.name
    TypeOfHeights[n] = g.typeOfLevel
    Heights[n]         = int(g.level)
    n = n + 1

# Examine the collected information
print(VariableName)
print(TypeOfHeights)
print(Heights)

# Close the .grb2 file
gfs.close()

#%%
####################################
### Sort and Get Unique Values #####
####################################

# Sort the variables so they have a logical order
VariableName  = np.sort(VariableName)
TypeOfHeights = np.sort(TypeOfHeights)
Heights       = np.sort(Heights)

# Only unique values are needed. Checking against repeats is simply redundant.
VarNameUniq      = np.unique(VariableName)
TypeOfHeightUniq = np.unique(TypeOfHeights)
HeightsUniq      = np.unique(Heights)

# Examine sorted, unique values
print(VarNameUniq)
print(TypeOfHeightUniq)
print(HeightsUniq)

#%%
#########################################
### Write Variable Names to a .txt ######
#########################################

# Write the variable names in a column format (easier for user to read) for
#   for user reference.
f = open('GFS_Variable_Names.txt', 'w')

for name in VarNameUniq:
    f.write(name + ',\n')

f.close()

# This file is for operational use in codes. The new line add an empty column
#   that makes the easier to read (for users) versions more difficult to use 
#   practically in codework.
f = open('GFS_VarNames.txt', 'w')

for name in VarNameUniq:
    if name == VarNameUniq[-1]:
        f.write(name)
    else:
        f.write(name + ',')

f.close()

#%%
#########################################
### Write Type of Heights to a .txt #####
#########################################

# Write the type of heights in a column format (easier for user to read) for
#   for user reference.
f = open('GFS_Type_Of_Heights.txt', 'w')

for typeOfHeight in TypeOfHeightUniq:
    f.write(typeOfHeight + ',\n')

f.close()

# This file is for operational use in codes. The new line add an empty column
#   that makes the easier to read (for users) versions more difficult to use 
#   practically in codework..
f = open('GFS_TypeOfHeights.txt', 'w')

for typeOfHeight in TypeOfHeightUniq:
    if typeOfHeight == TypeOfHeightUniq[-1]:
        f.write(typeOfHeight)
    else:
        f.write(typeOfHeight + ',')

f.close()

#%%
######################################
### Write the Heights to a .txt ######
######################################

# Write the heights in a column format (easier for user to read) for
#   for user reference.
f = open('GFS_Height_Values.txt', 'w')

for height in HeightsUniq:
    f.write(str(height) + ',\n')

f.close()

# This file is for operational use in codes. The new line add an empty column
#   that makes the easier to read (for users) versions more difficult to use 
#   practically in codework.
f = open('GFS_Heights.txt', 'w')

for height in HeightsUniq:
    if height == HeightsUniq[-1]:
        f.write(str(height))
    else:
        f.write(str(height) + ',')

f.close()

#%%
#######################################
### Write Forecast Hours to a .txt ####
#######################################

# Note, starting June 14, 2019 the GFS forecasts became a consistent
# 3 hour intervals from 000 to 384. Might modify this for those later.

# Initialize the forecast hours. It has 92 values (excluding 000).
FT = ['tmp'] * 92
n = 0

# Create the forecast hour strings for values below 240 hours
for i in range(3, 240+3, 3):
    if i < 10:
        FT[n] = '00' + str(i)
    elif (i < 100) & (i > 10):
        FT[n] = '0' + str(i)
    else:
        FT[n] = str(i)
    n = n + 1

# Create the forecast hour strings for values above 240
for i in range(240+12, 384+12, 12):
    FT[n] = str(i)
    n = n + 1

# Examine the forecast hours to make sure it is complete and correct
print(FT)

# Write the heights in a column format (easier for user to read) for
#   for user reference. This is not comma seperate for proper reading in the
#   bash script.
f = open('GFS_Forecast_Times.txt', 'w')

for ForecastTime in FT:
    f.write(ForecastTime + '\n')
    
f.close()

# This file is for operational use in codes. The new line add an empty column
#   that makes the easier to read (for users) versions more difficult to use 
#   practically in codework.
f = open('GFS_FH.txt', 'w')

for ForecastTime in FT:
    if ForecastTime == FT[-1]:
        f.write(ForecastTime)
    else:
        f.write(ForecastTime + ',')
    
f.close()

#####################
### End of Script ###
#####################