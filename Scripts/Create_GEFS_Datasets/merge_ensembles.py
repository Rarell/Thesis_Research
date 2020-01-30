#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 12:49:11 2020

@author: stuartedris
"""

#%%
#####################################
### Import some libraries ###########
#####################################

import os, sys, warnings
import numpy as np
from netCDF4 import Dataset
from glob import glob

#%% 
############################
### Main Function ##########
############################
def main():
    '''
    

    '''
    
    # Load arguements
    script = sys.argv[0]
    EMfile        = sys.argv[1]
    VarParameters = sys.argv[2]
    Parameters    = sys.argv[3]
    
    # Unpack the .txt file to obtain the model information
    VarName, TypeOfHeight, Height, VarSName = np.loadtxt(VarParameters,
                                                         usecols = (0,1,2,3),
                                                         dtype = str, delimiter = ',',
                                                         unpack = True)
    
    Year, Month, Day, ModelRun, Source = np.loadtxt(Parameters, 
                                                    usecols = (1,2,3,4,5), dtype = str, 
                                                    delimiter = ',', unpack = True)
    
    # Load all the ensemble members
    EnsembleData, FH, lat, lon, VarVD, mask, units = load_ensembles(VarSName, EM, 
                                                                    Year, Month, Day, ModelRun)
    
    # Merge the ensembles
    EnsembleMean, Spread = merge_ensembles()
    
    # Write the .nc files to contain the merged ensemble data
    write_nc()
    
    
#%%
# Function to load a single ensemble member
def load_ensemble(SName, file, EM):
    '''

    '''
    
    # Open the ensemble member file
    with Dataset(file, 'r') as nc:
        # For the complete dataset, only 1 lat/lon grid, mask, units, and
        #   valid date are needed. Load them only once.
        if str(EM) == '00':
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            
            FH = nc.variables['FH'][:]
            
            VD  = nc.variables['VD'][:]
            
            mask = nc.variables['mask'][:,:]
            
            units = nc.variables['units'][:]
        else:
            pass
        
        # Load the variable data. Note that it is in lat x lon x time format
        Var = nc.variables[str(SName)][:,:,:]
        
    #######################
    ### End of Function ###
    #######################
    
    if str(EM) == '00':
        return Var, FH, lat, lon, VD, mask, units
    else:
        return Var
    
#%%
# Fucntion to load all ensemble members
def load_ensembles(SName, EM, Year, Month, Day, ModelRun):
    '''

    '''
    
    # The path to the individual ensemble member data
    path = './Data/Ensemble_Members/'
    
    # Load a temporary file
    tmpFile = 'GEFS_' + str(SName) + '_' + str(Year) + str(Month) + str(Day) +\
        '_' + str(ModelRun) + '_00.nc'
    with Dataset(path +  tmpFile, 'r') as nc:
        tmp = nc.variables[str(SName)][:,:,:]
        
    # Collect dimensions
    EMlen = len(EM)
    J, I, T = tmp.shape
    
    # Initialize variables
    EnsData = np.ones((J, I, T, EMlen)) * np.nan
    FH = np.ones((T, )) * np.nan
    lat = np.ones((J, )) * np.nan
    lon = np.ones((I, )) * np.nan
    mask = np.ones((J, I)) * np.nan
    
    # Load all ensemble files
    TruncFileName = 'GEFS_' + str(SName) + '_' + str(Year) + str(Month) + str(Day) +\
        '_' + str(ModelRun)
    for n, file in enumerate(glob(path + TruncFileName + '_*.nc')):
        if file == path + tmpFile:
            EnsData[:,:,:,n], FH[:], lat[:], lon[:], VD, mask, units =\
            load_ensemble(SName, file, EM[n])
        else:
            EnsData[:,:,:,n] = load_ensemble(SName, file, EM[n])
            
    return EnsData, FH, lat, lon, VD, mask, units


