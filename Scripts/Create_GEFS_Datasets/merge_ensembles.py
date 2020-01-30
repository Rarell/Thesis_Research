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
    EM = np.loadtxt(EMfile, dtype = str, unpack = True)
    
    VarName, TypeOfHeight, Height, VarSName = np.loadtxt(VarParameters,
                                                         usecols = (0,1,2,3),
                                                         dtype = str, delimiter = ',',
                                                         unpack = True)
    
    Year, Month, Day, ModelRun, Source = np.loadtxt(Parameters, 
                                                    usecols = (1,2,3,4,5), dtype = str, 
                                                    delimiter = ',', unpack = True)
    
    # Load all the ensemble members
    EnsembleData, FH, lat, lon, VD, mask, units = load_ensembles(VarSName, EM, 
                                                                    Year, Month, Day, ModelRun)
    
    # Merge the ensembles
    EnsembleMean, Spread = merge_ensembles(EnsembleData)
    
    # Write the .nc files to contain the merged ensemble data
    write_nc(VarSName, VarName, EnsembleMean, FH, lat, lon, VD, mask, units, Year,
             Month, Day, ModelRun, TypeOfHeight, TYPE = 'mean')
    
    write_nc(VarSName, VarName, Spread, FH, lat, lon, VD, mask, units, Year,
             Month, Day, ModelRun, TypeOfHeight, TYPE = 'spread')
    
    
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

#%%
# Function to "merge" the ensembles (find their mean and spread)

def merge_ensembles(EnsData):
    '''
    
    '''
    
    # Compute the mean
    EnsMean = np.nanmean(EnsData[:,:,:,:], axi = -1)
    
    # Compute the spread (standard deviation w.r.t. the ensemble mean)
    Spread = np.nanstd(EnsData[:,:,:,:], axis = -1)
    
    # Perform bias correction
    
    return EnsMean, Spread


#%%
# Function to write the ensemble mean and spread information
def write_nc(VarSName, VarName, Var, FH, lat, lon, VD, mask, units,
             Year, Month, Day, ModelRun, TypeOfHeight, TYPE = 'mean'):
    '''

    '''    

    filename = 'GEFS_' + str(VarSName) + '_' + str(Year) + str(Month) +\
               str(Day) + '_' + str(ModelRun) + '_' + str(TYPE) + '.nc'
#    path = '/Users/Rarrell/Desktop/Thesis_Research/Data/'
    path = './Data/'
    
    # Collect the variable dimensions
    J, I , T = Var.shape
    
    # Begin writing the .nc file
    with Dataset(path + filename, 'w', format = 'NETCDF4') as nc:
        # Write the main description for the model variable data.
        nc.description = 'GEFS forecast data for ' + VarName + ' valid for ' +\
                         Day + '-' + Month + '-' + Year +\
                         '. This forecast is either for the ensemble mean or spread. The variable is ' +\
                         'given for all forecast hours.\n' +\
                         'Variable: ' + VarName +'(' + units + '). This variable ' +\
                         'is in the format of lat x lon x time/forecast hour.\n' +\
                         'mask: Land - sea mask. sea = 0, land = 1.\n' +\
                         'FH: Forecast hour (hour). The hour at which the ' +\
                         'forecast is made. E.g., FH = 3 is the 3 hour forecast ' +\
                         'from ' + Day + '-' + Month + '-' + Year + '.\n' +\
                         'VD: Date for which the GFS model run is valid.\n' +\
                         'units: The units of the variable.\n' +\
                         'lat: Latitude of the grid (vector form).\n' +\
                         'lon: Longitude of the grid (vector form).'
        
        # Create the spatial and temporal dimensions
        nc.createDimension('lat', size = J)
        nc.createDimension('lon', size = I)
        nc.createDimension('time', size = T)
        
        # Create the lat and lon variables
        nc.createVariable('lat', lat.dtype, ('lat', ))
        nc.createVariable('lon', lon.dtype, ('lon', ))
        
        nc.variables['lat'][:] = lat[:]
        nc.variables['lon'][:] = lon[:]
        
        # Create the valid date variable
        nc.createVariable('VD', np.str)
        nc.variables['VD'][0] = np.str(VD)
        
        # Create a variable for the forecast hour.
        nc.createVariable('FH', int, ('time', ))
        nc.variables['FH'][:] = FH[:]
        
        # Create the variable for the main variable
        nc.createVariable(str(VarSName), Var.dtype, ('lat', 'lon', 'time'))
        nc.setncatts({'long_name' : VarName, 'units' : units, 
                      'level_desc' : TypeOfHeight})
        nc.variables[str(VarSName)][:,:,:] = Var[:,:,:]
        
        # Create the units variable for the variable
        nc.createVariable('units', np.str)
        nc.variables['units'][0] = np.str(units)
        
        # Create the mask variable
        nc.createVariable('mask', mask.dtype, ('lat', 'lon'))
        nc.variables['mask'][:,:] = mask[:,:]
        
        #######################
        ### End of Function ###
        #######################
        
        
#%%
#########################################
### Call and Run the Main Function ######
#########################################
        
if __name__ == '__main__':
    main()
    
#####################
### End of Script ###
#####################