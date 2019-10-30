#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 14:57:15 2019

@author: Rarrell
"""

#%%
#####################################
### Import some libraries ###########
#####################################

import os, sys, warnings
import numpy as np
from netCDF4 import Dataset
from glob import glob

#%% Main function
def main():
    script = sys.argv[0]
    VarName      = sys.argv[1]
    VarSName     = sys.argv[2]
    TypeOfHeight = sys.argv[3]
    Parameters = sys.argv[4]
    
    Years, Months, Days, ModelRuns, Source = np.loadtxt(Parameters, dtype = str,
                                     delimiter = ' ', unpack = True)
    
    Var, VarFH, lat, lon, VarVD, mask, units = load_multiple_nc(VarSName, Years,
                                                                Months, Days,
                                                                ModelRuns)
    
    write_nc(VarSName, VarName, Var, VarFH, lat, lon, VarVD, mask, units,
             Years, Months, Days, TypeOfHeight)

#%%
##################################
### load_nc_full Function ########
##################################
    
def load_nc_full(VarSName, file):
    '''
    This function loads one of the temporary .nc files created in the grb_to_nc
      script. Some information is only loaded once as it is only needed once
      for the complete dataset.
    
    Inputs:
        VarSName - Variable short name
        file - The name of the .nc file
        VarFH - The forecast hour for the data
    '''
#    path = '/Users/Rarrell/Desktop/Thesis_Research/tmp/'
#    filename = 'GFS_' + str(VarSName) + '_' +\
#                  str(VarFH) + '.nc'
    
    # Open the temporary .nc file
    with Dataset(file, 'r') as nc:
        # Load the latitude and longitude
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        
        # Load the valid date
        VD  = nc.variables['VD'][:]
        
        # Load the mask
        mask = nc.variables['mask'][:,:]
        
        # Load the units
        units = nc.variables['units'][:]

        # Load the forecast hour
        FH = nc.variables['FH'][:]
        
        # Load the variable data. Note that it is in lat x lon format
        Var = nc.variables[str(VarSName)][:,:]
    
    #######################
    ### End of Function ###
    #######################
    
    return Var, FH, lat, lon, VD, mask, units
    
#%%
#####################################
### load_nc_partial Function ########
#####################################
    
def load_nc_partial(VarSName, file):
    '''
    This function loads one of the temporary .nc files created in the grb_to_nc
      script. Some information is only loaded once as it is only needed once
      for the complete dataset.
    
    Inputs:
        VarSName - Variable short name
        file - The name of the .nc file
        VarFH - The forecast hour for the data
    '''
#    path = '/Users/Rarrell/Desktop/Thesis_Research/tmp/'
#    filename = 'GFS_' + str(VarSName) + '_' +\
#                  str(VarFH) + '.nc'
    
    # Open the temporary .nc file
    with Dataset(file, 'r') as nc:
        # Load the forecast hour
        VD = nc.variables['VD'][:]
        
        # Load the variable data. Note that it is in lat x lon format
        Var = nc.variables[str(VarSName)][:,:]
    
    #######################
    ### End of Function ###
    #######################
    
    return Var, VD

#%%
####################################
### load_multiple_nc Function ######
####################################

def load_multiple_nc(VarSName, Years, Months, Days, ModelRuns):
    '''
    This function loads all the temporary .nc files created in the grb_to_nc
      script and places the information in a single set of variables.
      
    Inputs:
        VarSName - Variable short name
        FH - Vector of all forecast hours in the model run
    
    Outputs:
        Var - Variable data in a lat x lon x forecast hour/time format.
        VarFH - Vector of all forecast hours in the model run
        lat - Vector of all latitudes
        lon - Vector of all longitudes
        VarVD - Valid date for the model run
        mask - Gridded land-sea mask (lat x lon format)
        units - Variable units
    '''

#    path = '/Users/Rarrell/Desktop/Thesis_Research/tmp/'
    
    # Load a temporary file to collect the time and space dimensions
    path = './Data/tmp/'
    
    # Load a temporary file
    with Dataset(path + 'GFS_' + VarSName + '_' + str(Years[0]) + str(Months[0]) +\
                 str(Days[0]) + '_' + str(ModelRuns[0]) + '.nc', 'r') as nc:
        tmp = nc.variables[str(VarSName)][:,:]
    
    # Collect dimensions
    T = len(glob(path + 'GFS_' + str(VarSName) + '*.nc'))
    J, I = tmp.shape
    
    # Initialize variables
    VarVD = ['tmp'] * T
    lat = np.ones((J, )) * np.nan
    lon = np.ones((I, )) * np.nan
    Var = np.ones((J, I, T)) * np.nan
    mask = np.ones((J, I)) * np.nan
    
    # Load all temporary .nc files
    for n, file in enumerate(glob(path + 'GFS_' + str(VarSName) + '_*.nc')):
        if n == 0:
            Var[:,:,n], VarFH, lat[:], lon[:], VarVD[n], mask, units =\
            load_nc_full(VarSName, file)
        else:
            Var[:,:,n], VarVD[n] = load_nc_partial(VarSName, file)
    
    #######################
    ### End of Function ###
    #######################
    
    return Var, VarFH, lat, lon, VarVD, mask, units

#%%
#############################
### write_nc Function #######
#############################

def write_nc(VarSName, VarName, Var, VarFH, lat, lon, VarVD, mask, units,
             Years, Months, Days, TypeOfHeight):
    '''
    This function takes the complete model run data for a single variable
      collected in load_multiple_nc and writes it to a single .nc file for user
      use.
      
    Inputs:
        VarSName - Variable short name
        VarName  - Variable long name
        Var - Main variable data (lat x lon x forecast hour/tim format)
        VarFH - Vector of forecast hours
        lat - Vector of latitudes
        lon - Vector of longitudes
        VarVD - Variable valid date
        mask  - Gridded land-sea mask (lat x lon format)
        units - Units of the variable
        year  - Year the model run is valid for
        month - Month the model is valid for
        day   - Day the model is valid for
        ModelRun - The model run the model is valid for
        TypeOfHeight - The type of height the variable is located at
    '''
    
    # Create the .nc file name based on the date, model run, and variable short
    #   name.
    # Note that ed stands for end date
    filename = 'GFS_' + str(VarSName) + '_ed_' + str(Years[-1]) +\
               str(Months[-1]) + str(Days[-1]) + '_.nc'
#    path = '/Users/Rarrell/Desktop/Thesis_Research/Data/'
    path = './Data/Synthetic_Data/'
    
    # Collect the variable dimensions
    J, I , T = Var.shape
    
    # Begin writing the .nc file
    with Dataset(path + filename, 'w', format = 'NETCDF4') as nc:
        # Write the main description for the model variable data.
        nc.description = 'GFS forecast data for ' + VarName + ' valid from ' +\
                         str(VarVD[0]) + 'to' + str(VarVD[-1]) + '. The variable ' +\
                         'is given a minimum 4 week period and is composed of 3 ' +\
                         'hour forecasts across multiple model runs.\n' +\
                         'Variable: ' + VarName +'(' + units + '). This variable ' +\
                         'is in the format of lat x lon x time.\n' +\
                         'mask: Land - sea mask. sea = 0, land = 1.\n' +\
                         'FH: Forecast hour (hour). The hour at which the ' +\
                         'forecast is made.\n' +\
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
        nc.createVariable('FH', np.str)
        nc.variables['FH'][0] = np.str(VarFH)
        
        # Create a variable for the forecast hour.
        nc.createVariable('VD', np.str, ('time', ))
        for n in range(len(VarVD)):
            nc.variables['VD'][n] = np.str(VarVD[n])
        
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
