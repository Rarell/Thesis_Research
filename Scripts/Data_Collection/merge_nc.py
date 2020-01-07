#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:09:42 2019

@author: Rarrell

This is a script designed as a component of the extract_GFS_variable.sh program.
This is a script designed to take the temporary .nc files created in the
  grb_to_nc script and merge them into a single .nc file over all forecast hours.
  This is the final output (and goal) of the extract_GFS_variable.sh program
  for use by the user.
  
Arguements:
    parameters - The temporary .txt file from input_model_information. It 
                 contains request, year, month, day, source, var, TypeOfHeight,
                 height, model_run, and VarSName in that order.
    ForecastHour - A .txt file containing all forecast hours in the GFS model run.
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
    This is the main function. It loads the arguements, loads the temporary
      .nc files, and writes all the data in a single, combined .nc file.
    '''
    
    # Load arguements
    script = sys.argv[0]
    VarParameters = sys.argv[1]
    Parameters    = sys.argv[2]
    ForecastHour  = sys.argv[3]
    
    # Unpack the .txt file to obtain the model information
    VarName, TypeOfHeight, Height, VarSName = np.loadtxt(VarParameters,
                                                         usecols = (0,1,2,3),
                                                         dtype = str, delimiter = ',',
                                                         unpack = True)
    
    Year, Month, Day, ModelRun, Source = np.loadtxt(Parameters, 
                                                    usecols = (1,2,3,4,5), dtype = str, 
                                                    delimiter = ',', unpack = True)
    
    
    # Unpack the .txt file containing all the forecast hours
    FH = np.loadtxt(ForecastHour, dtype = str, delimiter = ',')
    
    # Load all the temporary .nc files and place them in a single set of variables
    Var, VarFH, lat, lon, VarVD, mask, units = load_multiple_nc(VarSName, FH)
    
    # Write a single .nc file containing the data for all forecast hours
    write_nc(VarSName, VarName, Var, VarFH, lat, lon, VarVD, mask, units,
             Year, Month, Day, ModelRun, TypeOfHeight)
    

#%%
#############################
### load_nc Function ########
#############################
    
def load_nc(VarSName, file, VarFH):
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
        # For the complete dataset, only 1 lat/lon grid, mask, units, and
        #   valid date are needed. Load them only once.
        if file == './Data/tmp/GFS_' + str(VarSName) + '_003.nc':
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            
            VD  = nc.variables['VD'][:]
            
            mask = nc.variables['mask'][:,:]
            
            units = nc.variables['units'][:]
            
        else:
            pass
        
        # Load the forecast hour
        FH = nc.variables['FH'][:]
        
        # Load the variable data. Note that it is in lat x lon format
        Var = nc.variables[str(VarSName)][:,:]
    
    #######################
    ### End of Function ###
    #######################
    
    if file == './Data/tmp/GFS_' + str(VarSName) + '_003.nc':
        return Var, FH, lat, lon, VD, mask, units
    else:
        return Var, FH
    

#%%
####################################
### load_multiple_nc Function ######
####################################

def load_multiple_nc(VarSName, FH):
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
    with Dataset(path + 'GFS_' + str(VarSName) + '_006.nc', 'r') as nc:
        tmp = nc.variables[str(VarSName)][:,:]
    
    # Collect dimensions
    T = len(glob(path + 'GFS_' + str(VarSName) + '_*.nc'))
    J, I = tmp.shape
    
    # Initialize variables
    VarFH = np.ones((T, )) * np.nan
    lat = np.ones((J, )) * np.nan
    lon = np.ones((I, )) * np.nan
    Var = np.ones((J, I, T)) * np.nan
    mask = np.ones((J, I)) * np.nan
    
    # Load all temporary .nc files
    for n, file in enumerate(glob(path + 'GFS_' + str(VarSName) + '_*.nc')):
        if file == path + 'GFS_' + str(VarSName) + '_003.nc':
            Var[:,:,n], VarFH[n], lat[:], lon[:], VarVD, mask, units =\
            load_nc(VarSName, file, FH[n])
        else:
            Var[:,:,n], VarFH[n] = load_nc(VarSName, file, FH[n])
    
    # Note because of a new nomeclature of the tmp files, glob loads things
    # in a random order. Sort everything to be in chronological order.
    VarFH2 = np.ones((T, )) * np.nan
    Var2   = np.ones((J, I, T)) * np.nan
    for n, i in enumerate(np.sort(VarFH)):
        ind = np.where(VarFH == i)[0]
        
        VarFH2[n]   = VarFH[ind[0]]
        Var2[:,:,n] = Var[:,:,ind[0]]
    
    #######################
    ### End of Function ###
    #######################
    
    return Var2, VarFH2, lat, lon, VarVD, mask, units

#%%
#############################
### write_nc Function #######
#############################

def write_nc(VarSName, VarName, Var, VarFH, lat, lon, VarVD, mask, units,
             year, month, day, ModelRun, TypeOfHeight):
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
    #   name
    filename = 'GFS_' + str(VarSName) + '_' + str(year) + str(month) +\
               str(day) + '_' + str(ModelRun) + '.nc'
#    path = '/Users/Rarrell/Desktop/Thesis_Research/Data/'
    path = './Data/'
    
    # Collect the variable dimensions
    J, I , T = Var.shape
    
    # Begin writing the .nc file
    with Dataset(path + filename, 'w', format = 'NETCDF4') as nc:
        # Write the main description for the model variable data.
        nc.description = 'GFS forecast data for ' + VarName + ' valid for ' +\
                         day + '-' + month + '-' + year + '. The variable is ' +\
                         'given for all forecast hours.\n' +\
                         'Variable: ' + VarName +'(' + units + '). This variable ' +\
                         'is in the format of lat x lon x time/forecast hour.\n' +\
                         'mask: Land - sea mask. sea = 0, land = 1.\n' +\
                         'FH: Forecast hour (hour). The hour at which the ' +\
                         'forecast is made. E.g., FH = 3 is the 3 hour forecast ' +\
                         'from ' + day + '-' + month + '-' + year + '.\n' +\
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
        nc.variables['VD'][0] = np.str(VarVD)
        
        # Create a variable for the forecast hour.
        nc.createVariable('FH', int, ('time', ))
        nc.variables['FH'][:] = VarFH[:]
        
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