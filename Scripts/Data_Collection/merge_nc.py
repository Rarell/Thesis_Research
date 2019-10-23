#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:09:42 2019

@author: Rarrell

This is a script designed to take the temporary .nc files created and merge
them into a single .nc file over all forecast hours.
"""

#%% Load modules
import os, sys, warnings
import numpy as np
from netCDF4 import Dataset
from glob import glob

#%% Define the main function
def main():
    script = sys.argv[0]
    parameters = sys.argv[1]
    ForecastHour = sys.argv[2]
    
    year, month, day, source, VarName, TypeOfHeight =\
        np.loadtxt(parameters, usecols = np.arange(1, 6+1), dtype = str, 
                   delimiter = ',')
    
    ModelRun, VarSName = np.loadtxt(parameters, usecols = np.arange(8, 9+1), 
                                    dtype = str, delimiter = ',')
    
    FH = np.loadtxt(ForecastHour, dtype = str, delimiter = ',')
    
    Var, VarFH, lat, lon, VarVD, mask, units = load_multiple_nc(VarSName, FH)
    
    write_nc(VarSName, VarName, Var, VarFH, lat, lon, VarVD, mask, units,
             year, month, day, ModelRun, TypeOfHeight)
    

#%% Define a function to load in netcdf files.
def load_nc(VarSName, file, VarFH):
#    path = '/Users/Rarrell/Desktop/Thesis_Research/tmp/'
#    filename = 'GFS_' + str(VarSName) + '_' +\
#                  str(VarFH) + '.nc'
    with Dataset(file, 'r') as nc:
        if str(VarFH) == '003':
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            
            VD  = nc.variables['VD'][:]
            
            mask = nc.variables['mask'][:,:]
            
            units = nc.variables['units'][:]
            
        else:
            pass
        
        FH = nc.variables['FH'][:]
    
        Var = nc.variables[str(VarSName)][:,:]
    
    if str(VarFH) == '003':
        return Var, FH, lat, lon, VD, mask, units
    else:
        return Var, FH
    

#%% Define a function to load the .nc files.
def load_multiple_nc(VarSName, FH):

#    path = '/Users/Rarrell/Desktop/Thesis_Research/tmp/'
    path = './Data/tmp/'
    # Load in a temporary file.
    with Dataset(path + 'GFS_' + str(VarSName) + '_6.nc', 'r') as nc:
        tmp = nc.variables[str(VarSName)][:,:]
    
    # Initialize some values.
    T = len(glob(path + 'GFS_' + str(VarSName) + '_*.nc'))
    J, I = tmp.shape
    
    VarFH = np.ones((T, )) * np.nan
    lat = np.ones((J, )) * np.nan
    lon = np.ones((I, )) * np.nan
    Var = np.ones((J, I, T)) * np.nan
    mask = np.ones((J, I)) * np.nan
    
    # Load in all data points
    for n, file in enumerate(glob(path + 'GFS_' + str(VarSName) + '_*.nc')):
        if FH[n] == '003':
            Var[:,:,n], VarFH[n], lat[:], lon[:], VarVD, mask, units =\
            load_nc(VarSName, file, FH[n])
        else:
            Var[:,:,n], VarFH[n] = load_nc(VarSName, file, FH[n])
    
    return Var, VarFH, lat, lon, VarVD, mask, units

#%% Define a function to write the combined data.
def write_nc(VarSName, VarName, Var, VarFH, lat, lon, VarVD, mask, units,
             year, month, day, ModelRun, TypeOfHeight):
    filename = 'GFS_' + str(VarSName) + '_' + str(year) + str(month) +\
               str(day) + '_' + str(ModelRun) + '.nc'
#    path = '/Users/Rarrell/Desktop/Thesis_Research/Data/'
    path = './Data/'
    
    J, I , T = Var.shape
    
    with Dataset(path + filename, 'w', format = 'NETCDF4') as nc:
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
        
#%% Call the main functions
        
if __name__ == '__main__':
    main()