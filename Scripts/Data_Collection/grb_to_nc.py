#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 15:48:20 2019

@author: Rarrell

This is a file designed to take in various inputs to create a url to download
a GFS grib file from the NCEI or NCEP, extract the desired variable at the 
desired height, and dump it into a temporary netcdf file.
"""

#%% Import needed modules
import sys, os, warnings
import wget
import numpy as np
import pygrib
from netCDF4 import Dataset
from glob import glob
from datetime import datetime


#%% Define the function that will extract the grb data.
def grb_extract_data(var, TypeOfHeight, height = 0, request = False, 
                     year = '2019', month = '06', day = '25', model_run = '00', 
                     forecast_hour = '000', source = 'NCEI'):
    if request is True:
        filename = 'gfs_3_' + year + month + day + '_' + '00' + model_run +\
                   '_' + forecast_hour + '.grb2'
        path = '/Users/Rarrell/Downloads/gfs_3_' + year + month + day +\
               model_run + '.g2/'
        grb_file = path + filename
    else:
        filename = 'gfs_3_' + year + month + day + '_' + '00' + model_run +\
                   forecast_hour + '.grb2'
        path = '/Users/Rarrell/Downloads/tmp/'
        if source is 'NCEI':
            url_start = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
            url_end   = year + month + '/' + year + month + day +'/' + filename
            url = url_start + url_end
        elif source is 'NCEP':
            urlstuff
        else:
            error message
        
        wget.download(url, path + filename)
    
    gfs = pygrib.open(path + filename)
    
    VarMsg = gfs.select(name = var, typeOfLevel = TypeOfHeight, level = height)
    
    VarData  = VarMsg[0].values
    VarUnits = VarMsg[0].units
    VarFH    = VarMsg[0].forecastTime
    VarVD    = VarMsg[0].validDate
    VarName  = VarMsg[0].name
    VarSName = varMsg[0].shortName
    
    lat, lon = VarMsg[0].latlons()
    
    MaskMsg = gfs.select(name = 'Land-sea mask',
                      typeOfLevel = 'surface')
    Mask    = MaskMsg[0].values
    
    JMask, IMask = np.where(Mask == 0)
    VarData[JMask, IMask] = 0
    return VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName, Mask


#%% Define the function that will write the NetCDF file.
    def create_tmp_nc(VarData, lat, lon, VarUnits, VarFH, VarVD, VarName,
                      VarSName, Mask, TypeOfHeight):
        filename = 'GFS_' + VarSName + '_' + str(VarVD) + '_' +\
                   str(VarFH) + '.nc'
        path = 'Users/Rarrell/Desktop/Thesis_Research/tmp/'
        
        J, I = VarData.shape
        
        with Dataset(path + filename, 'w', format = 'netCDF4') as nc:
            nc.description = "This is a temporary file will hold GFS data " +\
                             "for " + VarName + ". This will be for the " +\
                             str(VarFH) + " forecast hour. This file will" +\
                             " only exist until all the forecast hour files" +\
                             " are compressed together at the end of the " +\
                             "bash script."
            
            nc.createDimension('lat', size = J)
            nc.createDimension('lon', size = I)
            
            nc.createVariable('lat', lat.dtype, ('lat', ))
            nc.createVariable('lon', lon.dtype, ('lon', ))
            
            nc.variables['lat'][:] = lat[:,0]
            nc.variables['lon'][:] = lon[0,:]
            
            nc.createVariable('FH', int)
            nc.variables['FH'][:] = VarFH
            
            nc.createVariable('VD', np.str)
            nc.variables['VD'][0] = np.str(VarVD)
            
            nc.createVariable(VarSName,VarData.dtype, ('lat', 'lon'))
            nc.setncatts({'long_name' : VarName, 'units' : VarUnits,
                          'level_desc' : TypeOfHeight})
            nc.variables[VarSName][:,:] = VarData[:,:]
            
            nc.createVariable('units', np.str)
            nc.variables['units'][0] = VarUnits
            
            nc.createVariable('mask', Mask.dtype, ('lat', 'lon'))
            nc.variables['mask'][:,:] = Mask[:,:]
