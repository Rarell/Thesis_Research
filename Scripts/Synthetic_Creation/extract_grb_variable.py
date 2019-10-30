#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:07:36 2019

@author: Rarrell
"""

#%% Import libraries
import sys, os, warnings
import numpy as np
import pygrib
from netCDF4 import Dataset
from glob import glob

#%% Main function
def main():
    '''
    '''
    script = sys.argv[0]
    VarName      = sys.argv[1]
    TypeOfHeight = sys.argv[2]
    Height       = sys.argv[3]
    FH    = sys.argv[4]
    Year  = sys.argv[5]
    Month = sys.argv[6]
    Day   = sys.argv[7]
    ModelRun = sys.argv[8]
    Source   = sys.argv[9]
    
    VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName, Mask =\
        grbExtractData(VarName, TypeOfHeight, Height, FH, Year, Month, Day,
                       ModelRun, Source)
        
    create_tmp_nc(VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName,
                  Mask, Year, Month, Day, ModelRun, TypeOfHeight)
    
    
        
    
    
#%% grbExtractData function
def grbExtractData(VarName, TypeOfHeight, Height, FH, Year, Month, Day, 
                   ModelRun, Source):
    '''
    '''
    path = './Data/tmp/'
    
    if Source == 'Request':
        filename = 'gfs_3_' + str(Year) + str(Month) + str(Day) + '_' +\
                   str(ModelRun) +'00_0' + str(FH) + '.grb2'
        
        path = './Data/gfs_3_' + str(Year) + str(Month) + str(Day) +\
               str(ModelRun) + '.g2/'
        
        grb_file = path + filename
        
    elif Source == 'NCEI':
        filename = 'gfs_3_' + str(Year) + str(Month) + str(Day) + '_' +\
                   str(ModelRun) + '00_' + str(FH) + '.grb2'
        
        grb_file = path + filename
    
    else:
        filename = 'gfs.t' + str(ModelRun) + 'z.pgrb2.1p00.f' + str(FH)
        
        grb_file = path + filename
        
    gfs = pygrib.open(grb_file)
    
    VarMSG = gfs.select(name = VarName, typeOfLevel = TypeOfHeight, 
                        level = int(Height))
    
    # Collect variable data and information
    VarData  = VarMSG[0].values
    VarUnits = VarMSG[0].units
    VarFH    = VarMSG[0].forecastTime
    VarVD    = VarMSG[0].validDate
    VarName  = VarMSG[0].name
    VarSName = VarMSG[0].shortName
    
    # Collect lat, lon grid
    lat, lon = VarMSG[0].latlons()
    
    # Collect land-sea mask
    MaskMSG = gfs.select(name = 'Land-sea mask',
                      typeOfLevel = 'surface')
    Mask    = MaskMSG[0].values
    
    # Variable values are a mask by default. Convert sea values to 0 to make
    #   the variable a non-masked array.
    JMask, IMask = np.where(Mask == 0)
    VarData[JMask, IMask] = 0
    
    # Close the .grb2 file
    gfs.close()
    
    return VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName, Mask

#%% 
#################################
### create_tmp_nc Function ######
#################################
    
def create_tmp_nc(VarData, lat, lon, VarUnits, VarFH, VarVD, VarName,
                  VarSName, Mask, Year, Month, Day, ModelRun, TypeOfHeight):
    '''
    This function takes the variable data and information extracted from the
      GFS grib files (in grb_extract_data function) and puts it in a temporary
      netcdf file to hold that data.
      
    Inputs:
        VarData - The gridded variable data
        lat - The gridded latitude
        lon - The gridded longitude
        VarUnits - The variable units
        VarFH - The forecast hour for the variable
        VarVD - The valid date for the variable
        VarName - The variable long name
        VarSName - The variable short name
        Mask - The gridded land-sea mask
        TypeOfHeight - The type of height the variable is located at
    '''
    
    # Define the path and temporary file name
    filename = 'GFS_' + VarSName + '_' + str(Year) + str(Month) +\
               str(Day) + '_' + str(ModelRun) + '.nc'
    path = './Data/tmp/'
#    path = '/Users/Rarrell/Desktop/Thesis_Research/tmp/'
    
    # Collect the dimensions of the gridded data
    J, I = VarData.shape
    
    # Open the netcdf file and begin writing.
    with Dataset(path + filename, 'w', format = 'NETCDF4') as nc:
        # Write a short description. (This is more for completeness. Temporary
        #   files are deleted at the end of the extract_GFS_variable.sh program.
        #   So this will most likely not be seen.)
        nc.description = "This is a temporary file will hold GFS data " +\
                         "for " + VarName + ". This will be for the " +\
                         str(VarFH) + " forecast hour. This file will" +\
                         " only exist until all the forecast hour files" +\
                         " are compressed together at the end of the " +\
                         "bash script."
        
        # Create the latitude and longitude dimensions
        nc.createDimension('lat', size = J)
        nc.createDimension('lon', size = I)
        
        # Create the latitude and longitude variables and fill them with the
        #   (vector) data
        nc.createVariable('lat', lat.dtype, ('lat', ))
        nc.createVariable('lon', lon.dtype, ('lon', ))
        
        nc.variables['lat'][:] = lat[:,0]
        nc.variables['lon'][:] = lon[0,:]
        
        # Create and fill the forecast hour information. Note that the lack of
        #   of a specified dimension defaults the variable to a scalar.
        nc.createVariable('FH', int)
        nc.variables['FH'][:] = VarFH
        
        # Create and fill the valid date information. Note that the lack of
        #   of a specified dimension defaults the variable to a scalar.
        nc.createVariable('VD', np.str)
        nc.variables['VD'][0] = np.str(VarVD)
        
        # Create and fill the main variable. A short description is also added.
        nc.createVariable(VarSName,VarData.dtype, ('lat', 'lon'))
        nc.setncatts({'long_name' : VarName, 'units' : VarUnits,
                      'level_desc' : TypeOfHeight})
        nc.variables[VarSName][:,:] = VarData[:,:]
        
        # Create and fill the variable unit information. Note that the lack of
        #   of a specified dimension defaults the variable to a scalar.
        nc.createVariable('units', np.str)
        nc.variables['units'][0] = VarUnits
        
        # Create and fill the land-sea mask variable.
        nc.createVariable('mask', Mask.dtype, ('lat', 'lon'))
        nc.variables['mask'][:,:] = Mask[:,:]
        
        #######################
        ### End of Function ###
        #######################

#%% Call the main function
if __name__ == '__main__':
    main()