#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 13:07:36 2019

@author: Rarrell

This script is designed as part of the create_synthetic_dataset.sh program.
This script takes in date and variable information for a GFS grib file
  downloaded in the GFS_Download script, and opens that GFS file and extracts
  the sepcified variable and some associated information. That variable data
  and information is then dumped into a temporary netcdf file for later use in
  the create_sythetic_dataset.sh program.
  
Arguments:
    VarName      - The name of the variable to be extracted
    TypeOfHeight - The type of height the variable is located in
    Height       - The height level (unit determined by TypeOfHeight) the
                   variable is located at
    FH       - The forecast hour of the grib file
    Year     - The year the GFS run was made
    Month    - The month the GFS run was made
    Day      - The day the GFS run was made
    ModelRun - The model run of the GFS run
    Source   - The location of grib file came from (this effects the location
               and nomenclature of the grib file)
"""

#%% 
#####################################
### Import some libraries ###########
#####################################

import sys, os, warnings
import numpy as np
import pygrib
from netCDF4 import Dataset
from glob import glob

#%% 
############################
### Main Function ##########
############################

def main():
    '''
    This is the main function of the script. It collects the arguements and
      uses them in functions to collect the variable data and information and
      dumps them in a netcdf file.
    '''
    
    # Collect the arguments
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
    
    # Extract the specified data from the grib file
    VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName, Mask =\
        grbExtractData(VarName, TypeOfHeight, Height, FH, Year, Month, Day,
                       ModelRun, Source)
    
    # Dump the data into a netcdf file
    create_tmp_nc(VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName,
                  Mask, Year, Month, Day, ModelRun, TypeOfHeight)
    
    
        
    
#%%  
######################################
### grbExtractData Function ##########
######################################
    
def grbExtractData(VarName, TypeOfHeight, Height, FH, Year, Month, Day, 
                   ModelRun, Source):
    '''
    This function constructs a GFS filename given serveral parameters and
      collects a specified variable and associated information.
      
    Inputs:
        VarName      - String. The name of the variable to be extracted from 
                       the GFS file.
        TypeOfHeight - String. The type of height VarName is located with.
        Height       - String. The height VarName is located at.
        FH    - String. The forecast hour of the GFS file being examined.
        Year  - String. The year of the GFS file being examined.
        Month - String. The month of the GFS file being examined.
        Day    - String. The day of the GFS file being examined.
        ModelRun - String. The model run of the GFS file being examined.
        Source   - String. Where the GFS file came from. Defines how the filename
                   is constructed.
                   
    Outputs:
        VarData - Float array. The gridded data of VarName.
        lat - Float array. The gridded latitude associated with VarData.
        lon - Float array. The gridded longitude associated with VarData.
        VarUnits - String. The units of the VarData.
        VarFH - Int. The valid forecast hour of VarData.
        VarVD - String. The date for which VarData is valid.
        VarName  - String. The long name of VarData.
        VarSName - String. The short name of VarData.
        Mask - Float array. Land-sea mask of VarData.
    '''
    
    # Define the path where the data is located
    path = './Data/tmp/'
    
    # Construct the GFS file name based on where it was downloaded from
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
    
    # Open the grib file    
    gfs = pygrib.open(grb_file)
    
    # Collect the message of the specified variable
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
    
    # Close the grib file
    gfs.close()
    
    ###################
    # End of Function #
    ###################
    
    return VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName, Mask

#%% 
#################################
### create_tmp_nc Function ######
#################################
    
def create_tmp_nc(VarData, lat, lon, VarUnits, VarFH, VarVD, VarName,
                  VarSName, Mask, Year, Month, Day, ModelRun, TypeOfHeight):
    '''
    This function takes the variable data and information extracted from the
      GFS grib files (in grbExtractData function) and puts it in a temporary
      netcdf file to hold that data.
      
    Inputs:
        VarData - The gridded variable data
        lat - The gridded latitude
        lon - The gridded longitude
        VarUnits - The variable units
        VarFH - The forecast hour for the variable
        VarVD - The valid date for the variable
        VarName  - The variable long name
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
        #   files are deleted at the end of the create_synthetic_dataset.sh program.
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

#%%
#########################################
### Call and Run the Main Function ######
#########################################
   
if __name__ == '__main__':
    main()
    
#####################
### End of Script ###
#####################