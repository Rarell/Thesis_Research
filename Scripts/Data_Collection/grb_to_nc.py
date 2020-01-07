#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 15:48:20 2019

@author: Rarrell

This is a script designed as a component of the extract_GFS_variable.sh program.
This script uses the date and model run information obtained in the temporary
  file from the input_model_information script and forecast hour to create 
  the url the GFS file is located at (or the .g2 file if a data request was 
  made). If a url is constructed, the script then downloads the .grb2 file.
  
Once the .grb2 file is obtained, the script uses variable name, type of height,
  and height (also from the input_model_information script) to extract that
  variable, as well as the variable's short name, the lon/lat grid, and mask.
  
The variable, mask, lon/lat, forecast hour, valid date, and variable units to
  a temporary netcdf file.
  
Arguements:
    parameters - the temporary .txt file from input_model_information. It 
                 contains request, year, month, day, source, var, TypeOfHeight,
                 height, and model_run in that order.
    ForecastHour - the model forecast hour being examined.
"""

#%%
#####################################
### Import some libraries ###########
#####################################

import sys, os, warnings
import wget
import numpy as np
import pygrib
from netCDF4 import Dataset
from glob import glob
from datetime import datetime


#%%
############################
### Main Function ##########
############################
def main():
    '''
    This is the main function that loads the arguements and runs the functions
    to collect (if needed) the GFS .grb2 file, load the .grb2 file, extract
    the variable, and write it a temporary .nc file.
    '''
    
    # Load in the arguements
    script = sys.argv[0]
    ForecastHour = sys.argv[1]
    VarParameters = sys.argv[2]
    Parameters    = sys.argv[3]
    
    
    # Unpack the .txt file to obtain the model information
    Var, TypeOfHeight = np.loadtxt(VarParameters, usecols = (0, 1), 
                                       dtype = str, delimiter = ',', unpack = True)
    Height = np.loadtxt(VarParameters, usecols = 2, delimiter = ',')
    
    Request = np.loadtxt(Parameters, usecols = 0, delimiter = ',')
    Year, Month, Day, ModelRun, Source = np.loadtxt(Parameters, usecols = (1,2,3,4,5), 
                                                    dtype = str, delimiter = ',', unpack = True)

    # Download (if needed) the GFS .grb2 file, load it, and extract the
    #  the desired variable and associated information
    VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName, Mask =\
        grb_extract_data(var = Var, TypeOfHeight = TypeOfHeight, 
                         height = Height, year = Year, month = Month, day = Day,
                         model_run = ModelRun, forecast_hour = ForecastHour,
                         request = Request, source = Source)
    
    # Write the variable and information to a temporary .nc file
    create_tmp_nc(VarData = VarData, lat = lat, lon = lon, VarUnits = VarUnits,
                  VarFH = ForecastHour, VarVD = VarVD, VarName = VarName, 
                  VarSName = VarSName, Mask = Mask, TypeOfHeight = TypeOfHeight)
    
    # Append the variable short name to the temporary .txt file. Do only for
    #   for the first forecast hour so it does not keep appending.
    if str(ForecastHour) == '003':
        append_tmp(VarSName)
     

#%% Define the function that will extract the grb data.
#################################
### grb_extract_data Function ###
#################################
def grb_extract_data(var, TypeOfHeight, height, year, month, day, 
                     model_run = '00', forecast_hour = '003', 
                     request = False, source = 'NCEI'):
    '''
    This function uses the date and model run information to create the url the
      data file is located at (if request is false). The .grb2 file is then 
      downloaded. If request is true, it constructs the name of the .gr2 
      folder the GFS data is located in.
      
    The function then loads the .grb2 file and extracts the desired variable
      and associated information. The function finishes by removing the .grb2
      file, which is no longer needed after this function.
      
    Inputs:
        var - The name of the variable to be extracted
        TypeOfHeight - The type of height the variable is located at
        height - The height value the variable is located at (validity depends
                 on TypeOfHeight)
        year  - The year of the model run
        month - The month of the model run
        day   - The day of the model run
        model_run - The hour the model was initialized at
        forecast_hour - The forecast hour being considered
        request - Binary value. True if a data request was needed.
                  False otherwise.
        source  - The source of the .grb2 files if request is false (NCEI or 
                  NCEP).
    
    Outputs:
        VarData - The gridded variable data (in lat x lon)
        lat - Gridded latitude data
        lon - Gridded longtitude data
        VarUnits - The units of the variable
        VarFH - The forecast hour of the variable
        VarVD - The valid date for the variable
        VarName - The long name of the variable
        VarSName - The short name of the variable
        Mask - The gridded land-sea mask
    '''
    
    ##############################
    ### Collect the .grb2 File ###
    ##############################
    
    # If a data request was needed, construct location of the .g2 folder
    if int(request) == 1:
        filename = 'gfs_3_' + str(year) + str(month) + str(day) + '_' + '00' +\
        str(model_run) + '_' + str(forecast_hour) + '.grb2'
        path = './Data/tmp/gfs_3_' + str(year) + str(month) + str(day) +\
               str(model_run) + '.g2/'
#        path = '/Users/Rarrell/Downloads/gfs_3_' + str(year) + str(month) +\
#               str(day) + str(model_run) + '.g2/'
        grb_file = path + filename
        
    # If data has to be downloaded, construct the url (given a source) and
    #   download the .grb2 file
    else:
#        path = '/Users/Rarrell/Downloads/tmp/'
        path = './Data/tmp/'
        
        # Construct the url
        if source == np.str('NCEI'):
            filename = 'gfs_3_' + str(year) + str(month) + str(day) + '_' +\
                   str(model_run) + '00' + '_' + str(forecast_hour) + '.grb2'
            url_start = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
            url_end   = str(year) + str(month) + '/' + str(year) + str(month) +\
                        str(day) +'/' + filename
            url = url_start + url_end
        elif source == 'NCEP':
            filename  = 'gfs.t' + str(model_run) + 'z.pgrb2.1p00.f' +\
                        str(forecast_hour)
            url_start = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
            url_end   = 'gfs.' + str(year) + str(month) + str(day) + '/' +\
                        str(model_run) + '/' + filename
            url = url_start + url_end
        else:
            print('Error: The data should come from NCEI or NCEP.')
            sys.exit(1) # Exit the program if this happens.
        
        # Download the .grb2 data
        wget.download(url, path + filename)
        
        grb_file = path + filename
    
    #######################################
    ### Load the .grb2 and Extract Data ###
    #######################################
    
    # Open the .grb2 file
    gfs = pygrib.open(grb_file)
    
    # Find the message for the variable
    VarMsg = gfs.select(name = var, typeOfLevel = TypeOfHeight, level = height)
    
    # Collect variable data and information
    VarData  = VarMsg[0].values
    VarUnits = VarMsg[0].units
    VarFH    = VarMsg[0].forecastTime
    VarVD    = VarMsg[0].validDate
    VarName  = VarMsg[0].name
    VarSName = VarMsg[0].shortName
    
    # Collect lat, lon grid
    lat, lon = VarMsg[0].latlons()
    
    # Collect land-sea mask
    MaskMsg = gfs.select(name = 'Land-sea mask',
                      typeOfLevel = 'surface')
    Mask    = MaskMsg[0].values
    
    # Variable values are a mask by default. Convert sea values to 0 to make
    #   the variable a non-masked array.
    JMask, IMask = np.where(Mask == 0)
    VarData[JMask, IMask] = 0
    
    # Close the .grb2 file
    gfs.close()
    
    #######################
    ### End of Function ###
    #######################
    return VarData, lat, lon, VarUnits, VarFH, VarVD, VarName, VarSName, Mask


#%%
#################################
### create_tmp_nc Function ######
#################################
    
def create_tmp_nc(VarData, lat, lon, VarUnits, VarFH, VarVD, VarName,
                  VarSName, Mask, TypeOfHeight):
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
    filename = 'GFS_' + VarSName + '_' +\
                   str(VarFH) + '.nc'
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

#%%
############################
### append_tmp Function ####
############################

def append_tmp(VSName):
    '''
    This function appends the variable short name to the temporary .txt file
      containing the model information.
      
    Input:
        VSName - The variable short name to be appended to the .txt file
    '''
    
    # Define the path
    path = './'
    #path = '/Users/Rarrell/Desktop/Thesis_Research/'
    
    # Open the temporary file
    f = open(path + 'tmp_var.txt', 'a')
    
    # Append the short name
    f.write(',' + VSName)
    
    # Close the file
    f.close()
    
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