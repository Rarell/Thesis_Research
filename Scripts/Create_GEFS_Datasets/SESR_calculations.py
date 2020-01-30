#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:11:24 2019

@author: Rarrell

This script is designed as part of the create_synthetic_dataset.sh program.
This script takes in the latent heat (LE) and potential evaporation rate (PET)
  netcdf files made earlier in the create_sythetic_dataset program and calculate
  the standardized evaporative stress ratio (SESR). SESR is calculated on a
  daily basis using both daily averaged LE and PET, and daily summed LE and PET.
  The daily SESR data is placed in a netcdf file for use by the user.
"""

#%%
#####################################
### Import some libraries ###########
#####################################
 
import os, sys, warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.mpl.ticker as cticker
from scipy import stats
from netCDF4 import Dataset, num2date
from datetime import datetime, timedelta
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator
 
#%%
  # Main function, for later. For now, just define some raw variables that
  #   would be loaded in the main function
  
# SNamePET = 'pevpr'
# SNameLE  = 'lhtfl'
  
# Year = '2019'
# Month = '09'
# Day = '01'

#%%
  # Load data
def main():
    '''
    '''
    
    # Collect the arguments
    script = sys.argv[0]
    SNameLE  = sys.argv[1]
    SNamePET = sys.argv[2]
    Parameters = sys.argv[3]
    
    Year, Month, Day, ModelRun = np.loadtxt(Parameters, usecols = (1, 2, 3, 4), 
                                     delimiter = ',', dtype = str, unpack = True)
    
    LE  = load_nc(SNameLE, Year, Month, Day, ModelRun)
    PET = load_nc(SNamePET, Year, Month, Day, ModelRun)

    ESR = CalculateESR(LE['lhtfl'], PET['pevpr'], LE['ymd'])

    ESRM = load_nc_climatology('mean')
    ESRSTD = load_nc_climatology('std')

    SESR = CalculateSESR(ESR, ESRM['esrm'], ESRSTD['esrstd'], LE['lat'], LE['lon'],
                         ESRM['lat'], ESRM['lon'], LE['ymd'], ESRM['day'])

    write_nc('sesr', 'SESR', SESR, LE['FH'], ESRM['lat1d'], ESRM['lon1d'],
             np.unique(LE['ymd']), Year, Month, Day, ModelRun)


#%%
##############################
 ### load_nc function ########
##############################
 
def load_nc(VarSName, Year, Month, Day, ModelRun):
    '''
    This function loads a netcdf file given the main variable's short name.
      Note this function is created for the datasets created earlier in the 
      create_synthetic_dataset.sh program, as seen in filename structure, the
      variable names, and the date structure for the valid dates.
      
    Inputs:
        VarSName - The short name of the variable for the netcdf
        Year  - The ending year of variable dataset
        Month - The ending month of the variable dataset
        Day   - The ending day of the variable dataset
        
    Outputs:
        X - Dictonary containing information from the netcdf file.
            This information includes:
                VarSname - Variable data (lat x lon x time)
                lat and lon - latitude and longitude (both gridded)
                FH    - The forecast hour each variable was extracted at
                mask  - The land-sea mask of the variable
                units - The units of the variable
                date  - The full date at each time step
                ymd   - The date (just year, month, and day) for each time step
                year  - The year for each time step
                month - The month for each time step
                day   - The day for each time step
                hour  - The hour/model run for each time step
    '''
    
    # Define the path and filename
    path = './Data/'
    filename = 'GFS_' + str(VarSName) + '_' + str(Year) + str(Month) +\
               str(Day) + '_' + str(ModelRun) + '.nc'
    
    # Initialize the directory for the data
    X = {}
    
    # Give the time format for the valid dates
    DateFormat = '%Y-%m-%d %H:%M:%S'
    
    with Dataset(path + filename, 'r') as nc:
        # Load the latitude and longitude
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        
        lat, lon = np.meshgrid(lat, lon)
        
        X['lat'] = lat.T
        X['lon'] = lon.T
        
        # Collect the time data
        VD = nc.variables['VD'][:]
        FH = nc.variables['FH'][:]
        X['FH'] = FH
        
        time = datetime.strptime(VD, DateFormat)
        
        dFH = np.ones((len(FH))) * np.nan
        dFH[0] = 0
        for i in range(1, len(dFH)):
            dFH[i] = FH[i] - FH[0]
        
        date = np.asarray([time + timedelta(hours = dt) for dt in dFH])
        
        X['date']  = date
        X['ymd']   = np.asarray([datetime(d.year, d.month, d.day) for d in date])
        X['year']  = np.asarray([d.year for d in date])
        X['month'] = np.asarray([d.month for d in date])
        X['day']   = np.asarray([d.day for d in date])
        X['hour']  = np.asarray([d.hour for d in date])
        
        # Collect the variable information and associated units.
        X[str(VarSName)] = nc.variables[str(VarSName)][:,:,:] # lat x lon x time
        X['units'] = nc.variables['units'][:]
        
        # Collect the mask data
        X['mask'] = nc.variables['mask'][:,:]
        
        ###################
        # End of Function #
        ###################
        
    return X
        

#%%
  # Write a function load the climatology file

def load_nc_climatology(ClimateType = 'mean'):
    '''
    '''

    # Define a path and filename
    path = './Data/SESR_Climatology/'
    if ClimateType == 'mean':
        filename = 'esr_climatology_mean.nc'
        VarName  = 'esrm'
    else:
        filename = 'esr_climatology_std.nc'
        VarName  = 'esrstd'
    
    # Initialize the directory for the data
    X = {}
    
    DateFormat = '%m-%d'
    
    with Dataset(path + filename, 'r') as nc:
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        
        X['lat1d'] = lat
        X['lon1d'] = lon
        
        lat, lon = np.meshgrid(lat, lon)
        
        X['lat'] = lat.T
        X['lon'] = lon.T
        
        # Collect the time data
        time = nc.variables['date'][:]
        date = np.asarray([datetime.strptime(time[d], DateFormat) for d in range(len(time))])
        
        X['day'] = date
        
        # Collect the data
        X[VarName] = nc.variables[VarName][:,:,:]
        
    return X
  
#%%

#############################
### write_nc Function #######
#############################

def write_nc(VarSName, VarName, Var, VarFH, lat, lon, dates, year, month, day, ModelRun):
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
        Years  - Range of years the data is valid for
        Months - Month the months the data is valid for
        Days   - Day the days the data is valid for
        TypeOfHeight - The type of height the variable is located at
    '''
    
    # Create the .nc file name based on the date, model run, and variable short
    #   name.
    # Note that ed stands for end date
    filename = 'GFS_' + str(VarSName) + '_' + str(year) +\
               str(month) + str(day) + '_' + str(ModelRun) + '.nc'
#    path = '/Users/Rarrell/Desktop/Thesis_Research/Data/'
    path = './Data/'
    
    # Collect the variable dimensions
    J, I, T = Var.shape
    
    # Begin writing the .nc file
    with Dataset(path + filename, 'w', format = 'NETCDF4') as nc:
        # Write the main description for the model variable data.
        nc.description = 'GFS forecast data for ' + VarName + ' valid from ' +\
                         str(dates[0]) + ' to ' + str(dates[-1]) + '. The variable ' +\
                         'is given for at least a 4 week period and is composed of 3 ' +\
                         'hour forecasts across multiple model runs.\n' +\
                         'Variable: ' + VarName +' (unltess). This variable ' +\
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
        
        # Create the forecast hour variable
        # nc.createVariable('FH', int, ('time', ))
        # nc.variables['FH'][:] = VarFH[:]
        
        # Create a variable for the date.
        nc.createVariable('date', np.str, ('time', ))
        for n in range(len(dates)):
            nc.variables['date'][n] = np.str(dates[n].strftime('%Y-%m-%d'))
        
        # Create the variable for the main variable
        nc.createVariable(str(VarSName), Var.dtype, ('lat', 'lon', 'time'))
        nc.setncatts({'long_name' : VarName, 'units' : 'unitless', 
                      'level_desc' : 'surface'})
        nc.variables[str(VarSName)][:,:,:] = Var[:,:,:]

        
        #######################
        ### End of Function ###
        #######################



#%%
  # ESR calculations
def CalculateESR(LE, PET, dates, lower_cutoff = 0, upper_cutoff = 1):
    '''
    '''
    J, I, T = LE.shape
    
    ESRraw = LE/PET
    
    # Remove extranious values
    ESRraw = ESRraw.reshape(J*I, T, order = 'F')
    for t in range(T):
        ind = np.where( (ESRraw[:,t] < lower_cutoff) | (ESRraw[:,t] > upper_cutoff) )[0]
        ESRraw[ind,t] = np.nan
        
    ESRraw = ESRraw.reshape(J, I, T, order = 'F')
    
    # Calculate the daily average ESR
    T = len(np.unique(dates)) # Count the total number of days
    
    ESR = np.ones((J, I, T)) * np.nan
    
    for n, date in enumerate(np.unique(dates)):
        ind = np.where(dates == date)[0]
        ESR[:,:,n] = np.nanmean(ESRraw[:,:,ind], axis = -1)
        
    return ESR


#%%
  # Calculate sesr and delta sesr

def CalculateSESR(ESR, ESRmean, ESRstd, lat, lon, latreff, lonreff, dates, datereff):
    '''
    '''

    # Focus on the domain the climatology data is
    
    # Convert the ESR longitude to same method as the climatology
    lon_change = lon
    ind = np.where(lon[0,:] > 180)[0]
    lon_change[:,ind] = lon_change[:,ind] - 360
    
    latind = np.where( (lat[:,0] >= np.nanmin(latreff)) & 
                      (lat[:,0] <= np.nanmax(latreff)) )[0]
    lonind = np.where( (lon_change[0,:] >= np.nanmin(lonreff)) &
                      (lon_change[0,:] <= np.nanmax(lonreff)) )[0]
    
    ESRsubset = np.ones((latind.size, lonind.size, ESR.shape[-1])) * np.nan
    
    for n, j in enumerate(reversed(latind)):
        for m, i in enumerate(lonind):
            ESRsubset[n,m,:] = ESR[j,i,:]
    
    J, I, T = ESRsubset.shape
    
    # Calculate SESR
    SESR = np.ones((J, I, T)) * np.nan
    for n, date in enumerate(np.unique(dates)):
        day = datetime(1900, np.unique(dates)[n].month, np.unique(dates)[n].day)
        ind = np.where( datereff == day )[0]
        for i in range(I):
            for j in range(J):
                SESR[j,i,n] = (ESRsubset[j,i,n] - ESRmean[i,j,ind[0]])/ESRstd[i,j,ind[0]]
        
    return SESR


#%%
#########################################
### Call and Run the Main Function ######
#########################################
    
if __name__ == '__main__':
    main()
    
#####################
### End of Script ###
#####################
