#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 15:23:09 2019

@author: Rarrell

This is a script designed to test individual components of the grb_to_nc.py
script.
"""

#%% Create a testing point to test parts of the above function.

#%% Import needed modules
import sys, os, warnings
import wget
import numpy as np
import pygrib
from netCDF4 import Dataset
from glob import glob
from datetime import datetime, timedelta

#%%
# Create some example input
var = 'Potential evaporation rate'
typeOfHeight = 'surface'
height = 0
year = '2019'
month = '06'
day = '10'
model_run = '0000'
forecast_hour = '003'

#%%
# Create an example download

# Define the url. This is the https url from the NCDC. Everything before
# gfs-avn-hi/ is fixed. After it is /YYYYMM/YYYYMMDD/file_name
# The file name is in the format of:
#   gfs_3_YYYYMMDD_MRMR_FRR.grb2
#   gfs_3 is the GFS 003 model (1˚ domain)
#   MRMR is the model run (0000, 0600, 1200, and 1800)
#   FRR is the forecast hour (003, 006, 009, ..., 384)
#filename = 'gfs_3_20190610_0000_003.grb2' 
#url = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/201906/20190610/' + filename
path = '/Users/Rarrell/Downloads/tmp/'
filename = 'gfs_3_' + year + month + day + '_' +\
            model_run  + '_' + forecast_hour + '.grb2'

url_start = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
url_end   = year + month + '/' + year + month +\
            day + '/' + filename
url = url_start + url_end

print(url)

# Note when using wget.download. The first arguement is the url. The second is
# where to download to. This second step automatically appends a ./ to the path.
# To bypass this, use /Users/Username/... to tell it to start from the root
# directory.
# Download the GFS data into Downloads/tmp/ file with the filename of filename
wget.download(url, path + filename)    

#%%
# Note there are shortNames. Might look into this.
# Create an example retrevial of the grb data.

print(path + filename)
gfs = pygrib.open(path + filename)

#for g in gfs:
#    print(g.name, g.typeOfLevel)

# Collect the variable
var_msg = gfs.select(name = 'Potential evaporation rate',
                     typeOfLevel = 'surface')
var_data = var_msg[0].values
var_units = var_msg[0].units
var_forecastHour = var_msg[0].forecastTime
var_valDate = var_msg[0].validDate
var_name = var_msg[0].name
var_shortName = var_msg[0].shortName

# Collect the mask
mask_msg = gfs.select(name = 'Land-sea mask',
                      typeOfLevel = 'surface')
mask_data = mask_msg[0].values
# Collect gridded latitude and longitude data
lat, lon = var_msg[0].latlons()

gfs.close()

# Remove mask values. var_data will then automatically be turned into a float32
# datatype instead of a masked array.
j_mask, i_mask = np.where(mask_data == 0)
#var_data[j_mask, i_mask] = 0

#%%
# Write the extracted data into a netcdf file.
J, I = var_data.shape

with Dataset(path + 'tmp_fil.nc', 'w', format = 'NETCDF4') as nc:
    # Add a short description
    nc.description = "This is a test nc file. The variable will be called " + var_shortName
    
    # Create the dimensions for latitude and longitude
    nc.createDimension('lat', size = J)
    nc.createDimension('lon', size = I)
    
    # Create the latitude and longitude data
    nc.createVariable('lat', lat.dtype, ('lat', ))
    nc.createVariable('lon', lon.dtype, ('lon', ))
    nc.variables['lat'][:] = lat[:,0]
    nc.variables['lon'][:] = lon[0,:]
    
    # Add a variable for the valid date. dimensions is not specified so this is
    # a scalar. Note this is a datetime, and does not have a .dtype attribute.
    nc.createVariable('ValidDate', np.str)
    nc.variables['ValidDate'][0] = np.str(var_valDate)
    
    # Add a variable for the forecast hour. For the the tmp files, this is
    # also a scalar.
    nc.createVariable('ForecastHour', int)
    nc.variables['ForecastHour'][:] = var_forecastHour
    
    # Create the main variable
    nc.createVariable(var_shortName, var_data.dtype, ('lat', 'lon'))
    nc.setncatts({'long_name' : var_name, 'units' : var_units,\
                  'level_desc' : typeOfHeight})
    nc.variables[var_shortName][:,:] = var_data[:,:] # Data is in lat x lon
    
    # Place the units in their own variable
    nc.createVariable('units', np.str)
    nc.variables['units'][0] = var_units
    
    # Create a value for the land mask.
    nc.createVariable('mask', mask_data.dtype, ('lat', 'lon'))
    nc.variables['mask'][:,:] = mask_data[:,:]
    

#%%
# Open the created .nc file and compare with data to ensure everything was
# written correctly.
print(var_forecastHour.dtype)
print(Dataset('/Users/Rarrell/Downloads/tmp/tmp_fil.nc'))
#%%
with Dataset('/Users/Rarrell/Desktop/Thesis_Research/tmp/GFS_pevpr_3.nc', 'r') as nc:
    read_lat = nc.variables['lat'][:]
    read_lon = nc.variables['lon'][:]
    
    read_lat, read_lon = np.meshgrid(read_lat, read_lon)
    
    read_FH = nc.variables['FH'][:]
    read_VD = nc.variables['VD'][:]
    
    read_var = nc.variables[var_shortName][:,:]
    read_units = nc.variables['units'][:]
    
    read_mask = nc.variables['mask'][:,:]

# Note the date format is %Y-%m-%d %H:%M:%S
#print(lat)
print(read_lat.T == lat) # Verified as true. Note that the lon and lat grid is transposed.
#print(lon)
print(read_lon.T == lon) # Verified as true.
#print(var_valDate)
#print(datetime.strptime(read_VD, 'YYYY-MM-DD hh:mm:ss'))
print(var_valDate == datetime.strptime(read_VD, '%Y-%m-%d %H:%M:%S')) # Verified as true.
#print(var_forecastHour)
print(var_forecastHour == read_FH) # Verified as true.
#print(mask_data)
print(mask_data == read_mask) # Verified as true.
#print(var_data)
print(var_data == read_var) # Verified as true.
#print(var_units)
print(var_units == read_units) # Verified as true.

#%% Test the merge_nc.py output

with Dataset('Data/GFS_pevpr_20190610.nc', 'r') as nc:
    read_lat = nc.variables['lat'][:]
    read_lon = nc.variables['lon'][:]
    
    read_lat, read_lon = np.meshgrid(read_lat, read_lon)
    
    read_FH = nc.variables['FH'][:]
    read_VD = nc.variables['VD'][:]
    
    read_var = nc.variables['pevpr'][:,:,:]
    read_units = nc.variables['units'][:]
    read_mask = nc.variables['mask'][:,:]
    
print(read_lat.T == lat)
print(read_lon.T == lon)
print(read_VD)
print(read_VD == var_valDate)
print(read_FH[0] == var_forecastHour)
print(read_var[:,:,0] == var_data)
print(read_mask == mask_data)

#%%
  # Test for the SESR calculations script

path = './Data/'
filename = 'GFS_lhtfl_20190901_00.nc'

print(Dataset(path + filename, 'r'))

X = {}

with Dataset(path + filename, 'r') as nc:
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    
    FH = nc.variables['FH'][:]
    VD = nc.variables['VD'][:]
    
    LE = nc.variables['lhtfl'][:,:,:]
    
    
print(VD)
print(FH)

#%%
  # More SESR calculations testing

DateFormat = '%Y-%m-%d %H:%M:%S'

VD = datetime.strptime(VD, DateFormat)

initChange = timedelta(hours = 3)

trueVD = VD - initChange

print(trueVD)

#%%
  # Testing trying to turn FH into timedeltas
dFH = np.ones((len(FH))) * np.nan
dFH[0] = 0
for i in range(1, len(dFH)):
    dFH[i] = FH[i] - FH[0]
#dFH[1:] = FH2[1:] - FH2[:-1]


print(dFH)

time = VD

for dt in dFH:
    time = VD + timedelta(hours = dt)
    print(time)

time = VD
date = np.asarray([time + timedelta(hours = dt) for dt in dFH])
print(date)









