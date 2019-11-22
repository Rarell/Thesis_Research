#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 11:11:52 2019

@author: Rarrell

This script is designed to take the daily ESR NARR file and interpolate the 
  grid to a 1 degree by 1 degree or 0.5 degrees by 0.5 degrees grid, and then
  calculate the ESR and SESR climatologies, including the mean, standard
  deviation, and quantiles for ESR and SESR.
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
  # Examine the files

########
#### Hard coded path. Should change. ######
#########
path = '/Volumes/TestDrive/'

print(Dataset(path + 'esr_narr.nc', 'r'))
print(Dataset(path + 'lat_narr.nc', 'r'))
print(Dataset(path + 'lon_narr.nc', 'r'))
print(Dataset(path + 'sesr_narr.nc', 'r'))

#%%
  # Function to load the files
def load2Dnc(filename, SName, path = '/Volumes/TestDrive/'):
    '''
    '''
    
    with Dataset(path + filename, 'r') as nc:
        var = nc.variables[SName][:,:]
        
    return var

def load3Dnc(filename, SName, path = '/Volumes/TestDrive/'):
    '''
    '''
    
    with Dataset(path + filename, 'r') as nc:
        var = nc.variables[SName][:,:,:]
        
    return var

  
#%%
  # Interpolate the grid

def InterpolateToGFSGrid(var, lat, lon):
    '''
    
    Assumes var is in a time x lat x lon setup, and lat and lon are in a
    lat x lon grid.
    '''
    
    # Load a GFS grid for reference use
    path = './Data/Synthetic_Data/'
    filename = 'GFS_pevpr_ed_20190901.nc' # Could use changing later

    with Dataset(path + filename, 'r') as nc:
        lat_ref = nc.variables['lat'][:][::-1]
        lon_ref = nc.variables['lon'][:] - 180
    
    # Focus on the part of the grid that contains var
    lon_ind = np.where( (lon_ref[:] >= np.nanmin(lon)) & 
                        (lon_ref[:] <= np.nanmax(lon)) )[0]
    lat_ind = np.where( (lat_ref[:] >= np.nanmin(lat)) & 
                        (lat_ref[:] <= np.nanmax(lat)) )[0]
    
    lon_ref = lon_ref[lon_ind]
    lat_ref = lat_ref[lat_ind]
    
    # Initialize the interpolated variable
    T, J, I = var.shape
    
    # Reorder some variables for easier searching and referencing
    lon1d = lon.reshape(J*I, order = 'F')
    lat1d = lat.reshape(J*I, order = 'F')
    
    var2d = var.reshape(T, J*I, order = 'F')
    
    I = lon_ref.size
    J = lat_ref.size
    
    IntVar = np.ones((I, J, T)) * np.nan
    
    # Interpolate the data onto the GFS grid. Do this by averaging every
    #   data point located within a GFS grid.
    for i in range(lon_ref):
        for j in range(lat_ref):
            # Determine all points between the desired latitude and longitude
            #   seperately
            
            # Note fewer points are used in the final latitude/longitude since
            #   they go from that latitude to the grid end (less than 1 degree).
            #   Since these points should be in the Atlantic/Arctic sea, this
            #   should be unimportant though (loops go from south to north and
            #   west to east).
            if i == range(lon_ref)[-1]:
                iind = np.where( lon1d >= lon_ref[i] )[0]
            else:
                iind = np.where( (lon1d >= lon_ref[i]) & (lon1d <= lon_ref[i+1]) )[0]
            
            if j == range(lat_ref)[-1]:
                jind = np.where( lat1d >= lat_ref[j] )[0]
            else:
                jind = np.where( (lat1d >= lat_ref[j]) & (lat1d <= lat_ref[j+1]) )[0]
            
            # Determine the grid points between the desired latitude and
            #   longitude together
            ind = np.intersect1d(iind, jind)
            
            IntVar[i,j,:] = np.nanmean(var2d[:,ind], axis = -1)
            
    return IntVar, lat_ref, lon_ref

#%%
  # Calculate the climatological means and standard deviations
  
def CalculateClimatology(var):
    '''
    '''
    
    # Obtain the dimensions of the variable
    I, J, T = var.shape
    
    # Count the number of years
    NumYear = np.ceil(T/365)
    
    # Create a variable for each day, assumed starting at Jan 1 and no
    #   leap years (i.e., each year is only 365 days each)
    day = np.ones((T)) * np.nan
    
    n = 0
    for i in range(1, NumYear):
        if i >= NumYear:
            day[n:-1] = np.arange(1, len(day[n:-1]))
        else:
            day[n:n+365] = np.arange(1, 365+1)
        
        n = n + 365
    
    # Initialize the climatological mean and standard deviation variables
    ClimMean = np.ones((I, J, 365)) * np.nan
    ClimStd  = np.ones((I, J, 365)) * np.nan
    
    # Calculate the mean and standard deviation for each day and at each grid
    #   point
    for i in range(1, 365+1):
        ind = np.where(i == day)[0]
        ClimMean[:,:,i-1] = np.nanmean(var[:,:,ind], axis = -1)
        ClimStd[:,:,i-1]  = np.nanstd(var[:,:,ind], axis = -1)
    
    return ClimMean, ClimStd

#%%
  # Create functions to create datetime datasets

def DateRange(StartDate, EndDate):
    '''
    This function takes in two dates and outputs all the dates inbetween
    those two dates.
    
    Inputs:
        StartDate - A datetime. The starting date of the interval.
        EndDate - A datetime. The ending date of the interval.
        
    Outputs:
        All dates between StartDate and EndDate (inclusive)
    '''
    for n in range(int((EndDate - StartDate).days) + 1):
        yield StartDate + timedelta(n) 

#%%
  # Create a function to write a variable to a .nc file
  
def WriteNC(var, lat, lon, dates, filename = 'tmp.nc', VarName = 'tmp'):
    '''
    '''
    
    # Define the path
    path = './Data/'
    
    # Determine the spatial and temporal lengths
    I = lon.size
    J = lat.size
    T = len(dates)
    
    with Dataset(path + filename, 'w', format = 'NETCDF4') as nc:
        # Write a description for the .nc file
        nc.description = 'This is one of six .nc files containing synthetic ' +\
                         'data needed to calculate flash drought using SESR. ' +\
                         'These files contain the daily mean and standard ' +\
                         'deviation of ESR, daily SESR, daily change in SESR ' +\
                         'and daily meand and standard deviation of the change ' +\
                         'in SESR. Note all these variables are unitless. Also '+\
                         'note that all these variables are on a 1 degree '+\
                         'x 1 degree GFS grid. The variable this file contains ' +\
                         'is ' + str(VarName) + '.\n' +\
                         'Variable: ' + str(VarName) + '(unitless). This is the ' +\
                         'main variable for this file. It is in the format ' +\
                         'lon x lat x time.\n' +\
                         'lat: The latitude (vector form).\n' +\
                         'lon: The longitude (vector form).\n' +\
                         'date: List of dates starting from 01-01 to 12-31 ' +\
                         'for climatology variables (%m-%d format), or from ' +\
                         '01-01-1979 to 12-31-2016 for SESR/change in SESR ' +\
                         'datasets (%Y-%m-%d format). For change in SESR, date ' +\
                         'is the start date of the change.'
        
        # Create the spatial and temporal dimensions
        nc.createDimension('lat', size = J)
        nc.createDimension('lon', size = I)
        nc.createDimension('time', size = T)
        
        # Create the lat and lon variables
        nc.createVariable('lat', lat.dtype, ('lat', ))
        nc.createVariable('lon', lon.dtype, ('lon', ))
        
        nc.variables['lat'][:] = lat[:]
        nc.variables['lon'][:] = lon[:]
        
        # Create the date variable
        nc.createVariable('date', str, ('time', ))
        for n in range(len(dates)):
            nc.variables['date'][n] = np.str(dates[n])
            
        # Create the main variable
        nc.createVariable(VarName, var.dtype, ('lon', 'lat', 'time'))
        nc.variables[str(VarName)][:,:,:] = var[:,:,:]
        

#%%
  # Load the files
  
esr = load3Dnc('esr_narr.nc', 'esr') # Dataset is time x lat x lon
lat = load2Dnc('lat_narr.nc', 'lat') # Dataset is lat x lon
lon = load2Dnc('lon_narr.nc', 'lon') # Dataset is lat x lon

print(esr)
print(lat)
print(lon)

#%%
  # Turn positive lon values to negative and see if this improves the map.
  
for i in range(len(lon[:,0])):
    ind = np.where( lon[i,:] > 0 )[0]
    lon[i,ind] = -1*lon[i,ind]

#%%
  # Create a sample plot of ESR to see what the map/grid looks like

# Lonitude and latitude tick information
lat_int = 15
lon_int = 10

lat_label = np.arange(-90, 90, lat_int)
lon_label = np.arange(-180, 180, lon_int)

#lon_formatter = cticker.LongitudeFormatter()
#lat_formatter = cticker.LatitudeFormatter()

# Colorbar information
cmin = 0.0; cmax = 1.0; cint = 0.10
clevs = np.arange(cmin, cmax+cint, cint)
nlevs = len(clevs) - 1
cmap = plt.get_cmap(name = 'Reds', lut = nlevs)
#cmap  = plt.get_cmap(name = 'RdBu_r', lut = nlevs)

data_proj = ccrs.PlateCarree()
fig_proj  = ccrs.PlateCarree()

# Figure
fig = plt.figure(figsize = [12, 16])
ax = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()

ax.set_xticks(lon_label, crs = ccrs.PlateCarree())
ax.set_yticks(lat_label, crs = ccrs.PlateCarree())
ax.set_xticklabels(lon_label, fontsize = 16)
ax.set_yticklabels(lat_label, fontsize = 16)

ax.xaxis.tick_bottom()
ax.yaxis.tick_left()

#ax.xaxis.set_major_formatter(lon_formatter)
#ax.yaxis.set_major_formatter(lat_formatter)

cs = ax.contourf(lon, lat, esr[-1,:,:], levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)

plt.show(block = False)
  
#%%
  # Load one of the example datasets with the desired grid

########
#### Hard coded path. Should change. ######
#########
path = '/Users/Rarrell/Desktop/Thesis_Research/Data/Synthetic_Data/'
filename = 'GFS_pevpr_ed_20190901.nc' # Could use changing later

with Dataset(path + filename, 'r') as nc:
    lat_ref = nc.variables['lat'][:][::-1]
    lon_ref = nc.variables['lon'][:] - 180
    
    
print(lat_ref)
print(lon_ref)

#%%
  # Interpolate the NARR grid unto the GFS grid and test/examine it on a map
  
esr_GFS, lat_GFS, lon_GFS = InterpolateToGFSGrid(esr, lat, lon)

lon_grid, lat_grid = np.meshgrid(lon_GFS, lat_GFS)

# Reuse the same color map, labels, and projections as the first map on lines
#   194 - 212

data_proj = ccrs.PlateCarree()
fig_proj  = ccrs.PlateCarree()

# Figure
fig = plt.figure(figsize = [12, 16])
ax = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()

ax.set_xticks(lon_label, crs = ccrs.PlateCarree())
ax.set_yticks(lat_label, crs = ccrs.PlateCarree())
ax.set_xticklabels(lon_label, fontsize = 16)
ax.set_yticklabels(lat_label, fontsize = 16)

ax.xaxis.tick_bottom()
ax.yaxis.tick_left()

#ax.xaxis.set_major_formatter(lon_formatter)
#ax.yaxis.set_major_formatter(lat_formatter)

cs = ax.contourf(lon_grid, lat_grid, esr_GFS[:,:,-1], levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)

plt.show(block = False)  

#%%
  # Calculate the climatologies and test/examine by plotting a map for June 30
  
esr_mean, esr_std = CalculateClimatology(esr_GFS)

# Reuse the same color map, labels, and projections as the first map on lines
#   194 - 212

data_proj = ccrs.PlateCarree()
fig_proj  = ccrs.PlateCarree()

# Figure
fig = plt.figure(figsize = [12, 16])
ax = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()

ax.set_xticks(lon_label, crs = ccrs.PlateCarree())
ax.set_yticks(lat_label, crs = ccrs.PlateCarree())
ax.set_xticklabels(lon_label, fontsize = 16)
ax.set_yticklabels(lat_label, fontsize = 16)

ax.xaxis.tick_bottom()
ax.yaxis.tick_left()

#ax.xaxis.set_major_formatter(lon_formatter)
#ax.yaxis.set_major_formatter(lat_formatter)

cs = ax.contourf(lon_grid, lat_grid, esr_mean[:,:,181], levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)

plt.show(block = False)  
#%%
  # Calculate SESR
  
sesr = (esr_GFS - esr_mean)/esr_std
  
#%%
  # Calculate the change in SESR
  
delta_sesr = sesr[1:] - sesr[:-1]
  
#%% 
  # Calculate the climatology in the change in SESR
  
delta_sesr_mean, delta_sesr_std = CalculateClimatology(delta_sesr)

#%%
  # Write two sets of datetime variables. One for sesr and delta_sesr (every)
  #   date from the start point to the end point.
  # The other is for the climatologies and goes from Jan 01 to Dec 31
  
narr_start = datetime(1979, 1, 1)
narr_end   = datetime(2016, 12, 31)

narr_dates = DateRange(narr_start, narr_end)
narr_dates = narr_dates.strftime('%Y-%m-%d')


year_start = datetime(1980, 1, 1)
year_end   = datetime(1980, 12, 31)

year_dates = DateRange(year_start, year_end)
year_dates = year_dates.strftime('%m-%d')

#%%
  # Write a .nc file for a ESR climatology, change in SESR, and change in SESR
  #   climatology
  
# Write the ESR climatology files
WriteNC(esr_mean, lat_GFS, lon_GFS, year_dates, 'esr_climatology_mean.nc',
         'esrm')

WriteNC(esr_std, lat_GFS, lon_GFS, year_dates, 'esr_climatology_std.nc',
         'esrstd')

# Write the SESR file
WriteNC(sesr, lat_GFS, lon_GFS, narr_dates, 'GFS_grid_sesr.nc', 'sesr')

# Write the change in SESR file
WriteNC(delta_sesr, lat_GFS, lon_GFS, narr_dates[:-1], 'GFS_grid_delta_sesr.nc',
         'dsesr')

# Write the change in SESR climatology files
WriteNC(delta_sesr_mean, lat_GFS, lon_GFS, year_dates, 'delta_sesr_climatology_mean.nc',
         'dsesrm')

WriteNC(delta_sesr_std, lat_GFS, lon_GFS, year_dates, 'delta_sesr_climatology_std.nc',
         'dsesrsstd')
 