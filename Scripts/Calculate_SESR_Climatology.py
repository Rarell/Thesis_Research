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
from datetime import datetime
from matplotlib.dates import YearLocator, MonthLocator, DateFormatter
from matplotlib.ticker import MultipleLocator

#%%
  # Examine the files
  
path = '/Volumes/TestDrive/'

print(Dataset(path + 'esr_narr.nc', 'r'))
print(Dataset(path + 'lat_narr.nc', 'r'))
print(Dataset(path + 'lon_narr.nc', 'r'))
print(Dataset(path + 'sesr_narr.nc', 'r'))

#%%
  # Function to load the files
def load_2Dnc(filename, SName, path = '/Volumes/TestDrive/'):
    '''
    '''
    
    with Dataset(path + filename, 'r') as nc:
        var = nc.variables[SName][:,:]
        
    return var

def load_3Dnc(filename, SName, path = '/Volumes/TestDrive/'):
    '''
    '''
    
    with Dataset(path + filename, 'r') as nc:
        var = nc.variables[SName][:,:,:]
        
    return var

  
#%%
  # Load the files
  
esr = load_3Dnc('esr_narr.nc', 'esr') # Dataset is time x lat x lon
lat = load_2Dnc('lat_narr.nc', 'lat') # Dataset is lat x lon
lon = load_2Dnc('lon_narr.nc', 'lon') # Dataset is lat x lon

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


cmin = 0.0; cmax = 1.0; cint = 0.10
clevs = np.arange(cmin, cmax+cint, cint)
nlevs = len(clevs) - 1
cmap = plt.get_cmap(name = 'Reds', lut = nlevs)
#cmap  = plt.get_cmap(name = 'RdBu_r', lut = nlevs)

data_proj = ccrs.PlateCarree()
fig_proj  = ccrs.PlateCarree()

fig = plt.figure(figsize = [12, 16])
ax = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()

cs = ax.contourf(lon, lat, esr[-1,:,:], levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)

plt.show(block = False)
  
#%%
  # Load one of the example datasets with the desired grid
  
path = '/Users/Rarrell/Desktop/Thesis_Research/Data/Synthetic_Data/'
filename = 'GFS_pevpr_ed_20190901.nc' # Could use changing later

with Dataset(path + filename, 'r') as nc:
    lat_ref = nc.variables['lat'][:][::-1]
    lon_ref = nc.variables['lon'][:] - 180
    
    
print(lat_ref)
print(lon_ref)

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
  
path = '/Users/Rarrell/Desktop/Thesis_Research/Data/Synthetic_Data/'
filename = 'GFS_pevpr_ed_20190901.nc' # Could use changing later

with Dataset(path + filename, 'r') as nc:
    lat_ref = nc.variables['lat'][:][::-1]
    lon_ref = nc.variables['lon'][:] - 180
    
    
print(lat_ref)
print(lon_ref)

#%%
  # Interpolate the grid

def interpolate_to_GFS_grid(var, lat, lon):
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
    I = lon_ref.size
    J = lat_ref.size
    
    IntVar = np.ones((I, J, T)) * np.nan
    
    # Interpolate the data onto the GFS grid. Do this by averaging every
    #   data point located within a GFS grid.
    for i in range(lon_ref):
        for j in range(lat_ref):
            iind = np.where( (lon[0,:] >= lon_ref[i]) & (lon[0,:] <= lon_ref[i+1]) )[0]
            jind = np.where( (lat[:,0] >= lat_ref[j]) & (lat[:,0] <= lat_ref[j+1]) )[0]
            
            IntVar[i,j,:] = np.nanmean(var[:,jind,iind], axis = 0)
            
    return IntVar

#%%
  # Calculate the climatological means and standard deviations
  
def calculate_climatolgoy(var):
    '''
    '''
    
    # Obtain the dimensions of the variable
    I, J, T = var.shape
    
    # Count the number of years
    NumYear = np.ceil(T/365)
    
    # Create a variable for each day, assumed starting at Jan 1 and only no
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
        ClimMean[:,:,i-1] = np.nanmean(var[ind,:,:], axis = -1)
        ClimStd[:,:,i-1]  = np.nanstd(var[ind,:,:], axis = -1)
    
    return ClimMean, ClimStd

#%%
  # Calculate SESR
  
  
#%%
  # Calculate the change in SESR
  
#%% 
  # Calculate the climatology in the change in SESR
  
#%%
  # Create a function to write a variable to a .nc file
  
#%%
  # Write a .nc file for a ESR climatology, change in SESR, and change in SESR
  #   climatology
  