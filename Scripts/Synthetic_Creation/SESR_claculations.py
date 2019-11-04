#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 11:11:24 2019

@author: Rarrell
"""

#%%
 # Load some libraries
 
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
  # Main function, for later. For now, just define some raw variables that
  #   would be loaded in the main function
  
  SNamePET = 'pevpr'
  SNameLE  = 'lhtfl'
  
  Year = '2019'
  Month = '09'
  Day = '01'
  
#%%
  # Load nc function
def load_nc(VarSName, Year, Month, Day):
    '''
    '''
    
    # Define the path and filename
    path = './Data/Synthetic_Data/'
    filename = 'GFS_' + str(VarSName) + '_ed_' + str(Year) + str(Month) +\
               str(Day) + '.nc'
    
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
        time = nc.variables['VD'][:]
        date = np.asarray([datetime.strptime(time[d], DateFormat) for d in range(len(time))])
        
        X['date']  = date
        X['ymd']   = np.asarray([datetime(d.year, d.month, d.day) for d in date])
        X['year']  = np.asarray([d.year for d in date])
        X['month'] = np.asarray([d.month for d in date])
        X['day']   = np.asarray([d.day for d in date])
        X['hour']  = np.asarray([d.hour for d in date])
        
        # Collect for the forecast hour the variable is valid for
        X['FH'] = nc.variables['FH'][:]
        
        # Collect the variable information and associated units.
        X[str(VarSName)] = nc.variables[str(VarSName)][:,:,:] # lat x lon x time
        X['units'] = nc.variables['units'][:]
        
        # Collect the mask data
        X['mask'] = nc.variables['mask'][:,:]
        
    return X
        
        
  
#%%
  # Load data
  
LE = load_nc(SNameLE, Year, Month, Day)
PET = load_nc(SNamePET, Year, Month, Day)

print(LE['FH'], PET['FH'])
print(LE['ymd'], PET['ymd'])

#%%
  # Some LE calculations
J, I, T = LE['lhtfl'].shape
T_days = len(np.unique(LE['ymd'])) # Count the total number of days
LE['sum'] = np.ones((J, I, T_days)) * np.nan
LE['avg'] = np.ones((J, I, T_days)) * np.nan
  
for n, date in enumerate(np.unique(LE['ymd'])):
    ind = np.where(LE['ymd'] == date)[0]
    LE['sum'][:,:,n] = np.nansum(LE['lhtfl'][:,:,ind], axis = -1)
    LE['avg'][:,:,n] = np.nanmean(LE['lhtfl'][:,:,ind], axis = -1)


print(LE['lhtfl'])
print(LE['sum'])
print(LE['avg'])

#%%
  # Some PET calculations
J, I, T = PET['pevpr'].shape
T_days = len(np.unique(PET['ymd'])) # Count the total number of days
PET['sum'] = np.ones((J, I, T_days)) * np.nan
PET['avg'] = np.ones((J, I, T_days)) * np.nan
  
for n, date in enumerate(np.unique(PET['ymd'])):
    ind = np.where(PET['ymd'] == date)[0]
    PET['sum'][:,:,n] = np.nansum(PET['pevpr'][:,:,ind], axis = -1)
    PET['avg'][:,:,n] = np.nanmean(PET['pevpr'][:,:,ind], axis = -1)


print(PET['pevpr'])
print(PET['sum'])
print(PET['avg'])

#%%
  # ESR calculations

ESR = LE['lhtfl']/PET['pevpr']
ESR_sum = LE['sum']/PET['sum']
ESR_avg = LE['avg']/PET['avg']

#%%
  # Plot ESR

cmin = 0.0; cmax = 1.0; cint = 0.1
clevs = np.arange(cmin, cmax+cint, cint)
nlevs = len(clevs) - 1
cmap  = plt.get_cmap(name = 'Reds', lut = nlevs)

data_proj = ccrs.PlateCarree()
fig_proj  = ccrs.PlateCarree()

fig = plt.figure(figsize = [12, 16])
ax = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()

cs = ax.contourf(LE['lon'], LE['lat'], ESR_sum[:,:,-1], levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

plt.show(block = False)

  