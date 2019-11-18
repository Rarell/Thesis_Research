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
##############################
 ### load_nc function ########
##############################
 
def load_nc(VarSName, Year, Month, Day):
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
        
        ###################
        # End of Function #
        ###################
        
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

# Calculate the daily averages and daily sums
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

# Calculate the daily averages and daily sums
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

ESRsum = np.ones((J, I, T_days)) * np.nan
ESRavg = np.ones((J, I, T_days)) * np.nan

# Calculate the daily averages and daily sums
for n, date in enumerate(np.unique(PET['ymd'])):
    ind = np.where(PET['ymd'] == date)[0]
    ESRsum[:,:,n] = np.nansum(ESR[:,:,ind], axis = -1)
    ESRavg[:,:,n] = np.nanmean(ESR[:,:,ind], axis = -1)
    
print(ESRsum == ESR_sum)
print(ESRavg == ESR_avg)

#%%
  # Plot ESR/make and example plot

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

cs = ax.contourf(LE['lon'], LE['lat'], ESRavg[:,:,-1], levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

plt.show(block = False)

#%%
  # Plot some histograms

ESR_avec = ESR_avg.reshape(I*J*T_days, order = 'F')
ESR_svec = ESR_sum.reshape(I*J*T_days, order = 'F')

ind = np.where( ((ESR_avec[:] > 10) | (ESR_avec[:] < -10)) )[0]
ESR_avec[ind] = np.nan

ind = np.where( ((ESR_svec[:] > 10) | (ESR_svec[:] < -10)) )[0]
ESR_svec[ind] = np.nan
  
fig, axes = plt.subplots(figsize = [12, 14], ncols = 1, nrows = 2, sharex = True)
ax1 = axes[0]; ax2 = axes[1]

plt.subplots_adjust(left = 0.1, right = 0.9, bottom = 0.1, top = 0.9,
                    wspace = 0.10, hspace = 0.10)


ax1.set_title('Daily average ESR histogram')
navg = ax1.hist(ESR_avec[:], bins = 100, density = False, color = 'r', alpha = 0.7)


ax2.set_title('Daily sum ESR histogram')
nsum = ax2.hist(ESR_svec[:], bins = 100, density = False, color = 'b', alpha = 0.7)



plt.show(block = False)



#%%
neg_avg = np.where( navg[1] < 0 )[0]
neg_sum = np.where( nsum[1] < 0 )[0]


print(sum(navg[0][neg_avg]/sum(navg[0])))
print(sum(nsum[0][neg_sum]/sum(nsum[0])))

#%%
  # Try some fresh plotting
  
ESR_avg = ESR_avec.reshape(J, I, T_days, order = 'F')
ESR_sum = ESR_svec.reshape(J, I, T_days, order = 'F')

cmin = 0.0; cmax = 1.0; cint = 0.10
clevs = np.arange(cmin, cmax+cint, cint)
nlevs = len(clevs) - 1
cmap = plt.get_cmap(name = 'Reds', lut = nlevs)
#cmap  = plt.get_cmap(name = 'RdBu_r', lut = nlevs)

data_proj = ccrs.PlateCarree()
fig_proj  = ccrs.PlateCarree()

fig = plt.figure(figsize = [12, 16])
ax = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.set_title('One day sum of ESR calculated from GFS data.\n Valid for ' +\
             str(LE['date'][-1]), size = 18)

ax.coastlines()

cs = ax.contourf(LE['lon'], LE['lat'], ESR_sum[:,:,-1], levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.355, 0.02, 0.29])
cbar = fig.colorbar(cs, cax = cbax)

plt.show(block = False)