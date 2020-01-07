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
# cell 1
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
# cell 2
  # Examine the files

########
#### Hard coded path. Should change. ######
#########
path = '/Volumes/My Book/'

print(Dataset(path + 'esr_narr.nc', 'r'))
print(Dataset(path + 'lat_narr.nc', 'r'))
print(Dataset(path + 'lon_narr.nc', 'r'))
print(Dataset(path + 'sesr_narr.nc', 'r'))

#%% 
# cell 3
  # Function to load the files
def load2Dnc(filename, SName, path = '/Volumes/My Book/'):
    '''
    '''
    
    with Dataset(path + filename, 'r') as nc:
        var = nc.variables[SName][:,:]
        
    return var

def load3Dnc(filename, SName, path = '/Volumes/My Book/'):
    '''
    '''
    
    with Dataset(path + filename, 'r') as nc:
        var = nc.variables[SName][:,:,:]
        
    return var

  
#%% 
# cell 4
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
    for i in range(I):
        for j in range(J):
            # Determine all points between the desired latitude and longitude
            #   seperately
            
            # Note fewer points are used in the final latitude/longitude since
            #   they go from that latitude to the grid end (less than 1 degree).
            #   Since these points should be in the Atlantic/Arctic sea, this
            #   should be unimportant though (loops go from south to north and
            #   west to east).
            if i == range(I)[-1]:
                iind = np.where( lon1d >= lon_ref[i] )[0]
            else:
                iind = np.where( (lon1d >= lon_ref[i]) & (lon1d <= lon_ref[i+1]) )[0]
            
            if j == range(J)[-1]:
                jind = np.where( lat1d >= lat_ref[j] )[0]
            else:
                jind = np.where( (lat1d >= lat_ref[j]) & (lat1d <= lat_ref[j+1]) )[0]
            
            # Determine the grid points between the desired latitude and
            #   longitude together
            ind = np.intersect1d(iind, jind)
            
            IntVar[i,j,:] = np.nanmean(var2d[:,ind], axis = -1)
            
    return IntVar, lat_ref, lon_ref

#%% 
# cell 5
  # Calculate the climatological means and standard deviations
  
def CalculateClimatology(var, pentad = False):
    '''
    '''
    
    # Obtain the dimensions of the variable
    I, J, T = var.shape
    
    # Count the number of years
    if pentad is True:
        yearLen = int(365/5)
    else:
        yearLen = int(365)
        
    NumYear = int(np.ceil(T/yearLen))
    
    # Create a variable for each day, assumed starting at Jan 1 and no
    #   leap years (i.e., each year is only 365 days each)
    day = np.ones((T)) * np.nan
    
    n = 0
    for i in range(1, NumYear+1):
        if i >= NumYear:
            day[n:T+1] = np.arange(1, len(day[n:T+1])+1)
        else:
            day[n:n+yearLen] = np.arange(1, yearLen+1)
        
        n = n + yearLen
    
    # Initialize the climatological mean and standard deviation variables
    ClimMean = np.ones((I, J, yearLen)) * np.nan
    ClimStd  = np.ones((I, J, yearLen)) * np.nan
    
    # Calculate the mean and standard deviation for each day and at each grid
    #   point
    for i in range(1, yearLen+1):
        ind = np.where(i == day)[0]
        ClimMean[:,:,i-1] = np.nanmean(var[:,:,ind], axis = -1)
        ClimStd[:,:,i-1]  = np.nanstd(var[:,:,ind], axis = -1)
    
    return ClimMean, ClimStd

#%% 
# cell 6
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
# cell 7
  # Create a function to write a variable to a .nc file
  
def WriteNC(var, lat, lon, dates, filename = 'tmp.nc', VarName = 'tmp'):
    '''
    '''
    
    # Define the path
    path = './Data/SESR_Climatology/'
    
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
                         'and daily mean and standard deviation of the change ' +\
                         'in SESR. Note all these variables are unitless. Also '+\
                         'note that all these variables are on a 1 degree '+\
                         'x 1 degree GFS grid. The variable this file contains ' +\
                         'is ' + str(VarName) + '.\n' +\
                         'Variable: ' + str(VarName) + ' (unitless). This is the ' +\
                         'main variable for this file. It is in the format ' +\
                         'lon x lat x time.\n' +\
                         'lat: The latitude (vector form).\n' +\
                         'lon: The longitude (vector form).\n' +\
                         'date: List of dates starting from 01-01 to 12-31 ' +\
                         'for climatology variables (%m-%d format), or from ' +\
                         '01-01-1979 to 12-31-2016 for SESR/change in SESR ' +\
                         'datasets (%Y-%m-%d format). Leap year additions are excluded.' +\
                         'For change in SESR, date ' +\
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
# cell 8
  # Load the files
  
esr = load3Dnc('esr_narr.nc', 'esr') # Dataset is time x lat x lon
lat = load2Dnc('lat_narr.nc', 'lat') # Dataset is lat x lon
lon = load2Dnc('lon_narr.nc', 'lon') # Dataset is lat x lon

print(esr)
print(lat)
print(lon)

#%%
# cell 9
  # Turn positive lon values to negative and see if this improves the map.
  
for i in range(len(lon[:,0])):
    ind = np.where( lon[i,:] > 0 )[0]
    lon[i,ind] = -1*lon[i,ind]

#%%
# cell 10
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
# cell 11
  # Load one of the example datasets with the desired grid

########
#### Hard coded path. Should change. ######
#########
path = './Data/Synthetic_Data/'
filename = 'GFS_pevpr_ed_20190901.nc' # Could use changing later

with Dataset(path + filename, 'r') as nc:
    lat_ref = nc.variables['lat'][:][::-1]
    lon_ref = nc.variables['lon'][:] - 180
    
    
print(lat_ref)
print(lon_ref)

#%%
# cell 12
  # Interpolate the NARR grid unto the GFS grid and test/examine it on a map
  
esr_GFS, lat_GFS, lon_GFS = InterpolateToGFSGrid(esr, lat, lon)

lon_grid, lat_grid = np.meshgrid(lon_GFS, lat_GFS)

# Reuse the same color map, labels, and projections as the first map on lines
#   274 - 292

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

cs = ax.contourf(lon_grid, lat_grid, esr_GFS[:,:,-1].T, levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)

plt.show(block = False)  

#%%
# cell 13
  # Calculate the climatologies and test/examine by plotting a map for June 30
  
esr_mean, esr_std = CalculateClimatology(esr_GFS)

# Reuse the same color map, labels, and projections as the first map on lines
#   274 - 292

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

cs = ax.contourf(lon_grid, lat_grid, esr_mean[:,:,181].T, levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)

plt.show(block = False)  
#%%
# cell 14
  # Calculate SESR
sesr = np.ones((esr_GFS.shape)) * np.nan
climate_days = np.arange(1, 365+1)
NumYear = int(np.ceil(sesr.shape[-1]/365))
    
# Create a variable for each day, assumed starting at Jan 1 and no
#   leap years (i.e., each year is only 365 days each)
days = np.ones((sesr.shape[-1])) * np.nan
T = len(days)

n = 0
for i in range(1, NumYear+1):
    if i >= NumYear:
        days[n:T+1] = np.arange(1, len(days[n:T+1])+1)
    else:
        days[n:n+365] = np.arange(1, 365+1)
                                 
    n = n + 365

for day in climate_days:
    print('Working on day %i' %day)
    ind = np.where( days == day )[0]
    for t in ind:
        sesr[:,:,t] = (esr_GFS[:,:,t] - esr_mean[:,:,day-1])/esr_std[:,:,day-1]


#%%
# cell 15
  # Make an example plot of SESR

# Colorbar information
cmin = -3.0; cmax = 3.0; cint = 0.20
clevs = np.arange(cmin, cmax+cint, cint)
nlevs = len(clevs) - 1
# cmap = plt.get_cmap(name = 'Reds', lut = nlevs)
cmap  = plt.get_cmap(name = 'RdBu_r', lut = nlevs)

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

cs = ax.contourf(lon_grid, lat_grid, sesr[:,:,-1].T, levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)
#%%
# cell 16
  # Calculate the change in SESR
  
delta_sesr = sesr[:,:,1:] - sesr[:,:,:-1]
  
#%%
# cell 17
  # Calculate the climatology in the change in SESR
  
delta_sesr_mean, delta_sesr_std = CalculateClimatology(delta_sesr)

#%%
# cell 18
  # Write two sets of datetime variables. One for sesr and delta_sesr (every)
  #   date from the start point to the end point.
  # The other is for the climatologies and goes from Jan 01 to Dec 31
  
narr_start = datetime(1979, 1, 1)
narr_end   = datetime(2016, 12, 31)

narr_dates_gen = DateRange(narr_start, narr_end)
narr_dates = ['tmp'] * sesr.shape[-1]
n = 0
for date in narr_dates_gen:
    if date.strftime('%m-%d') == '02-29': # Exclude leap years
        pass
    else:
        narr_dates[n] = date.strftime('%Y-%m-%d')
        n = n + 1


year_start = datetime(1981, 1, 1)
year_end   = datetime(1981, 12, 31)

year_dates_gen = DateRange(year_start, year_end)
year_dates = ['tmp'] * 365
n = 0
for day in year_dates_gen:
    year_dates[n] = day.strftime('%m-%d')
    n = n + 1

#%%
# cell 19
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

#%%
# cell 20
  # Print the datasets and see if everything seems to have been made properly

path = './Data/SESR_Climatology/'
print(Dataset(path + 'esr_climatology_mean.nc', 'r'))
print(Dataset(path + 'esr_climatology_std.nc', 'r'))
print(Dataset(path + 'GFS_grid_sesr.nc', 'r'))
print(Dataset(path + 'GFS_grid_delta_sesr.nc', 'r'))
print(Dataset(path + 'delta_sesr_climatology_mean.nc', 'r'))
print(Dataset(path + 'delta_sesr_climatology_std.nc', 'r'))


#%%
# cell 21
  # Test the calculated SESR to the already given SESR pentad

# Find and load the given sesr pentad data
sesrPentad = load3Dnc('sesr_narr.nc', 'sesr')

# Convert to GFS grid
sesrPentadGFS, GFS_lat, GFS_lon = InterpolateToGFSGrid(sesrPentad, lat, lon)

# Plot the data to ensure it was interpolated correctly
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

cs = ax.contourf(lon_grid, lat_grid, sesrPentadGFS[:,:,-1].T, levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)

#%%
# cell 22
  # Create the date data

pentad_dates = ['tmp'] * sesrPentadGFS.shape[-1]
n = 0
m = 5
narr_dates_gen = DateRange(narr_start, narr_end) # Note that this has to be redone for every loop
for date in narr_dates_gen:
    if date.strftime('%m-%d') == '02-29': # Exclude leap years
        pass
    elif m < 4:
        m = m + 1
    else:
        pentad_dates[n] = date.strftime('%Y-%m-%d')
        n = n + 1
        m = 0

# Print the pentad dates to ensure they came out correctly.
print(pentad_dates)
    
#%%
# cell 23
  # Create a plot for Iowa from May to August (the output of pentad
  # and daily can be compared to Fig. 6a in Christian et al. 2019 to
  # check for correctness)

dateFMT = DateFormatter('%m-%d')

# Convert the strings to datetimes
pentad_datetime = np.asarray([datetime.strptime(d, '%Y-%m-%d') for d in pentad_dates])
narr_datetime   = np.asarray([datetime.strptime(d, '%Y-%m-%d') for d in narr_dates])

# Find the dates corresponding to the summer of 2012
pentadind = np.where( (pentad_datetime >= datetime(2012, 5, 1)) &
                     (pentad_datetime <= datetime(2012, 8, 31)) )[0]
narrdatesind = np.where( (narr_datetime >= datetime(2012, 5, 1)) &
                        (narr_datetime <= datetime(2012, 8, 31)) )[0]

# Find the coordinates for Iowa
lonind = np.where( (GFS_lon >= -96) & (GFS_lon <= -91) )[0]
latind = np.where( (GFS_lat >= 41) & (GFS_lat <= 43) )[0]


# Create some temporary variables that has latitudinally averaged values
narrtmp   = np.nanmean(sesr[lonind,:,:], axis = 0)
pentadtmp = np.nanmean(sesrPentadGFS[lonind,:,:], axis = 0)

# Create the variable that will be plotted
narrsesr   = np.nanmean(narrtmp[latind,:], axis = 0)
pentadsesr = np.nanmean(pentadtmp[latind,:], axis = 0)

# Create the plot
fig = plt.figure(figsize = [18, 12])
ax  = fig.add_subplot(1, 1, 1)

ax.set_title('Average SESR for Iowa during the summer of 2012', size = 25)

ax.plot(narr_datetime[narrdatesind], narrsesr[narrdatesind], 'b-', label = 'Daily SESR')
ax.plot(pentad_datetime[pentadind], pentadsesr[pentadind], 'r-.', label = 'Pentad SESR')

ax.legend(fontsize = 22, shadow = True)

ax.set_ylim([-3.5, 1])
ax.set_yticks(np.arange(-3, 1+1, 1))
ax.xaxis.set_major_formatter(dateFMT)

ax.set_ylabel('SESR (unitless)', size = 22)
ax.set_xlabel('Time', size = 22)

for i in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
    i.set_size(22)

for i in ax.xaxis.get_ticklabels():
    i.set_rotation(0)

path = './Figures/Tests/'
savename = 'Daily_SESR_vs_Pentad_SESR.png'
plt.savefig(path + savename, bbox_inches = 'tight')
plt.show(block = False)


#%%
# cell 24
  # Calculate pentads to compare to the given values

pentadTest = np.ones((sesrPentadGFS.shape)) * np.nan

n = 0
for i in range(0, sesr.shape[-1], 5):
    pentadTest[:,:,n] = np.nanmean(sesr[:,:,i:i+5], axis = -1) 
    n = n + 1
    
# Plot the test pentad to see how it appears
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

cs = ax.contourf(lon_grid, lat_grid, pentadTest[:,:,-1].T, levels = clevs, cmap = cmap, 
                  transform = data_proj, extend = 'both', zorder = 1)

cbax = fig.add_axes([0.92, 0.325, 0.02, 0.35])
cbar = fig.colorbar(cs, cax = cbax)

ax.set_extent([np.nanmin(lon), np.nanmax(lon), np.nanmin(lat), np.nanmax(lat)], 
                crs = fig_proj)


#%%
# cell 25
  # Plot the test pentad to see how it compares to the given


# Create some temporary variables that has latitudinally averaged values
testtmp   = np.nanmean(pentadTest[lonind,:,:], axis = 0)
pentadtmp = np.nanmean(sesrPentadGFS[lonind,:,:], axis = 0)

# Create the variable that will be plotted
testsesr   = np.nanmean(testtmp[latind,:], axis = 0)
pentadsesr = np.nanmean(pentadtmp[latind,:], axis = 0)

# Create the plot
fig = plt.figure(figsize = [18, 12])
ax  = fig.add_subplot(1, 1, 1)

ax.set_title('Average SESR for Iowa during the summer of 2012, calculated pentads from daily SESR', size = 18)

ax.plot(pentad_datetime[pentadind], testsesr[pentadind], 'b-', label = 'Calculated Pentad SESR')
ax.plot(pentad_datetime[pentadind], pentadsesr[pentadind], 'r-.', label = 'Pentad SESR')

ax.legend(fontsize = 18, shadow = True)

ax.set_ylim([-3.5, 1])

ax.set_ylabel('SESR (unitless)', size = 18)
ax.set_xlabel('Time', size = 18)

for i in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
    i.set_size(18)
    
for i in ax.xaxis.get_ticklabels():
    i.set_rotation(45)

path = './Figures/Tests/'
savename = 'Calculated_Pentad_SESR_vs_Pentad_SESR.png'
plt.savefig(path + savename, bbox_inches = 'tight')
    
plt.show(block = False)

#%%
# cell 26
  # There pentads will be used in criteria 3 and 4, so calculated dsesr for pentads, and climatologies
dsesr_pentad = pentadTest[:,:,1:] - pentadTest[:,:,:-1]

dsesr_pentad_mean, dsesr_pentad_std = CalculateClimatology(dsesr_pentad, pentad = True)

# Create the pentad climatology dates
pentad_year = ['tmp'] * int(365/5)
n = 0
m = 5
pentad_year_gen = DateRange(datetime(1900, 1, 1), datetime(1900, 12, 31)) # Note that this has to be redone for every loop
for date in pentad_year_gen:
    if date.strftime('%m-%d') == '02-29': # Exclude leap years
        pass
    elif m < 4:
        m = m + 1
    else:
        pentad_year[n] = date.strftime('%m-%d')
        n = n + 1
        m = 0


#%%
# cell 27
  # Write the pentad sesr, dsesr, and dsesr climatologies to .nc files

# Write the SESR pentad file
WriteNC(pentadTest, lat_GFS, lon_GFS, pentad_dates, 'GFS_grid_sesr_pentad.nc',
         'sesr')

# Write the Delta SESR pentad file
WriteNC(dsesr_pentad, lat_GFS, lon_GFS, pentad_dates[:-1], 'GFS_grid_delta_sesr_pentad.nc',
         'dsesr')

# Write the change in SESR pentad climatology files
WriteNC(dsesr_pentad_mean, lat_GFS, lon_GFS, pentad_year, 'delta_sesr_climatology_mean_pentad.nc',
         'dsesrm')

WriteNC(dsesr_pentad_std, lat_GFS, lon_GFS, pentad_year, 'delta_sesr_climatology_std_pentad.nc',
         'dsesrsstd')

#%%
# cell 28
  # Calculate pentads in esr, then sesr to see if this changes things
ESRpentadTest = np.ones((sesrPentadGFS.shape)) * np.nan

n = 0
for i in range(0, esr_GFS.shape[-1], 5):
    ESRpentadTest[:,:,n] = np.nanmean(esr_GFS[:,:,i:i+5], axis = -1) 
    n = n + 1


ESRpentadm, ESRpentadstd = CalculateClimatology(ESRpentadTest, pentad = True)


# Calculate SESR
sesrTest = np.ones((ESRpentadTest.shape)) * np.nan
climate_days = np.arange(1, 365/5+1)
yearLen = int(365/5)
NumYear = int(np.ceil(sesrTest.shape[-1]/yearLen))
    
# Create a variable for each day, assumed starting at Jan 1 and no
#   leap years (i.e., each year is only 365 days each)
daysTest = np.ones((sesrTest.shape[-1])) * np.nan
T = len(daysTest)

n = 0
for i in range(1, NumYear+1):
    if i >= NumYear:
        daysTest[n:T+1] = np.arange(1, len(daysTest[n:T+1])+1)
    else:
        daysTest[n:n+yearLen] = np.arange(1, yearLen+1)
                                 
    n = n + yearLen

for day in climate_days:
    print('Working on day %i' %day)
    ind = np.where( daysTest == day )[0]
    for t in ind:
        sesrTest[:,:,t] = (ESRpentadTest[:,:,t] - ESRpentadm[:,:,int(day-1)])/ESRpentadstd[:,:,int(day-1)]
        
#%%
# cell 29
 # plot the test SESR with given sesr pentad

# Create some temporary variables that has latitudinally averaged values
testtmp   = np.nanmean(sesrTest[lonind,:,:], axis = 0)
pentadtmp = np.nanmean(sesrPentadGFS[lonind,:,:], axis = 0)

# Create the variable that will be plotted
testsesr   = np.nanmean(testtmp[latind,:], axis = 0)
pentadsesr = np.nanmean(pentadtmp[latind,:], axis = 0)

# Create the plot
fig = plt.figure(figsize = [18, 12])
ax  = fig.add_subplot(1, 1, 1)

ax.set_title('Average SESR for Iowa during the summer of 2012', size = 18)

ax.plot(pentad_datetime[pentadind], testsesr[pentadind], 'b-', label = 'Daily SESR')
ax.plot(pentad_datetime[pentadind], pentadsesr[pentadind], 'r-.', label = 'Pentad SESR')

ax.legend(fontsize = 18, shadow = True)

ax.set_ylim([-3.5, 1])

ax.set_ylabel('SESR (unitless)', size = 18)
ax.set_xlabel('Time', size = 18)

for i in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
    i.set_size(18)
    
plt.show(block = False)













