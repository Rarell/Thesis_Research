#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 18 16:58:03 2019

@author: stuartedris

This is a script designed to intake the four week synthetic data of SESR, the
  GFS forecast of SESR, and climatology variables to identify and forecast
  flash drought (identification is at t = 0 in the GFS forecast, and the
  and the forecast is t > 0 in the GFS forecast). The method used for this is
  outlined in Christian et al. 2019
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
  # Name some user defined variables to change for the GFS forecast date and model run

Year  = '2012'
Month = '07'
Day   = '20'
ModelRun = '00'

esrmFN       = 'esr_climatology_mean.nc'
esrstdFN     = 'esr_climatology_std.nc'
sesrClimFN   = 'GFS_grid_sesr.nc'
sesrPentadFN = 'GFS_grid_sesr_pentad.nc'
dsesrClimFN  = 'GFS_grid_delta_sesr_pentad.nc'
dsesrmFN     = 'delta_sesr_climatology_mean_pentad.nc'
dsesrstdFN   = 'delta_sesr_climatology_std_pentad.nc'

syntheticFN = 'GFS_sesr_ed_' + Year + Month + Day + '.nc'

forecastFN  = 'GFS_sesr_' + Year + Month + Day + '_' + ModelRun + '.nc'

esrmSName      = 'esrm'
esrstdSName    = 'esrstd'
sesrClimSName  = 'sesr'
dsesrClimSName = 'dsesr'
dsesrmSName    = 'dsesrm'
dsesrstdSName  = 'dsesrsstd'

syntheticSName = 'sesr'

forecastSName  = 'sesr'

OutPath = './Figures/'

#%%
# cell 3
  # Examine the .nc files

# Examine climatology values
print(Dataset('./Data/SESR_Climatology/esr_climatology_mean.nc', 'r'))
print(Dataset('./Data/SESR_Climatology/esr_climatology_std.nc', 'r'))
print(Dataset('./Data/SESR_Climatology/GFS_grid_sesr.nc', 'r'))
print(Dataset('./Data/SESR_Climatology/GFS_grid_delta_sesr.nc', 'r'))
print(Dataset('./Data/SESR_Climatology/delta_sesr_climatology_mean.nc', 'r'))
print(Dataset('./Data/SESR_Climatology/delta_sesr_climatology_std.nc', 'r'))

# Examine the synthetic data
print(Dataset('./Data/Synthetic_Data/GFS_sesr_ed_20190901.nc', 'r'))

# Examine the forecast data
print(Dataset('./Data/GFS_sesr_20190901_00.nc', 'r'))


#%%
# cell 4
  # Create a function to import the climatology and annual data
def load_climatology(SName, file, path = './Data/SESR_Climatology/'):
    '''

    '''
    
    X = {}
    DateFormat = '%m-%d'
    
    with Dataset(path + file, 'r') as nc:
        # Load the grid
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        
        lat, lon = np.meshgrid(lat, lon)
        
        X['lat'] = lat
        X['lon'] = lon
        
        # Collect the time information
        time = nc.variables['date'][:]
        dates = np.asarray([datetime.strptime(time[d], DateFormat) for d in range(len(time))])
        
        X['date'] = dates
        X['month'] = np.asarray([d.month for d in dates])
        X['day']   = np.asarray([d.day for d in dates])
        X['ymd']   = np.asarray([datetime(d.year, d.month, d.day) for d in dates])
        
        # Collect the data itself
        X[str(SName)] = nc.variables[str(SName)][:,:,:]
        
    return X

#%%
# cell 5
  # Create a function to import the climatology and annual data
def load_full_climatology(SName, file, path = './Data/SESR_Climatology/'):
    '''

    '''
    
    X = {}
    DateFormat = '%Y-%m-%d'
    
    with Dataset(path + file, 'r') as nc:
        # Load the grid
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        
        lat, lon = np.meshgrid(lat, lon)
        
        X['lat'] = lat
        X['lon'] = lon
        
        # Collect the time information
        time = nc.variables['date'][:]
        dates = np.asarray([datetime.strptime(time[d], DateFormat) for d in range(len(time))])
        
        X['date'] = dates
        X['year']  = np.asarray([d.year for d in dates])
        X['month'] = np.asarray([d.month for d in dates])
        X['day']   = np.asarray([d.day for d in dates])
        X['ymd']   = np.asarray([datetime(d.year, d.month, d.day) for d in dates])
        
        # Collect the data itself
        X[str(SName)] = nc.variables[str(SName)][:,:,:]
        
    return X
        
        
#%%
# cell 6
  # Create a function to import the synthetic data
def load_synthetic(SName, file, path = './Data/Synthetic_Data/'):
    '''
    
    '''
    
    X = {}
    DateFormat = '%Y-%m-%d'
    
    with Dataset(path + file, 'r') as nc:
        # Load the grid
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        
        lat, lon = np.meshgrid(lat, lon)
        
        X['lat'] = lat
        X['lon'] = lon
        
        # Collect the time information
        time = nc.variables['date'][:]
        dates = np.asarray([datetime.strptime(time[d], DateFormat) for d in range(len(time))])
        
        X['date'] = dates
        X['year']  = np.asarray([d.year for d in dates])
        X['month'] = np.asarray([d.month for d in dates])
        X['day']   = np.asarray([d.day for d in dates])
        X['ymd']   = np.asarray([datetime(d.year, d.month, d.day) for d in dates])
        
        # Load the forecast hour that the synthetic data was created at
        X['FH'] = nc.variables['FH'][:]
        
        # Load the actual data
        X[str(SName)] = nc.variables[str(SName)][:,:,:]
        
    return X

#%%
# cell 7
  # Create a function to load the forecast data
def load_forecast(SName, file, path = './Data/'):
    '''
    '''
    
    X = {}
    DateFormat = '%Y-%m-%d'
    
    with Dataset(path + file, 'r') as nc:
        # Load the grid
        lat = nc.variables['lat'][:]
        lon = nc.variables['lon'][:]
        
        lat, lon = np.meshgrid(lat, lon)
        
        X['lat'] = lat
        X['lon'] = lon
        
        # Collect the time information
        time = nc.variables['date'][:]
        dates = np.asarray([datetime.strptime(time[d], DateFormat) for d in range(len(time))])
        
        X['date'] = dates
        X['year']  = np.asarray([d.year for d in dates])
        X['month'] = np.asarray([d.month for d in dates])
        X['day']   = np.asarray([d.day for d in dates])
        X['ymd']   = np.asarray([datetime(d.year, d.month, d.day) for d in dates])
        
        # Load the forecast data
        X[str(SName)] = nc.variables[str(SName)][:,:,:]
        
    return X

#%%
# cell 8
  # Load the data

esrm      = load_climatology('esrm', 'esr_climatology_mean.nc')
Climsesr  = load_full_climatology(sesrClimSName, sesrClimFN)
Climdsesr = load_full_climatology(dsesrClimSName, dsesrClimFN)
Pensesr   = load_full_climatology(sesrClimSName, sesrPentadFN)
dsesrm    = load_climatology(dsesrmSName, dsesrmFN)
dsesrstd  = load_climatology(dsesrstdSName, dsesrstdFN)

Syn  = load_synthetic(syntheticSName, syntheticFN)

Forc = load_forecast(forecastSName, forecastFN)


#%%
# cell 9
  # Create a map of the last time in the synthetic data and first in forecast to
  #   ensure they are the same

# Colorbar information
cmin = -3.0; cmax = 3.0; cint = 0.2
clevs = np.arange(cmin, cmax + cint, cint)
nlevs = len(clevs) - 1
cmap = plt.get_cmap(name = 'RdBu_r', lut = nlevs)

# Projection informatino
data_proj = ccrs.PlateCarree()
fig_proj  = ccrs.PlateCarree()

# Synthetic data figure
fig = plt.figure(figsize = [12, 18])
ax  = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()

cs = ax.contourf(Syn['lon'], Syn['lat'], Syn['sesr'][:,:,-1].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'both')

cbax = fig.add_axes([0.12, 0.35, 0.78, 0.02])
cbar = fig.colorbar(cs, cax = cbax, orientation = 'horizontal')

plt.show()


# Forecast data figure
fig = plt.figure(figsize = [12, 18])
ax  = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()

cs = ax.contourf(Forc['lon'], Forc['lat'], Forc['sesr'][:,:,0].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'both')

cbax = fig.add_axes([0.12, 0.35, 0.78, 0.02])
cbar = fig.colorbar(cs, cax = cbax, orientation = 'horizontal')

plt.show()


#%%
# cell 10
  # Create a new time variable. This variable is 0 for forecast hour 0/realtime
  #  analysis. It is > 0 for forecast, and < 0 for synthetic data.

# Initialize the variable

T = len(Syn['date']) - 1 + len(Forc['date'])

counter = np.ones((T, )) * np.nan

for n, i in enumerate(range(-1 * (len(Syn['date']) - 1), len(Forc['date']))):
    counter[n] = i
    
# Check the counter to ensure the calculations are correct
print(counter)

#%%
# cell 11
  # Combine the synthetic and forecast datasets into a single set.

###################
### Important Note ###  
###################
# Forecast time = 0 is taken as real time here, because the synthetic data ends
#   at 00Z, so its last day only has one datapoint for its mean, and the forecast
#   has all of them (it starts at 00Z) for its mean. This needs to be 
#   adjustable for if the synthetic ends on 18Z (has a full day mean), and forecast
#   starts at 18Z (only one or two values in its mean), or any value inbetween


# Initialize the combined variables
J, I, T = Syn['sesr'].shape
T = counter.size

sesr = np.ones((J, I, T)) * np.nan

# Create a combined time variable
time = np.concatenate((Syn['date'][:-1], Forc['date']), axis = 0)

# Create a combined SESR variable

ind = np.where(counter < 0)[0]
sesr[:,:,ind] = Syn['sesr'][:,:,:-1]

ind = np.where(counter >= 0)[0]
sesr[:,:,ind] = Forc['sesr'][:,:,:]

# Check the new variables
print(time)
print(sesr)


#%%
# cell 12
  # Create a plot of SESR for GFS and NARR for comparison
  # Note: This focuses on Iowa, 2012 a it has a reference figure in Fig 6a, Christien et al. 2019
    
lonmin = -96
lonmax = -91
latmin = 41
latmax = 43

dateFMT = DateFormatter('%m-%d')

# Find the year corresponding to the GFS interval
narrdatesind = np.where( (Climsesr['ymd'] >= time[0]) &
                        (Climsesr['ymd'] <= time[-1]) )[0]

# Find the coordinates for the given location
lonind = np.where( (Syn['lon'][:,0] >= lonmin) & (Syn['lon'][:,0] <= lonmax) )[0]
latind = np.where( (Syn['lat'][0,:] >= latmin) & (Syn['lat'][0,:] <= latmax) )[0]

# Create some temporary variables that has latitudinally averaged values
narrtmp = np.nanmean(Climsesr['sesr'][lonind,:,:], axis = 0)
gfstmp  = np.nanmean(sesr[latind,:,:], axis = 0)

# Create the variables that will be plotted
narrsesr = np.nanmean(narrtmp[latind,:], axis = 0)
gfssesr  = np.nanmean(gfstmp[lonind,:], axis = 0)

# Create the plot
fig = plt.figure(figsize = [18, 12])
ax  = fig.add_subplot(1, 1, 1)

ax.set_title('NARR and GFS SESR time series for Iowa in June and July, 2012', size = 25)

ax.plot(time, gfssesr, 'b-', label = 'GFS')
ax.plot(Climsesr['ymd'][narrdatesind], narrsesr[narrdatesind], 'r-.', label = 'NARR')
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
    
savename = 'SESR_timeseries_for_NARR_GFS.png'
plt.savefig('./Figures/' + savename, bbox_inches = 'tight')
    
plt.show(block = False)

#%%
# cell 13
# Standardize the narr dsesr.
ClimSdsesr = np.ones((I, J, Climdsesr['dsesr'].shape[-1])) * np.nan
for n, t in enumerate(Climdsesr['date']):
    date = datetime(1900, Climdsesr['date'][n].month, Climdsesr['date'][n].day)
    ind  = np.where( dsesrm['date'] == date )[0]
    ClimSdsesr[:,:,n] = (Climdsesr['dsesr'][:,:,n] - dsesrm['dsesrm'][:,:,ind[0]])/dsesrstd['dsesrsstd'][:,:,ind[0]]

#%%
# cell 14
  # Use the same algorithm to identify flash drought for the NARR to test it and ensure it works.

# Initialize the criteria variables for flash drought identification
I, J, T   = Climsesr['sesr'].shape
T = 365

narrcrit1 = np.ones((J, I, T)) * np.nan
narrcrit2 = np.ones((J, I, T)) * np.nan
narrcrit3 = np.ones((J, I, T)) * np.nan
narrcrit4 = np.ones((J, I, T)) * np.nan

narr_time = np.asarray([datetime(2012, d.month, d.day) for d in esrm['ymd']]) # Create a time variable for a full year

#%%
# cell 15
  # Find criteria for the NARR

StartDate  = Climdsesr['ymd'][-1] # Initialize start_date so it always has some value
MinChange  = timedelta(days = 30)
mdDates    = np.asarray([datetime(1900, d.month, d.day) for d in Climdsesr['date']]) # Create a month/day array with all the months and days in the pentad data
CritPercentile = 40
NumExceptions  = 1
count = -99
FiveDays = timedelta(days = 5)
OneDay   = timedelta(days = 1)
ZeroDays = timedelta(days = 0)
criteria3 = 0 # Start by assuming criterias 3 and 4 are false
criteria4 = 0
# m = 0
# n = 0

for j in range(J):
    print('Currently %4.2f %% done.' %(j/J * 100))
    for i in range(I):
        StartDate = Climdsesr['ymd'][-1] # Reset the start date so it does not carry over to the next grid point.
        count = -99 # Reset the counter
        criteria3 = 0 # Reset criteria 3 and 4 at the start of each grid point.
        criteria4 = 0
        for t, date in enumerate(narr_time):
            # Find all days in the full dataset equal to the current day
            ind = np.where( (Climsesr['month'] == narr_time[t].month) & 
                           (Climsesr['day'] == narr_time[t].day) )[0]
            
            # Find all days in the pentad dataset equal to the current day
            penind = np.where( (Climdsesr['month'] == narr_time[t].month) & 
                              (Climdsesr['day'] == narr_time[t].day) )[0]
            
            # Find the current date indice
            tind = np.where(narr_time[t] == Climsesr['ymd'])[0]
            
            # Find the current pentad indice
            DateDelta = date - Pensesr['ymd']
            if (np.mod(date.year, 4) == 0) & (date.month == 3) & (date.day == 1): # Exclude leap years. Feburary 29 can add an extra day, making time delta = 5 days at March 1.
                DateDelta = DateDelta - OneDay
            else:
                DateDelta = DateDelta
                
            pentind   = np.where( (DateDelta < FiveDays) & (DateDelta >= ZeroDays) )[0]
#            print(Climdsesr['date'][ind])
            
            
            ## Calculate the 20th quantile for criteria 2 ##
            percent20 = np.nanpercentile(Climsesr['sesr'][i,j,ind], 20)
            
            
            ## Calculate the quantile for criteria 3 that location and time #
            Crit3Quant = np.nanpercentile(ClimSdsesr[i,j,penind], CritPercentile)
            
            # Determine if this loop is the start of a flash drought
            Crit3Part2 = (criteria3 == 0)
            
            
            ## Determine the mean change in SESR between start and end dates ##
            StartInd   = np.where(Climdsesr['ymd'] == StartDate)[0]
            MeanChange = np.nanmean(ClimSdsesr[i,j,StartInd[0]:pentind[0]]) # Note this exclude the current pentad, which is still being tested
                                                                         # Also note if StardInd > pentind (as in the initialized case), the slice returns an empty array.
            
            # Determine the mean change in SESR between start and end date for all years
            tmpDateStart = datetime(1900, Climdsesr['month'][StartInd], Climdsesr['day'][StartInd])
            tmpDateEnd   = datetime(1900, Climdsesr['month'][pentind], Climdsesr['day'][pentind])
            Crit4Ind = np.where( (mdDates >= tmpDateStart) & (mdDates <= tmpDateEnd) )[0] # This could cause problems for winter droughts (month 12 turns back 1), but should be fine since they are unimportant.
            
            # Calculate the 25th percentile for the mean changes
            percent25 = np.nanpercentile(ClimSdsesr[i,j,Crit4Ind], 25)
            
            ### Determine criteria 3 for the given pentad ###
            if (ClimSdsesr[i,j,pentind] <= Crit3Quant) & (Crit3Part2 == 1):
                criteria3 = 1
                count = 0
                StartDate = date
            elif (ClimSdsesr[i,j,pentind] <= Crit3Quant):
                criteria3 = 1
            elif (ClimSdsesr[i,j,pentind] > Crit3Quant) & (abs(count) < NumExceptions):
                criteria3 = 1
                count = count + 1
            elif (ClimSdsesr[i,j,pentind] > Crit3Quant) & (abs(count) >= NumExceptions): # This should catch all points when criteria 3 fails
                criteria3 = 0
                count = -99
            else:
                criteria3 = criteria3 # Do nothing. Between pentads, Crit3Quant is nan (penind is empty), so if block should come to here.
                                      #   Pentad dates mark the beginning of the pentad, so between pentad dates (i.e., within the pentad
                                      #   in which criteria 3 is true or false) leave criteria 3 unchanged as it remains the same within the pentad.
            
            # Assign criteria 3 for each day
            if criteria3 == 1:
                narrcrit3[j,i,t] = 1
            else:
                narrcrit3[j,i,t] = 0
               
            
            ### Determine criteria 4 for the given pentad ###
            if (MeanChange <= percent25):# | (narrcrit3[j,i,t] == 0):
                criteria4 = 1
            elif (MeanChange > percent25):
                criteria4 = 0
            else:
                criteria4 = criteria4
            
            # Assign criteria 4 for each day
            if criteria4 == 1:
                narrcrit4[j,i,t] = 1
            else:
                narrcrit4[j,i,t] = 0
                
                
            #### Determine Criteria 1 ###
            # Five days are added here because, while date marks the beginning of Delta SESR, the pentad at the end of Delta SESR
            #   must be included as well, which is five days. 
            if (( (date - StartDate) + FiveDays ) < MinChange):# | (narrcrit3[j,i,t] == 0) ):
                narrcrit1[j,i,t] = 0
            else:
                narrcrit1[j,i,t] = 1
                
                
            ### Determine criteria 2 ###
            if Climsesr['sesr'][i,j,tind] < percent20:
                narrcrit2[j,i,t] = 1
            else:
                narrcrit2[j,i,t] = 0
                
            # if (j == latind[0]) & (i == lonind[0]):
            #     print(date - StartDate, ClimSdsesr[i,j,pentind], Crit3Quant, date)
            #     print(count, criteria3)
                # print(ClimSdsesr[i,j,pentind], Crit3Quant, date)
                # print(Climdsesr['dsesr'][i,j,pentind])
        
                        
print(narrcrit1)
print('\n')
print(narrcrit2)
print('\n')                                         
print(narrcrit3)
print('\n')
print(narrcrit4)


#%%
# cell 23
  # Make a time series plot of criteria 3 for the NARR to see how it looks. To start, focus
  # focus it on 34 N and 100 W or 42 N and 94 W

  # This is used to help fine tune criteria 3 for daily SESR instead of pentads

StartDate = datetime(2012, 5, 1)
EndDate   = datetime(2012, 8, 31)
FocusLon = -94
FocusLat = 42
lonind = np.where(Syn['lon'][:,0] == FocusLon)[0]
latind = np.where(Syn['lat'][0,:] == FocusLat)[0]
dateind = np.where( (narr_time >= StartDate) & (narr_time <= EndDate) )[0]

tmpcrit3 = np.nanmean(narrcrit3[latind,:,:], axis = 0)
plotcrit3 = np.nanmean(tmpcrit3[lonind,:], axis = 0)

fig = plt.figure(figsize = [18, 12])
ax  = fig.add_subplot(1, 1, 1)

ax.set_title('Criteria 3 time series for %4.0f N and %4.0f E' %(FocusLat, FocusLon), size = 18)

#ax.plot(time, crit3[latind[0],lonind[0],:], 'k')
ax.plot(Climsesr['ymd'][dateind], plotcrit3[dateind], 'k')

ax.set_ylim([-0.1, 1.1])

ax.set_ylabel('Criteria 3 (0 Fail, 1 Succeed)', size = 18)
ax.set_xlabel('Time', size = 18)

for i in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
    i.set_size(18)

for i in ax.xaxis.get_ticklabels():
    i.set_rotation(45)

plt.show(block = False)



#%%
# cell 25
  # Initialize and identify the flash drought with the NARR data

narrFD = np.ones((J, I, T)) * np.nan

for t in range(T):
    for j in range(J):
        for i in range(I):
            FlashDrought = ((narrcrit1[j,i,t] == 1) & (narrcrit2[j,i,t] == 1) & 
                            (narrcrit3[j,i,t] == 1) & (narrcrit4[j,i,t] == 1))
            # Determine if all criteria have tested true for a given grid point
            if FlashDrought == 1:
                narrFD[j,i,t] = 1
            elif  (t != 0) & (narrFD[j,i,t-1] != 0):
                # narrFD[j,i,t] = 2
                narrFD[j,i,t] = 1
            else:
                narrFD[j,i,t] = 0
                
            if ((narr_time[t] >= datetime(2012, 7, 15)) & (narr_time[t] <= datetime(2012, 7, 25)) & 
                (FlashDrought == 1)):
                narrFD[j,i,t] = 2
                
print(narrFD)

#%%
# cell 26
  # Plot the flash drought for the NARR

# Use this to select the desired date to examine
PlotDate = datetime(2012, 7, 20)
dateind  = np.where(narr_time == PlotDate)[0]

# Color information
cmin = 0; cmax = 2; cint = 0.5
clevs = np.arange(cmin, cmax + cint, cint)
nlevs = len(clevs)
cmap = plt.get_cmap(name = 'binary', lut = nlevs)


fig = plt.figure(figsize = [12, 18])
ax  = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)

cs = ax.contourf(Syn['lon'], Syn['lat'], narrFD[:,:,dateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax.set_extent([-105, -85, 30, 50])

# cbax = fig.add_axes([0.12, 0.35, 0.78, 0.02])
# cbar = fig.colorbar(cs, cax = cbax, orientation = 'horizontal')

plt.show()

#%%
  # Calculate the SESR pentads for the GFS data.

# Find the nearest labeled pentad to the start of the investigated time
DeltaDate = time[0] - Pensesr['ymd']
if (np.mod(time[0].year, 4) == 0) & (time[0].month == 3) & (time[0].day == 1): # Exclude leap years. Feburary 29 can add an extra day, making time delta = 5 days at March 1.
    DeltaDate = DeltaDate - timedelta(days = 1)
else:
    DateDelta = DateDelta
    
ptind = np.where( (DeltaDate < timedelta(days = 5)) & (DeltaDate >= timedelta(days = 0)) )[0]

m = DeltaDate[ptind[0]].days # Initialize m so that the calculate pentads line up with climatology.
                             #   E.g., the nearet pentad is March 1, and the beginning of the dataset is March 4,
                             #   then m is initialized to 3 so that remaining pentad start dates line up with the
                             #   calculated climatology. Note this means the synthetic data should be at least 34 days
                             #   so that the "true start" of the dataset includes a full pentad. dsesr[0] and dsesr[-1]
                             #   are also problematic, since they may not contain a full five days. But dsesr[0] is not
                             #   the "true start" so is not as important. dsesr[-1] is the end of the model run and is
                             #   already part of "fantasy land" as a result, and is thus less important.

# Initialize the sesr pentad data
J, I, T = sesr.shape

sesr_pentad = np.ones((J, I, int(T/5+2))) * np.nan # At most, sesr_pentad can have T/5 (pentad) + 2 time points 
                                                   #   (one for beginning and end which may not include a full pentad).
                                                   #   This may have some unused time points.
time_pentad_str = ['tmp'] * int(T/5+2) # Initialize a time variable to help keep track of the pentad data.

n = 0
for i in range(0, sesr.shape[-1], 5):
    if i == 0: # for the first time point, average from m to i+5 to ensure the remaining pentads correspond to climatology pentads
        sesr_pentad[:,:,n] = np.nanmean(sesr[:,:,m:i+5], axis = -1)
        time_pentad_str[n] = time[i].strftime('%Y-%m-%d')
        n = n + 1
        # print(time[i])
    elif (sesr.shape[-1] - i) < 5:
        sesr_pentad[:,:,n] = np.nanmean(sesr[:,:,i:-1], axis = -1) 
        time_pentad_str[n] = time[i].strftime('%Y-%m-%d')
        # print(time[i])
    else:
        sesr_pentad[:,:,n] = np.nanmean(sesr[:,:,i:i+5], axis = -1) 
        time_pentad_str[n] = time[i].strftime('%Y-%m-%d')
        n = n + 1
        # print(time[i])

# Remove any unused time points that may exist
n = 0
for t in range(sesr_pentad.shape[-1]): 
    if np.isnan(np.nanmean(sesr_pentad[:,:,n])): # If the entire grid is nan (i.e., unused time point) np.nanmean returns nan)
        sesr_pentad = np.delete(sesr_pentad, obj = n, axis = -1)
        n = n - 1
    else:
        pass
    
    n = n + 1

# Find and remove unused time points in the time variable
n = 0
for t in range(len(time_pentad_str)):
    if time_pentad_str[n] == 'tmp':
        time_pentad_str.remove(time_pentad_str[n])
        n = n - 1
    else:
        pass
    
    n = n + 1

# Convert the time variable to a datetime
time_pentad = np.asarray([datetime.strptime(d, '%Y-%m-%d') for d in time_pentad_str])

print(sesr_pentad)
print(time_pentad)

#%%
# cell 12
  # Calculate delta sesr

J, I, T = sesr_pentad.shape

dsesr = np.ones((J, I, T-1)) * np.nan

dsesr = sesr_pentad[:,:,1:] - sesr_pentad[:,:,:-1]

#%%
# cell 13
  # Standardize delta sesr
Sdsesr = np.ones((J, I, T-1)) * np.nan

for n, t in enumerate(time_pentad[:-1]):
    date = datetime(1900, time_pentad[n].month, time_pentad[n].day)
    ind  = np.where( dsesrm['date'] == date )[0]
    Sdsesr[:,:,n] = (dsesr[:,:,n] - dsesrm['dsesrm'][:,:,ind[0]].T)/dsesrstd['dsesrsstd'][:,:,ind[0]].T


#%%
# cell 15
 # Initialize the criteria variables for flash drought identification

J, I, T = sesr.shape

crit1 = np.ones((J, I, T)) * np.nan
crit2 = np.ones((J, I, T)) * np.nan
crit3 = np.ones((J, I, T)) * np.nan
crit4 = np.ones((J, I, T)) * np.nan

#%%
# cell 23
  # Perform Flash Drought analysis.

############
### Important Note ###
############
# This method uses a pixel by pixel check to account for heteorgenuity, but
#  there should be better ways to do this, and this should be updated when
#  they are known.

# Criteria 1:
#   The flash drought must be at least 30 days in length

# Criteria 2:
#   The SESR value on the final day of the drought must be below the
#   20th percentile for that location and day.

# Criteria 3:
#   a) The change in SESR (delta SESR) must be below the 40th percentile
#      for that location and day (start day of the change).
#   b) The change in SESR (delta SESR) is allowed to be above the the 
#      40th percentile once (1) and only once.

# Criteria 4:
#   The mean change in SESR between start and end dates must be below
#   the 25th percentile of mean delta SESRs for that date range and location

# For purposes of a GFS analysis, each day in the GFS forecast and realtime
#   analysis are treated as a 'final day' for the drought. Hence, this loop
#   is for all times when counter >= 0.

#tmp = np.ones((time.size)) * np.nan

StartDate  = time[-1] # Initialize start_date so it always has some value
MinChange  = timedelta(days = 30)
mdDates    = np.asarray([datetime(1900, d.month, d.day) for d in Climdsesr['date']]) # Create a month/day array with all the months and days in the pentad data
CritPercentile = 40
NumExceptions  = 1
count = -99
FiveDays = timedelta(days = 5)
OneDay   = timedelta(days = 1)
ZeroDays = timedelta(days = 0)
criteria3 = 0 # Start by assuming criterias 3 and 4 are false
criteria4 = 0

for j in range(J):
    print('Currently %4.2f %% done.' %(j/J * 100))
    for i in range(I):
        StartDate = time[-1] # Reset the start date so it does not carry over to the next grid point.
        count = -99 # Reset the counter
        criteria3 = 0 # Reset criteria 3 and 4 at the start of each grid point.
        criteria4 = 0
        for t, date in enumerate(time):
            # Find all days in the full dataset equal to the current day
            ind = np.where( (Climsesr['month'] == time[t].month) & 
                           (Climsesr['day'] == time[t].day) )[0]
            
            # Find all days in the pentad dataset equal to the current day
            penind = np.where( (Climdsesr['month'] == time[t].month) & 
                              (Climdsesr['day'] == time[t].day) )[0]
            
            # Find the current pentad indice for the GFS data
            DateDelta = date - time_pentad
            if (np.mod(date.year, 4) == 0) & (date.month == 3) & (date.day == 1): # Exclude leap years. Feburary 29 can add an extra day, making time delta = 5 days at March 1.
                DateDelta = DateDelta - OneDay
            else:
                DateDelta = DateDelta
                
            ptind   = np.where( (DateDelta < FiveDays) & (DateDelta >= ZeroDays) )[0]
            
            ## Calculate the 20th quantile for criteria 2 ##
            percent20 = np.nanpercentile(Climsesr['sesr'][i,j,ind], 20)
            
            
            ## Calculate the quantile for criteria 3 that location and time #
            Crit3Quant = np.nanpercentile(ClimSdsesr[i,j,penind], CritPercentile)
            
            # Determine if this loop is the start of a flash drought
            Crit3Part2 = (criteria3 == 0)
            
            
            ## Determine the mean change in SESR between start and end dates ##
            StartInd   = np.where(time == StartDate)[0]
            MeanChange = np.nanmean(Sdsesr[j,i,StartInd[0]:ptind[0]]) # Note this exclude the current pentad, which is still being tested
                                                                      # Also note if StardInd > pentind (as in the initialized case), the slice returns an empty array.
            
            # Determine the mean change in SESR between start and end date for all years
            tmpDateStart = datetime(1900, StartDate.month, StartDate.day)
            tmpDateEnd   = datetime(1900, time_pentad[ptind[0]].month, time_pentad[ptind[0]].day)
            Crit4Ind = np.where( (mdDates >= tmpDateStart) & (mdDates < tmpDateEnd) )[0] # This could cause problems for winter droughts (month 12 turns back 1), but should be fine since they are unimportant.
            
            # Calculate the 25th percentile for the mean changes
            percent25 = np.nanpercentile(ClimSdsesr[i,j,Crit4Ind], 25)
            
            
             ### Determine criteria 3 for the given pentad ###
            if time_pentad[ptind[0]] == time_pentad[-1]:
                criteria3 = 0 # The final pentad does not have another pentad to compare to, so no Delta SESR. Thus assume criteria3 = 0
            else:
                if (Sdsesr[j,i,ptind] <= Crit3Quant) & (Crit3Part2 == 1):
                    criteria3 = 1
                    count = 0
                    StartDate = date
                elif (Sdsesr[j,i,ptind] <= Crit3Quant):
                    criteria3 = 1
                elif (Sdsesr[j,i,ptind] > Crit3Quant) & (abs(count) < NumExceptions):
                    criteria3 = 1
                    count = count + 1
                elif (Sdsesr[j,i,ptind] > Crit3Quant) & (abs(count) >= NumExceptions): # This should catch all points when criteria 3 fails
                    criteria3 = 0
                    count = -99
                else:
                    criteria3 = criteria3 # Do nothing. Between pentads, Crit3Quant is nan (penind is empty), so if block should come to here.
                                          #   Pentad dates mark the beginning of the pentad, so between pentad dates (i.e., within the pentad
                                          #   in which criteria 3 is true or false) leave criteria 3 unchanged as it remains the same within the pentad.
            
            # Assign criteria 3 for each day
            if criteria3 == 1:
                crit3[j,i,t] = 1
            else:
                crit3[j,i,t] = 0
               
            
            ### Determine criteria 4 for the given pentad ###
            if (MeanChange <= percent25):# | (narrcrit3[j,i,t] == 0):
                criteria4 = 1
            elif (MeanChange > percent25):
                criteria4 = 0
            else:
                criteria4 = criteria4
            
            # Assign criteria 4 for each day
            if criteria4 == 1:
                crit4[j,i,t] = 1
            else:
                crit4[j,i,t] = 0
            
            
            #### Determine Criteria 1 ###
            # Five days are added here because, while date marks the beginning of Delta SESR, the pentad at the end of Delta SESR
            #   must be included as well, which is five days. 
            if (( (date - StartDate) + FiveDays ) < MinChange):# | (narrcrit3[j,i,t] == 0) ):
                crit1[j,i,t] = 0
            else:
                crit1[j,i,t] = 1
                
                
            ### Determine criteria 2 ###
            if sesr[j,i,t] < percent20:
                crit2[j,i,t] = 1
            else:
                crit2[j,i,t] = 0


        
                        
print(crit1)
print('\n')
print(crit2)
print('\n')                                       
print(crit3)
print('\n')
print(crit4)


#%%
# cell 17
  # Make a time series plot of criteria 3 to see how it looks. To start, focus
  # focus it on 34 N and 100 W or 42 N and 94 W

  # This is used to help fine tune criteria 3 for daily SESR instead of pentads

FocusLon = -94
FocusLat = 42
lonind = np.where(Syn['lon'][:,0] == FocusLon)[0]
latind = np.where(Syn['lat'][0,:] == FocusLat)[0]

tmp = np.nanmean(crit3[latind,:,:], axis = 0)

fig = plt.figure(figsize = [18, 12])
ax  = fig.add_subplot(1, 1, 1)

ax.set_title('Criteria 3 time series for %4.0f N and %4.0f E' %(FocusLat, FocusLon), size = 18)

#ax.plot(time, crit3[latind[0],lonind[0],:], 'k')
ax.plot(time, np.nanmean(tmp[lonind,:], axis = 0), 'k')

ax.set_ylim([-0.1, 1.1])

ax.set_ylabel('Criteria 3 (0 Fail, 1 Succeed)', size = 18)
ax.set_xlabel('Time', size = 18)

for i in ax.xaxis.get_ticklabels() + ax.yaxis.get_ticklabels():
    i.set_size(18)

plt.show(block = False)

    

#%%
# cell 19
  # Flash drought occurs where all four criteria are true. Determine this.

FD = np.ones((J, I, T)) * np.nan

for t in range(T):
    for j in range(J):
        for i in range(I):
            
            # Determine if all criteria have tested true for a given grid point
            if ((crit1[j,i,t] == 1) & (crit2[j,i,t] == 1) & 
                (crit3[j,i,t] == 1) & (crit4[j,i,t] == 1)):
                FD[j,i,t] = 1
            elif  (t != 0) & (FD[j,i,t-1] != 0):
                # FD[j,i,t] = 2
                FD[j,i,t] = 1
            else:
                FD[j,i,t] = 0
                
print(FD)


#%%
# cell 20
  # Print a flash drought map and see how it comes out

# Use this to select the desired date to examine
PlotDate = datetime(2012, 7, 20)
dateind  = np.where(time == PlotDate)[0]

# Color information
cmin = 0; cmax = 2; cint = 0.5
clevs = np.arange(cmin, cmax + cint, cint)
nlevs = len(clevs)
cmap = plt.get_cmap(name = 'binary', lut = nlevs)


fig = plt.figure(figsize = [12, 18])
ax  = fig.add_subplot(1, 1, 1, projection = fig_proj)

ax.coastlines()
ax.add_feature(cfeature.BORDERS)
ax.add_feature(cfeature.STATES)

cs = ax.contourf(Syn['lon'], Syn['lat'], FD[:,:,dateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax.set_extent([-105, -85, 30, 50])

# cbax = fig.add_axes([0.12, 0.35, 0.78, 0.02])
# cbar = fig.colorbar(cs, cax = cbax, orientation = 'horizontal')

plt.show()


#%%
  # The next few cells are specifically for creating two figures to be used on the poster for this project

# Lat/Lon tick information
LatInt = 5
LonInt = 10

LatLabel = np.arange(5, 80, LatInt)
LonLabel = np.arange(-170, -70, LonInt)

LatFormatter = cticker.LatitudeFormatter()
LonFormatter = cticker.LongitudeFormatter()

# Use this to select the desired date to examine
PlotDate = datetime(2012, 7, 20)
dateind  = np.where(time == PlotDate)[0]

# Color information
cmin = 0; cmax = 2; cint = 0.5
clevs = np.arange(cmin, cmax + cint, cint)
nlevs = len(clevs)
cmap = mcolors.LinearSegmentedColormap.from_list("", ["white", "black", "red"], 3)

fig = plt.figure(figsize = [12, 18], frameon = True)
fig.suptitle('Flash Drought Criteria for ' + time[dateind[0]].strftime('%Y/%m/%d') + ' (0 hour forecast for GFS)', y = 0.6, size = 20)

### Begin the plots for the GFS data ###
# Set the first part of the figure
ax1 = fig.add_subplot(1, 4, 1, projection = fig_proj)

ax1.set_title('Criteria 1', size = 18)

ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.STATES)

ax1.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax1.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax1.set_yticklabels(LatLabel, fontsize = 14)

ax1.xaxis.tick_bottom()
ax1.yaxis.tick_left()

ax1.xaxis.set_major_formatter(LonFormatter)
ax1.yaxis.set_major_formatter(LatFormatter)

ax1.set_xticklabels(' ')

cs = ax1.contourf(Syn['lon'], Syn['lat'], crit1[:,:,dateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

# cbax = fig.add_axes([0.12, 0.35, 0.78, 0.02])
# cbar = fig.colorbar(cs, cax = cbax, orientation = 'horizontal')

ax1.set_extent([-105, -85, 30, 50])

ax1.set_ylabel('GFS', size = 18)


# Set the second part of the figure
ax2 = fig.add_subplot(1, 4, 2, projection = fig_proj)

ax2.set_title('Criteria 2', size = 18)

ax2.coastlines()
ax2.add_feature(cfeature.BORDERS)
ax2.add_feature(cfeature.STATES)

ax2.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax2.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax2.xaxis.tick_bottom()
ax2.yaxis.tick_left()

ax2.xaxis.set_major_formatter(LonFormatter)
ax2.yaxis.set_major_formatter(LatFormatter)

ax2.set_xticklabels(' ')
ax2.set_yticklabels(' ')

cs = ax2.contourf(Syn['lon'], Syn['lat'], crit2[:,:,dateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax2.set_extent([-105, -85, 30, 50])


# Set the third part of the figure
ax3 = fig.add_subplot(1, 4, 3, projection = fig_proj)

ax3.set_title('Criteria 3', size = 18)

ax3.coastlines()
ax3.add_feature(cfeature.BORDERS)
ax3.add_feature(cfeature.STATES)

ax3.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax3.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax3.set_yticklabels(LatLabel, fontsize = 14)

ax3.xaxis.tick_bottom()
ax3.yaxis.tick_left()

ax3.xaxis.set_major_formatter(LonFormatter)
ax3.yaxis.set_major_formatter(LatFormatter)

ax3.set_xticklabels(' ')
ax3.set_yticklabels(' ')

cs = ax3.contourf(Syn['lon'], Syn['lat'], crit3[:,:,dateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax3.set_extent([-105, -85, 30, 50])


# Set the fourth part of the figure
ax4 = fig.add_subplot(1, 4, 4, projection = fig_proj)

ax4.set_title('Criteria 4', size = 18)

ax4.coastlines()
ax4.add_feature(cfeature.BORDERS)
ax4.add_feature(cfeature.STATES)

ax4.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax4.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax4.xaxis.tick_bottom()
ax4.yaxis.tick_right()

ax4.set_yticklabels(LatLabel, fontsize = 14)

ax4.xaxis.set_major_formatter(LonFormatter)
ax4.yaxis.set_major_formatter(LatFormatter)

ax4.set_xticklabels(' ')

cs = ax4.contourf(Syn['lon'], Syn['lat'], crit4[:,:,dateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax4.set_extent([-105, -85, 30, 50])

savename = 'GFS_criteria_20120720.png'
plt.savefig('./Figures/' + savename, bbox_inches = 'tight')

plt.show(block = False)

#%%
  # The bottom half of the figure

PlotDate = datetime(2012, 7, 20)
narrdateind  = np.where(narr_time == PlotDate)[0]

fig = plt.figure(figsize = [12, 18])
### Begin the plots with the NARR data ###
# Begin the first part
ax1 = fig.add_subplot(1, 4, 1, projection = fig_proj)

ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.STATES)

ax1.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax1.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax1.set_yticklabels(LatLabel, fontsize = 14)
ax1.set_xticklabels(LonLabel, fontsize = 14)

ax1.xaxis.tick_bottom()
ax1.yaxis.tick_left()

ax1.xaxis.set_major_formatter(LonFormatter)
ax1.yaxis.set_major_formatter(LatFormatter)

cs = ax1.contourf(Syn['lon'], Syn['lat'], narrcrit1[:,:,narrdateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

# cbax = fig.add_axes([0.12, 0.35, 0.78, 0.02])
# cbar = fig.colorbar(cs, cax = cbax, orientation = 'horizontal')

ax1.set_extent([-105, -85, 30, 50])

ax1.set_ylabel('NARR', size = 18)

# Begin the second part
ax2 = fig.add_subplot(1, 4, 2, projection = fig_proj)

ax2.coastlines()
ax2.add_feature(cfeature.BORDERS)
ax2.add_feature(cfeature.STATES)

ax2.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax2.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax2.set_yticklabels(LatLabel, fontsize = 14)
ax2.set_xticklabels(LonLabel, fontsize = 14)

ax2.xaxis.tick_bottom()
ax2.yaxis.tick_left()

ax2.xaxis.set_major_formatter(LonFormatter)
ax2.yaxis.set_major_formatter(LatFormatter)

ax2.set_yticklabels(' ')

cs = ax2.contourf(Syn['lon'], Syn['lat'], narrcrit2[:,:,narrdateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

# cbax = fig.add_axes([0.12, 0.35, 0.78, 0.02])
# cbar = fig.colorbar(cs, cax = cbax, orientation = 'horizontal')

ax2.set_extent([-105, -85, 30, 50])


# Set the third part of the figure
ax3 = fig.add_subplot(1, 4, 3, projection = fig_proj)

ax3.coastlines()
ax3.add_feature(cfeature.BORDERS)
ax3.add_feature(cfeature.STATES)

ax3.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax3.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax3.set_xticklabels(LonLabel, fontsize = 14)

ax3.xaxis.tick_bottom()
ax3.yaxis.tick_left()

ax3.xaxis.set_major_formatter(LonFormatter)
ax3.yaxis.set_major_formatter(LatFormatter)

ax3.set_yticklabels(' ')

cs = ax3.contourf(Syn['lon'], Syn['lat'], narrcrit3[:,:,narrdateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax3.set_extent([-105, -85, 30, 50])


# Set the fourth part of the figure
ax4 = fig.add_subplot(1, 4, 4, projection = fig_proj)

ax4.coastlines()
ax4.add_feature(cfeature.BORDERS)
ax4.add_feature(cfeature.STATES)

ax4.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax4.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax4.xaxis.tick_bottom()
ax4.yaxis.tick_right()

ax4.set_yticklabels(LatLabel, fontsize = 14)
ax4.set_xticklabels(LonLabel, fontsize = 14)

ax4.xaxis.set_major_formatter(LonFormatter)
ax4.yaxis.set_major_formatter(LatFormatter)

cs = ax4.contourf(Syn['lon'], Syn['lat'], narrcrit4[:,:,narrdateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax4.set_extent([-105, -85, 30, 50])


savename = 'NARR_criteria_20120720.png'
plt.savefig('./Figures/' + savename, bbox_inches = 'tight')

plt.show(block = False)


#%%
  # Next, create a side by side figure of the flash drought identified.

fig = plt.figure(figsize = [12, 18])
fig.suptitle('Flash Drought Idenified ' + time[dateind[0]].strftime('%Y/%m/%d') + ' (0 hour forecast for GFS)', y = 0.66, size = 20)

### Begin the plots for the GFS data ###
# Set the GFS part of the figure
ax1 = fig.add_subplot(1, 2, 1, projection = fig_proj)

ax1.set_title('GFS', size = 18)

ax1.coastlines()
ax1.add_feature(cfeature.BORDERS)
ax1.add_feature(cfeature.STATES)

ax1.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax1.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax1.set_yticklabels(LatLabel, fontsize = 14)
ax1.set_xticklabels(LonLabel, fontsize = 14)

ax1.xaxis.tick_bottom()
ax1.yaxis.tick_left()

ax1.xaxis.set_major_formatter(LonFormatter)
ax1.yaxis.set_major_formatter(LatFormatter)

cs = ax1.contourf(Syn['lon'], Syn['lat'], FD[:,:,dateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax1.set_extent([-105, -85, 30, 50])


# Set the NARR part of the figure
ax2 = fig.add_subplot(1, 2, 2, projection = fig_proj)

ax2.set_title('NARR', size = 18)

ax2.coastlines()
ax2.add_feature(cfeature.BORDERS)
ax2.add_feature(cfeature.STATES)

ax2.set_xticks(LonLabel, crs = ccrs.PlateCarree())
ax2.set_yticks(LatLabel, crs = ccrs.PlateCarree())

ax2.set_yticklabels(LatLabel, fontsize = 14)
ax2.set_xticklabels(LonLabel, fontsize = 14)

ax2.xaxis.tick_bottom()
ax2.yaxis.tick_right()

ax2.xaxis.set_major_formatter(LonFormatter)
ax2.yaxis.set_major_formatter(LatFormatter)

cs = ax2.contourf(Syn['lon'], Syn['lat'], narrFD[:,:,narrdateind[0]].T, levels = clevs,
                 cmap = cmap, transform = data_proj, extend = 'max')

ax2.set_extent([-105, -85, 30, 50])

savename = 'Flash_Drought_Identified_20120720.png'
plt.savefig('./Figures/' + savename, bbox_inches = 'tight')

plt.show(block = False)



