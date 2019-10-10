#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 16:09:42 2019

@author: Rarrell

This is a script designed to take the temporary .nc files created and merge
them into a single .nc file over all forecast hours.
"""

#%% Load modules
import os, sys, warnings
import numpy as np
from netCDF4 import Dataset
from glob import glob

#%% Define a function to load in netcdf files.
def load_nc(VarSName, VarFH):
    path = '/Users/Rarrell/Desktop/Thesis_Research/tmp/'
    filename = 'GFS_' + VarSName + '_' +\
                   str(VarFH) + '.nc'
    with Dataset(path + filename, 'r') as nc:
        if str(VarFH) == '003':
            lat = nc.variables['lat'][:]
            lon = nc.variables['lon'][:]
            
            VD  = nc.variables['VD'][:]
            
            mask = nc.variables['mask'][:,:]
            
            units = nc.variables['units'][:]
            
        else:
            continue
        
    FH = nc.variables['FH'][:]
    
    Var = nc.variables[VarSName][:,:]
    
    return Var, FH, lat, lon, VD, mask, units if str(VarFH) == '003' else return Var, FH