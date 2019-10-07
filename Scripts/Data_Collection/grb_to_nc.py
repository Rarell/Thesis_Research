#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 15:48:20 2019

@author: Rarrell

This is a file designed to take in various inputs to create a url to download
a GFS grib file from the NCEI or NCEP, extract the desired variable at the 
desired height, and dump it into a temporary netcdf file.
"""

#%% Import needed modules
import sys, os, warnings
import numpy as np
import pygrib
from glob import glob


#%% Define the function that will perform the work.
def grb_to_nc(var, TypeOfHeight, height = None, request = False, year = 2019, 
              month = 06, day = 25, model_run = 00, forecast_hour = 00, 
              url_start = ''):
    print('This is the start of the program')
    
#%% Create a testing point to test parts of the above function.
    