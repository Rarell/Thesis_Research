#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:15:07 2019

@author: Rarrell
"""

#%% Import some libraries
import sys, os, warnings
import numpy as np
from datetime import datetime, timedelta


#%% Create an error class for inpurt errors.
class InputError(Exception):
    """This is raised when there is an error in the user input"""
    pass

#%% Define a function to collect the date
def collect_date():
    # Prompt the user for the desired year.
    year = input('Enter the year of desired model run (between 2005 and present):  ')
    # Check the input value
    try:
        int(year)
    except (InputError, ValueError):
        raise InputError('The year must be a numeric integer.')
            
    try:
        if ( (datetime(int(year), 1, 1) >= datetime(2005, 1, 1)) & 
            (datetime(int(year), 1, 1) <= datetime.now()) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The year must be between 2005 and the present year.')


    # Prompt the user for the desired month.
    month = input('Please enter the month of the desired model run ' +\
                  '(between 1 and 12):  ')
    # Check the input value
    try:
        int(month)
    except (InputError, ValueError):
        raise InputError('The month must be a numeric integer.')
    
    try:
        if ( (int(month) >= 1) & (int(month) <= 12) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The month must be between 1 and 12.')
    
    # Prompt the user for the desired day
    day = input('Enter the day of the model run (between 1 and 31):  ')
    # Check the input value
    try:
        int(day)
    except (InputError, ValueError):
        raise InputError('The day must be a numeric whole number.')
    
    try:
        if int(day) < 1:
            raise InputError('The day must be a positive whole number.')
        else:
            pass
    
        if int(month) == 2:
            if int(day) > 28:
                raise InputError('There are no more than 28 days in February.')
            else:
                pass
        elif ( (int(month) == 4) | (int(month) == 6) | (int(month) == 9) |
                (int(month) == 11) ):
            if int(day) > 30:
                raise InputError('There are no more than 30 days in' +\
                                 'April, June, September, and November.')
            else:
                pass
        else:
            if int(day) > 31:
                raise InputError('There are no more than 31 days in Janruary, ' +\
                                 'March, May, July, August, October, and December.')
            else:
                pass
    except InputError as ie:
        raise InputError(ie)
    
    return year, month, day

#%% Create a function to determine the data source (NCEI or NCEP)
def determine_source(year, month, day):
    
    # Initialize some check variables.
    Request_page = 'https://www.ncdc.noaa.gov/has/HAS.FileAppRouter?' +\
                   'datasetname=GFS3&subqueryby=STATION&applname=&outdest=FILE'
    NCEI_source  = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
    NCEP_source  = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
    
    DeltaYear = timedelta(days = 365)
    DeltaWeek = timedelta(days = 7)
    
    date = datetime(int(year), int(month), int(day))
    
    DeltaTime = datetime.now() - date
    
    # Check if a data request is needed
    if DeltaTime > DeltaYear:
        # Ensure the user knows to make a request for old data.
        print('The requested date is over 1 year old. A data request is ' +\
              'required for this. This can be made at: \n' +\
              Request_page + '\n' +\
              'Please ensure the request is made, the ' +\
              '.g2 file is downloaded and unzipped in the Downloads ' +\
              'directory. Failing to do so can result in an error and ' +\
              'termination of the program.')
        
        request = 1
        source = None
    else:
        request = 0
        
        # Check if the data should come from the NCEI or NCEP
        if DeltaTime > DeltaWeek:
            print('Data will be downloaded from the NCEI database. This is ' +\
                  'located at: \n' + NCEI_source)
            source = 'NCEI'
        else:
            print('New data will be downloaded from the NCEP database. ' +\
                  'This is located at: \n' + NCEP_source)
            source = 'NCEP'

    return source, request

#%% Create a new function to prompt for variable information.

# Prompt for the variable name
Var = input('Please enter the variable to be extract. ' +\
            'See GFS_Variable_Names.txt for a list of options.'  )
# Check the input is valid


TypeOfHeight = input('Please enter the type of height the variable is located' +\
                     'at. See GFS_Type_Of_Heights.txt for a list of options.'  )
# Check the input is valid

# Prompt the user for the height
height = input('Type in the height the variable is located at. The units ' +\
               '(e.g. Pa, km, PVU, etc.) is determined by the type of height.' +\
               'If the type of height is the surface, then the height is 0.'  )
# Check the input is valid





