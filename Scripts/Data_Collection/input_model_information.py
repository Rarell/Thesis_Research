#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:15:07 2019

@author: Rarrell
"""

#%% Import some libraries
import sys, os, warnings
import numpy as np
import pygrib
from datetime import datetime, timedelta

#%% Create the main function
def main():
    year, month, day, ModelRun = collect_date()
    
    source, request = determine_source(year, month, day, ModelRun)
    
    Var, TypeOfHeight, height = input_model_variables()
    
    write_tmp(year, month, day, source, request, Var, TypeOfHeight, height,
              ModelRun)


#%% Create an error class for inpurt errors.
class InputError(Exception):
    """This is raised when there is an error in the user input"""
    pass

#%% Create a function to count the characters in a string.
def CharCount(string):
    count = 0
    for c in string:
        count = count + 1
    
    return count
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
        
    month_count = CharCount(month)
    if (int(month) < 10) & (month_count < 2):
        month = '0' + month
    else:
        pass
    
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
        
    day_count = CharCount(day)
    if (int(day) < 10) & (day_count < 2):
        day = '0' + day
    else:
        pass
        
    # Input the user for the model run
    ModelRun = input('Input the desired model (either 00, 06, 12, or 18):  ')
    # Check the input
    try:
        if ((ModelRun == '00') | (ModelRun == '06') | (ModelRun == '12') | 
                (ModelRun == '18')):
            pass
        else:
            raise InputError('There only 00Z, 06Z, 12Z, and 18Z model runs. ' +\
                            'Please try again with 00, 06, 12, or 18.')
    
    except InputError as ie:
        raise InputError(ie)
              
                
    
    return year, month, day, ModelRun

#%% Create a function to determine the data source (NCEI or NCEP)
def determine_source(year, month, day, ModelRun):
    
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
              '.g2 file is downloaded and unzipped in the Data ' +\
              'directory.')
        
        request = 1
        source = None
        
        try:
            filename = 'gfs_3_' + str(year) + str(month) + str(day) + '_' +\
                       '00' + str(ModelRun) + '_003.grb2'
            path = './Data/gfs_3_' + str(year) + str(month) + str(day) +\
                   str(ModelRun) + '.g2/'
#           path = '/Users/Rarrell/Downloads/gfs_3_' + str(year) + str(month) +\
#                  str(day) + str(model_run) + '.g2/'
            grb_file = path + filename
            
            grb = pygrib.open(grb_file)
            grb.close()
        except OSError:
            raise OSError('.grb2 folder not found in the Data folder.' +\
                          'Please make to make the request at: \n' +\
                          Request_page + '\n' +\
                          'and download the data to the Data folder.')
        
        # Use try here to see if the file exists. Prompt the user to make the
        # data request if not and exit the program.
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
def input_model_variables():
    
    # Collect the valid input options
    path = '/Users/Rarrell/Desktop/Thesis_Research/Scripts/Data_Collection/txt_files/'
    ValidVar = np.loadtxt(path + 'GFS_VarNames.txt', dtype = str, delimiter = ',')
    ValidTypeOfHeights = np.loadtxt(path + 'GFS_TypeOfHeights.txt', dtype = str,
                                    delimiter = ',')
    ValidHeights = np.loadtxt(path + 'GFS_Heights.txt', dtype = str, delimiter = ',')
    
    
    # Prompt for the variable name
    Var = input('Please enter the variable to be extract. ' +\
                'See GFS_Variable_Names.txt for a list of options. \n'  )
    # Check the input is valid
    try:
        n = 0
        for ValVar in ValidVar:
            if Var == ValVar:
                n = n + 1
            else: 
                n = n
        
        if n == 1:
            pass
        else:
            raise InputError(Var + ' is not a valid name. ' +\
                             'See GFS_Variable_Names.txt for a list of valid ' +\
                             'options.')
    
    except InputError as ie:
        raise InputError(ie)


    # Prompt the user for the type of height
    TypeOfHeight = input('Please enter the type of height the variable is located' +\
                         'at. See GFS_Type_Of_Heights.txt for a list of options. \n'  )
    # Check the input is valid
    try:
        n = 0
        for ValToH in ValidTypeOfHeights:
            if TypeOfHeight == ValToH:
                n = n + 1
            else:
                n = n
        
        if n == 1:
            pass
        else:
            raise InputError(TypeOfHeight + ' is not a valid type of height. ' +\
                             'See GFS_Type_Of_Heights.txt for a list of valid ' +\
                             'options.')
    
    except InputError as ie:
        raise InputError(ie)
    
    # Prompt the user for the height
    height = input('Type in the height the variable is located at. The units ' +\
                   '(e.g. Pa, km, PVU, etc.) is determined by the type of height.' +\
                   'If the type of height is the surface, then the height is 0. \n'  )
    # Check the input is valid
    try:
        if (TypeOfHeight == 'surface') & (height != '0'):
            raise InputError('The height must be 0 if the TypeOfHeight is surface.')
        else:
            pass  
        
        n = 0
        for ValHeight in ValidHeights:
            if (height + '.0') == ValHeight:
                n = n + 1
            else:
                n = n
        
        if n == 1:
            pass
        else:
            raise InputError(height + ' is not a valid height. See ' +\
                             'GFS_Height_Values.txt for a list of valid options.')
    
    except InputError as ie:
        raise InputError(ie)
    
    return Var, TypeOfHeight, height

#%% Define a function to append the tmp file with the model information.
def write_tmp(year, month, day, source, request, Var, TypeOfHeight, height,
              ModelRun):
#    path = '/Users/Rarrell/Desktop/Thesis_Research/'
    path = './'
    f = open(path + 'tmp.txt', 'a')
    
    f.write(str(request) + ',')
    f.write(year + ',')
    f.write(month + ',')
    f.write(day + ',')
    f.write(str(source) + ',')
    f.write(Var + ',')
    f.write(TypeOfHeight + ',')
    f.write(height + ',')
    f.write(ModelRun)
    
    f.close()


#%% Call the main function
if __name__ == '__main__':
    main()
