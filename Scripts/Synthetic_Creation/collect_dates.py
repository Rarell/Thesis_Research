#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:31:38 2019

@author: Rarrell

This script is designed as part of the create_syntetic_dataset.sh program.
This script will prompt the user for a start and end date of a range of dates
  that is a minimum of 30 days. Inputs are checked for validity, and then the
  source of the GFS data is determined. The information is then appended to an
  empty temporary file in the main program.
"""

#%% 
#####################################
### Import some libraries ###########
#####################################

import sys, os, warnings
import numpy as np
import pygrib
from datetime import datetime, timedelta

#%% 
############################
### Main Function ##########
############################

def main():
    '''
    This is the main function of the script. It prompts the user for the start
      and end dates, determines the source of the data, and appends it to a
      temporary file.
    '''
    
    # Prompt the user for the start and end dates of the range and check the
    #   validity of the inputs and output all dates for all model runs
    Years, Months, Days, ModelRuns = CollectDates()
    
    # Determine what source the data is coming from
    Sources = DetermineSources(Years, Months, Days, ModelRuns)
    
    # Append all the information into a temporary file.
    AppendTMP(Years, Months, Days, ModelRuns, Sources)


#%% 
##########################################
### Error Class for Inpurt Errors. #######
##########################################
class InputError(Exception):
    """This is raised when there is an error in the user input"""
    pass

#%%
#############################
### DateRange Function ######
#############################
    
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
################################
### CollectDates Function ######
################################
def CollectDates():
    '''
    This function prompts the user for the start date and end date of 
      date range, minimum of 30 days, and generates all dates between 
      (inclusive) for all model runs (00Z, 06Z, 12Z, and 18Z).
      
    Outputs:
        Years     - List all year(s) for every date entry
        Months    - List of all months in the date interval
        Days      - List all days in the date interval
        ModelRuns - List of all model runs for the date interval.
    '''
    
    #####################################################
    # Prompt the user for the year of the starting date #
    #####################################################
    StartYear = input('Enter the starting year of dataset (between 2005 and present):  ')
    
    # Check that the year is a number
    try:
        int(StartYear)
    except (InputError, ValueError):
        raise InputError('The year must be a numeric integer.')
    
    # Check that the year is inside the valid year range (2005 - present)       
    try:
        if ( (datetime(int(StartYear), 1, 1) >= datetime(2005, 1, 1)) & 
            (datetime(int(StartYear), 1, 1) <= datetime.now()) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The year must be between 2005 and the present year.')
    
    ######################################################
    # Prompt the user for the month of the starting date #
    ######################################################
    StartMonth = input('Enter the starting month of the dataset ' +\
                  '(between 1 and 12):  ')
    
    # Check the input month is an integer number
    try:
        int(StartMonth)
    except (InputError, ValueError):
        raise InputError('The month must be a numeric integer.')
    
    # Check that the month is between 1 and 12
    try:
        if ( (int(StartMonth) >= 1) & (int(StartMonth) <= 12) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The month must be between 1 and 12.')
    
    ####################################################
    # Prompt the user for the day of the starting date #
    ####################################################
    StartDay = input('Enter the starting day of the dataset (between 1 and 31):  ')
    
    # Check that the input day is an integer number
    try:
        int(StartDay)
    except (InputError, ValueError):
        raise InputError('The day must be a numeric whole number.')
    
    # Check that the day is positive, and within a valid range
    #  (1 to 28 for February, 1 to 30 for appropriate months, and 1 to 31
    #  for the remaining months)
    try:
        if int(StartDay) < 1:
            raise InputError('The day must be a positive whole number.')
        else:
            pass
    
        if int(StartMonth) == 2:
            if int(StartDay) > 28:
                raise InputError('There are no more than 28 days in February.')
            else:
                pass
        elif ( (int(StartMonth) == 4) | (int(StartMonth) == 6) | 
                (int(StartMonth) == 9) | (int(StartMonth) == 11) ):
            if int(StartDay) > 30:
                raise InputError('There are no more than 30 days in' +\
                                 'April, June, September, and November.')
            else:
                pass
        else:
            if int(StartDay) > 31:
                raise InputError('There are no more than 31 days in Janruary, ' +\
                                 'March, May, July, August, October, and December.')
            else:
                pass
    except InputError as ie:
        raise InputError(ie)
    
    
    ###################################################
    # Prompt the user for the year of the ending date #
    ###################################################
    EndYear = input('Enter the ending year of dataset (between 2005 and present):  ')
    
    # Check that the year is a number
    try:
        int(EndYear)
    except (InputError, ValueError):
        raise InputError('The year must be a numeric integer.')
    
    # Check that the year is inside the valid year range (2005 - present)
    try:
        if ( (datetime(int(EndYear), 1, 1) >= datetime(2005, 1, 1)) & 
            (datetime(int(EndYear), 1, 1) <= datetime.now()) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The year must be between 2005 and the present year.')
    
    ####################################################
    # Prompt the user for the month of the ending date #
    ####################################################
    EndMonth = input('Enter the ending month of the dataset ' +\
                  '(between 1 and 12):  ')
    
    # Check the input month is an integer number 
    try:
        int(EndMonth)
    except (InputError, ValueError):
        raise InputError('The month must be a numeric integer.')
    
    # Check that the month is between 1 and 12
    try:
        if ( (int(EndMonth) >= 1) & (int(EndMonth) <= 12) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The month must be between 1 and 12.')
    
    ##################################################
    # Prompt the user for the day of the ending date #
    ##################################################
    EndDay = input('Enter the ending day of the dataset (between 1 and 31):  ')
    
    # Check that the input day is an integer number
    try:
        int(EndDay)
    except (InputError, ValueError):
        raise InputError('The day must be a numeric whole number.')
    
    # Check that the day is positive, and within a valid range
    #  (1 to 28 for February, 1 to 30 for appropriate months, and 1 to 31
    #  for the remaining months)
    try:
        if int(EndDay) < 1:
            raise InputError('The day must be a positive whole number.')
        else:
            pass
    
        if int(EndMonth) == 2:
            if int(EndDay) > 28:
                raise InputError('There are no more than 28 days in February.')
            else:
                pass
        elif ( (int(EndMonth) == 4) | (int(EndMonth) == 6) | 
                (int(EndMonth) == 9) | (int(EndMonth) == 11) ):
            if int(EndDay) > 30:
                raise InputError('There are no more than 30 days in' +\
                                 'April, June, September, and November.')
            else:
                pass
        else:
            if int(EndDay) > 31:
                raise InputError('There are no more than 31 days in Janruary, ' +\
                                 'March, May, July, August, October, and December.')
            else:
                pass
    except InputError as ie:
        raise InputError(ie)
    
    # Determine how many days lies between the start and end dates
    StartDate = datetime(int(StartYear), int(StartMonth), int(StartDay))
    EndDate = datetime(int(EndYear), int(EndMonth), int(EndDay))
    
    MinChange = timedelta(days = 30)
    TimeDelta = EndDate - StartDate
    
    # Check that the start and end dates are at least 30 days apart and the 
    #   end date is after the start dates
    try:
        if StartDate > EndDate:
            raise InputError('The starting date must be before the ending date.')
        else:
            pass    

        if MinChange > TimeDelta:
            raise InputError('The start and end dates must be at least 30 days apart.')
        else:
            pass
    except InputError as ie:
        raise InputError(ie)
    
    # Prompt the user for the starting/ending model run (e.g., if the start/last 
    #   should be 00Z, 12Z, etc.)
    ModelRun = input('Select a model run to start at (00, 06, 12, or 18)  ')
    
    # Check that the input is one of the four valid inputs
    try:
        if ((ModelRun == '00') | (ModelRun == '06') | (ModelRun == '12') | 
                (ModelRun == '18')):
            pass
        else:
            raise InputError('There only 00Z, 06Z, 12Z, and 18Z model runs. ' +\
                            'Please try again with 00, 06, 12, or 18.')
    
    except InputError as ie:
        raise InputError(ie)
    
    # Initialize some inputs
    Years     = ['tmp'] * (4*int(TimeDelta.days)+1)
    Months    = ['tmp'] * (4*int(TimeDelta.days)+1)
    Days      = ['tmp'] * (4*int(TimeDelta.days)+1)
    ModelRuns = ['tmp'] * (4*int(TimeDelta.days)+1)
    
    n = 0
    m = 0
    
    # Determine all the dates between the start and end dates (inclusive).
    #   Each date should be repeated four times (one for each model run), save
    #   the start and end dates, which start/end on a specific model run
    for dt in DateRange(StartDate, EndDate):
        # Add an entry on each date for each model run
        for run in ['00', '06', '12', '18']:
            # Ensure it starts and ends on the specified on model run
            if ((int(ModelRun) > int(run)) & (n == 0)):
                continue
            elif ((datetime(dt.year, dt.month, dt.day) == EndDate) &\
                  (int(ModelRun) < int(run))):
                continue
            else:
                pass
            
            # Reset the model run counter for each new day (starts at 00Z)
            if run == '00':
                m = 0
            
            Years[n+m]  = str(dt.year)
        
            if dt.month < 10:
                Months[n+m] = '0' + str(dt.month)
            else:
                Months[n+m] = str(dt.month)
            
            if dt.day < 10:
                Days[n+m] = '0' + str(dt.day)
            else:
                Days[n+m]   = str(dt.day)
            
            ModelRuns[n+m] = run
            m = m + 1

        n = n + m
    
    ###################
    # End of Function #
    ###################
    
    return Years, Months, Days, ModelRuns

#%% 
####################################
### DetermineSources Function ######
####################################
    
def DetermineSources(Years, Months, Days, ModelRuns):
    '''
    This function determines the source of the GFS data (data request needed,
      NCEI, or NCEP) given a set of dates.
    
    Inputs:
        Years     - List. List of year(s) for each date whose source is to be determined
        Months    - List. List of month(s) for each date whose souce is to be determined
        Days      - List. List of days for each date whose source is to be determined
        ModelRuns - List. List of model runs corresponding to each date
        
    Outputs:
        Sources - List. List of sources of GFS data corresponding to each input date  
    '''
    
    # Create some initial values
    Request_page = 'https://www.ncdc.noaa.gov/has/HAS.FileAppRouter?' +\
                   'datasetname=GFS3&subqueryby=STATION&applname=&outdest=FILE'
    NCEI_source  = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
    NCEP_source  = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
    
    DeltaWeek = timedelta(days = 7)
    DeltaYear = timedelta(days = 365)
    
    Sources = ['tmp'] * len(Years)
    
    # Inform the user where the data will be coming from
    print('If any dates are older than a year, then a data request is needed.' +\
          'The data request can be made at: \n' +\
          Request_page + '\n' +\
          'Please make the request for all dates older than a year and for ' +\
          'each model run, and make sure the file is located in the Data folder.\n')
    
    print('If any of the dates are older than a week, but newer than the a ' +\
          'than a year, the data will be downloaded from the NCEI from: \n' +\
          NCEI_source + '.\n')
    
    print('If any of the dates are newer than a week, it will be downloaded ' +\
          'from the NCEP page at: \n' +\
          NCEP_source + '.\n')
    
    # Determine where the GFS data is coming from for each given date
    for n in range(len(Years)):
        Date = datetime(int(Years[n]), int(Months[n]), int(Days[n]))
        TimeDelta = datetime.now() - Date
        
        if TimeDelta > DeltaYear:
            Sources[n] = 'Request'
            
            # Check that the .g2 folder exists in the Data folder.
            try:
                filename = 'gfs_3_' + str(Years[n]) + str(Months[n]) +\
                           str(Days[n]) + '_' + '00' + str(ModelRuns[n]) + '_003.grb2'
                path = './Data/gfs_3_' + str(Years[n]) + str(Months[n]) +\
                       str(Days[n]) + str(ModelRuns[n]) + '.g2/'
#               path = '/Users/Rarrell/Downloads/gfs_3_' + str(year) + str(month) +\
#                      str(day) + str(model_run) + '.g2/'
                grb_file = path + filename
                
                grb = pygrib.open(grb_file)
                grb.close()
            except OSError:
                raise OSError('.grb2 folder not found in the Data folder.' +\
                              'Please make to make the request at: \n' +\
                              Request_page + '\n' +\
                              'and download the data to the Data folder.')
            
        elif ((TimeDelta < DeltaYear) & (TimeDelta > DeltaWeek)):
            Sources[n] = 'NCEI'
        else:
            Sources[n] = 'NCEP'
    
    ###################
    # End of Function #
    ###################
    return Sources

#%%
#############################
### AppendTMP Function ######
#############################
    
def AppendTMP(Years, Months, Days, ModelRuns, Sources):
    '''
    This function appends the information from the previous two functions onto
      a temporary file.
      
    Inputs:
        Years     - List. List of year(s) for each date 
        Months    - List. List of month(s) for each date 
        Days      - List. List of days for each date 
        ModelRuns - List. List of model runs corresponding to each date
        Sources   - List. List of sources corresponding to each date
    '''
    
    # Define the location of the temporary file
    path = './'
    
    # Open the file for appending
    f = open(path + 'tmp_list.txt', 'a')
    
    # Write the information to it
    for n in range(len(Years)):
        f.write(Years[n] + ' ' + Months[n] + ' ' + Days[n] + ' ' +\
                ModelRuns[n] + ' ' + Sources[n] + '\n')
    
    # Close the function
    f.close()
    
    ###################
    # End of Function #
    ###################
    
#%%
#########################################
### Call and Run the Main Function ######
#########################################
    
if __name__ == '__main__':
    main()
    
#####################
### End of Script ###
#####################