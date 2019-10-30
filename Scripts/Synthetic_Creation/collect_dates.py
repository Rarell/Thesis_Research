#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:31:38 2019

@author: Rarrell
"""

#%% Import libraries
import sys, os, warnings
import numpy as np
import pygrib
from datetime import datetime, timedelta

#%% Main function
def main():
    '''
    '''
    
    Years, Months, Days, ModelRuns = CollectDates()
    
    Sources = DetermineSources(Years, Months, Days, ModelRuns)
    
    AppendTMP(Years, Months, Days, ModelRuns, Sources)


#%% 
##########################################
### Error Class for Inpurt Errors. #######
##########################################
class InputError(Exception):
    """This is raised when there is an error in the user input"""
    pass

#%% DateRange function
def DateRange(StartDate, EndDate):
    for n in range(int((EndDate - StartDate).days) + 1):
        yield StartDate + timedelta(n)


#%% Collect dates
def CollectDates():
    '''
    '''
    
    StartYear = input('Enter the starting year of dataset (between 2005 and present):  ')
    
    # Check that the year is inside the valid year range (2005 - present)
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
    
    
    EndYear = input('Enter the ending year of dataset (between 2005 and present):  ')
    
    # Check that the year is inside the valid year range (2005 - present)
    try:
        if ( (datetime(int(EndYear), 1, 1) >= datetime(2005, 1, 1)) & 
            (datetime(int(EndYear), 1, 1) <= datetime.now()) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The year must be between 2005 and the present year.')
    
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
        
    StartDate = datetime(int(StartYear), int(StartMonth), int(StartDay))
    EndDate = datetime(int(EndYear), int(EndMonth), int(EndDay))
    
    MinChange = timedelta(days = 30)
    TimeDelta = EndDate - StartDate
    
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
        
    Years     = ['tmp'] * (4*int(TimeDelta.days)+1)
    Months    = ['tmp'] * (4*int(TimeDelta.days)+1)
    Days      = ['tmp'] * (4*int(TimeDelta.days)+1)
    ModelRuns = ['tmp'] * (4*int(TimeDelta.days)+1)
    
    n = 0
    m = 0
    for dt in DateRange(StartDate, EndDate):

        for run in ['00', '06', '12', '18']:
            if ((int(ModelRun) > int(run)) & (n == 0)):
                continue
            elif ((datetime(dt.year, dt.month, dt.day) == EndDate) &\
                  (int(ModelRun) < int(run))):
                continue
            else:
                pass
            
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
            
    return Years, Months, Days, ModelRuns

#%% Determine sources
def DetermineSources(Years, Months, Days, ModelRuns):
    '''
    '''
    Request_page = 'https://www.ncdc.noaa.gov/has/HAS.FileAppRouter?' +\
                   'datasetname=GFS3&subqueryby=STATION&applname=&outdest=FILE'
    NCEI_source  = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
    NCEP_source  = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
    
    DeltaWeek = timedelta(days = 7)
    DeltaYear = timedelta(days = 365)
    
    Sources = ['tmp'] * len(Years)
    
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
    
    return Sources

#%% Append the tmp_list.txt file
def AppendTMP(Years, Months, Days, ModelRuns, Sources):
    '''
    '''
    path = './'
    f = open(path + 'tmp_list.txt', 'a')
    
    for n in range(len(Years)):
        f.write(Years[n] + ' ' + Months[n] + ' ' + Days[n] + ' ' +\
                ModelRuns[n] + ' ' + Sources[n] + '\n')
        
    f.close()
    
#%% Call the main function
if __name__ == '__main__':
    main()