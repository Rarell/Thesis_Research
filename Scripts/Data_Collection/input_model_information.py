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
    return something