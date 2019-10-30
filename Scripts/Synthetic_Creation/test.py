#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 13:05:37 2019

@author: Rarrell
"""

#%% Import libraries
import sys, os, warnings
import numpy as np
import pygrib
from datetime import datetime, timedelta


#%% Test the creation of multiple days and years

YearStart = input("Enter start year  ")
MonthStart = input("Enter start month  ")
DayStart = input("Enter start day  ")


YearEnd = input("Enter end year  ")
MonthEnd = input("Enter the end month  ")
DayEnd = input("Enter the end day  ")

MinChange = timedelta(days = 30)

DateStart = datetime(int(YearStart), int(MonthStart), int(DayStart))
DateEnd = datetime(int(YearEnd), int(MonthEnd), int(DayEnd))

TimeDelta = DateEnd - DateStart

def daterange(StartDate, EndDate):
    for n in range(int((EndDate - StartDate).days) + 1):
        yield StartDate + timedelta(n)

#%%
for dt in daterange(DateStart, DateEnd):
    year  = dt.year
    month = dt.month
    day   = dt.day
    
    year  = str(year)
    month = str(month)
    day   = str(day)
    
    print(year, month, day)
    

