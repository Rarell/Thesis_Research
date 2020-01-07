#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 15:30:53 2019

@author: stuartedris
"""

#%%
#####################################
### Import some libraries ###########
#####################################

import sys, os, warnings
import numpy as np
from datetime import datetime


#%%
  # Main function
def main():
    '''
    '''
    
    # Collect the arguments
    script = sys.argv[0]
    Parameters = sys.argv[1]
    
    # Unpack the date information
    Year, Month, Day = np.loadtxt(Parameters, usecols = (1,2,3), dtype = str, 
                                  delimiter = ',', unpack = True)
    
    # Initialize the cutoff dates for when the GFS changes grids:
    #   Pre 2019 (assumed date, not truly accurate), GFS grid changes after hour 192 (8 days)
    #   Pre July 2019, GFS grid changes after hour 240 (10 days)
    #   Post July, 2019, the GFS grid does not change
    EightDayDate = datetime(2019, 1, 1)
    TenDayDate   = datetime(2019, 7, 1)
    
    # Determine the date of this GFS run
    CurrentDate = datetime(int(Year), int(Month), int(Day))
    
    # Create the forecast hours dataset depending on date
    if CurrentDate < EightDayDate:
        GridChangeHours = 192
        FH = ObtainForecastHours(GridChangeHours)
    elif (CurrentDate < TenDayDate) & (CurrentDate >= EightDayDate):
        GridChangeHours = 240
        FH = ObtainForecastHours(GridChangeHours)
    else:
        GridChangeHours = 384
        FH = ObtainForecastHours(GridChangeHours)
        
    # Write the forecast hours dataset to iterable and comma seperate files
    write_txt(FH)

#%%
  # Function for writing the forecast hour variables
def ObtainForecastHours(EndPoint):
    '''
    
    '''
    
    # Determine the length of the forecast hour variable (3 hour resolutioin)
    VarLen = EndPoint/3
    
    # Initialize the forecast hour variable
    FH = ['tmp'] * int(VarLen)
    n = 0
    
    # Create the forecast hour strings for values below EndPoint hours
    for i in range(3, EndPoint+3, 3):
        if i < 10:
            FH[n] = '00' + str(i)
        elif (i < 100) & (i > 10):
            FH[n] = '0' + str(i)
        else:
            FH[n] = str(i)
        n = n + 1
        
    return FH

#%%
  # Function to write the .txt files

def write_txt(FH):
    '''
    
    '''
    
    # Define the path and filename for the iteratable forecast hours
    path = './Scripts/Data_Collection/'
    filename = 'GFS_Forecast_Times.txt'
    
    # Write the iterable forecast hours
    f = open(path + filename, 'w')
    for ForecastHour in FH:
        f.write(ForecastHour + '\n')
        
    f.close()
    
    # Define the path and filename for the comma seperated forecast hours
    path = './Scripts/Data_Collection/txt_files/'
    filename = 'GFS_FH.txt'
    
    # Write the comma seperate forecast hours
    f = open(path + filename, 'w')
    for ForecastHour in FH:
        if ForecastHour == FH[-1]:
            f.write(ForecastHour)
        else:
            f.write(ForecastHour + ',')
        
    f.close()


#%%
#########################################
### Call and Run the Main Function ######
#########################################
    
if __name__ == '__main__':
    main()
    
#####################
### End of Script ###
#####################