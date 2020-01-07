#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 14 14:15:07 2019

@author: Rarrell

This is a script designed as a component of the extract_GFS_variable.sh program.
This script prompts the user for information on the GFS file and variable to be
  extracted. It will ask the user for the date, determine where the data comes 
  from depending on the date, and determine if .g2 file exists in the Data 
  folder if it a data request is needed
This script then prompts the user for the variable name to be extracted and
  what height it is located at.
This script then concludes by appending this information to a temporary .txt
  file for later use in the extract_GFS_variable.sh
"""

#%% 
#####################################
### Import some libraries ###########
#####################################
import sys, os, warnings
import numpy as np
import pygrib
import wget
import tarfile
from datetime import datetime, timedelta

#%% 
############################
### Main Function ##########
############################

def main():
    '''
    This is the main function that calls the others to obtain the GFS
      information write to a temporary .txt file.
    '''
    
    script = sys.argv[0]
    parameter = sys.argv[1]
    
    # Collect the date and model run of the GFS data being examined.
    year, month, day, ModelRun = collect_date()
    
    # Determine where the data comes from (NCEI or NCEP). If a data request
    #   is needed, determine if the .g2 folder exists in the Data folder.
    source, request = determine_source(year, month, day, ModelRun)
    
    if parameter == '-v':
        # Collect variable information of the variable to be extracted.
        Var, TypeOfHeight, height = input_model_variables()
    
        # Write the information into the temporary .txt file
        write_tmp(year, month, day, source, request, ModelRun,
                  Var, TypeOfHeight, height)
    else:
        write_tmp(year, month, day, source, request,  ModelRun)


#%% 
##########################################
### Error Class for Inpurt Errors. #######
##########################################
class InputError(Exception):
    """This is raised when there is an error in the user input"""
    pass

#%% 
###############################
### CharCount Function ########
###############################
    
def CharCount(string):
    '''
    Function designed to count the characters in a given string.
    
    Input:
        string - the string whose characters will be counted
    
    Output:
        count - the number of characters in string
    '''
    
    # Initialize count
    count = 0
    
    # Increment count by 1 for every character in string
    for c in string:
        count = count + 1
    
    return count
#%% 
################################
### collect_date Function ######
################################
    
def collect_date():
    '''
    This function is designed to prompt the user for the date and model run of
      of the GFS file to be examined. The inputs are also checked to ensure
      they are valid. If not, an InputError is raised.
      
    Outputs:
        year  - year of the GFS file being examined. Given by user.
        month - month of the GFS file being examined. Given by user.
        day   - day of the GFS file being examined. Given by the user.
        ModelRun - the model run (00Z, 06Z, 12Z, or 18Z) of the GFS file being
                     examined. Given by the user.
    '''
    
    ##############################################
    #### Prompt the user for the desired year. ###
    ##############################################
    year = input('Enter the year of desired model run (between 2005 and present):  ')
    
    # Check that the year is an integer number
    try:
        int(year)
    except (InputError, ValueError):
        raise InputError('The year must be a numeric integer.')
    
    # Check that the year is inside the valid year range (2005 - present)       
    try:
        if ( (datetime(int(year), 1, 1) >= datetime(2005, 1, 1)) & 
            (datetime(int(year), 1, 1) <= datetime.now()) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The year must be between 2005 and the present year.')

    ##############################################
    ### Prompt the user for the desired month. ###
    ##############################################
    month = input('Please enter the month of the desired model run ' +\
                  '(between 1 and 12):  ')
    
    # Check the input month is an integer number
    try:
        int(month)
    except (InputError, ValueError):
        raise InputError('The month must be a numeric integer.')
    
    # Check that the month is between 1 and 12
    try:
        if ( (int(month) >= 1) & (int(month) <= 12) ):
            pass
        else:
            raise InputError
    except InputError:
        raise InputError('The month must be between 1 and 12.')
    
    # Check if month values below 10 are entered as X or 0X
    # If the former, convert it to 0X.
    month_count = CharCount(month)
    if (int(month) < 10) & (month_count < 2):
        month = '0' + month
    else:
        pass
    
    ###########################################
    ### Prompt the user for the desired day ###
    ###########################################
    day = input('Enter the day of the model run (between 1 and 31):  ')
    
    # Check that the input day is an integer number
    try:
        int(day)
    except (InputError, ValueError):
        raise InputError('The day must be a numeric whole number.')
    
    # Check that the day is positive, and within a valid range
    #  (1 to 28 for February, 1 to 30 for appropriate months, and 1 to 31
    #  for the remaining months)
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
    
    # Check if the day input is enter as X or 0X.
    # If it is the former, convert to 0X.   
    day_count = CharCount(day)
    if (int(day) < 10) & (day_count < 2):
        day = '0' + day
    else:
        pass
    
    #########################################
    ### Prompt the user for the model run ###
    #########################################
    ModelRun = input('Input the desired model (either 00, 06, 12, or 18):  ')
    
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
              
    #######################
    ### End of Function ###
    #######################
    return year, month, day, ModelRun

#%% 
####################################
### determine_source Function ######
####################################
    
def determine_source(year, month, day, ModelRun):
    '''
    This function is designed to determine where the GFS data is obtained from
      given the date that the user input. There are three different sources:
          NCEP (for data less than 1 week old)
          NCEI (for data between 1 week to 1 year old)
          Data request from the NCEI (for data older than 1 year)
      If a data request is needed, the script also checks that the .g2 folder
      exists in the Data folder, as it should for the rest of the
      extract_GFS_variable.sh program.
    
    Outputs:
        source - The source data will come from (NCEI, or NCEP)
        request - Binary number if a data request is needed
    '''
    
    # Initialize some variables. This are the base urls of the GFS data and
    #   does not change. This will be printed for the user so s/he knows where
    #   the data is coming from.
    Request_page = 'https://www.ncdc.noaa.gov/has/HAS.FileAppRouter?' +\
                   'datasetname=GFS3&subqueryby=STATION&applname=&outdest=FILE'
    NCEI_source  = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
    NCEP_source  = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
    
    # Initialize change in time variables equivalent to 1 week and 1 year
    DeltaYear = timedelta(days = 365)
    DeltaWeek = timedelta(days = 7)
    
    # Construct the date from the user's given year, month, and day.
    date = datetime(int(year), int(month), int(day))
    
    # Determine how old the date is
    DeltaTime = datetime.now() - date
    
    # Check if the date is 1 week or older and if it is older than 1 year.
    #   If it is newer than 1 week, the source is NCEP.
    #   If it is newer than 1 year, but older than 1 week, the source is NCEI.
    #   If it is older than 1 year, a data request is needed. Check that the
    #   .g2 folder exists in the Data folder.
    if DeltaTime > DeltaYear:
        # Ensure the user knows to make a request for old data.
        print('The requested date is over 1 year old. A data request is ' +\
              'required for this. This can be made at: \n' +\
              Request_page + '\n' +\
              'Please ensure the request is made, then input the request ID.' +\
              'here.\n')
        
        RequestID = input( )
        request = 1
        source = None
        
        # Collect the requested datafile.
        print('Collecting the requested data.')
        
        path = './Data/tmp/'
        filename  = 'gfs_3_' + str(year) + str(month) + str(day) + str(ModelRun) +\
            '.g2.tar'
        url_start = 'https://www1.ncdc.noaa.gov/pub/has/model/'
        url = url_start + str(RequestID) + '/' + filename
        wget.download(url, path + filename)
        
        print('\n Unpacking the .tar file.')
        tf = tarfile.open(path + filename)
        
        directoryname = 'gfs_3_' + str(year) + str(month) + str(day) + str(ModelRun) +\
            '.g2/'
        tf.extractall(path = './Data/tmp/' + directoryname)
        
        # Check that the .g2 folder exists in the Data folder.
#         try:
#             filename = 'gfs_3_' + str(year) + str(month) + str(day) + '_' +\
#                         '00' + str(ModelRun) + '_003.grb2'
#             path = './Data/gfs_3_' + str(year) + str(month) + str(day) +\
#                     str(ModelRun) + '.g2/'
# #           path = '/Users/Rarrell/Downloads/gfs_3_' + str(year) + str(month) +\
# #                  str(day) + str(model_run) + '.g2/'
#             grb_file = path + filename
            
#             grb = pygrib.open(grb_file)
#             grb.close()
#         except OSError:
#             raise OSError('.grb2 folder not found in the Data folder.' +\
#                           'Please make to make the request at: \n' +\
#                           Request_page + '\n' +\
#                           'and download the data to the Data folder.')
    
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

    #######################
    ### End of Function ###
    #######################
    return source, request

#%% 
########################################
### input_model_variables Function #####
########################################
    
def input_model_variables():
    '''
    This function prompts the user for the variable name and height of the
      variable to be extracted from the GFS file. The function also checks
      that the input names and heights are valid/exist in the GFS grib file.
      
    Outputs:
        Var - name of the variable to be extracted
        TypeOfHeight - Type of the height of where the variable is located
                       (e.g., surface, isobaricInhPa, etc.)
        height - The height Var is located at. A valid value of height depends
                 on TypeOfHeight (e.g., if TypeOfHeight = surface, then 
                 height = 0. If TypeOfHeight = isobaricInhPa, then height can
                 have many values.)
    '''
    
    # Collect the valid input options
#    path = '/Users/Rarrell/Desktop/Thesis_Research/Scripts/Data_Collection/txt_files/'
    path = './Scripts/Data_Collection/txt_files/'
    ValidVar = np.loadtxt(path + 'GFS_VarNames.txt', dtype = str, delimiter = ',')
    ValidTypeOfHeights = np.loadtxt(path + 'GFS_TypeOfHeights.txt', dtype = str,
                                    delimiter = ',')
    ValidHeights = np.loadtxt(path + 'GFS_Heights.txt', dtype = str, delimiter = ',')
    
    #####################################
    ### Prompt for the variable name ####
    #####################################
    Var = input('Please enter the variable to be extract. ' +\
                'See GFS_Variable_Names.txt for a list of options. \n'  )
    
    # Check the input name is one of the valid variable names
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

    ##############################################
    ### Prompt the user for the type of height ###
    ##############################################
    TypeOfHeight = input('Please enter the type of height the variable is located' +\
                         'at. See GFS_Type_Of_Heights.txt for a list of options. \n'  )

    # Check the input type of height is one of the valid type of heights
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
    
    ######################################
    ### Prompt the user for the height ###
    ######################################
    height = input('Type in the height the variable is located at. The units ' +\
                   '(e.g. Pa, km, PVU, etc.) is determined by the type of height.' +\
                   'If the type of height is the surface, then the height is 0. \n'  )
    
    # Check the input height is one of the valid heights. This will remove non-integers.
    # Check that the height is 0 if the TypeOfHeight is surface
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
    
    #######################
    ### End of Function ###
    #######################
    return Var, TypeOfHeight, height

#%% 
############################
### write_tmp Function #####
############################
def write_tmp(year, month, day, source, request, ModelRun, 
              Var = None, TypeOfHeight = None, height = None):
    '''
    This function takes all the GFS file and variable information obtained in
      collect_data, determine_source, and input_model_variable and writes it to
      a temporary file (called tmp.txt) that was created at the start of the
      extract_GFS_variable.sh program for later use in that program.
    '''
#    path = '/Users/Rarrell/Desktop/Thesis_Research/'
    # Define the path and open the .txt file. Not this is designed to append,
    #   not to write a new file.
    path = './'
    f = open(path + 'tmp.txt', 'a')
    
    # Write the variable information. For values that are not a string
    #   already, turn them into a string. The file is comma seperated in a
    #   single row (columns are read differently in the np.loadtxt used in
    #   later scripts).
    f.write(str(request) + ',')
    f.write(year + ',')
    f.write(month + ',')
    f.write(day + ',')
    f.write(ModelRun + ',')
    f.write(str(source) + ',')
    
    # Close the .txt file
    f.close()
    
    # Append to the variable tmp file if they are entered.
    if (Var is None) & (TypeOfHeight is None) & (height is None):
        pass
    else:
        f = open(path + 'tmp_var.txt', 'a')
        f.write(str(Var) + ',')
        f.write(str(TypeOfHeight) + ',')
        f.write(str(height))
        
        f.close()
    
    #######################
    ### End of Function ###
    #######################

#%% 
#########################################
### Call and Run the Main Function ######
#########################################
    
if __name__ == '__main__':
    main()

#####################
### End of Script ###
#####################