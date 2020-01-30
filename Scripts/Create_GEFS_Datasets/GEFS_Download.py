#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:15:01 2019

@author: Rarrell

This script is designed as part of the create_synthetic_dataset.sh program.
This script uses the given arguements to construct the url in which the GFS
  data is located. It then downloads data for use later in program.
  
Arguements:
    script - The name of the python script
    FH - String. The forecast hour (in a 00X or 0XX, etc. form) of the GFS data.
    Year     - String. The year for which the GFS run was made.
    Month    - String. The month for which the GFS run was made.
    Day      - String. The day for which the GFS run was made.
    ModelRun - String. The model run of the GFS (00, 06, 12, or 18).
    Source   - The source of the GFS data (whether it is from the NCEI, or NCEP
               url or if a data request is needed).
"""

#%% 
#####################################
### Import some libraries ###########
#####################################

import sys, os, warnings
import wget
import numpy as np
import pygrib
import tarfile


#%% 
############################
### Main Function ##########
############################

def main():
    '''
    This is the main function of the script. This is where the arguements are
      collected and the other function is called.
    '''
    
    # Collect the arguements
    script = sys.argv[0]
    FH     = sys.argv[1]
    EM     = sys.argv[2]
    Parameters = sys.argv[3]
    
    # Unpack the parameters
    Request, Year, Month, Day, ModelRun, Source, RequestID = np.loadtxt(Parameters,
                                                                        dtype = str, delimieter = ',', unpack = True)
    
    # Construct the url and download the data.
    DownloadData(Year, Month, Day, ModelRun, FH, Source, EM, Request, RequestID)

#%%
####################################
### DownloadData Function ##########
####################################

def DownloadData(Year, Month, Day, ModelRun, FH, Source, EM, Request, RequestID):
    '''
    This function constructs the url for a GFS data file given the date, model
      run, forecast hour, and location of the data, and downloads the data to
      a temporary file.
      
    Inputs:
        Year - String. The year in which the GFS model was made.
        Month - String. The year in which the GFS model was made.
        Day - String. The year in which the GFS model was made.
        ModelRun - String. The model run of the GFS model (00, 06, 12, or 18)
        FH - String. The forecast hour of the GFS data being collected.
        Source - String. The location of the GFS data (either NCEI, NCEP, or
                 a data request).
    '''
    
    # Define the path to which the data will be downloaded
    path = './Data/tmp/'
    
    # Construct the url based on where the data will come from. If a data
    #   request is needed, check that the file exists.
    if int(Request) == 1:
        filename  = 'gens_3_' + str(Year) + str(Month) + str(Day) + str(ModelRun) +\
            '_' + str(EM) + '.g2.tar'
        url_start = 'https://www1.ncdc.noaa.gov/pub/has/'
        url = url_start + str(RequestID) + '/' + filename
        wget.download(url, path + filename)
        tf = tarfile.open(path + filename)
        
        directoryname = 'gens_3_' + str(Year) + str(Month) + str(Day) + str(ModelRun) +\
            '_' + str(EM) + '.g2/'
        tf.extractall(path = './Data/tmp/' + directoryname)
        # Check that the .g2 folder exists in the Data folder.
#         try:
#             filename = 'gfs_3_' + str(Year) + str(Month) +\
#                        str(Day) + '_' + '00' + str(ModelRun) + '_003.grb2'
#             path = './Data/gfs_3_' + str(Year) + str(Month) +\
#                    str(Day) + str(ModelRun) + '.g2/'
# #           path = '/Users/Rarrell/Downloads/gfs_3_' + str(year) + str(month) +\
# #                  str(day) + str(model_run) + '.g2/'
#             grb_file = path + filename
                
#             grb = pygrib.open(grb_file)
#             grb.close()
#         except OSError:
#             raise OSError('.grb2 folder not found in the Data folder.' +\
#                           'Please make to make the request at (see above ' +\
#                           'instructions) and download the data to the ' +\
#                           'Data folder.')

    elif Source == 'NCEI':
        # Construct the NCEI url
        filename = 'gens-b_3_' + str(Year) + str(Month) + str(Day) + '_' +\
                   str(ModelRun) + '00' + '_' + str(FH) + '_' + str(EM) + '.grb2'
        url_start = 'https://nomads.ncdc.noaa.gov/data/gens/'
        url_end   = str(Year) + str(Month) + '/' + str(Year) + str(Month) +\
                    str(Day) +'/' + filename
        url = url_start + url_end
        
        # Download the data
        wget.download(url, path + filename)
        
    else:
        # Construct the NCEP url
        filename  = 'gep' + str(EM) + '.t' + str(ModelRun) + 'z.pgrb2f' +\
                        str(FH)
        url_start = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gens/prod/'
        url_end   = 'gefs.' + str(Year) + str(Month) + str(Day) + '/' +\
                    str(ModelRun) + '/pgrb2/' + filename
        url = url_start + url_end
        
        # Download the data
        wget.download(url, path + filename)
    
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
