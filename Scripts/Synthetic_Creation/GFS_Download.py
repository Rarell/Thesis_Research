#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 12:15:01 2019

@author: Rarrell
"""

#%% Import libraries
import sys, os, warnings
import wget
import numpy as np
import pygrib


#%% Main function
def main():
    script = sys.argv[0]
    FH     = sys.argv[1]
    Year   = sys.argv[2]
    Month  = sys.argv[3]
    Day    = sys.argv[4]
    ModelRun = sys.argv[5]
    Source   = sys.argv[6]
    
    DownloadData(Year, Month, Day, ModelRun, FH, Source)

#%% Download data function
def DownloadData(Year, Month, Day, ModelRun, FH, Source):
    '''
    '''
    path = './Data/tmp/'
    
    if Source == 'Request':
        # Check that the .g2 folder exists in the Data folder.
        try:
            filename = 'gfs_3_' + str(Year) + str(Month) +\
                       str(Day) + '_' + '00' + str(ModelRun) + '_003.grb2'
            path = './Data/gfs_3_' + str(Year) + str(Month) +\
                   str(Day) + str(ModelRun) + '.g2/'
#           path = '/Users/Rarrell/Downloads/gfs_3_' + str(year) + str(month) +\
#                  str(day) + str(model_run) + '.g2/'
            grb_file = path + filename
                
            grb = pygrib.open(grb_file)
            grb.close()
        except OSError:
            raise OSError('.grb2 folder not found in the Data folder.' +\
                          'Please make to make the request at (see above ' +\
                          'instructions) and download the data to the ' +\
                          'Data folder.')
                              
    elif Source == 'NCEI':
        filename = 'gfs_3_' + str(Year) + str(Month) + str(Day) + '_' +\
                   str(ModelRun) + '00' + '_' + str(FH) + '.grb2'
        url_start = 'https://nomads.ncdc.noaa.gov/data/gfs-avn-hi/'
        url_end   = str(Year) + str(Month) + '/' + str(Year) + str(Month) +\
                    str(Day) +'/' + filename
        url = url_start + url_end
        
        wget.download(url, path + filename)
        
    else:
        filename  = 'gfs.t' + str(ModelRun) + 'z.pgrb2.1p00.f' +\
                        str(FH)
        url_start = 'https://www.ftp.ncep.noaa.gov/data/nccf/com/gfs/prod/'
        url_end   = 'gfs.' + str(Year) + str(Month) + str(Day) + '/' +\
                    str(ModelRun) + '/' + filename
        url = url_start + url_end
        
        wget.download(url, path + filename)
        

#%% Call main function
if __name__ == '__main__':
    main()