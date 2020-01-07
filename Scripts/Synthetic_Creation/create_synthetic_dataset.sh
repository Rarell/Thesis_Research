# This is a shell script designed to take a range of dates (minimum 30 days), defined by
#   the user, and collect data from GFS outputs and calculated drought indices. For each
#   date, the 3 hour forecast is taken (some variables are not calculated for the initial
#   conditions) for each model in a given day and a variable(s) is extracted from the
#   GFS file. The full file is kept placed in the Data/Synthetic_Data folder for user
#   reference. The data is used to then calculate drought indices (the exact index(ices)
#   is given in the arguments for this script). The drought indice is also placed in the
#   Data/Synthetic_Data folder.
# 
# Note this script assumes that it is located in a folder called "Thesis_Research/" and 
#   the a directory map the same as the one shown in directory_map.pdf.
#   (see https://github.com/Rarell/Thesis_Research for more detail).
#
# Drought indices used as arguments must be used in their acronym form (e.g., SESR, CTP,
#   HI, etc.).
#
# Current drought indices that this script will calculate and are valid arguments:
#   SESR - The standardized evaporative stress ratio. The SESR data is calculated on a
#          daily time scale and it is calculated via the daily sum and average latent
#          heat and potential evaporation separately.


# Create some initial temporary files
touch ./tmp_list.txt
mkdir ./Data/tmp

# Prompt the user for the start and end dates, and collect all the dates and sources
python ./Scripts/Synthetic_Creation/collect_dates.py

# Check that the tmp file is filled with values. If not, then there was an error.
#   In the case of an error, stop the program and reference the python error message.
if [ -s ./tmp_list.txt ]
then
    :
else
    echo 'An error occurred. See the above message.'
    exit 1
fi

# Inform the user what variables will be extracted and will be used in the calculations.
#   The variables depends on what arguments are used in the inputs.
echo 'The follow variables will be collected:'
if [ "$@" = 'SESR' ]
then
	echo 'Latent heat net flux'
	echo 'Potential evaporation rate'
elif [ "$@" = 'HI' ]
then
	echo 'Working on it.'
fi

# Loop over each date and model run and collect the required variable, store the variable
#   in a temporary netcdf file, and delete the memory expensive grib files.
while read tmp
    do
    
    # Inform the user what day and model run is being worked on
    echo " "
    echo $tmp
    
    # Download the grib file. Note many surface variables are not calculated until the
    #   3 hour forecast, so the three hour forecast is downloaded. Note some variables are
    #   a 0 - 3 hour average (e.g., latent heat flux). The forecast hour for these variables
    #   are listed as 0.
    python ./Scripts/Synthetic_Creation/GFS_Download.py 003 $tmp
#    python ./Scripts/Synthetic_Creation/GFS_Download.py 006 $tmp
    
    # From the grib files, extract the required variables needed for a specific index, 
    #   and place them in a temporary netcdf file in the Data/tmp/ folder.
    if [ "$@" = 'SESR' ]
    then
    	python ./Scripts/Synthetic_Creation/extract_grb_variable.py 'Latent heat net flux' 'surface' '0' '003' $tmp
    	python ./Scripts/Synthetic_Creation/extract_grb_variable.py 'Potential evaporation rate' 'surface' '0' '003' $tmp
    fi
    
    # Since the required data is stored elsewhere, remove the grib files as they take a
    #   lot of memory.
    # Check if a .g2 file exists (it should for a data request), then remove the whole,
    #   unnecessary .g2 file. If not, there was no data request, and remove the .grb2 file.
    if [ -d ./Data/tmp/*.g2 ]
    then
    	rm -r ./Data/tmp/*.g2/
    	rm -r ./Data/tmp/*.g2.tar
    else
    	rm ./Data/tmp/*.grb2
    fi
done < ./tmp_list.txt

# Add new line after the downloading data is done
echo " "

# Merge the extracted data into a single netcdf file for all time periods, then perform
#   the index calculations. The resulting data is placed in the Data/Synthetic_Data/ folder.
if [ "$@" = 'SESR' ]
then
	echo 'Merging LE and PET .nc files'
	python ./Scripts/Synthetic_Creation/merge_nc.py 'Latent heat net flux' 'lhtfl' 'surface' tmp_list.txt
	python ./Scripts/Synthetic_Creation/merge_nc.py 'Potential evaporation rate' 'pevpr' 'surface' tmp_list.txt
	
	echo 'Calculating SESR'
	# Calculate SESR
	python ./Scripts/Synthetic_Creation/SESR_calculations.py lhtfl pevpr tmp_list.txt
fi

# Remove tmp files/folders
echo ' '
echo 'Removing temporary files and folders.'
rm tmp_list.txt
rm -r Data/tmp/

# End of program
echo 'Done'

# Short names:
#   Potential evaporation rate: pevpr
#   Latent heat net flux: lhtfl
    
