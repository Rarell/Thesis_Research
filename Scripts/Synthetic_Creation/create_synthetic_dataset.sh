# Documentation

# Create some initial temporary files
touch ./tmp_list.txt
mkdir ./Data/tmp

# Collect all the dates and sources
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

# Inform the user what variables will be extracted and will be used in the calculations
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
    #   a 0 - 3 hour average (e.g., latent heat flux).
    python ./Scripts/Synthetic_Creation/GFS_Download.py 003 $tmp
#    python ./Scripts/Synthetic_Creation/GFS_Download.py 006 $tmp
    
    # From the grib files, extract the required variables needed for a specific indice, 
    #   and place them in a temporary netcdf file in the Data/tmp/ folder.
    if [ "$@" = 'SESR' ]
    then
    	python ./Scripts/Synthetic_Creation/extract_grb_variable.py 'Latent heat net flux' 'surface' '0' '003' $tmp
    	python ./Scripts/Synthetic_Creation/extract_grb_variable.py 'Potential evaporation rate' 'surface' '0' '003' $tmp
    fi
    
    # Since the required data is stored elsewhere, remove the grib files as they take a
    #   lot of memory.
    rm ./Data/tmp/*.grb2
done < ./tmp_list.txt

# Add new line after the downloaded data is done
echo " "

# Merge the extracted data into a single netcdf file for all time periods, then perform
#   the indice calculations. The resulting data is placed in the Data/Synthetic_Data/ folder.
if [ "$@" = 'SESR' ]
then
	echo 'Merging LE and PET .nc files'
	python ./Scripts/Synthetic_Creation/merge_nc.py 'Latent heat net flux' 'lhtfl' 'surface' tmp_list.txt
	python ./Scripts/Synthetic_Creation/merge_nc.py 'Potential evaporation rate' 'pevpr' 'surface' tmp_list.txt
	
	echo 'Calculating SESR'
	# Calculate SESR and create figures
fi

# Remove tmp files/folders
echo ' '
echo 'Removing temporary files and folders.'
rm tmp_list.txt
rm -r Data/tmp/

# Done
echo 'Done'

# Short names:
#   Potential evaporation rate: pevpr
#   Latent heat net flux: lhtfl
    
