# This is a shell script designed, for a user defined date and model run, download a GFS 
#   file, extract a single variable and place it into a .nc file for all forecast hours
#   in the model run for user use. This program does not accept any arguments. User
#   defined values are input during the function. This script assumes that it is located 
#   in a folder called "Thesis_Research/" and the a directory map the same as the one 
#   shown in directory_map.pdf.
#   (see https://github.com/Rarell/Thesis_Research for more detail).


# Create some initial, temporary files/folders.
touch ./tmp.txt
mkdir ./Data/tmp

# Collect information for model date, and variable name and location
echo 'Collecting date and variable information'
python ./Scripts/Data_Collection/input_model_information.py

# If an error occured in the script, end the program. The tmp file will be empty if there
#   was an error. Reference the python error message to determine the error.
if [ -s ./tmp.txt ]
then
    :
else
    echo 'An error occured. See the above message.'
    exit 1
fi

# Loop over all the forecast hours. For each hour, extract the variable and place in a
#   temporary file.
while read FH
    do
    echo " "
    echo "Extracting data from the $FH Forecast Hour."
    python Scripts/Data_Collection/grb_to_nc.py tmp.txt $FH
done < ./Scripts/Data_Collection/GFS_Forecast_Times.txt

# Merge the temporary .nc files
echo 'Merging .nc files.'
python Scripts/Data_Collection/merge_nc.py tmp.txt Scripts/Data_Collection/txt_files/GFS_FH.txt

# Remove temporary files/folders
echo ' '
echo 'Removing temporary files and folders.'
rm tmp.txt
rm -r Data/tmp/

# End of program
echo 'Done'