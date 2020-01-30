# This is a shell script designed, for a user defined date and model run, download a GFS 
#   file, extract a single variable and place it into a .nc file for all forecast hours
#   in the model run for user use. This program does not accept any arguments. User
#   defined values are input during the function. This script assumes that it is located 
#   in a folder called "Thesis_Research/" and the a directory map the same as the one 
#   shown in directory_map.pdf.
#   (see https://github.com/Rarell/Thesis_Research for more detail).


# Create some initial, temporary files/folders.
touch ./tmp.txt
touch ./tmp_var.txt
mkdir ./Data/tmp

# Collect information for model date, and variable name and location
echo 'Collecting date and variable information'
python ./Scripts/Create_GEFS_Datasets/input_model_information.py "$1"

# If an error occurred in the script, end the program. The tmp file will be empty if there
#   was an error. Reference the python error message to determine the error.
if [ -s ./tmp.txt ]
then
    :
else
    echo 'An error occurred. See the above message.'
    exit 1
fi

# Create the forecast hours that will be looped over depending on the date
python ./Scripts/Create_GEFS_Datasets/write_FH.py tmp.txt

# If the user did not ask for a specific variable, but an indice then inform the user what
#   variable(s) will be extracted
if [ "$1" = '-v' ]
then
	:
elif [ "$@" = 'SESR' ]
then
	echo 'The following variables will be collected'
	echo 'Latent heat surface flux'
	echo 'Potential evapotranspiration'
elif [ "$@" = 'HI' ]
then
	echo 'Working on it.'
fi

# Loop over all the forecast hours. For each hour, extract the variable and place in a
#   temporary file.
while read FH
    do
    echo " "
    echo "Extracting data from the $FH Forecast Hour."
    
    for EM in ./Scripts/Create_GEFS_Datasets/Ensemble_Members.txt # Loop through each ensemble member
    	do
    	
    	# Download the data
    	python ./Scripts/Create_GEFS_Datasets/GFS_Download.py $FH $EM tmp.txt
    	
    	# If the user asked for a specific variable, then tmp_var is not empty, and collect that
    	#   variable. If an indice was asked for, tmp_var is empty, and collect those variables.
    	if [ -s ./tmp_var.txt ]
    	then
    		python ./Scripts/Create_GEFS_Datasets/grb_to_nc.py $FH $EM tmp_var.txt tmp.txt
    	elif [ "$@" = 'SESR' ]
    	then
    		echo 'Latent heat net flux','surface','0' > tmp_var.txt
    		python ./Scripts/Create_GEFS_Datasets/grb_to_nc.py $FH $EM tmp_var.txt tmp.txt
    		> tmp_var.txt
    	
    		echo 'Potential evaporation rate','surface','0' > tmp_var.txt
    		python ./Scripts/Create_GEFS_Datasets/grb_to_nc.py $FH $EM tmp_var.txt tmp.txt
    		> tmp_var.txt
    	elif ["$@" = 'HI' ]
    	then
    		echo 'Something happens'
    	fi
    
    	if [ -d ./Data/tmp/*.g2 ]
		then
			:
		else
			rm ./Data/tmp/*.grb2
		fi
    done	
#    python Scripts/Data_Collection/grb_to_nc.py tmp.txt $FH
done < ./Scripts/Create_GEFS_Datasets/GFS_Forecast_Times.txt

# Merge the temporary .nc files into full forecast files for each ensemble member
for EM in ./Scripts/Create_GEFS_Datasets/Ensemble_Members.txt
	do
	echo " "
	echo "Merging .nc files for ensemble member $EM."
	
	# If the user asked for a specific variable, then tmp_var is not empty, and collect that
    #   variable. If an indice was asked for, tmp_var is empty, and collect those variables.
    if [ -s ./tmp_var.txt ]
	then
    	python ./Scripts/Create_GEFS_Datasets/merge_nc.py $EM tmp_var.txt tmp.txt ./Scripts/Create_GEFS_Datasets/GFS_FH.txt
    elif [ "$@" = 'SESR' ]
	then
    	echo 'Latent heat net flux','surface','0' > tmp_var.txt
    	python ./Scripts/Create_GEFS_Datasets/merge_nc.py $EM tmp_var.txt tmp.txt ./Scripts/Create_GEFS_Datasets/GFS_FH.txt
		> tmp_var.txt
    	
    	echo 'Potential evaporation rate','surface','0' > tmp_var.txt
    	python ./Scripts/Create_GEFS_Datasets/merge_nc.py $EM tmp_var.txt tmp.txt ./Scripts/Create_GEFS_Datasets/GFS_FH.txt
		> tmp_var.txt
    elif ["$@" = 'HI' ]
	then
    	echo 'Something happens'
    fi
	
done
	
echo " " # Start a new line

# Merge the ensemble members
echo 'Merging ensemble members.'
if [ -s ./tmp_var.txt ]
then
	python ./Scripts/Create_GEFS_Datasets/merge_ensembles.py ./Scripts/Create_GEFS_Datasets/Ensemble_Members.txt tmp_var.txt tmp.txt

elif [ "$@" = 'SESR' ]
then
	echo 'Latent heat net flux','surface','0','lhtfl' > tmp_var.txt
	python ./Scripts/Create_GEFS_Datasets/merge_ensembles.py ./Scripts/Create_GEFS_Datasets/Ensemble_Members.txt tmp_var.txt tmp.txt
	> tmp_var.txt
	
	echo 'Potential evaporation rate','surface','0','pevpr' > tmp_var.txt
	python ./Scripts/Create_GEFS_Datasets/merge_ensembles.py ./Scripts/Create_GEFS_Datasets/Ensemble_Members.txt tmp_var.txt tmp.txt
	
	echo 'Calculating SESR'
	python ./Scripts/Data_Collection/SESR_calculations.py 'lhtfl' 'pevpr' tmp.txt
elif [ "$@" = 'HI' ]
then
	echo 'More stuff happens'
fi

# Remove temporary files/folders
echo ' '
echo 'Removing temporary files and folders.'
rm ./tmp.txt
rm ./tmp_var.txt
rm -r ./Data/tmp/

# Check if a .g2 file exits (as it should for a data request) and remove at the end of the program
if [ -d ./Data/tmp/*.g2 ]
then
    rm -r ./Data/tmp/*.g2/
    rm -r ./Data/tmp/*.g2.tar
else
	:
fi

# End of program
echo 'Done'