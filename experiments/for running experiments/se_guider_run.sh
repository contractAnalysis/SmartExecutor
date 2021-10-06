#!/bin/bash

index=$1
data_file=$2
cli_timeout=$3
create_timeout=$4
execution_timeout=$5
control_level=$6

#cli_timeout=960
#create_timeout=10
#execution_timeout=900


result_folder_prefix='results_se_guider_'
contract_folder='__contracts_1818'


rm -rf ${result_folder_prefix}${index}_cl${control_level}
mkdir ${result_folder_prefix}${index}_cl${control_level}

exec < ${data_file} || exit 1
 # read header # read (and ignore) the first line
while IFS="," read file_name solc_version contract_nam 
do
	 echo ""
	 contract_name=${contract_nam:0:-1}
	 echo "++++ ${file_name}  :  ${solc_version}  :  ${contract_name} ++++" |& tee -a ./${result_folder_prefix}${index}_cl${control_level}/"${file_name}__${contract_name}.txt" 
	 solc-select use "$solc_version" 	 
	

	 start=$SECONDS
	 
	timeout ${cli_timeout} myth analyze $(echo ./${contract_folder}/$file_name:${contract_name}) -fdg -cl ${control_level} --create-timeout $create_timeout --execution-timeout $execution_timeout |& tee -a ./${result_folder_prefix}${index}_cl${control_level}/"${file_name}__${contract_name}.txt"


	end=$SECONDS
	runtime=$((end-start))
	echo "time_used: "${runtime}" seconds"	|& tee -a ./${result_folder_prefix}${index}_cl${control_level}/"${file_name}__${contract_name}.txt"
	
	# save the time for each contract
	echo "#@contract_info_time" |& tee -a ./${result_folder_prefix}${index}_cl${control_level}/"${file_name}__${contract_name}.txt"
	echo ${file_name[$i]}:${solc_version}:${contract_name}:${runtime}:${cli_timeout}:${create_timeout}:${execution_timeout} |& tee -a ./${result_folder_prefix}${index}_cl${control_level}/"${file_name}__${contract_name}.txt"

	# one example of directing terminal ouput to file
	# echo "Solc version:" $solc_version >> ./analysis_result/"${file_name}__${contract_name}.txt" 2>&1
done

