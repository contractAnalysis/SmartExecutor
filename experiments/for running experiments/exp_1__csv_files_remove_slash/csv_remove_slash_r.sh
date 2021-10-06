#!/bin/bash



#num=$1
num=20

for j in `seq 0 ${num} `
do

	echo "token_contracts_info_part_"${j}".csv"
	sed 's/\r$//' "contracts_info_part_"${j}".csv"  > "contracts_info_part_"${j}"_r.csv"



done



