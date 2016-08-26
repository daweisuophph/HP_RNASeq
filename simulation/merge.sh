#!/bin/bash
test_folder_name=reads_repeated_exp

COND=( N T )

for exp in {1..10}
do
	test_name=${test_folder_name}/exp${exp}
	for c in "${COND[@]}"
	do
		for i in {1..5}
		do
			READ_DIR=./${test_name}/${i}${c}/
			rm -rf ${READ_DIR}/paired_1.fa ${READ_DIR}/paired_2.fa
			cp ${READ_DIR}/pairedDE_1.fa ${READ_DIR}/paired_1.fa
			cat ${READ_DIR}/pairedNONDE_1.fa >> ${READ_DIR}/paired_1.fa
			cp ${READ_DIR}/pairedDE_2.fa ${READ_DIR}/paired_2.fa
			cat ${READ_DIR}/pairedNONDE_2.fa >> ${READ_DIR}/paired_2.fa
		done
	done
done
