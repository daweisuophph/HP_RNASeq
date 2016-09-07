#!/bin/bash
test_folder_name=reads_repeated_exp

READLEN=100
NREAD=300000
PAIREDEND="200,20"
ERRORRATE=0

REFERENCE=./1.fa

POSBIAS=./input/sampleposbias.txt
READERR=./input/samplereaderror.txt

for exp in {1..10}
do
	test_name=${test_folder_name}/exp${exp}
	for i in {1..5}
	do
		BED=./DE/sample.bed

		FASTAFILE=./${test_name}/${i}N/pairedDE
		RANDEXPLV=./${test_name}/${i}N/explvprofileDE.txt
		CMD1="gensimreads.py -s ${exp}${i}1 -e $RANDEXPLV  -n $NREAD -b $POSBIAS -l $READLEN -p $PAIREDEND  $BED "
		CMD2="getseqfrombed.py -s ${exp}${i}2 -b $READERR -f A -r $ERRORRATE -l $READLEN -  $REFERENCE"
		echo "$CMD1  | $CMD2 | splitfasta.py -o $FASTAFILE" > scripts/gen_reads_${exp}_${i}NDE.sh

		FASTAFILE=./${test_name}/${i}T/pairedDE
		RANDEXPLV=./${test_name}/${i}T/explvprofileDE.txt
		CMD1="gensimreads.py -s ${exp}${i}3 -e $RANDEXPLV  -n $NREAD -b $POSBIAS -l $READLEN -p $PAIREDEND  $BED "
		CMD2="getseqfrombed.py -s ${exp}${i}4 -b $READERR -f A -r $ERRORRATE -l $READLEN -  $REFERENCE"
		echo "$CMD1  | $CMD2 | splitfasta.py -o $FASTAFILE" > scripts/gen_reads_${exp}_${i}TDE.sh
	done

	for i in {1..5}
	do
		BED=./NON_DE/sample.bed

		FASTAFILE=./${test_name}/${i}N/pairedNONDE
		RANDEXPLV=./${test_name}/${i}N/explvprofileNONDE.txt
		CMD1="gensimreads.py -s ${exp}${i}5 -e $RANDEXPLV  -n $NREAD -b $POSBIAS -l $READLEN -p $PAIREDEND  $BED "
	   CMD2="getseqfrombed.py -s ${exp}${i}6 -b $READERR -f A -r $ERRORRATE -l $READLEN -  $REFERENCE"
		echo "$CMD1  | $CMD2 | splitfasta.py -o $FASTAFILE" > scripts/gen_reads_${exp}_${i}NNONDE.sh

		FASTAFILE=./${test_name}/${i}T/pairedNONDE
		RANDEXPLV=./${test_name}/${i}T/explvprofileNONDE.txt
		CMD1="gensimreads.py -s ${exp}${i}7 -e $RANDEXPLV  -n $NREAD -b $POSBIAS -l $READLEN -p $PAIREDEND  $BED "
		CMD2="getseqfrombed.py -s ${exp}${i}8 -b $READERR -f A -r $ERRORRATE -l $READLEN -  $REFERENCE"
		echo "$CMD1  | $CMD2 | splitfasta.py -o $FASTAFILE" > scripts/gen_reads_${exp}_${i}TNONDE.sh
	done
done
