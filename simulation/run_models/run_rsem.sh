#!/bin/bash

RSEM_ARGS="-paired-end --no-qualities --no-bam-output --bowtie2 -p 1"
REFERENCE="reference"
READ_FOLDER=reads
OUTPUT_MATRIX=isoform.matrix
OUTPUT_FILE=isoform.result
CONDITIONS=( "N" "T" )

for i in {1..5}
do
   for c in ${CONDITIONS[@]}
   do
      CMD="rsem-calculate-expression ${RSEM_ARGS} ${READ_FOLDER}/${i}${c}/paired_1.fa ${RSEM_ARGS}/${i}${c}/paried_2.fa ${REFERENCE} ${i}${c}"
      echo $CMD
   done
done

for i in {1..5}
do
   CMD="rsem-generate-data-matrix"
   for c in ${CONDITIONS[@]}
   do
      CMD="${CMD} ${i}${c}.isoforms.results"
   done
   CMD="${CMD} > $OUTPUT_MATRIX"
   echo ${CMD}
done

CMD="rsem-run-ebseq $OUTPUT_MATRIX 5,5 $OUTPUT_FILE"
echo ${CMD}
