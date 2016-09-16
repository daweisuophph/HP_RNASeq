#!/bin/bash

CONDITIONS=( "N" "T" )
MISO_ARGS="--read-len 100 --paired-end 200 20"
OUTPUT_DIR=output
BAM_DIR=bams

for c in ${CONDITIONS[@]}
do
   for i in {1..5}
   do
      CMD="miso --run indexed ${BAM_DIR}/$i$c/accepted_hits.bam --output-dir ${OUTPUT_DIR}/$i$c/ ${MISO_ARGS}"
      echo ${CMD}
   done
done
