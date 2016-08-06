READ_DIR="./read/"
REPLICATES=( "11" "22" "28" "30" )
CONDITION=( "N" "T" )
FILENAME="RUM.sorted.bam"
OUTPUT="readCounts"
for c in "${CONDITION[@]}"
do
	rm -f ${OUTPUT}${c}.txt
	for r in "${REPLICATES[@]}"
	do
		../third_party/samtools/samtools idxstats ${READ_DIR}${r}${c}/${FILENAME} | \
		awk ' BEGIN {sum = 0;} { sum = sum + $3; } END { print sum } ' >> ${OUTPUT}${c}.txt
	done
done
