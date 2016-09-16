DEISOM=HP_RNASeq/run
READ_DIR=bams
GFF_DIR=indexed
OUTPUT_DIR=output
DEISOM_ARGS="--trunk-size 5 --gff-dir ${GFF_DIR} --min-read 0 --read-len 100 --paired-end 200.0 20.0 --out-iter 5 --in-iter 500 --human-readable"
CONDITIONS=( N T )

for c in ${CONDITIONS[@]}
do
   BAMS=
   for i in {1..5}
   do
      BAMS=${BAMS}${READ_DIR}/${i}${c}/accepted_hits.bam,
   done

   CMD="${DEISOM} --bams ${BAMS} $DEISOM_ARGS --output ${OUTPUT_DIR}/${c}"
   echo ${CMD}
done
