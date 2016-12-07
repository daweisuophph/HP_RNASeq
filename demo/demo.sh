DIR=..

echo "Counting reads..."
./count_reads.sh

echo ""
echo "Run N group..."
${DIR}/run \
   --gff ENSG00000100065.gff3 \
   --bams read/11N/RUM.sorted.bam,read/22N/RUM.sorted.bam,read/28N/RUM.sorted.bam,read/30N/RUM.sorted.bam \
   --read-len 90 \
   --paired-end 178.0 56.0 \
   --out-iter 5 \
   --in-iter 5000 \
   --read-count readCountsN.txt \
   --output ./output/N/ \
   --human-readable

echo ""
echo "Run T group..."
${DIR}/run \
   --gff ENSG00000100065.gff3 \
   --bams read/11T/RUM.sorted.bam,read/22T/RUM.sorted.bam,read/28T/RUM.sorted.bam,read/30T/RUM.sorted.bam \
   --read-len 90 \
   --paired-end 178.0 56.0 \
   --out-iter 5 \
   --in-iter 5000 \
   --read-count readCountsT.txt \
   --output ./output/T/ \
   --human-readable


#echo "Compute KL..."
#${DIR}/kl ./indexed ./output/T/ ./output/N/ > kl.txt
