DIR=..

echo "Creating indexed gff files..."
${DIR}/deisomIndexGFF ENSG00000100065.gff3 ./indexed/

echo "Counting reads..."
./count_reads.sh

echo ""
echo "Creating scripts for N group..."
${DIR}/deisomSplit --path ${DIR}/deisomRun \
		--trunk-size 1\
		--gff-dir indexed\
		--bams read/11N/RUM.sorted.bam,read/22N/RUM.sorted.bam,read/28N/RUM.sorted.bam,read/30N/RUM.sorted.bam \
		--read-len 90 \
		--paired-end 178.0 56.0 \
		--out-iter 5 \
		--in-iter 5000 \
		--read-count readCountsN.txt \
		--output ./output/N/ \
      --human-readable

echo ""
echo "Creating scripts for T group..."
${DIR}/deisomSplit --path ${DIR}/deisomRun \
		--trunk-size 1\
		--gff-dir indexed\
		--bams read/11T/RUM.sorted.bam,read/22T/RUM.sorted.bam,read/28T/RUM.sorted.bam,read/30T/RUM.sorted.bam \
		--read-len 90 \
		--paired-end 178.0 56.0 \
		--out-iter 5 \
		--in-iter 5000 \
		--read-count readCountsT.txt \
		--output ./output/T/ \
      --human-readable

echo ""
echo "Running scripts for N group..."
chmod 777 ./output/N/cluster_scripts/task0.sh
./output/N/cluster_scripts/task0.sh

echo ""
echo "Creating scripts for T group..."
chmod 777 ./output/T/cluster_scripts/task0.sh
./output/T/cluster_scripts/task0.sh

echo "Compute KL..."
${DIR}/deisomKL ./indexed ./output/T/ ./output/N/ > kl.txt
