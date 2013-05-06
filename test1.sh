./run --gene-ids ENSG00000100065 \
		--gffs data/chr22/ENSG00000100065.gff3\
		--bams read/11N/RUM.sorted.bam,read/22N/RUM.sorted.bam,read/28N/RUM.sorted.bam,read/30N/RUM.sorted.bam \
		--min-read 1 \
		--read-len 90 \
		--paired-end 178.0 56.0 \
		--out-iter 5 \
		--in-iter 5000 \
		--output ./ \
#		--bams  /Volumes/LENOVO_USB_HDD/data/11N/RUM.sorted.bam,/Volumes/LENOVO_USB_HDD/data/22N/RUM.sorted.bam,/Volumes/LENOVO_USB_HDD/data/28N/RUM.sorted.bam,/Volumes/LENOVO_USB_HDD/data/30N/RUM.sorted.bam\
#--bams  read/11N/RUM.sorted.bam \
