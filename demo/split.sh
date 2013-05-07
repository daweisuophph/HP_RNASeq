./split	--path /Users/pengh/Desktop/HP_RNAseq/run \
		--trunk-size 10 \
		--gff-dir indexed\
		--bams /Users/pengh/Desktop/HP_RNAseq/demo/read/11N/RUM.sorted.bam,/Users/pengh/Desktop/HP_RNAseq/demo/read/22N/RUM.sorted.bam,/Users/pengh/Desktop/HP_RNAseq/demo/read/28N/RUM.sorted.bam,/Users/pengh/Desktop/HP_RNAseq/demo/read/30N/RUM.sorted.bam \
		--min-read 1 \
		--read-len 90 \
		--paired-end 178.0 56.0 \
		--out-iter 20 \
		--in-iter 5000 \
		--output ~/Desktop/HP_RNAseq/output \
#--bams  read/11N/RUM.sorted.bam \
#		--bams  /Volumes/LENOVO_USB_HDD/data/11N/RUM.sorted.bam,/Volumes/LENOVO_USB_HDD/data/22N/RUM.sorted.bam,/Volumes/LENOVO_USB_HDD/data/28N/RUM.sorted.bam,/Volumes/LENOVO_USB_HDD/data/30N/RUM.sorted.bam\
