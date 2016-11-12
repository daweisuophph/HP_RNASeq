Author: Hao Peng (pengh@purdue.edu)
Date: Nov. 12, 2016
Version: 1.2v

Prerequisites:
1. boost v1.53.0+ 
	http://www.boost.org/
	For old versions, please remove -lboost_system in makefile
	Remember to set ${BOOST_HOME}
2. samtools v0.1.19+
	Old versions may apply. not tested.
	Remember to set ${SAMTOOLS_HOME}

Notes: Install boost v1.53 on MAX OS:
install_name_tool -id /boost/lib/libboost_system.dylib libboost_system.dylib
install_name_tool -change libboost_system.dylib @executable_path/../Frameworks/libboost_system.dylib libboost_filesystem.dylib

How to install:
1. set makefile
2. make

How to run:
1. To index reference:
indexGFF [gff3 file] [indexed reference dir]

2. To split jobs:
split	--path [path to HP_RNASeq]/run 
	--gff-dir [indexed reference directory] 
	--bams [bam 1],[bam 2],...,[bam M] 
	--read-len [read length]
	--paired-end [mean and std for paired-end reads]
	--out-iter [number of iterations for updating alpha]
	--in-iter [number of iterations for updating variational parameters]
	--output [output folder]
	
Other Options for split:
--trunk-size [trunk size]
--human-readable  (default is binary)

3. To evaluate of our method, we can use:
kl [indexed GFF] [output of group 1] [output of group 2]  
for general evaluations (either paired or unpaired between groups)
or
kl --beta [indexed GFF] [output of group 1] [output of group 2] 
for paired evalutions only

Please check demo folder for an example run of our model on a small data.
In simulation folder, we keep the codes for generating synthetic data in our paper.
