Author: Hao Peng (pengh@purdue.edu)
Date: Aug 26, 2016
Version: 1.1v

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


Please check demo folder for an example run of our model on a small data.
In simulation folder, we keep the codes for generating synthetic data in our paper.

For evaluation of our method, we can use 
kl [indexed GFF] [output of group 1] [output of group 2]  
for general evaluations (either paired or unpaired between groups)
or
kl --beta [indexed GFF] [output of group 1] [output of group 2] 
for paired evalutions only
