BOOST_HOME=/Users/pengh/Documents/boost/
BOOST_INCLUDE=-I${BOOST_HOME}/include
BOOST_LIB=-L${BOOST_HOME}/lib
BOOST_OPTION=-lboost_system -lboost_filesystem

SAMTOOLS_HOME=/Users/pengh/Documents/samtools/
SAMTOOLS_INCLUDE=-I${SAMTOOLS_HOME}
SAMTOOLS_LIB=-L${SAMTOOLS_HOME}/ -L${SAMTOOLS_HOME}/bcftools
SAMTOOLS_OPTION=-lbam -lbcf -lcurses -lz -lpthread -lm

INCLUDE=-I/Users/pengh/Documents/boost/include 
LIB=-L/Users/pengh/Documents/boost/lib -Lsamtools -Lsamtools/bcftools
OPTIONS=-lboost_system -lboost_filesystem -lbam -lbcf  -lcurses  -lz

all: run indexGFF estimate split kl

run: HP_Model.o HP_RunTask.o HP_Gene.o HP_Gff.o HP_Read.o HP_Param.o asa121.o
	g++ -O2 -o run HP_Model.o HP_RunTask.o HP_Gene.o HP_Gff.o HP_Read.o HP_Param.o asa121.o ${BOOST_LIB} ${SAMTOOLS_LIB} ${BOOST_OPTION} ${SAMTOOLS_OPTION}

indexGFF: HP_IndexGff.o HP_Gene.o HP_Gff.o 
	g++  -O2 -o indexGFF HP_IndexGff.o HP_Gene.o HP_Gff.o ${BOOST_LIB} ${SAMTOOLS_LIB} ${BOOST_OPTION} ${SAMTOOLS_OPTION}  

split: HP_SplitTasks.o HP_Param.o HP_Gff.o
	g++ -O2 -o split HP_SplitTasks.o HP_Param.o HP_Gff.o HP_Gene.o ${BOOST_LIB} ${BOOST_OPTION}

estimate: HP_Estimate.o HP_Read.o HP_Gene.o HP_Gff.o
	g++ -O2 -o estimate HP_Estimate.o HP_Read.o HP_Gene.o HP_Gff.o ${SAMTOOLS_LIB} ${SAMTOOLS_OPTION}

kl: HP_ComputeKL.o HP_Gff.o HP_Gene.o
	g++ -O2 -o kl HP_ComputeKL.o HP_Gff.o HP_Gene.o ${BOOST_LIB} ${BOOST_OPTION}

HP_Model.o: HP_Model.cpp HP_Model.h
	g++ -O2 -c -o HP_Model.o HP_Model.cpp ${SAMTOOLS_INCLUDE} ${BOOST_INCLUDE}

HP_RunTask.o: HP_RunTask.cpp HP_Model.h
	g++ -O2 -c -o HP_RunTask.o HP_RunTask.cpp ${BOOST_INCLUDE}

HP_IndexGff.o: HP_IndexGff.cpp HP_Gene.h HP_Gff.h
	g++ -O2 -c -o HP_IndexGff.o HP_IndexGff.cpp ${BOOST_INCLUDE}

HP_Gene.o: HP_Gene.cpp HP_Gene.h
	g++ -O2 -c -o HP_Gene.o HP_Gene.cpp

HP_Param.o: HP_Param.cpp HP_Param.h
	g++ -O2 -c -o HP_Param.o HP_Param.cpp

HP_Gff.o: HP_Gff.cpp HP_Gff.h
	g++ -O2 -c -o HP_Gff.o HP_Gff.cpp

HP_Read.o: HP_Read.cpp HP_Read.h
	g++ -O2 -c -o HP_Read.o HP_Read.cpp ${BOOST_INCLUDE}

HP_Estimate.o: HP_Estimate.cpp
	g++ -O2 -c -o HP_Estimate.o HP_Estimate.cpp ${SAMTOOLS_INCLUDE} ${BOOST_INCLUDE}

HP_SplitTasks.o: HP_SplitTasks.cpp
	g++ -O2 -c -o HP_SplitTasks.o HP_SplitTasks.cpp ${BOOST_INCLUDE}

HP_ComputeKL.o: HP_ComputeKL.cpp
	g++ -O2 -c -o HP_ComputeKL.o HP_ComputeKL.cpp ${BOOST_INCLUDE}

asa121.o: asa121.h asa121.cpp
	g++ -O2 -c -o asa121.o asa121.cpp

clean:
	rm -f *.o
	rm -f run indexGFF estimate split
