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

all: run indexGFF estimate 

run: HP_Model.o HP_RunTask.o HP_Gene.o HP_Gff.o HP_Read.o
	g++ -o run HP_Model.o HP_RunTask.o HP_Gene.o HP_Gff.o HP_Read.o ${BOOST_LIB} ${SAMTOOLS_LIB} ${BOOST_OPTION} ${SAMTOOLS_OPTION}

indexGFF: HP_IndexGff.o HP_Gene.o HP_Gff.o 
	g++  -o indexGFF HP_IndexGff.o HP_Gene.o HP_Gff.o ${BOOST_LIB} ${SAMTOOLS_LIB} ${BOOST_OPTION} ${SAMTOOLS_OPTION} 

estimate: HP_Estimate.o HP_Read.o HP_Gene.o HP_Gff.o
	g++ -o estimate HP_Estimate.o HP_Read.o HP_Gene.o HP_Gff.o ${SAMTOOLS_LIB} ${SAMTOOLS_OPTION}

HP_Model.o: HP_Model.cpp HP_Model.h
	g++ -c -o HP_Model.o HP_Model.cpp ${SAMTOOLS_INCLUDE} ${BOOST_INCLUDE}

HP_RunTask.o: HP_RunTask.cpp HP_Model.h
	g++ -c -o HP_RunTask.o HP_RunTask.cpp ${BOOST_INCLUDE}

HP_IndexGff.o: HP_IndexGff.cpp HP_Gene.h HP_Gff.h
	g++ -c -o HP_IndexGff.o HP_IndexGff.cpp ${BOOST_INCLUDE}

HP_Gene.o: HP_Gene.cpp HP_Gene.h
	g++ -c -o HP_Gene.o HP_Gene.cpp

HP_Gff.o: HP_Gff.cpp HP_Gff.h
	g++ -c -o HP_Gff.o HP_Gff.cpp

HP_Read.o: HP_Read.cpp HP_Read.h
	g++ -c -o HP_Read.o HP_Read.cpp

HP_Estimate.o: HP_Estimate.cpp
	g++ -c -o HP_Estimate.o HP_Estimate.cpp ${SAMTOOLS_INCLUDE}

clean:
	rm -f *.o
	rm run
	rm indexGFF
