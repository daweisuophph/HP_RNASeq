#BOOST_HOME=
BOOST_INCLUDE=-I${BOOST_HOME}/include
BOOST_LIB=-L${BOOST_HOME}/lib
#BOOST_OPTION=${BOOST_HOME}/lib/libboost_system.a ${BOOST_HOME}/lib/libboost_filesystem.a
BOOST_OPTION=-lboost_system -lboost_filesystem

SAMTOOLS_HOME=third_party/samtools
SAMTOOLS_INCLUDE=-I${SAMTOOLS_HOME}
SAMTOOLS_LIB=-L${SAMTOOLS_HOME}/ #-L${SAMTOOLS_HOME}/bcftools
SAMTOOLS_OPTION=-lbam  -lcurses -lz -lpthread -lm

HTSLIB_INCLUDE=-I${SAMTOOLS_HOME}/htslib-1.3.1
HTSLIB_LIB=-L${SAMTOOLS_HOME}/htslib-1.3.1/
HTSLIB_OPTION=-lhts

LBFGS_INCLUDE=-Ithird_party/liblbfgs/include
LBFGS_LA=third_party/liblbfgs/lib/liblbfgs.la
LIBTOOL=third_party/liblbfgs/libtool
LINK_OPTION=--mode=link --tag CXX

CC=g++
CC_OPTION=-O2

INCLUDE=${BOOST_INCLUDE} ${SAMTOOLS_INCLUDE} ${HTSLIB_INCLUDE} ${LBFGS_INCLUDE}
LIB=${BOOST_LIB} ${SAMTOOLS_LIB} ${HTSLIB_LIB} ${BOOST_OPTION} ${SAMTOOLS_OPTION} ${HTSLIB_OPTION}


all: run indexGFF estimate split kl

run: HP_Model.o HP_RunTask.o HP_Gene.o HP_Gff.o HP_Read.o HP_Param.o asa121.o
	${LIBTOOL} ${LINK_OPTION} ${CC} ${CC_OPTION} -o run HP_Model.o HP_RunTask.o HP_Gene.o HP_Gff.o HP_Read.o HP_Param.o asa121.o ${LBFGS_LA} ${LIB}

indexGFF: HP_IndexGff.o HP_Gene.o HP_Gff.o 
	${CC}  ${CC_OPTION} -o indexGFF HP_IndexGff.o HP_Gene.o HP_Gff.o ${LIB}

split: HP_SplitTasks.o HP_Param.o HP_Gff.o
	${CC} ${CC_OPTION} -o split HP_SplitTasks.o HP_Param.o HP_Gff.o HP_Gene.o ${LIB}

estimate: HP_Estimate.o HP_Read.o HP_Gene.o HP_Gff.o
	${CC} ${CC_OPTION} -o estimate HP_Estimate.o HP_Read.o HP_Gene.o HP_Gff.o ${LIB}

kl: HP_ComputeKL.o HP_Gff.o HP_Gene.o
	${CC} ${CC_OPTION} -o kl HP_ComputeKL.o HP_Gff.o HP_Gene.o ${LIB}

HP_Model.o: HP_Model.cpp HP_Model.h
	${CC} ${CC_OPTION} -c -o HP_Model.o HP_Model.cpp ${INCLUDE}

HP_RunTask.o: HP_RunTask.cpp HP_Model.h
	${CC} ${CC_OPTION} -c -o HP_RunTask.o HP_RunTask.cpp ${INCLUDE}

HP_IndexGff.o: HP_IndexGff.cpp HP_Gene.h HP_Gff.h
	${CC} ${CC_OPTION} -c -o HP_IndexGff.o HP_IndexGff.cpp ${BOOST_INCLUDE}

HP_Gene.o: HP_Gene.cpp HP_Gene.h
	${CC} ${CC_OPTION} -c -o HP_Gene.o HP_Gene.cpp

HP_Param.o: HP_Param.cpp HP_Param.h
	${CC} ${CC_OPTION} -c -o HP_Param.o HP_Param.cpp

HP_Gff.o: HP_Gff.cpp HP_Gff.h
	${CC} ${CC_OPTION} -c -o HP_Gff.o HP_Gff.cpp

HP_Read.o: HP_Read.cpp HP_Read.h
	${CC} ${CC_OPTION} -c -o HP_Read.o HP_Read.cpp ${BOOST_INCLUDE}

HP_Estimate.o: HP_Estimate.cpp
	${CC} ${CC_OPTION} -c -o HP_Estimate.o HP_Estimate.cpp ${INCLUDE}

HP_SplitTasks.o: HP_SplitTasks.cpp
	${CC} ${CC_OPTION} -c -o HP_SplitTasks.o HP_SplitTasks.cpp ${BOOST_INCLUDE}

HP_ComputeKL.o: HP_ComputeKL.cpp
	${CC} ${CC_OPTION} -c -o HP_ComputeKL.o HP_ComputeKL.cpp ${BOOST_INCLUDE}

asa121.o: asa121.h asa121.cpp
	${CC} ${CC_OPTION} -c -o asa121.o asa121.cpp

clean:
	rm -f *.o
	rm -f run indexGFF estimate split kl
