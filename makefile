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


all: deisomRun deisomIndexGFF deisomEstimate deisomSplit deisomKL

deisomRun: DEIsoM_Model.o DEIsoM_RunTask.o DEIsoM_Gene.o DEIsoM_Gff.o DEIsoM_Read.o DEIsoM_Param.o asa121.o
	${LIBTOOL} ${LINK_OPTION} ${CC} ${CC_OPTION} -o deisomRun DEIsoM_Model.o DEIsoM_RunTask.o DEIsoM_Gene.o DEIsoM_Gff.o DEIsoM_Read.o DEIsoM_Param.o asa121.o ${LBFGS_LA} ${LIB}

deisomIndexGFF: DEIsoM_IndexGff.o DEIsoM_Gene.o DEIsoM_Gff.o 
	${CC}  ${CC_OPTION} -o deisomIndexGFF DEIsoM_IndexGff.o DEIsoM_Gene.o DEIsoM_Gff.o ${LIB}

deisomSplit: DEIsoM_SplitTasks.o DEIsoM_Param.o DEIsoM_Gff.o
	${CC} ${CC_OPTION} -o deisomSplit DEIsoM_SplitTasks.o DEIsoM_Param.o DEIsoM_Gff.o DEIsoM_Gene.o ${LIB}

deisomEstimate: DEIsoM_Estimate.o DEIsoM_Read.o DEIsoM_Gene.o DEIsoM_Gff.o
	${CC} ${CC_OPTION} -o deisomEstimate DEIsoM_Estimate.o DEIsoM_Read.o DEIsoM_Gene.o DEIsoM_Gff.o ${LIB}

deisomKL: DEIsoM_ComputeKL.o DEIsoM_Gff.o DEIsoM_Gene.o
	${CC} ${CC_OPTION} -o deisomKL DEIsoM_ComputeKL.o DEIsoM_Gff.o DEIsoM_Gene.o ${LIB}

DEIsoM_Model.o: DEIsoM_Model.cpp DEIsoM_Model.h
	${CC} ${CC_OPTION} -c -o DEIsoM_Model.o DEIsoM_Model.cpp ${INCLUDE}

DEIsoM_RunTask.o: DEIsoM_RunTask.cpp DEIsoM_Model.h
	${CC} ${CC_OPTION} -c -o DEIsoM_RunTask.o DEIsoM_RunTask.cpp ${INCLUDE}

DEIsoM_IndexGff.o: DEIsoM_IndexGff.cpp DEIsoM_Gene.h DEIsoM_Gff.h
	${CC} ${CC_OPTION} -c -o DEIsoM_IndexGff.o DEIsoM_IndexGff.cpp ${BOOST_INCLUDE}

DEIsoM_Gene.o: DEIsoM_Gene.cpp DEIsoM_Gene.h
	${CC} ${CC_OPTION} -c -o DEIsoM_Gene.o DEIsoM_Gene.cpp

DEIsoM_Param.o: DEIsoM_Param.cpp DEIsoM_Param.h
	${CC} ${CC_OPTION} -c -o DEIsoM_Param.o DEIsoM_Param.cpp

DEIsoM_Gff.o: DEIsoM_Gff.cpp DEIsoM_Gff.h
	${CC} ${CC_OPTION} -c -o DEIsoM_Gff.o DEIsoM_Gff.cpp

DEIsoM_Read.o: DEIsoM_Read.cpp DEIsoM_Read.h
	${CC} ${CC_OPTION} -c -o DEIsoM_Read.o DEIsoM_Read.cpp ${BOOST_INCLUDE}

DEIsoM_Estimate.o: DEIsoM_Estimate.cpp
	${CC} ${CC_OPTION} -c -o DEIsoM_Estimate.o DEIsoM_Estimate.cpp ${INCLUDE}

DEIsoM_SplitTasks.o: DEIsoM_SplitTasks.cpp
	${CC} ${CC_OPTION} -c -o DEIsoM_SplitTasks.o DEIsoM_SplitTasks.cpp ${BOOST_INCLUDE}

DEIsoM_ComputeKL.o: DEIsoM_ComputeKL.cpp
	${CC} ${CC_OPTION} -c -o DEIsoM_ComputeKL.o DEIsoM_ComputeKL.cpp ${BOOST_INCLUDE}

asa121.o: asa121.h asa121.cpp
	${CC} ${CC_OPTION} -c -o asa121.o asa121.cpp

clean:
	rm -f *.o
	rm -f deisomRun deisomIndexGFF deisomEstimate deisomSplit desiomKL
