#include "HP_Param.h"
#include <sstream>

using namespace std;
HP_Param::HP_Param() {
	isSingleEnd = true;
	meanInsertedLen = 0;
	stdInsertedLen = 0;
	minRead = 1;
	readLen = 0;
	geneID = string();
	gff = string();
	bams = list<string>();
	overhangLen = 1;
	outputDir = string();
	numInIters = 5;
	numOutIters = 5000;
}

// to string
string HP_Param::toString() const {
	stringstream sstm;
	sstm << "Gene ID: " << geneID << endl;
	sstm << "GFF filename: " << gff << endl;
	sstm << "Number of subjects: " << bams.size() << endl;
	sstm << "BAM filename: ";
	for (list<string>::const_iterator ii = bams.begin(); ii != bams.end(); ii++)
		sstm << (*ii) << " ";
	sstm << endl;
	sstm << "Output directory: " << outputDir << endl;
	sstm << "Minimum number of reads: " << minRead << endl;
	sstm << "Read length: " << readLen << endl;
	sstm << "Overhang length: " << overhangLen << endl;
	if (isSingleEnd)
		sstm << "Single end" << endl;
	else
		sstm << "Pair end with mean " << meanInsertedLen
		<< " and std " << stdInsertedLen << endl;
	sstm << "Inner iteration limit: " << numInIters << endl;
	sstm << "Outer iteration limit: " << numOutIters << endl;
	return sstm.str();
}