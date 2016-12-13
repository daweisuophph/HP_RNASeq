/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#include "DEIsoM_Param.h"
#include <sstream>

using namespace std;
DEIsoM_Param::DEIsoM_Param() {
	isSingleEnd = true;
	meanInsertedLen = 0;
	stdInsertedLen = 0;
	bfgsUpdate = true;
	minRead = 1;
	readLen = 0;
	geneID = string();
	gff = string();
	bams = vector<string>();
	overhangLen = 1;
	outputDir = string();
	numInIters = 5;
	numOutIters = 5000;
	outputBinary = true;
}

// to string
string DEIsoM_Param::toString() const {
	stringstream sstm;
	sstm << "Gene ID: " << geneID << endl;
	sstm << "GFF filename: " << gff << endl;
	sstm << "Number of subjects: " << bams.size() << endl;
	sstm << "BAM filename: ";
   for (int i = 0; i < bams.size(); i++)
      sstm << bams[i] << " ";
	sstm << endl;
	sstm << "Output directory: " << outputDir << endl;
   if (readCountFile.size()) sstm << "Read count file: " << readCountFile << endl;
   sstm << "Update method: " << (bfgsUpdate? "BFGS": "Netwon's") << endl;
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
	if (outputBinary) {
		sstm << "Output in binary format." << endl;
	} else {
		sstm << "Output in human readable format." << endl;
	}
	return sstm.str();
}
