/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#ifndef _DEIsoM_PARAM
#define _DEIsoM_PARAM
#include <string>
#include <vector>

using namespace std;

class DEIsoM_Param {
public:
	bool isSingleEnd; // single end or pair end
	double meanInsertedLen; // mean of inserted Length
	double stdInsertedLen; //  standard deviation of inserted length
	bool bfgsUpdate; // whether to use bfgs (default) or newton's method to update alpha
	int minRead; // minimum number of reads
	int readLen; // read length
	int overhangLen; // Length of overhang constraints imposed on junctions.
	string geneID; // Gene ID
	string gff; //GFF filename
	vector<string> bams;
	string readCountFile; 
	string outputDir; //output directory
	int numOutIters; // number of interations to update alpha
	int numInIters; // number of interations to update hidden variables
	bool outputBinary; //whether to output in binary format
	
	DEIsoM_Param();
	string toString() const;
};

#endif
