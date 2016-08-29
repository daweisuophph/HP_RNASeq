/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: Aug 29, 2016
 Version: 1.1v
 */
#ifndef _HP_PARAM
#define _HP_PARAM
#include <string>
#include <list>

using namespace std;

class HP_Param {
public:
	bool isSingleEnd; // single end or pair end
	double meanInsertedLen; // mean of inserted Length
	double stdInsertedLen; //  standard deviation of inserted length
	bool bfgsUpdate; // whether to use bfgs (default) or newton's method to update alpha
	int minRead; // minimum number of reads
	int readLen; // read length
	int overhangLen; // Length of overhang constraints imposed on junctions.
   double initAlpha; // initial value of alpha
	string geneID; // Gene ID
	string gff; //GFF filename
	list<string> bams;
	string readCountFile; 
	string outputDir; //output directory
	int numOutIters; // number of interations to update alpha
	int numInIters; // number of interations to update hidden variables
	bool outputBinary; //whether to output in binary format
	
	HP_Param();
	string toString() const;
};

#endif
