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
	int minRead; // minimum number of reads
	int readLen; // read length
	int overhangLen; // Length of overhang constraints imposed on junctions.
	string geneID; // Gene ID
	string gff; //GFF filename
	list<string> bams;
	string outputDir; //output directory
	int numOutIters; // number of interations to update alpha
	int numInIters; // number of interations to update hidden variables
	
	HP_Param();
	string toString() const;
};

#endif