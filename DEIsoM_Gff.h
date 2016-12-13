/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */

#ifndef _DEIsoM_GFF
#define _DEIsoM_GFF

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "DEIsoM_Gene.h"

using namespace std;

class DEIsoM_Gff {
private:
	unordered_map <string, vector<DEIsoM_Exon> > exonByParent;
	
	void readGFFV3(ifstream &ifs);
	void readRecordV3(char **fields);
public:
   vector<DEIsoM_MRNA> mRNAs;
	DEIsoM_Gff(string gffFile);
};
#endif
