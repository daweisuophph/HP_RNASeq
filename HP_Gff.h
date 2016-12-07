/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */

#ifndef _HP_GFF
#define _HP_GFF

#include <fstream>
#include <string>
#include <unordered_map>
#include <vector>
#include "HP_Gene.h"

using namespace std;

class HP_Gff {
private:
	unordered_map <string, vector<HP_Exon> > exonByParent;
	
	void readGFFV3(ifstream &ifs);
	void readRecordV3(char **fields);
public:
   vector<HP_MRNA> mRNAs;
	HP_Gff(string gffFile);
};
#endif
