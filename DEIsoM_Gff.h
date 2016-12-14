/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */

#ifndef _DEIsoM_GFF
#define _DEIsoM_GFF

#include <fstream>
#include <string>
#include <map>
#include <list>
#include "DEIsoM_Gene.h"

using namespace std;

class DEIsoM_Gff {
private:
	map <string, list<DEIsoM_MRNA> > mRNAByParent;
	map <string, list<DEIsoM_Exon> > exonByParent;
	/* ignore all cdss */
	//map <string, list<CDS> > cdssByParent;
	
	void readGFFV3(ifstream &ifs);
	void readRecordV3(char **fields);
public:
	list <DEIsoM_Gene> genes;
	DEIsoM_Gff(string gffFile);
	void getGenesBySeqid(map<string, list<DEIsoM_Gene> > &genesBySeqid);
};
#endif
