/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */

#ifndef _HP_GFF
#define _HP_GFF

#include <fstream>
#include <string>
#include <map>
#include <list>
#include "HP_Gene.h"

using namespace std;

class HP_Gff {
private:
	map <string, list<HP_MRNA> > mRNAByParent;
	map <string, list<HP_Exon> > exonByParent;
	/* ignore all cdss */
	//map <string, list<CDS> > cdssByParent;
	
	void readGFFV3(ifstream &ifs);
	void readRecordV3(char **fields);
public:
	list <HP_Gene> genes;
	HP_Gff(string gffFile);
	void getGenesBySeqid(map<string, list<HP_Gene> > &genesBySeqid);
};
#endif
