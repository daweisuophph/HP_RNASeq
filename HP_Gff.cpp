/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */
#include "HP_Gff.h"
#include <iostream>
#include <cstring>
#include <cstdlib>

#define MAX_LINE_SIZE 1000

using namespace std;

HP_Gff::HP_Gff(string gffFile) {
	ifstream ifs(gffFile.c_str());
	if (!ifs) {
		cerr << "Error: the GFF file \"" << gffFile << "\" cannot be opened." << endl;
		exit(1);
	}
	char *line = new char[MAX_LINE_SIZE];
	ifs.getline(line, MAX_LINE_SIZE);
	if (ifs.good()) {
		char *dummy = new char[MAX_LINE_SIZE];
		int versionNum = 0;
		if (sscanf(line, "%s %d", dummy, &versionNum) != EOF && versionNum == 3){
		//if (strcmp(line, "##gff-version 3") == 0 || strcmp(line, "##gff-version\t3") == 0) {
			delete dummy;
			readGFFV3(ifs);
		} else {
			cerr << "Error: the version \"" << line << "\" is not supported" << endl;
			delete dummy;
			delete line;
			exit(1);
		}
	}
	delete line;
	ifs.close();
}

static bool compareByStart(HP_Exon e1, HP_Exon e2) {
	return (e1.start < e2.start);
}

void HP_Gff::readGFFV3(ifstream &ifs) {
	int lineNumber = 1;
	char *line = new char[MAX_LINE_SIZE];
	while (!ifs.eof()) {
		ifs.getline(line, MAX_LINE_SIZE);
		lineNumber++;
		// skip empty lines and comments
		if (line[0] != 0 && line[0] != '#') {
			char *fields[10];
			int i = 0;
			fields[i] = strtok(line, "\t");
			while ((fields[++i] = strtok(NULL, "\t")) && i < 9);
			if (i != 9 || fields[9] != 0) {
				cerr << "Error: invalid number of fields (should be 9) at line " << lineNumber << "." << endl;
				exit(1);
			}
			readRecordV3(fields);
			if (lineNumber % 100000 == 0) {
				cerr << lineNumber << " lines loaded."<< endl;
			}
		}
	}
	delete line;
	
	
	//build heriachy
	for (list<HP_Gene>::iterator ig = genes.begin();
		 ig != genes.end(); ig++) {
		if (mRNAByParent.find(ig->ID) != mRNAByParent.end()) {
			ig->mRNAs = mRNAByParent[ig->ID];
			for (list<HP_MRNA>::iterator im = ig->mRNAs.begin();
				 im != ig->mRNAs.end(); im++) {
				if (exonByParent.find(im->ID) != exonByParent.end()) {
					im->exons  = exonByParent[im->ID];
					im->exons.sort(compareByStart);
					/*
					for (list<Exon>::iterator ie = im->exons.begin();
						 ie != im->exons.end(); ie++) {
						if (cdssByParent.find(ie->ID) != cdssByParent.end()) {
							ie->cdss = cdssByParent[ie->ID];
						}
					}
					 */
				}
				/*
				if (cdssByParent.find(im->ID) != cdssByParent.end()) {
					im->cdss = cdssByParent[im->ID];
				}
				 */
			}							
		}
	}
}

void HP_Gff::readRecordV3(char **fields) {
	HP_Record *r = 0;
	if (strcmp(fields[2], "gene") == 0) {
		r = new HP_Gene();
	} else if (strcmp(fields[2], "transcript") == 0 ||
				strcmp(fields[2], "aberrant_processed_transcript") == 0 ||
				strcmp(fields[2], "nc_primary_transcript") == 0 ||
				strcmp(fields[2], "NMD_transcript_variant") == 0 ||
				strcmp(fields[2], "processed_transcript") == 0 ||
				strcmp(fields[2], "pseudogenic_transcript") == 0 ||
				strcmp(fields[2], "RNA") == 0 ||
				strcmp(fields[2], "rRNA") == 0 ||
				strcmp(fields[2], "snoRNA") == 0 ||
				strcmp(fields[2], "snRNA") == 0 ||
				strcmp(fields[2], "miRNA") == 0 ||
				strcmp(fields[2], "lincRNA") == 0 ||
			   strcmp(fields[2], "mRNA") == 0) {
		r = new HP_MRNA();
	} else if (strcmp(fields[2], "exon") == 0) {
		r = new HP_Exon();
	} else if (strcmp(fields[2], "CDS") == 0) {
		r = new HP_CDS();
	} else {
		return;
	}
	/*
	if (fields[0][0] != 'c' ||
		fields[0][1] != 'h' ||
		fields[0][2] != 'r') {
		r->seqid = "chr" + string(fields[0]);
	} else {
		r->seqid = string(fields[0]);
	}
	*/
	r->seqid = string(fields[0]);
	r->source = string(fields[1]);
	r->type = string(fields[2]);
	r->start = atoi(fields[3]);
	r->end = atoi(fields[4]);
	if (strcmp(fields[5], ".") == 0) {
		r->score = -1;
	} else {
		r->score = atof(fields[5]);
	}
	r->strand = fields[6][0];
	if (strcmp(fields[7], ".") == 0) {
		r->phase = -1;
	} else {
		r->phase = atoi(fields[7]);
	}
	char *attribute = strtok(fields[8], ";");
	do {
		string att(attribute);
		int i = 0;
		while (i < att.size() && att[i] != '=') i++;
		if (attribute[i] == '=') {
			string tag = att.substr(0, i);
			string value = att.substr(i+1, att.size());
			if (tag.compare("ID") == 0) {
				r->ID = value;
			} else if (tag.compare("Name") == 0) {
				r->name = value;
			} else if (tag.compare("Parent") == 0) {
				r->parent = value;
			}
		}
	} while ((attribute = strtok(NULL, ";")) != NULL);
	
	if (r->type[0] == 'g') {
		genes.push_back(*(HP_Gene *)r);
	} else if (r->type[0] == 'e') {
		if (exonByParent.find(r->parent) == exonByParent.end()) {
			exonByParent[r->parent] = list<HP_Exon>();
		}
		exonByParent[r->parent].push_back(*(HP_Exon *)r);
	} else if (r->type[0] == 'C') {
		/*
		if (cdssByParent.find(r->parent) == cdssByParent.end()) {
			cdssByParent[r->parent] = list<CDS>();
		}
		cdssByParent[r->parent].push_back(*(CDS *)r);
		 */
	} else {
		if (mRNAByParent.find(r->parent) == mRNAByParent.end()) {
			mRNAByParent[r->parent] = list<HP_MRNA>();
		}
		mRNAByParent[r->parent].push_back(*(HP_MRNA *)r);
	}
	delete r;
}

void HP_Gff::getGenesBySeqid(map<string, list<HP_Gene> > &genesBySeqid) {
	for (list<HP_Gene>::iterator ig = genes.begin(); ig != genes.end(); ig++) {
		if (genesBySeqid.find(ig->seqid) == genesBySeqid.end()) {
			genesBySeqid[ig->seqid] = list<HP_Gene>();
		}
		genesBySeqid[ig->seqid].push_back(*ig);
	}
}

