/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#include "DEIsoM_Read.h"
#include "DEIsoM_Gff.h"
#include "DEIsoM_Gene.h"
#include "sam.h"
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <cmath>
#include <cstdlib>
#include <sstream>

using namespace std;

static int fetch_func(const bam1_t *b, void *data) {
    map<string, list<DEIsoM_Read> > &readsByName = *(map<string, list<DEIsoM_Read> >*)data;
	
	DEIsoM_Read read;
	read.pos = b->core.pos+1; // change from 0 based to 1 based
	read.name = bam1_qname(b);
	read.len = b->core.l_qseq;
	read.cigar = vector<int32_t>(b->core.n_cigar);
	uint32_t *cigar = bam1_cigar(b);
	for (int i = 0; i < b->core.n_cigar; i++ ) {
		read.cigar[i] = cigar[i];
	}
	if (readsByName.find(read.name) == readsByName.end()) {
		readsByName[read.name] = list<DEIsoM_Read>();
	}
	readsByName[read.name].push_back(read);
	
    return 0;
}

int main(int argc, char **argv) {
	if (argc < 2) {
		cerr << "Usage: " << argv[0] << " <GFF> <bam 1> <bam 2> ..." << endl;
		exit(1);
	}
	DEIsoM_Gff gff = DEIsoM_Gff(string(argv[1]));
	
	list<int> insertedLengths;
	for (int i = 2; i < argc; i++) {
		samfile_t *in = samopen(argv[i], "rb", 0);
		if (in == 0) {
			cerr << "Error: cannot read BAM file \"" << argv[i] << "\"." << endl;
			exit(1);
		}
		
		int count = 0;
		for (list<DEIsoM_Gene>::iterator ii = gff.genes.begin();
			 ii != gff.genes.end(); ii++, count++) {
			if (count % 10000 == 0) {
				cerr << "finished " << count << " genes." << endl;
			}
			map<string, list<DEIsoM_Read> >readsByName;
			int beg, end;
			ii->getBounds(beg, end);
			stringstream sstm;
			sstm << ii->seqid << ":" << beg << "-" << end;
			int newBeg, newEnd, ref;
			bam_parse_region(in->header, sstm.str().c_str(),
							 &ref, &newBeg, &newEnd);
			//cout << ref << " " << newBeg << " " << newEnd << endl;
			if (ref < 0) {
				//cerr << "Error: region \"" << sstm.str() << "\" is invalid." << endl;
				continue;
			}
			bam_index_t *idx;
			bam_plbuf_t *buf;
			idx = bam_index_load(argv[i]); // load BAM index
			if (idx == 0) {
				cerr << "Error: BAM indexing file is not available." << endl;
				exit(1);
			}
			bam_fetch(in->x.bam, idx, ref, newBeg, newEnd, &readsByName, fetch_func);
			bam_index_destroy(idx);
			
			for (list<DEIsoM_MRNA>::iterator ik = ii->mRNAs.begin();
				 ik != ii->mRNAs.end(); ik++) {
				for (map<string, list<DEIsoM_Read> >::iterator ij = readsByName.begin();
					 ij != readsByName.end(); ij++) {
					if (ij->second.size() == 2) {
						list<DEIsoM_Read>::iterator read1 = ij->second.begin();
						list<DEIsoM_Read>::iterator read2 = ij->second.begin();
						read2++;
						if (read1->doesAlignTo(*ik) &&
							read2->doesAlignTo(*ik)) {
							int iLen = DEIsoM_GetInsertedLength(*read1, *read2, *ik);
							insertedLengths.push_back(iLen);
						}
					}
				}
			}
		}
		samclose(in);
	}
	double sum = 0;
	for (list<int>::iterator ii = insertedLengths.begin();
		 ii != insertedLengths.end(); ii++) {
		sum += *ii;
	}
	double mean = sum / insertedLengths.size();
	double var = 0;
	for (list<int>::iterator ii = insertedLengths.begin();
		 ii != insertedLengths.end(); ii++) {
		var += ((*ii) - mean)*((*ii) - mean);
	}
	double stderivation = sqrt(var/insertedLengths.size());
	cout << "Inserted length" << endl;
	cout << "Mean: " << mean << " std: " << stderivation << endl;
	return 0;
}