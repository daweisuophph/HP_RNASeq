#include <cstdlib>
#include <sstream>
#include <limits>
#include <cmath>
#include "HP_Model.h"
#include "HP_Gff.h"
#include <boost/math/distributions/normal.hpp>
#include "sam.h"

using namespace std;
using namespace boost::math;

HP_Param::HP_Param() {
	isSingleEnd = true;
	meanInsertedLen = 0;
	stdInsertedLen = 0;
	minRead = 1;
	readLen = 0;
	geneID = string();
	gff = string();
	bams = list<string>();
	overhangLen = 1;
	outputDir = string();
	numInIters = 5;
	numOutIters = 5000;
}

// to string
string HP_Param::toString() const {
	stringstream sstm;
	sstm << "Gene ID: " << geneID << endl;
	sstm << "GFF filename: " << gff << endl;
	sstm << "Number of subjects: " << bams.size() << endl;
	sstm << "BAM filename: ";
	for (list<string>::const_iterator ii = bams.begin(); ii != bams.end(); ii++)
		sstm << (*ii) << " ";
	sstm << endl;
	sstm << "Output directory: " << outputDir << endl;
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
	return sstm.str();
}

HP_Model::HP_Model(const HP_Param &param) {
	this->param = param;
	readsByName = vector<map<string, list<HP_Read> > >(param.bams.size());
}

// load gene from indexed GFF file
void HP_Model::loadGene() {
	HP_Gff gff(param.gff);
	if (gff.genes.size() != 1) {
		cerr << "Erorr: number of genes is not exactly one in GFF file \""
			<< param.gff << "\"." << endl;
		exit(1);
	}
	gene = *gff.genes.begin();
}

struct _HP_TMPSTRUCT {
	HP_Model *model;
	int index;
};

// callback for bam_fetch()
static int fetch_func(const bam1_t *b, void *data) {
	HP_Read read;
	read.pos = b->core.pos+1; // change from 0 based to 1 based
	read.name = bam1_qname(b);
	read.len = b->core.l_qseq;
	read.cigar = vector<int32_t>(b->core.n_cigar);
	uint32_t *cigar = bam1_cigar(b);
	for (int i = 0; i < b->core.n_cigar; i++ ) {
		read.cigar[i] = cigar[i];
	}
	//cout << read.toString() << endl;
	struct _HP_TMPSTRUCT *tmp = (struct _HP_TMPSTRUCT *) data;
	tmp->model->addRead(read, tmp->index);
    return 0;
}

// load reads that intersect with the Gene.
void HP_Model::loadReads() {
	int beg, end;
	gene.getBounds(beg, end);
	stringstream sstm;
	sstm << gene.seqid << ":" << beg << "-" << end;
	//cerr << "Range: " << sstm.str() << endl;
	int ind = 0;
	for (list<string>::iterator ii = param.bams.begin();
		 ii != param.bams.end(); ii++, ind++) {
		samfile_t *in = samopen(ii->c_str(), "rb", 0);
		if (in == 0) {
			cerr << "Error: cannot read BAM file \"" << (*ii) << "\"." << endl;
			exit(1);
		}
		int newBeg, newEnd, ref;
		bam_parse_region(in->header, sstm.str().c_str(),
						 &ref, &newBeg, &newEnd);
		//cout << ref << " " << newBeg << " " << newEnd << endl;
		if (ref < 0) {
			cerr << "Warning: region \"" << sstm.str() << "\" is invalid." << endl;
			return;
		}
		bam_index_t *idx;
        bam_plbuf_t *buf;
		idx = bam_index_load(ii->c_str()); // load BAM index
        if (idx == 0) {
            cerr << "Error: BAM indexing file is not available." << endl;
            exit(1);
        }
		struct _HP_TMPSTRUCT tmp;
		tmp.model = this;
		tmp.index = ind;
		bam_fetch(in->x.bam, idx, ref, newBeg, newEnd, &tmp, fetch_func);
        bam_index_destroy(idx);
		samclose(in);
	}
}

void HP_Model::computeLogScore() {
	vector<int> isoLengths(gene.mRNAs.size());
	int i = 0;
	for (list<HP_MRNA>::iterator ii = gene.mRNAs.begin();
		 ii != gene.mRNAs.end(); ii++, i++){
		isoLengths[i] = ii->getLength();
	}
	
	normal_distribution<double> pInsertLen(param.meanInsertedLen,
								   param.stdInsertedLen);
	
	numSubs = alignmentsBySub.size();
	numReadsBySub = vector<int>(numSubs);
	logScore = vector<vector<vector<double> > >(numSubs);
	
	list<list<vector<bool> > >::iterator iSub1 = alignmentsBySub.begin();
	list<list<vector<int> > >::iterator iSub2 = insertedLensBySub.begin();
	for (int subInd = 0; subInd < numSubs; subInd++, iSub1++, iSub2++) {
		numReadsBySub[subInd] = iSub1->size();
		logScore[subInd] = vector<vector<double> >(numReadsBySub[subInd]);
		list<vector<bool> >::iterator iRead1 = iSub1->begin();
		list<vector<int> >::iterator iRead2 = iSub2->begin();
		for (int readInd = 0; readInd < numReadsBySub[subInd]; readInd++, iRead1++, iRead2++) {
			numIsos = iRead1->size();
			logScore[subInd][readInd] = vector<double>(numIsos);
			for (int isoInd = 0; isoInd < numIsos; isoInd++) {
				if ((*iRead1)[isoInd]) {
					double logProbFrags = 0;
					if (!param.isSingleEnd) {
						logProbFrags = log(pdf(pInsertLen, (*iRead2)[isoInd]));
					}
					int numReadsPossible = isoLengths[isoInd]-param.readLen+1;
					logScore[subInd][readInd][isoInd] = log(1.0) - log(numReadsPossible) + logProbFrags;
					
				} else {
					logScore[subInd][readInd][isoInd] = -numeric_limits<double>::max();
				}
			}
		}
	}
}

void HP_Model::initVariables() {
	alpha = vector<double>(numIsos);
	for (int i = 0; i < numIsos; i++) {
		alpha[i] = 1.0;
	}
	betas = vector<vector<double> >(numSubs);
	for (int i = 0; i < numSubs; i++) {
		betas[i] = vector<double>(numIsos);
	}
	rs = vector<vector<vector<double> > >(numSubs);
	for (int i = 0; i < numSubs; i++) {
		rs[i] = vector<vector<double> >(numReadsBySub[i]);
		for (int j = 0; j < numReadsBySub[i]; j++) {
			rs[i][j] = vector<double>(numIsos);
		}
	}
}


// load data and initialize variables
void HP_Model::preprocessing() {
	cerr << "Loading gene..." << endl;
	loadGene();
	if (gene.mRNAs.size() <= 1) {
		cerr << "Warning: number of isoforms is less than 2. Skipping..." << endl;
		return;
	}
	cerr << "Loading read..." << endl;
	loadReads();
	cerr << "Computing alignments..." << endl;
	computeAlignments();
	if (alignmentsBySub.size() == 0) {
		cerr << "Warning: numbers of aligned reads for all subjects are below the threshold: " << param.minRead << ". skippping ..." << endl;
		return;
	}
	cerr << "Computing log score..." << endl;
	computeLogScore();
	cerr << "Initializing variables..." << endl;
	initVariables();
	cerr << endl;
}

void HP_Model::addRead(const HP_Read &read, int ind) {
	if (readsByName[ind].find(read.name) == readsByName[ind].end()) {
		list<HP_Read> &reads = readsByName[ind][read.name] = list<HP_Read>();
		reads.push_back(read);
	} else {
		readsByName[ind][read.name].push_back(read);
	}
}


void HP_Model::computeAlignments() {
	for (int ind = 0; ind < readsByName.size(); ind++) {
		list<vector<bool> >  alignments;
		list<vector<int> > insertedLens;
		for (map<string, list<HP_Read> >::iterator ii = readsByName[ind].begin();
			 ii != readsByName[ind].end(); ii++) {
			if (param.isSingleEnd) {
				if (ii->second.size() != 1) {
					continue;
				}
			} else {
				if (ii->second.size() != 2) {
					continue;
				}
			}
			vector<bool> alignment = vector <bool>(gene.mRNAs.size());
			vector<int> insertedLen = vector <int>(gene.mRNAs.size());
			int isoformIndex = 0;
			bool allMismatch = true;
			for (list<HP_MRNA>::iterator ij = gene.mRNAs.begin();
				 ij != gene.mRNAs.end(); ij++, isoformIndex++) {
				alignment[isoformIndex] = true;
				if (param.isSingleEnd) {
					list<HP_Read>::iterator read = ii->second.begin();
					if (!read->doesAlignTo(*ij)) {
						alignment[isoformIndex] = false;
					}
					if (!read->isOverhangOK(param.overhangLen)) {
						alignment[isoformIndex] = false;
					}
				} else {
					list<HP_Read>::iterator read1 = ii->second.begin();
					list<HP_Read>::iterator read2 = read1;
					read2++;
					if (!read1->doesAlignTo(*ij) ||
						!read2->doesAlignTo(*ij)) {
						alignment[isoformIndex] = false;
					}
					if (!read1->isOverhangOK(param.overhangLen)||
						!read2->isOverhangOK(param.overhangLen)) {
						alignment[isoformIndex] = false;
					}
					if (alignment[isoformIndex]) {
						int iLen = HP_GetInsertedLength(*read1, *read2, *ij);
						insertedLen[isoformIndex] = iLen;
					} else {
						insertedLen[isoformIndex] = 0;
					}
				}
				if (alignment[isoformIndex]) {
					allMismatch = false;
				}
			}
			
			/*
			 cout << ii->first << ": ";
			 for (int i = 0; i < alignment.size(); i++) {
			 cout << alignment[i] << " ";
			 }
			 cout << endl;
			 */
			/* remove reads that does not match any isoform */
			if (!allMismatch) {
				alignments.push_back(alignment);
				insertedLens.push_back(insertedLen);
			}
		}
		cerr << "Number of effective reads for subject " << (ind+1) << ": "
			<< alignments.size() << endl;
		if (alignments.size() >= param.minRead) {
			alignmentsBySub.push_back(alignments);
			insertedLensBySub.push_back(insertedLens);
		}
	}
	cerr << "Number of effective subjects: " << alignmentsBySub.size() << endl;
}

void HP_Model::performBVI() {
	
}



