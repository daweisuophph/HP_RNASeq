#include <cstdlib>
#include <sstream>
#include <fstream>
#include <limits>
#include <cmath>
#include "HP_Model.h"
#include "HP_Gff.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/filesystem.hpp>
#include "sam.h"
#include "asa121.h"

using namespace std;
using namespace boost::math;

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
	betasBySub = vector<vector<double> >(numSubs);
	for (int i = 0; i < numSubs; i++) {
		betasBySub[i] = vector<double>(numIsos);
	}
	rsBySub = vector<vector<vector<double> > >(numSubs);
	for (int i = 0; i < numSubs; i++) {
		rsBySub[i] = vector<vector<double> >(numReadsBySub[i]);
		for (int j = 0; j < numReadsBySub[i]; j++) {
			rsBySub[i][j] = vector<double>(numIsos);
		}
	}
	weights = vector<double>(numIsos);
	weightedScore = vector<double>(numIsos);
	ss = vector<double>(numIsos);
	q = vector<double>(numIsos);
	g = vector<double>(numIsos);
}


// load data and initialize variables
bool HP_Model::preprocessing() {
	cerr << "Loading gene..." << endl;
	loadGene();
	if (gene.mRNAs.size() <= 1) {
		cerr << "Warning: number of isoforms is less than 2. Skipping..." << endl;
		return false;
	}
	cerr << "Loading read..." << endl;
	loadReads();
	cerr << "Computing alignments..." << endl;
	computeAlignments();
	if (alignmentsBySub.size() == 0) {
		cerr << "Warning: numbers of aligned reads for all subjects are below the threshold: " << param.minRead << ". skippping ..." << endl;
		return false;
	}
	cerr << "Computing log score..." << endl;
	computeLogScore();
	cerr << "Initializing variables..." << endl;
	initVariables();
	// create directory to store output
	param.outputDir += "/" + gene.seqid;
	try {
		boost::filesystem::create_directories(param.outputDir);
	} catch (const boost::filesystem::filesystem_error& e) {
		cerr << "Warning: cannot create dicrectory \"" << param.outputDir << "\"."<< endl;
		return false;
	}
	cerr << endl;
	return true;
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

static inline bool allClose(const vector<double> &a, const vector<double> &b, double thres) {
	if (a.size() != b.size()) return false;
	for (int i = 0; i < a.size(); i++) {
		if (abs(a[i]-b[i])>thres) {
			return false;
		}
	}
	return true;
}

static inline double logSumExp(vector<double> &a) {
	double maxA = -numeric_limits<double>::max();
	for (int i = 0 ; i < a.size(); i++ ){
		if (a[i] > maxA) {
			maxA = a[i];
		}
	}
	double sum = 0;
	for (int i = 0; i < a.size(); i++) {
		sum += exp(a[i]-maxA);
	}
	return maxA + log(sum);
	
}

void HP_Model::updateRs(int subInd) {
	vector<double> &betas = betasBySub[subInd];
	vector<vector<double> > &rs = rsBySub[subInd];
	double sumBetas = 0.0;
	for (int l = 0; l < numIsos; l++) {
		sumBetas += betas[l];
	}
	double digammaSumBetas = digamma(sumBetas);
	for (int k = 0; k < numIsos; k++) {
		weights[k] = digamma(betas[k]) - digammaSumBetas;
	}
	for (int n = 0; n < numReadsBySub[subInd]; n++) {
		for (int k = 0; k < numIsos; k++) {
			weightedScore[k] = logScore[subInd][n][k] + weights[k];
		}
		double normalizer = logSumExp(weightedScore);
		for (int k = 0; k < numIsos; k++) {
			rs[n][k] = exp(weightedScore[k] - normalizer);
		}
	}
}

void HP_Model::updateBetas(int subInd) {
	vector<double> &betas = betasBySub[subInd];
	vector<vector<double> > &rs = rsBySub[subInd];
	for (int k = 0; k < numIsos; k++) {
		betas[k] = alpha[k];
		for (int n = 0; n < numReadsBySub[subInd]; n++) {
			betas[k] += rs[n][k];
		}
	}
}


void HP_Model::updateAlpha() {
	for (int iter = 0; iter < MAX_NEWTON_ITER; iter++) {
		vector<double> oldAlpha = alpha;
		// compute ss
		for (int k = 0; k < numIsos; k++) {
			ss[k] = 0.0;
		}
		for (int m = 0; m < numSubs; m++) {
			double sumBetas = 0.0;
			for (int k = 0; k < numIsos; k++) {
				sumBetas += betasBySub[m][k];
			}
			double digammaSumBetas = digamma(sumBetas);
			for (int k = 0; k < numIsos; k++) {
				ss[k] += digamma(betasBySub[m][k])-digammaSumBetas;
			}
		}
		// compute q
		int fault;
		for (int k = 0; k < numIsos; k++) {
			q[k] = -numSubs*trigamma(alpha[k], &fault);
		}
		//compute g
		double sumAlpha = 0.0;
		for (int k = 0; k < numIsos; k++) {
			sumAlpha += alpha[k];
		}
		double digammaSumAlpha = digamma(sumAlpha);
		for (int k = 0; k < numIsos; k++) {
			g[k] = numSubs*(digammaSumAlpha-digamma(alpha[k]))+ss[k];
		}
		//compute b
		double z = numSubs*trigamma(sumAlpha, &fault);
		double b = 0.0;
		for (int k = 0; k < numIsos; k++) {
			b += g[k]/q[k];
		}
		double c = 0.0;
		for (int k = 0; k < numIsos; k++) {
			c += 1/q[k];
		}
		b /= 1/z+c;
		//finally update alpha
		for (int k = 0; k < numIsos; k++) {
			alpha[k] -= (g[k]-b)/q[k];
			if (alpha[k] < 1e-15) {
				cerr << "Warning: alpha too small." << endl;
				alpha[k] = 1e-10;
			}
		}
		if (allClose(oldAlpha, alpha, 1e-10)) {
			cerr << "Newtons converge in " << iter << " iterations." << endl;
			break;
		}
	}
}

void HP_Model::performBVI() {
	isFinished = false;
	for (int currOutIter = 0; currOutIter < param.numOutIters; currOutIter++) {
		vector<double> oldAlpha = alpha;
		cerr << "Outer iter " << currOutIter << "..." << endl;
		for (int subInd = 0; subInd < numSubs; subInd++) {
			// init betas for this sub
			for (int isoInd = 0; isoInd < numIsos; isoInd++) {
				betasBySub[subInd][isoInd] = alpha[isoInd] + numReadsBySub[subInd]/(double) numIsos;
			}
			for (int currInIter = 0; currInIter < param.numInIters; currInIter++) {
				vector<double> oldBetas = betasBySub[subInd];
				updateRs(subInd);
				updateBetas(subInd);
				if (allClose(betasBySub[subInd], oldBetas, 1e-15)) {
					cerr << "Subject " << (subInd+1) << " converged in " << (currInIter+1) << " iterations."<< endl;
					break;
				}
			}
		}
		save();
		if (numSubs == 1) {
			cerr << "Only one subject. Do not update alpha." << endl;
			break;
		}
		updateAlpha();
		if (allClose(oldAlpha, alpha, 1e-6)) {
			cerr << "Alpha converge in " << currOutIter << " iterations." << endl;
			break;
		}
	}
	isFinished = true;
	save();
	cout << "Finished succussfull." << endl;
}

void HP_Model::save() {
	ofstream of((param.outputDir+"/"+gene.ID).c_str());
	if (of) {
		of << "# gene" << endl << gene.ID << endl;
		of << "# numIsos" << endl << numIsos << endl;
		of << "# isoforms" << endl;
		for (list<HP_MRNA>::iterator ii = gene.mRNAs.begin();
			 ii != gene.mRNAs.end(); ii++) {
			of << ii->ID << " ";
		}
		of << endl;
		of << "# isFinished" << endl << isFinished << endl;
		of << "# numSubs" << endl << numSubs << endl;
		of << "# alpha" << endl;
		for (int k = 0; k < numIsos; k++ ){
			of << alpha[k] << " ";
		}
		of << endl;
		of << "# betas" << endl;
		for (int m = 0; m < numSubs; m++) {
			for (int k = 0; k < numIsos; k++) {
				of << betasBySub[m][k] << " ";
			}
			of << endl;
		}
		of << "# rs" << endl;
		for (int m = 0; m < numSubs; m++) {
			for (int n = 0; n < numReadsBySub[m]; n++) {
				for (int k = 0; k < numIsos; k++) {
					of << rsBySub[m][n][k] << ",";
				}
				of << ";";
			}
			of << endl;
		}
		of << "# done";
		of.close();
	}
}



