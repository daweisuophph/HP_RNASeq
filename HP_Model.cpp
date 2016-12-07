/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */
#include <cstdlib>
#include <sstream>
#include <limits>
#include <thread>
#include <cmath>
#include <omp.h>
#include "HP_Model.h"
#include "HP_Gff.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/filesystem.hpp>
#include "sam.h"
#include "asa121.h"

//const int M_PI = acos(-1);

using namespace std;
using namespace boost::math;

HP_Model::HP_Model(const HP_Param &param) {
	this->param = param;
	readsByName = vector<unordered_map<string, vector<HP_Read> > >(param.bams.size());
	outputBuffer = new char[MAX_BUFFER_SIZE];
}

HP_Model::~HP_Model() {
	delete outputBuffer;
}

// load mRNAs from indexed GFF file
void HP_Model::loadMRNAs() {
	HP_Gff gff(param.gff);
   mRNAs = gff.mRNAs;
	if (mRNAs.size() <= 0) {
		cerr << "Erorr: number of mRNAs should be positive in\""
			<< param.gff << "\"." << endl;
		exit(1);
	}
	numIsos = mRNAs.size();
}

struct _HP_TMPSTRUCT {
	HP_Model *model;
	int index;
};

// load reads that intersect with the Gene.
void HP_Model::loadReads() {
   for (int subInd = 0; subInd < param.bams.size(); subInd++) {
		samfile_t *in = samopen(param.bams[subInd].c_str(), "rb", 0);
		if (in == 0) {
			cerr << "Error: cannot read BAM file \"" << 
           param.bams[subInd] << "\"." << endl;
			exit(1);
		}

      bam1_t* b=bam_init1();
      while (bam_read1(in->x.bam, b) >= 0) {
         HP_Read read;
         read.pos = b->core.pos+1; // change from 0 based to 1 based
         read.name = bam1_qname(b);
         read.len = b->core.l_qseq;
         read.cigar = vector<int32_t>(b->core.n_cigar);
         uint32_t *cigar = bam1_cigar(b);
         for (int i = 0; i < b->core.n_cigar; i++ ) {
            read.cigar[i] = cigar[i];
         }
	      readsByName[subInd][read.name].push_back(read);
      }

      bam_destroy1(b);
		samclose(in);
	}
}

// load the total read count of each subject
void HP_Model::loadReadCount() {
	ifstream ifs(param.readCountFile.c_str());
	if (ifs.is_open()) {
		while (ifs.good()) {
			int count;
			ifs >> count;
			totalNumReadsBySub.push_back(count);
		}
	} else {
		if (param.readCountFile.size()) {
			cerr << "Warning: file " << param.readCountFile << " not openned" << endl;
		}
		for (int i = 0; i < numSubs; i++) {
			totalNumReadsBySub.push_back(1);
		}
	}
}

void HP_Model::computeLogScore() {

	vector<int> isoLengths(numIsos);
	for (int i = 0; i < mRNAs.size(); i++) {
      isoLengths[i] = mRNAs[i].getLength(); 
   }
	
	numReadsPossible.resize(numIsos);
	for (int isoInd = 0; isoInd < numReadsPossible.size(); isoInd++) {
		numReadsPossible[isoInd] = max(0, isoLengths[isoInd] - param.readLen) + 1; 
	}

	normal_distribution<double> pInsertLen;
	
	if (!param.isSingleEnd) {
		pInsertLen = normal_distribution<double>(param.meanInsertedLen,
												 param.stdInsertedLen);
	}

	numSubs = alignmentsBySub.size();
	numReadsBySub = vector<int>(numSubs);
	logScore = vector<vector<vector<double> > >(numSubs);

	for (int subInd = 0; subInd < numSubs; subInd++) {
		numReadsBySub[subInd] = alignmentsBySub[subInd].size();
		logScore[subInd] = vector<vector<double> >(numReadsBySub[subInd]);
		for (int readInd = 0; readInd < numReadsBySub[subInd]; readInd++) {
			logScore[subInd][readInd] = vector<double>(numIsos);
			for (int isoInd = 0; isoInd < numIsos; isoInd++) {
				if (alignmentsBySub[subInd][readInd][isoInd]) {
					double logProbFrags = 0;
					if (!param.isSingleEnd) {
                  double tmp = insertedLensBySub[subInd][readInd][isoInd] - (double) param.meanInsertedLen;
						logProbFrags = -0.5*log(2*M_PI)-log((double)param.stdInsertedLen)-0.5*tmp*tmp/param.stdInsertedLen/param.stdInsertedLen;
					}
					logScore[subInd][readInd][isoInd] = log(1.0) - log(numReadsPossible[isoInd]) + logProbFrags;
					
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
	weightsBySub = vector<vector<double> > (numSubs, vector<double>(numIsos));
	weightedScoreBySub = vector<vector<double> > (numSubs, vector<double>(numIsos));
	ss = vector<double>(numIsos);
	q = vector<double>(numIsos);
	g = vector<double>(numIsos);
}


// load data and initialize variables
bool HP_Model::preprocessing() {
   startTime = clock();
	cerr << "Loading mRNAs..." << endl;
	loadMRNAs();
	if (mRNAs.size() <= 1) {
		cerr << "Warning: number of isoforms is less than 2. Skipping..." << endl;
		return false;
	}
	cerr << "Loading read..." << endl;
	loadReads();
	cerr << "Computing alignments..." << endl;
	computeAlignments();
	/*
	if (alignmentsBySub.size() == 0) {
		cerr << "Warning: numbers of aligned reads for all subjects are below the threshold: " << param.minRead << ". skippping ..." << endl;
		return false;
	}
	*/
	cerr << "Computing log score..." << endl;
	computeLogScore();
	cerr << "Loading read count..." << endl;
	loadReadCount();
	cerr << "Initializing variables..." << endl;
	initVariables();
	// create directory to store output
	try {
		boost::filesystem::create_directories(param.outputDir);
	} catch (const boost::filesystem::filesystem_error& e) {
		cerr << "Warning: cannot create dicrectory \"" << param.outputDir << "\"."<< endl;
		return false;
	}
	cerr << endl;
	return true;
}

void HP_Model::computeAlignments() {
	for (int ind = 0; ind < readsByName.size(); ind++) {
		vector<vector<bool> >  alignments;
		vector<vector<int> > insertedLens;
		for (unordered_map<string, vector<HP_Read> >::iterator ii = readsByName[ind].begin(); ii != readsByName[ind].end(); ii++) {
			if (param.isSingleEnd) {
				if (ii->second.size() != 1) {
					continue;
				}
			} else {
				if (ii->second.size() != 2) {
					continue;
				}
			}
			vector<bool> alignment = vector <bool>(mRNAs.size());
			vector<int> insertedLen = vector <int>(mRNAs.size());
			int isoformIndex = 0;
			bool allMismatch = true;
			for (vector<HP_MRNA>::iterator ij = mRNAs.begin();
				 ij != mRNAs.end(); ij++, isoformIndex++) {
				alignment[isoformIndex] = true;
				if (param.isSingleEnd) {
					vector<HP_Read>::iterator read = ii->second.begin();
					if (!read->doesAlignTo(*ij)) {
						alignment[isoformIndex] = false;
					}
					if (!read->isOverhangOK(param.overhangLen)) {
						alignment[isoformIndex] = false;
					}
				} else {
					vector<HP_Read>::iterator read1 = ii->second.begin();
					vector<HP_Read>::iterator read2 = read1;
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
		if (alignments.size() < param.minRead) {
			cerr << "num of reads aligned is less than the limit" << endl;
			exit(1);
		}
		alignmentsBySub.push_back(alignments);
		insertedLensBySub.push_back(insertedLens);
	}
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
		weightsBySub[subInd][k] = digamma(betas[k]) - digammaSumBetas;
	}
	for (int n = 0; n < numReadsBySub[subInd]; n++) {
		for (int k = 0; k < numIsos; k++) {
			weightedScoreBySub[subInd][k] = logScore[subInd][n][k] + weightsBySub[subInd][k];
		}
		double normalizer = logSumExp(weightedScoreBySub[subInd]);
		for (int k = 0; k < numIsos; k++) {
			rs[n][k] = exp(weightedScoreBySub[subInd][k] - normalizer);
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

lbfgsfloatval_t HP_Model::_evaluate(
		void *instance,
		const lbfgsfloatval_t *x,
		lbfgsfloatval_t *g,
		const int n,
		const lbfgsfloatval_t step) {
	return reinterpret_cast<HP_Model*>(instance)->evaluate(x, g, n, step);
}

lbfgsfloatval_t HP_Model::evaluate(
		const lbfgsfloatval_t *logAlpha,
		lbfgsfloatval_t *g,
		const int numIsos,
		const lbfgsfloatval_t step
		)
{
	lbfgsfloatval_t fx = 0.0;
	vector<double> alpha(numIsos);
	double sumAlphas = 0.0;
	double sumLogGammaAlphas = 0.0;
	for (int k = 0; k < numIsos; k++) {
		alpha[k] = exp(logAlpha[k]);
		sumAlphas += alpha[k];
		sumLogGammaAlphas += lgamma(alpha[k]);
	}
	double digammaSumAlphas = digamma(sumAlphas);

	// compute fx
	fx = numSubs * (sumLogGammaAlphas - lgamma(sumAlphas));

	for (int m = 0; m < numSubs; m++) {
		double sumBetas = 0.0;
		for (int k = 0; k < numIsos; k++) {
			sumBetas += betasBySub[m][k];
		}
		double digammaSumBetas = digamma(sumBetas);

		for (int k = 0; k < numIsos; k++) {
			fx += (alpha[k] - 1) * (digammaSumBetas - digamma(betasBySub[m][k]));
		}
	}

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

	// compute g
	for (int k = 0; k < numIsos; k++) {
		g[k] = (numSubs * (digamma(alpha[k]) - digammaSumAlphas) - ss[k]) * alpha[k];
	}

	return fx;
}

int HP_Model::_progress(
		void *instance,
		const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls) {
	return reinterpret_cast<HP_Model*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
}

int HP_Model::progress(
		const lbfgsfloatval_t *x,
		const lbfgsfloatval_t *g,
		const lbfgsfloatval_t fx,
		const lbfgsfloatval_t xnorm,
		const lbfgsfloatval_t gnorm,
		const lbfgsfloatval_t step,
		int n,
		int k,
		int ls
		)
{
   /*
	cerr << "Iteration " << k << ": ";
	cerr << "fx = " << fx <<  " x =";
	
	for (int i = 0; i < n; i++) {
		cerr << " " << x[i];
	}
	cerr << endl;
   */
	return 0;
}

void HP_Model::updateAlpha() {
   if (param.bfgsUpdate) {
      updateAlphaBFGS();
   } else {
      updateAlphaNewton();
   }
}

void HP_Model::updateAlphaNewton() {
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

void HP_Model::updateAlphaBFGS() {
	lbfgsfloatval_t fx;
	lbfgsfloatval_t *logAlpha = lbfgs_malloc(numIsos);

	for (int k = 0; k < numIsos; k++) {
		logAlpha[k] = log(alpha[k]);
	}

	int ret =  lbfgs(numIsos, logAlpha, &fx, _evaluate, _progress, this, NULL);
	
	for (int k = 0; k < numIsos; k++) {
		alpha[k] = exp(logAlpha[k]);
	}

	lbfgs_free(logAlpha);
}


double HP_Model::computeLogVBLBBySub(int m) {
	double logVBLB = 0.0;
	double sumAlphas = 0.0;
	for (int k = 0; k < numIsos; k++)
		sumAlphas += alpha[k];
	double logGammaSumAlphas = lgamma(sumAlphas);
	double sumLogGammaAlphas = 0.0;
	for (int k = 0; k < numIsos; k++)
		sumLogGammaAlphas += lgamma(alpha[k]);

	logVBLB += logGammaSumAlphas - sumLogGammaAlphas;
	double sumBetas = 0.0;
	for (int k = 0; k < numIsos; k++)
		sumBetas += betasBySub[m][k];
	double digammaSumBetas = digamma(sumBetas);
	for (int k = 0; k < numIsos; k++)
		logVBLB += (alpha[k] - 1.0) * 
			(digamma(betasBySub[m][k]) - digammaSumBetas);
	
	double logGammaSumBeta = lgamma(sumBetas);
	logVBLB -= logGammaSumBeta;
	for (int k = 0;k < numIsos; k++) 
		logVBLB += lgamma(betasBySub[m][k]) - (betasBySub[m][k] - 1.0) *
			(digamma(betasBySub[m][k]) - digammaSumBetas);

	for (int n = 0; n < numReadsBySub[m]; n++) {
		for (int k = 0; k < numIsos; k++) {
			if (rsBySub[m][n][k] >= 1e-50) {
				logVBLB += rsBySub[m][n][k] *
					( (digamma(betasBySub[m][k]) - digammaSumBetas)
					  + logScore[m][n][k]
					  - log(rsBySub[m][n][k]));
			}
		}
	}

	return logVBLB;
}

double HP_Model::computeLogVBLB() {
	double logVBLB = 0.0;
	for (int m = 0; m < numSubs; m++) {
		logVBLB += computeLogVBLBBySub(m);
	}

	return logVBLB;
}

void HP_Model::performBVI() {
	isFinished = false;
	for (int currOutIter = 0; currOutIter < param.numOutIters; currOutIter++) {
		if (currOutIter) {
			save(); // save intermediate steps
		}
		vector<double> oldAlpha = alpha;
		cerr << "Outer iter " << currOutIter << "..." << endl;
         
      vector<thread> threads;
      #pragma omp parallel for
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

		cerr << "VBLB: " << computeLogVBLB() << endl;
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
	cerr << "Finished succussfull." << endl;
}

void HP_Model::saveHumanReadable(ofstream &of) {
	of << "# isFinished" << endl << isFinished << endl;
	of << "# numIsos" << endl << numIsos << endl;
	of << "# isoforms" << endl;
	for (vector<HP_MRNA>::iterator ii = mRNAs.begin();
		 ii != mRNAs.end(); ii++) {
		of << ii->ID << " ";
	}
	of << endl;
	of << "# numSubs" << endl << numSubs << endl;
	of << "# totalNumReadsBySub" << endl;
	for (int m = 0; m < numSubs; m++ ){
		of << totalNumReadsBySub[m] << " ";
	}
	of << endl;
	of <<"# numReadsBySub" << endl;
	for (int m = 0; m < numSubs; m++ ){
		of << numReadsBySub[m] << " ";
	}
	of <<endl;
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
   /*
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
	of << "# log score" << endl;
	for (int m = 0; m < numSubs; m++) {
		for (int n = 0; n < numReadsBySub[m]; n++) {
			for (int k = 0; k < numIsos; k++) {
				of << logScore[m][n][k] << ",";
			}
			of << ";";
		}
		of << endl;
	}
   */
	of << "# done";
}


void HP_Model::saveFPKM(ofstream &of) {
   of << "# time" << endl << (clock()-startTime)/(double)CLOCKS_PER_SEC << endl; 
	of << "# isFinished" << endl << isFinished << endl;
	of << "# numIsos" << endl << numIsos << endl;
	of << "# isoforms" << endl;
	for (vector<HP_MRNA>::iterator ii = mRNAs.begin();
		 ii != mRNAs.end(); ii++) {
		of << ii->ID << " ";
	}
	of << endl;
	of << "# numSubs" << endl << numSubs << endl;
	of << "# totalNumReadsBySub" << endl;
	for (int m = 0; m < numSubs; m++ ){
		of << totalNumReadsBySub[m] << " ";
	}
	of << endl;
	of << "# numReadsBySub" << endl;
	for (int m = 0; m < numSubs; m++ ){
		of << numReadsBySub[m] << " ";
	}
	of <<endl;
	of << "# FPKM" << endl;
	for (int m = 0; m < numSubs; m++) {
		double betaSum = 0.0;
		for (int k = 0; k < numIsos; k++) {
			betaSum += betasBySub[m][k];
		}
		for (int k = 0; k < numIsos; k++) {
			of << betasBySub[m][k] / betaSum * numReadsBySub[m] / numReadsPossible[k] / totalNumReadsBySub[m] * 1.0e9 << " ";
		}
		of << endl;
	}
	of << "# done";
}


void HP_Model::saveBinary(ofstream &of) {
	int currSize = 0;
	// isFinished
	if (currSize + sizeof(int) >= MAX_BUFFER_SIZE) {
		of.write(outputBuffer, currSize);
		currSize = 0;
	}
	*(bool *) (outputBuffer+currSize) = isFinished;
	currSize += sizeof(bool);
	
	// numOfIsos
	if (currSize + sizeof(int) >= MAX_BUFFER_SIZE) {
		of.write(outputBuffer, currSize);
		currSize = 0;
	}
	*(int *) (outputBuffer+currSize) = numIsos;
	currSize += sizeof(int);
	
	// numOfSubs
	if (currSize + sizeof(int) >= MAX_BUFFER_SIZE) {
		of.write(outputBuffer, currSize);
		currSize = 0;
	}
	*(int *) (outputBuffer+currSize) = numSubs;
	currSize += sizeof(int);
	
	/*
	// numReads
	for (int m = 0; m < numSubs; m++) {
		if (currSize + sizeof(int) >= MAX_BUFFER_SIZE) {
			of.write(outputBuffer, currSize);
			currSize = 0;
		}
		*(int *) (outputBuffer+currSize) = numReadsBySub[m];
		currSize += sizeof(int);
	}
	 */
	// alpha
	for (int k = 0; k < numIsos; k++) {
		if (currSize + sizeof(double) >= MAX_BUFFER_SIZE) {
			of.write(outputBuffer, currSize);
			currSize = 0;
		}
		*(double *) (outputBuffer+currSize) = alpha[k];
		currSize += sizeof(double);
	}
	
	// betas
	for (int m = 0; m < numSubs; m++) {
		for (int k = 0; k < numIsos; k++) {
			if (currSize + sizeof(double) >= MAX_BUFFER_SIZE) {
				of.write(outputBuffer, currSize);
				currSize = 0;
			}
			*(double *) (outputBuffer+currSize) = betasBySub[m][k];
			currSize += sizeof(double);
		}
	}
	/*
	// rs
	for (int m = 0; m < numSubs; m++) {
		for (int n = 0; n < numReadsBySub[m]; n++) {
			for (int k = 0; k < numIsos; k++) {
				if (currSize + sizeof(double) >= MAX_BUFFER_SIZE) {
					of.write(outputBuffer, currSize);
					currSize = 0;
				}
				*(double *) (outputBuffer+currSize) = rsBySub[m][n][k];
				currSize += sizeof(double);
			}
		}
	}
	// log score
	for (int m = 0; m < numSubs; m++) {
		for (int n = 0; n < numReadsBySub[m]; n++) {
			for (int k = 0; k < numIsos; k++) {
				if (currSize + sizeof(double) >= MAX_BUFFER_SIZE) {
					of.write(outputBuffer, currSize);
					currSize = 0;
				}
				*(double *) (outputBuffer+currSize) = logScore[m][n][k];
				currSize += sizeof(double);
			}
		}
	}
	 */
	if (currSize) {
		of.write(outputBuffer, currSize);
	}
}


void HP_Model::save() {
	ofstream of((param.outputDir+"/output").c_str());
	if (of) {
		if (param.outputBinary) {
			saveBinary(of);
		} else {
			saveHumanReadable(of);
		}
		of.close();
	}
   if (param.readCountFile.size() > 0) {
      ofstream ofFPKM((param.outputDir+"/output.fpkm").c_str());
      if (ofFPKM) {
         saveFPKM(ofFPKM);
      }
   }
}



