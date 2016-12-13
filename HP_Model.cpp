/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */
#include <cstdlib>
#include <sstream>
#include <limits>
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
   lbfgs_parameter_init(&lbfgsParam);
   lbfgsParam.max_step = 100;
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
  
	numReadsPossible.resize(numIsos);
	for (int isoInd = 0; isoInd < numIsos; isoInd++) {
		numReadsPossible[isoInd] = max(0, 
            mRNAs[isoInd].getLength() - param.readLen) + 1; 
	}

   maxEnd.resize(numIsos);
   for (int isoInd = 0; isoInd < numIsos; isoInd++) {
      maxEnd[isoInd] = max(mRNAs[isoInd].end, isoInd == 0? -1: maxEnd[isoInd-1]);
   }
}

struct _HP_TMPSTRUCT {
	HP_Model *model;
	int index;
};

HP_Read loadRead(bam1_t * b) {
   HP_Read read;
   read.pos = b->core.pos+1; // change from 0 based to 1 based
   read.name = bam1_qname(b);
   read.len = b->core.l_qseq;
   read.cigar = vector<int32_t>(b->core.n_cigar);
   uint8_t * hi = bam_aux_get_core(b, "HI");
   if (hi[0] == 'C' || hi[0] == 'c') {
      read.hi = *(int8_t *) (hi + 1);
   } else if (hi[0] == 's' || hi[0] == 'S') {
      read.hi = *(int16_t *) (hi + 2);
   } else {
      read.hi = *(int32_t *) (hi + 4);
   }
   uint32_t *cigar = bam1_cigar(b);
   for (int i = 0; i < b->core.n_cigar; i++ ) {
      read.cigar[i] = cigar[i];
   }
   read.pnext = b->core.mpos+1;
   return read;
}

// load reads and compute log score for each isoforms
void HP_Model::loadReadsAndComputeLogScore() {
	numSubs = param.bams.size();
   numReadsBySub.resize(numSubs);

   #pragma omp parallel for
   for (int subInd = 0; subInd < numSubs; subInd++) {
      char buffer[sizeof(int)*2+sizeof(double)];

		samfile_t *in = samopen(param.bams[subInd].c_str(), "rb", 0);
		if (in == 0) {
			cerr << "Error: cannot read BAM file \"" << 
           param.bams[subInd] << "\"." << endl;
			exit(1);
		}

	   ofstream of((param.outputDir+"/tmp"+to_string(static_cast<long long>(subInd))).c_str(), ifstream::binary);
      if (!of) {
         cerr << "Error: cannot create tmp file" << endl;
         exit(1);
      }

	   normal_distribution<double> pInsertLen = normal_distribution<double>(param.meanInsertedLen, param.stdInsertedLen);

      // for pair-end reads
      HP_Read read1, read2;
      bam1_t* b = bam_init1();
      int numReads = 0;
      int ret = 0;

      unordered_map<int32_t, vector<HP_Read> > readsByHI;
      while (ret >= 0) { 
         ret = bam_read1(in->x.bam, b);
         read1 = loadRead(b);
         if (ret < 0 || read1.name != read2.name) {
            unordered_map<int, double> scores;

            for (unordered_map<int32_t, vector<HP_Read> >::iterator ii = readsByHI.begin(); ii != readsByHI.end(); ii++) {
               if (ii->second.size() != 2) continue; // only consider paired
               
               HP_MRNA ref;
               ref.start = min(ii->second[0].pos, ii->second[1].pos);

               vector<HP_MRNA>::iterator jj = upper_bound(mRNAs.begin(), mRNAs.end(), ref);
               while (jj != mRNAs.begin()) {
                  jj--;

                  int isoInd = jj - mRNAs.begin();
                  if (maxEnd[isoInd] < ref.start) break;

                  if (!ii->second[0].doesAlignTo(*jj)) continue;
                  if (!ii->second[1].doesAlignTo(*jj)) continue;
                  if (!ii->second[0].isOverhangOK(param.overhangLen)) continue;
                  if (!ii->second[1].isOverhangOK(param.overhangLen)) continue;

                  int iLen = HP_GetInsertedLength(ii->second[0], ii->second[1], *jj);
                  double tmp = (iLen - (double) param.meanInsertedLen) /
                     param.stdInsertedLen;
                  double logScore = -0.5 * log(2 * M_PI)
                     - log((double)param.stdInsertedLen) - 0.5 * tmp * tmp
                     + log(1.0) - log(numReadsPossible[jj - mRNAs.begin()]);

                  // use the maximum score`
                  if (scores.find(isoInd) == scores.end()) {
                     scores[isoInd] = logScore; 
                  } else {
                     scores[isoInd] = max(scores[isoInd], logScore);
                  } 
               }
            }
            *((int*) buffer) = scores.size();
            of.write(buffer, sizeof(int));
            for (unordered_map<int, double>::iterator ii = scores.begin();  ii != scores.end(); ii++) {
               *((int*) buffer) = ii->first;
               of.write(buffer, sizeof(int));
               *((double*) buffer) = ii->second;
               of.write(buffer, sizeof(double));
            }
            numReads++;
            readsByHI.clear();
         }
         readsByHI[read1.hi].push_back(read1);
         read2 = read1;
      }
      bam_destroy1(b);
      samclose(in);
		of.close();
      numReadsBySub[subInd] = numReads;
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
      // init betas for this sub
      for (int isoInd = 0; isoInd < numIsos; isoInd++) {
         betasBySub[i][isoInd] = alpha[isoInd] 
            + numReadsBySub[i]/(double) numIsos;
      }
	}

	weightsBySub = vector<vector<double> > (numSubs, vector<double>(numIsos));
   partBoundBySub = vector<double> (numSubs);

	ss = vector<double>(numIsos);
	q = vector<double>(numIsos);
	g = vector<double>(numIsos);
}

// load data and initialize variables
bool HP_Model::preprocessing() {
   startTime = clock();
	// create directory to store output
	try {
		boost::filesystem::create_directories(param.outputDir);
	} catch (const boost::filesystem::filesystem_error& e) {
		cerr << "Warning: cannot create dicrectory \"" << param.outputDir << "\"."<< endl;
		return false;
	}
	cerr << "Loading mRNAs..." << endl;
	loadMRNAs();
	if (mRNAs.size() <= 1) {
		cerr << "Warning: number of isoforms is less than 2. Skipping..." << endl;
		return false;
	}
	cerr << "Loading read and computing log score..." << endl;
	loadReadsAndComputeLogScore();
	cerr << "Initializing variables..." << endl;
	initVariables();
	cerr << endl;
	return true;
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

	int ret =  lbfgs(numIsos, logAlpha, &fx, _evaluate, _progress, this, &lbfgsParam);
	
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

   logVBLB += partBoundBySub[m];

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

      #pragma omp parallel for
      for (int subInd = 0; subInd < numSubs; subInd++) {
         vector<double> &betas = betasBySub[subInd];
         char buffer[sizeof(int)*2+sizeof(double)];

         // inner loop
         for (int currInIter = 0; currInIter < param.numInIters; currInIter++) {
            vector<double> oldBetas(betasBySub[subInd]);

            double sumBetas = 0.0;
            for (int k = 0; k < numIsos; k++) {
               sumBetas += betas[k];
            }
            double digammaSumBetas = digamma(sumBetas);
            for (int k = 0; k < numIsos; k++) {
               weightsBySub[subInd][k] = digamma(betas[k]) - digammaSumBetas;
            }
            for (int k = 0; k < numIsos; k++) {
               betas[k] = alpha[k];
            }
            partBoundBySub[subInd] = 0.0; 

            ifstream is((param.outputDir+"/tmp"+to_string(static_cast<long long>(subInd))).c_str(), ifstream::binary);

            if (!is) {
               cerr << "Error: cannot read tmp file" << endl;
               exit(1);
            }

            while (is.good()) {
               is.read(buffer, sizeof(int));
               int len = *(int *) buffer;
               vector<int> isoInds;
               vector<double> weightedScores;
               vector<double> scores;

               for (int l = 0; l < len && is.good(); l++) {
                  is.read(buffer, sizeof(int) + sizeof(double));
                  int k = *(int *) buffer; 
                  double logScore = *(double *) (buffer + sizeof(int));
                  isoInds.push_back(k);
                  scores.push_back(logScore);
                  weightedScores.push_back(logScore + weightsBySub[subInd][k]);
               }
               double normalizer = logSumExp(weightedScores);

               for (int l = 0; l < isoInds.size(); l++) {
                  double r = exp(weightedScores[l] - normalizer);
                  partBoundBySub[subInd] += r * 
                     ((digamma(betas[isoInds[l]]) - digammaSumBetas)
                     + scores[l] - log(r));
                  betas[isoInds[l]] += r;
               }
            }
            is.close();

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
			of << betasBySub[m][k] / betaSum / numReadsPossible[k] * 1.0e9 << " ";
		}
		of << endl;
	}
	of << "# done";
}


void HP_Model::saveBinary(ofstream &of) {
   char buffer[100];
	// isFinished
	*(bool *) buffer = isFinished;
   of.write(buffer, sizeof(bool));
	
	// numOfIsos
	*(int *) buffer = numIsos;
	of.write(buffer, sizeof(int));
	
	// numOfSubs
	*(int *) buffer = numSubs;
   of.write(buffer, sizeof(int));
	
	// alpha
	for (int k = 0; k < numIsos; k++) {
		*(double *) buffer = alpha[k];
      of.write(buffer, sizeof(double));
	}
	
	// betas
	for (int m = 0; m < numSubs; m++) {
		for (int k = 0; k < numIsos; k++) {
			*(double *) buffer = betasBySub[m][k];
         of.write(buffer, sizeof(double));
		}
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
   ofstream ofFPKM((param.outputDir+"/output.fpkm").c_str());
   if (ofFPKM) {
      saveFPKM(ofFPKM);
   }
}
