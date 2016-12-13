/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#ifndef _DEIsoM_MODEL
#define _DEIsoM_MODEL

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include<time.h>
#include <lbfgs.h>
#include "DEIsoM_Read.h"
#include "DEIsoM_Gene.h"
#include "DEIsoM_Param.h"

#define MAX_NEWTON_ITER 1000
using namespace std;

class DEIsoM_Model {
private:
	DEIsoM_Param param;

   // mRNAs
   vector<DEIsoM_MRNA> mRNAs;
   // max end before a given index
   vector<int> maxEnd;
	//sub->[readName, reads]
	vector<unordered_map<string, vector<DEIsoM_Read> > > readsByName;
	//sub->read->isoform->insertedLength
	vector<vector<vector<int> > > insertedLensBySub;
	//sub->read->isoform->match
	vector<vector<vector <bool> > > alignmentsBySub;
	//num of subjects
	int numSubs;
	//num of isoforms
	int numIsos;
	//num of reads by subject
	vector<int> numReadsBySub;
	//num of reads possible
	vector<int> numReadsPossible; 
	//isoform->alpha
	vector<double> alpha;
	//sub->isoform->beta
	vector<vector<double> >betasBySub;
	//isoform->(digamma(beta_k) - diagmma(sum(beta)))
	vector<vector<double> > weightsBySub;
	//isoform->sum(digamma(beta_k)-digamma(sum(beta)))
	vector<double> ss;
	//isoform->diagnal_Hessian
	vector<double> q;
	//isoform->gradient
	vector<double> g;
   // part of the lower bound ralted to reads
   vector<double> partBoundBySub;
   // param for lbfgs
   lbfgs_parameter_t lbfgsParam;

   // start time
   clock_t startTime;

	bool isFinished;

	static lbfgsfloatval_t _evaluate(
			void *instance,
			const lbfgsfloatval_t *x,
			lbfgsfloatval_t *g,
			const int n,
			const lbfgsfloatval_t step);

	lbfgsfloatval_t evaluate(
			const lbfgsfloatval_t *x,
			lbfgsfloatval_t *g,
			const int n,
			const lbfgsfloatval_t step);

	static int _progress(
			void *instance,
			const lbfgsfloatval_t *x,
			const lbfgsfloatval_t *g,
			const lbfgsfloatval_t fx,
			const lbfgsfloatval_t xnorm,
			const lbfgsfloatval_t gnorm,
			const lbfgsfloatval_t step,
			int n,
			int k,
			int ls);

	int progress(
			const lbfgsfloatval_t *x,
			const lbfgsfloatval_t *g,
			const lbfgsfloatval_t fx,
			const lbfgsfloatval_t xnorm,
			const lbfgsfloatval_t gnorm,
			const lbfgsfloatval_t step,
			int n,
			int k,
			int ls);

	void loadMRNAs();
   void loadReadsAndComputeLogScore();
	void initVariables();
	void updateRs(int subInd);
	void updateBetas(int subInd);
	void updateAlpha();
	void updateAlphaBFGS();
	void updateAlphaNewton();
	void saveHumanReadable(ofstream &of);
	void saveBinary(ofstream &of);
	void saveFPKM(ofstream &of);

public:
	DEIsoM_Model(const DEIsoM_Param &param);
	double computeLogVBLBBySub(int subInd);
	double computeLogVBLB();
	bool preprocessing();
	//perform Bayesian variational inference
	void performBVI();
	void save();
};
#endif
