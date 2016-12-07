/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#ifndef _HP_MODEL
#define _HP_MODEL

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <unordered_map>

#include<time.h>
#include <lbfgs.h>
#include "HP_Read.h"
#include "HP_Gene.h"
#include "HP_Param.h"

#define MAX_NEWTON_ITER 1000
#define MAX_BUFFER_SIZE 4096
using namespace std;

class HP_Model {
public:
	HP_Param param;

   // mRNAs
   vector<HP_MRNA> mRNAs;
	//sub->[readName, reads]
	vector<unordered_map<string, vector<HP_Read> > > readsByName;
	//sub->read->isoform->insertedLength
	vector<vector<vector<int> > > insertedLensBySub;
	//sub->read->isoform->match
	vector<vector<vector <bool> > > alignmentsBySub;
	//log score ln P(R,lamda|isoform) or ln P(R|isoform)
	//sub->read->isoform->score
	vector<vector<vector<double> > > logScore;
	//num of subjects
	int numSubs;
	//num of isoforms
	int numIsos;
	//total num of reads by subject (including reads mapped on other genes)
	vector<int> totalNumReadsBySub;
	//num of reads by subject
	vector<int> numReadsBySub;
	//num of reads possible
	vector<int> numReadsPossible; 
	//isoform->alpha
	vector<double> alpha;
	//sub->isoform->beta
	vector<vector<double> >betasBySub;
	//sub->read->isoform->r
	vector<vector<vector<double> > > rsBySub;
	//isoform->(digamma(beta_k) - diagmma(sum(beta)))
	vector<vector<double> > weightsBySub;
	//isoform->(log_score_k + weights_k)
	vector<vector<double> > weightedScoreBySub;
	//isoform->sum(digamma(beta_k)-digamma(sum(beta)))
	vector<double> ss;
	//isoform->diagnal_Hessian
	vector<double> q;
	//isoform->gradient
	vector<double> g;

   // start time
   clock_t startTime;

	bool isFinished;
	char *outputBuffer;

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
	void loadReads();
	void loadReadCount();
	void computeAlignments();
	void computeLogScore();
	void initVariables();
	void updateRs(int subInd);
	void updateBetas(int subInd);
	void updateAlpha();
	void updateAlphaBFGS();
	void updateAlphaNewton();
	void saveHumanReadable(ofstream &of);
	void saveBinary(ofstream &of);
	void saveFPKM(ofstream &of);

	void addRead(const HP_Read &read, int index);
	HP_Model(const HP_Param &param);
	~HP_Model();
	double computeLogVBLBBySub(int subInd);
	double computeLogVBLB();
	bool preprocessing();
	//perform Bayesian variational inference
	void performBVI();
	void save();
};
#endif
