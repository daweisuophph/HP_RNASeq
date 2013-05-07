#include <boost/filesystem.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <string>
#include <iostream>
#include <list>
#include <vector>
#include "HP_Gff.h"

using namespace std;
using namespace boost::math;
namespace fs = boost::filesystem;


void tranverseDir(fs::path dir, list<HP_Gene> &genes) {
	if (fs::exists(dir)) {
		if (fs::is_directory(dir)) {
			fs::directory_iterator end_iter;
			for ( fs::directory_iterator dir_itr(dir);
				 dir_itr != end_iter;
				 ++dir_itr ) {
				tranverseDir(dir_itr->path(), genes);
			}
		} else {
			HP_Gff gff(dir.string());
			if (gff.genes.size() == 1) {
				genes.push_back(*gff.genes.begin());
			}
		}
	}
}

bool loadBetas(ifstream &ifs, vector<vector<double> > &betas) {
	bool isFinished;
	int numOfIsos, numOfSubs;
	int numOfReads;
	ifs.read((char *) &isFinished, sizeof(bool));
	if (*(char*)&isFinished == '#') {
		cerr << "Warning: human readable file not supported" << endl;
		return false;
	}
	if (!ifs) return false;
	ifs.read((char *) &numOfIsos, sizeof(int));
	if (!ifs) return isFinished;
	ifs.read((char *) &numOfSubs, sizeof(int));
	if (!ifs) return isFinished;
	/*
	for (int i = 0; i < numOfReads; i++ ) {
		ifs.read((char *) &numOfReads, sizeof(int));
		if (!ifs) return isFinished;
	}
	 */
	double *buf = new double[numOfSubs*numOfIsos];
	ifs.read((char *) buf, sizeof(double)*numOfSubs*numOfIsos);
	if (!ifs.good()) {
		delete buf;
		return isFinished;
	}
	betas = vector<vector<double> >(numOfSubs);
	for (int i = 0; i < numOfSubs; i++) {
		betas[i] = vector<double>(numOfIsos);
		for (int j = 0 ; j < numOfIsos; j++) {
			betas[i][j] = buf[i*numOfIsos+j];
		}
	}
	delete buf;	
	
	return isFinished;
}

double kl(vector<vector<double> > &b1, vector<vector<double> > &b2) {
	if (b1.size() != b2.size()) {
		return -1;
	}
	double kl = 0;
	for (int i = 0; i < b1.size(); i++) {
		
		double a1 = 0;
		double a2 = 0;
		if (b1[i].size() != b2[i].size()) {
			return -1;
		}
		for (int k = 0; k < b1[i].size(); k++) {
			a1 += b1[i][k];
			a2 += b2[i][k];
		}
		// kl(q(b1)||q(b2)) + kl(q(b2)||q(b1))
		for (int k = 0; k < b1[i].size(); k++) {
			kl += (b1[i][k]-b2[i][k])*(digamma(b1[i][k])-digamma(a1));
			kl += (b2[i][k]-b1[i][k])*(digamma(b2[i][k])-digamma(a2));
		}
	}
	return kl;
}


int main(int argc, char **argv) {
	if (argc != 4) {
		cerr << "Usage: " << argv[0] << " <indexed GFF> <group 1> <group 2>" << endl;
	}
	
	fs::path dir(fs::initial_path<fs::path>());
	dir = fs::system_complete(argv[1]);
	
	list<HP_Gene> genes;
	tranverseDir(dir, genes);
	
	
	for (list<HP_Gene>::iterator ii = genes.begin();
		 ii != genes.end(); ii++) {
		cout << ii->seqid << " " << ii->ID << " ";
		string path1 = string(argv[2]) + "/" + ii->seqid + "/" + ii->ID;
		string path2 = string(argv[3]) + "/" + ii->seqid + "/" + ii->ID;
		ifstream ifs1(path1.c_str());
		ifstream ifs2(path2.c_str());
		vector<vector<double> > b1, b2;
		if (ifs1) {
			bool isFinished = loadBetas(ifs1, b1);
			if (isFinished) {
				cout << " 1";
			} else {
				cout << " 2";
			}
		} else {
			cout << " 0";
		}
		if (ifs2) {
			bool isFinished  = loadBetas(ifs2, b2);
			if (isFinished) {
				cout << " 1";
			} else {
				cout << " 2";
			}
		} else {
			cout << " 0";
		}
		if (ifs1 && ifs2) {
			cout << " " << kl(b1, b2);
		} else {
			cout << " -1";
		}
		cout << endl;
	}
	return 0;
}