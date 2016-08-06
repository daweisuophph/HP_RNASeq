/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
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

bool loadBetasBinary(ifstream &ifs, vector<vector<double> > &betas) {
	bool isFinished = false;
	int numOfIsos = 0, numOfSubs = 0;
	ifs.read((char *) &isFinished, sizeof(bool));
	if (!ifs) return false;
	ifs.read((char *) &numOfIsos, sizeof(int));
	if (!ifs) return isFinished;
	ifs.read((char *) &numOfSubs, sizeof(int));
	if (!ifs) return isFinished;

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

bool loadBetasHuman(ifstream &ifs, vector<vector<double> > &betas) {
	bool isFinished = false;
	int numOfIsos = 0, numOfSubs = 0;
   int num;
   bool flag = false;

   char *buf = new char[1000];
   while (ifs.good()) {
      ifs.getline(buf, 1000);
      if (strcmp(buf, "# isFinished") == 0) {
         ifs >> num;
         if (num == 1) isFinished = true;
      } else if (strcmp(buf, "# numIsos") == 0) {
         if (!ifs.good()) return false;
         ifs >> numOfIsos;
      } else if (strcmp(buf, "# numSubs") == 0) {
         if (!ifs.good()) return false;
         ifs >> numOfSubs;
      } else if (strcmp(buf, "# betas") == 0) {
         betas = vector<vector<double> >(numOfSubs, vector<double>(numOfIsos));
         for (int i = 0; i < numOfSubs; i++) {
            for (int j = 0; j < numOfIsos; j++) {
               if (!ifs.good()) return false;
               ifs >> betas[i][j];
            }
         }
         flag = true;
         break;
      }
   }

   if (!flag) return false;

	return isFinished;
}

bool loadBetas(ifstream &ifs, vector<vector<double> > &betas) {
   char c;
	ifs.read((char *) &c, sizeof(char));
   ifs.putback(c);
   if (c == '#') {
      return loadBetasHuman(ifs, betas);
   } else {
      return loadBetasBinary(ifs, betas);
   }
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
		//cerr << i << ": " << endl;
		//	cerr << "kl " << kl << endl;
		//cerr << "a1 : " << a1 << " a2 : " << a2 << endl;
		for (int k = 0; k < b1[i].size(); k++) {
			//cerr << "b1 : " << b1[i][k] << " b2 : " << b2[i][k] << endl;
			kl += (b1[i][k]-b2[i][k])*(digamma(b1[i][k])-digamma(a1));
			kl += (b2[i][k]-b1[i][k])*(digamma(b2[i][k])-digamma(a2));
			//cerr << (b1[i][k]-b2[i][k])*(digamma(b1[i][k])-digamma(a1)) << " ";
			//cerr << (b2[i][k]-b1[i][k])*(digamma(b2[i][k])-digamma(a2));
			//cerr << endl;
			//cerr << "kl " << kl << endl;
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
      bool loaded1 = false, loaded2 = false;
      if (ifs1) {
			bool isFinished = loadBetas(ifs1, b1);
			if (isFinished) {
            loaded1 = true;
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
            loaded2 = true;
				cout << " 1";
			} else {
				cout << " 2";
			}
		} else {
			cout << " 0";
		}
		if (loaded1 && loaded2) {
			cout << " " << kl(b1, b2);
		} else {
			cout << " -1";
		}
		cout << endl;
	}
	return 0;
}
