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
#include "DEIsoM_Gff.h"

using namespace std;
using namespace boost::math;
namespace fs = boost::filesystem;


void tranverseDir(fs::path dir, list<DEIsoM_Gene> &genes) {
	if (fs::exists(dir)) {
		if (fs::is_directory(dir)) {
			fs::directory_iterator end_iter;
			for ( fs::directory_iterator dir_itr(dir);
				 dir_itr != end_iter;
				 ++dir_itr ) {
				tranverseDir(dir_itr->path(), genes);
			}
		} else {
			DEIsoM_Gff gff(dir.string());
			if (gff.genes.size() == 1) {
				genes.push_back(*gff.genes.begin());
			}
		}
	}
}

bool loadBinary(ifstream &ifs, 
      vector<double> &alpha, 
      vector<vector<double> > &betas) {
	bool isFinished = false;
	int numOfIsos = 0, numOfSubs = 0;
	ifs.read((char *) &isFinished, sizeof(bool));
	if (!ifs) return false;
	ifs.read((char *) &numOfIsos, sizeof(int));
	if (!ifs) return isFinished;
	ifs.read((char *) &numOfSubs, sizeof(int));
	if (!ifs) return isFinished;

   int size = numOfIsos + numOfSubs * numOfIsos;

	double *buf = new double[size];
	ifs.read((char *) buf, sizeof(double)*size);
	if (!ifs.good()) {
		delete buf;
		return isFinished;
	}

   alpha = vector<double> (numOfIsos);
   for (int i = 0; i < numOfIsos; i++) {
      alpha[i] = buf[i];
   }

	betas = vector<vector<double> >(numOfSubs);
	for (int i = 0; i < numOfSubs; i++) {
		betas[i] = vector<double>(numOfIsos);
		for (int j = 0 ; j < numOfIsos; j++) {
			betas[i][j] = buf[i*numOfIsos+j + numOfIsos];
		}
	}
	delete buf;

	return isFinished;
}

bool loadHuman(ifstream &ifs, 
      vector<double> &alpha,
      vector<vector<double> > &betas) {
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
      } else if (strcmp(buf, "# alpha") == 0) {
			alpha = vector<double>(numOfIsos);
			for (int i = 0; i < numOfIsos; i++) {
				if (!ifs.good()) return false;
				ifs >> alpha[i];
			}
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

bool loadModel(ifstream &ifs, 
      vector<double> &alpha,
      vector<vector<double> > &betas) {
   char c;
	ifs.read((char *) &c, sizeof(char));
   ifs.putback(c);
   if (c == '#') {
      return loadHuman(ifs, alpha, betas);
   } else {
      return loadBinary(ifs, alpha, betas);
   }
}

double betaKL(vector<vector<double> > &b1, vector<vector<double> > &b2) {
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

double alphaKL(vector<double> &a1, vector<double> &a2) {
	if (a1.size() != a2.size()) {
		return -1;
	}
	double kl = 0;
	double s1 = 0;
	double s2 = 0;

	for (int k = 0; k < a1.size(); k++) {
		s1 += a1[k];
		s2 += a2[k];
	}
	// kl(q(a1)||q(a2)) + kl(q(a2)||q(a1))
	for (int k = 0; k < a1.size(); k++) {
      kl += (a1[k]/s1) * log((a1[k]/s1) / (a2[k]/s2));
		kl += (a2[k]/s2) * log((a2[k]/s2) / (a1[k]/s1));
	   //kl += (a1[k]-a2[k])*(digamma(a1[k])-digamma(s1));
		//kl += (a2[k]-a1[k])*(digamma(a2[k])-digamma(s2));
	}
	return kl;
}


int main(int argc, char **argv) {
	if (argc != 4 && argc != 5) {
		cerr << "Usage: " << argv[0] << " [--beta] <indexed GFF> <group 1> <group 2>" << endl;
      return 1;
	}

   bool computeBetaKL = false;
   char *indexDir, *outputDir1, *outputDir2;

   if (argc == 5) {
      if (strcmp(argv[1], "--beta") != 0) {
         cerr << "Unknown argument: " <<  argv[1] << endl;
         cerr << "Usage: " << argv[0] << " [--beta] <indexed GFF> <group 1> <group 2>" << endl;
         return 1;
      }
      computeBetaKL = true;
      indexDir = argv[2];
      outputDir1 = argv[3];
      outputDir2 = argv[4];
   } else {
      indexDir = argv[1];
      outputDir1 = argv[2];
      outputDir2 = argv[3];
   }
	
	fs::path dir(fs::initial_path<fs::path>());
	dir = fs::system_complete(indexDir);
	
	list<DEIsoM_Gene> genes;
	tranverseDir(dir, genes);
	
	
	for (list<DEIsoM_Gene>::iterator ii = genes.begin();
		 ii != genes.end(); ii++) {
		cout << ii->seqid << " " << ii->ID << " ";
		string path1 = string(outputDir1) + "/" + ii->seqid + "/" + ii->ID;
		string path2 = string(outputDir2) + "/" + ii->seqid + "/" + ii->ID;
		ifstream ifs1(path1.c_str());
		ifstream ifs2(path2.c_str());
      vector<double> a1, a2; //alphas
		vector<vector<double> > b1, b2; //betas
      bool loaded1 = false, loaded2 = false;
      if (ifs1) {
			bool isFinished = loadModel(ifs1, a1, b1);
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
			bool isFinished  = loadModel(ifs2, a2, b2);
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
         if (computeBetaKL) {
            cout << " " << betaKL(b1, b2);
         } else {
            cout << " " << alphaKL(a1, a2);
         }
		} else {
			cout << " -1";
		}
		cout << endl;
	}
	return 0;
}
