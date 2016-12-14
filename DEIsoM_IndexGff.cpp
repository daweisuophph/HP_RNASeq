/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: June 22, 2016
 Version: 1.1v
 */
#include <string>
#include <cstdlib>
#include <iostream>
#include "DEIsoM_Gff.h"
#include "DEIsoM_Gene.h"
#include <boost/filesystem.hpp>

using namespace std;


void parseArg(int argc, char **argv, string &gffFile, string &outputDir) {
	bool errorFlag = false;
	if (argc != 3)
		errorFlag = true;
	else {
		gffFile = string(argv[1]);
		if (!boost::filesystem::exists(argv[1])) {
			errorFlag = true;
			cerr << "Error: cannot access GFF file \"" << argv[1] << "\"." << endl;
		}
		outputDir = string(argv[2]);
		try {
			boost::filesystem::create_directories(outputDir);
		} catch (const boost::filesystem::filesystem_error& e) {
			errorFlag = true;
			cerr << "Error: cannot create directory \"" << argv[2] << "\"." << endl;
		}
	}
	
	if (errorFlag) {
		cerr << "Usage: " << argv[0] << " <GFF filename> <output directory>" << endl;
		exit(1);
	}
}


int main (int argc, char **argv) {
	string gffFile;
	string outputDir;
	parseArg(argc, argv, gffFile, outputDir);
	cerr << "Loading GFF file" << endl;
	DEIsoM_Gff gff(gffFile);
	cerr << "Finished loading" << endl;
	/*
	// print all genes
	for (list<Gene>::iterator ii = gff.genes.begin(); ii != gff.genes.end(); ii++) {
		cout << (*ii).toString() << endl;
	}
	*/
	map<string, list<DEIsoM_Gene> > genesBySeqid;
	gff.getGenesBySeqid(genesBySeqid);
	cerr << "Output indexed GFF files" << endl;
	for (map<string, list<DEIsoM_Gene> >::iterator ii = genesBySeqid.begin();
		 ii != genesBySeqid.end(); ii++) {
		string chr = ii->first;
		try {
			boost::filesystem::create_directories(outputDir+"/"+chr);
		} catch (const boost::filesystem::filesystem_error& e) {
			cerr << "Error: cannot create dicrectory \"" <<
				outputDir << "/" << chr << "\"."<< endl;
			exit(1);
		}
		for (list<DEIsoM_Gene>::iterator ij = ii->second.begin();
			 ij != ii->second.end(); ij++) {
			ofstream geneFile((outputDir+"/"+chr+"/"+ij->ID+".gff3").c_str());
			if (geneFile) {
				geneFile << "##gff-version 3" << endl;
				geneFile << ij->toString();
				geneFile.close();
			} else {
				cerr << "Warning: cannot write to file \"" <<
					outputDir << "/" << chr << "/" << ij->ID << "\"." << endl;
			}
		}
	}
	
	return 0;
}
