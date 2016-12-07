/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
/*
Run the analysis for one gene
*/

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <boost/filesystem.hpp>
#include "HP_Model.h"

using namespace std;

HP_Param parseArg(int argc, char **argv)  {
	HP_Param param;
	bool errorFlag = false;
	// parse options
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--gff")==0 && i+1<argc) {
         param.gff = string(argv[++i]);
         if (!boost::filesystem::exists(param.gff)) {
            errorFlag = true;
            cerr << "Error: cannot access GFF file \"" << param.gff << "\"." << endl;
         }
		} else if (strcmp(argv[i], "--bams")==0 && i+1<argc) {
			char *filename = strtok(argv[++i], ",");
			do {
				param.bams.push_back(filename);
				if (!boost::filesystem::exists(filename)) {
					errorFlag = true;
					cerr << "Error: cannot access BAM file \"" << filename << "\"." << endl;
				}
			} while ((filename = strtok(NULL, ",")) != NULL);
		} else if (strcmp(argv[i], "--output")==0 && i+1<argc) {
			param.outputDir = string(argv[++i]);
			try {
				boost::filesystem::create_directories(param.outputDir);
			} catch (const boost::filesystem::filesystem_error& e) {
				errorFlag = true;
				cerr << "Error: cannot write to directory \"" << argv[i] << "\"." << endl;
			}
		} else if (strcmp(argv[i], "--min-read")==0 && i+1<argc) {
			param.minRead = atoi(argv[++i]);
			if (param.minRead < 0) {
				errorFlag = true;
				cerr << "Error: min # of reads must be non-negative" << endl;
			}
         //cerr << "Warning: --min-read not supported yet" << endl;
         //++i; //ignore this option (not supported yet)
		} else if (strcmp(argv[i], "--read-len")==0 && i+1<argc) {
			param.readLen = atoi(argv[++i]);
			if (param.readLen == 0) { 
				errorFlag = true;
				cerr << "Error: invalid read length \"" << argv[i] << "\"." << endl;
			}
		} else if (strcmp(argv[i], "--overhang-len")==0 && i+1<argc) {
			param.overhangLen = atoi(argv[++i]);
			if (param.overhangLen == 0) { 
				errorFlag = true;
				cerr << "Error: invalid overhang length \"" << argv[i] << "\"." << endl;
			}
		} else if (strcmp(argv[i], "--paired-end") == 0 && i+2<argc) {
			param.isSingleEnd = 0;
			param.meanInsertedLen = atof(argv[++i]); // now it does not check for error
			param.stdInsertedLen = atof(argv[++i]); // not it does not check for error
		} else if (strcmp(argv[i], "--in-iter")==0 && i+1<argc) {
			param.numInIters = atoi(argv[++i]);
			if (param.numInIters == 0) { 
				errorFlag = true;
				cerr << "Error: invalid number of inner iterations \"" << argv[i] << "\"." << endl;
			}
		} else if (strcmp(argv[i], "--out-iter")==0 && i+1<argc) {
			param.numOutIters = atoi(argv[++i]);
			if (param.numOutIters == 0) { 
				errorFlag = true;
				cerr << "Error: invalid number of outer iterations \"" << argv[i] << "\"." << endl;
			}
		} else if (strcmp(argv[i], "--read-count")==0 && i+1<argc) {
			param.readCountFile = string(argv[++i]);
		} else if (strcmp(argv[i], "--human-readable")==0) {
			param.outputBinary = false;
		} else if (strcmp(argv[i], "--newton-update")==0) {
			param.bfgsUpdate = false;
		} else {
			cerr << "Error: cannot recognize option " << i << " : \"" << argv[i] << "\"." << endl;
			errorFlag = true;
		}
	}

	if (param.readLen == 0) {
		errorFlag = true;
		cerr << "Error: read length is required" << endl;
	}
	if (param.bams.empty()) {
		errorFlag = true;
		cerr << "Error: BAM filenames are required" << endl;
	}
	if (param.outputDir.empty()) {
		errorFlag = true;
		cerr << "Error: Output directory is required" << endl;
	}
	// print usage
	if (errorFlag) {
		cerr << "Usage: " << argv[0] << " <options>" << endl;
		cerr << "Options:" << endl;
		cerr << "--gff <GFF filename>" << endl;
		cerr << "--bams <BAM filename 1>,<BAM filename 2>,..." << endl;	
		cerr << "--output <output directory>" << endl;
		cerr << "--min-read <minimum number of reads>" << endl;
		cerr << "--read-len <read length>" << endl;
		cerr << "--overhang-len <overhang length>" << endl;	
		cerr << "--paired-end <mean> <std>" << endl;	
		cerr << "--in-iter <# of inner iters>" << endl;	
		cerr << "--out-iter <# of outer iters>" << endl;
		cerr << "--read-count <filename of the read count file>" << endl;
		cerr << "--human-readable" << endl;
		cerr << "--newton-update" << endl;
		
		exit(1);
	}
	
	return param;
}

int main(int argc, char **argv) {
	HP_Param param = parseArg(argc, argv);
   HP_Model model(param);
   if (model.preprocessing()) {
      model.performBVI();
      model.save();
   }
	return 0;
}
