/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
/*
 Split tasks
 */

#include <iostream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <list>
#include <sstream>
#include <boost/filesystem.hpp>
#include "DEIsoM_Param.h"
#include "DEIsoM_Gff.h"

using namespace std;
namespace fs = boost::filesystem;

typedef struct _Temp_Param {
	DEIsoM_Param param;
	int trunkSize;
	string mainPath;
} Temp_Param;

Temp_Param parseArg(int argc, char **argv)  {
	Temp_Param tmpParam;
	tmpParam.trunkSize = 1;
	DEIsoM_Param &param = tmpParam.param;
	
	bool errorFlag = false;
	// parse options
	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "--path")==0 && i+1<argc) {
			tmpParam.mainPath = string(argv[++i]);
			if (!boost::filesystem::exists(argv[i])){
				errorFlag = true;
				cerr << "Error: cannot find program \"" << argv[i] << "\"." << endl;
			}
		} else if (strcmp(argv[i], "--gff-dir")==0 && i+1<argc) {
			param.gff = string(argv[++i]);
			if (!boost::filesystem::exists(param.gff)){
				errorFlag = true;
				cerr << "Error: cannot access indexed GFF folder \"" << param.gff << "\"." << endl;
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
		} else if (strcmp(argv[i], "--trunk-size")==0 && i+1<argc) {
			tmpParam.trunkSize = atoi(argv[++i]);
			if (tmpParam.trunkSize == 0) {
				errorFlag = true;
				cerr << "Error: invalid trunk size \"" << argv[i] << "\"." << endl;
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
	
	if (tmpParam.mainPath.empty()) {
		errorFlag = true;
		cerr << "Error: main path is required" << endl;
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

	if (errorFlag) {
	// print usage
		cerr << "Usage: " << argv[0] << " <options>" << endl;
		cerr << "Options:" << endl;
		cerr << "--path <path to program>" << endl;
		cerr << "--gff-dir <indexed Gff folder>" << endl;
		cerr << "--bams <BAM filename 1>,<BAM filename 2>,..." << endl;
		cerr << "--output <output directory>" << endl;
		cerr << "--min-read <minimum number of reads>" << endl;
		cerr << "--read-len <read length>" << endl;
		cerr << "--overhang-len <overhang length>" << endl;
		cerr << "--trunk-size <trunk size>" << endl;
		cerr << "--paired-end <mean> <std>" << endl;
		cerr << "--in-iter <# of inner iters>" << endl;
		cerr << "--out-iter <# of outer iters>" << endl;
		cerr << "--read-count <read count filename>" << endl;
		cerr << "--human-readable" << endl;
		cerr << "--newton-update" << endl;
		
		exit(1);
	}
	
	return tmpParam;
}

typedef struct _Task {
	string geneID;
	string gff;
} Task;

void tranverseDir(fs::path dir, list<Task> &tasks) {
	if (fs::exists(dir)) {
		if (fs::is_directory(dir)) {
			fs::directory_iterator end_iter;
			for ( fs::directory_iterator dir_itr(dir);
				 dir_itr != end_iter;
				 ++dir_itr ) {
				tranverseDir(dir_itr->path(), tasks);
			}
		} else {
			DEIsoM_Gff gff(dir.string());
			if (gff.genes.size() == 1) {
				Task task;
				task.geneID = gff.genes.begin()->ID;
				task.gff = dir.string();
				tasks.push_back(task);
			}
		}
	}
}

int main(int argc, char **argv) {
	Temp_Param tmpParam = parseArg(argc, argv);
	fs::path dir(fs::initial_path<fs::path>());
	dir = fs::system_complete(tmpParam.param.gff);
	list<Task> tasks;
	tranverseDir(dir, tasks);
	
	list<Task>::iterator ii = tasks.begin();
	string &outputDir = tmpParam.param.outputDir;
	try {
		boost::filesystem::create_directories(outputDir+"/cluster_scripts/");
	} catch (const boost::filesystem::filesystem_error& e) {
		cerr << "Error: cannot create dicrectory \"" <<
		outputDir << outputDir << "/cluster_scripts/"<< endl;
		exit(1);
	}
	int taskID = 0;
	
	while (ii != tasks.end()) {
		stringstream sstm1, sstm2;
		sstm1 << ii->geneID;
		sstm2 << ii->gff;
		ii++;
		for (int i = 1; i < tmpParam.trunkSize && ii != tasks.end();
				i++, ii++) {
			sstm1 << "," << ii->geneID;
			sstm2 << "," << ii->gff;
		}
		stringstream sstm3;
		sstm3 << outputDir << "/cluster_scripts/task" << taskID << ".sh";
		ofstream script(sstm3.str().c_str());
		if (script) {
			script << tmpParam.mainPath;
			script << " --gene-ids " << sstm1.str();
			script << " --gffs " << sstm2.str();
			script << " --bams ";
			for (list<string>::iterator ij = tmpParam.param.bams.begin();
				 ij != tmpParam.param.bams.end(); ij++ ){
				if (ij == tmpParam.param.bams.begin()) {
					script << *ij;
				} else {
					script << "," << *ij;
				}
			}
			script << " --output " << tmpParam.param.outputDir;
			script << " --min-read " << tmpParam.param.minRead;
			script << " --read-len " << tmpParam.param.readLen;
			script << " --overhang-len " << tmpParam.param.overhangLen;
			if (!tmpParam.param.isSingleEnd) {
				script << " --paired-end " << tmpParam.param.meanInsertedLen
					<< " " << tmpParam.param.stdInsertedLen;
			}
			script << " --in-iter " << tmpParam.param.numInIters;
			script << " --out-iter " << tmpParam.param.numOutIters;
			if (tmpParam.param.readCountFile.size()) {
				script << " --read-count " << tmpParam.param.readCountFile;
			}
			if (!tmpParam.param.outputBinary) {
				script << " --human-readable";
			}
			if (!tmpParam.param.bfgsUpdate) {
				script << " --newton-update";
			}
			script.close();
		}
		taskID++;
	}
	return 0;
}
