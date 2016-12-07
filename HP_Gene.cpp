/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#include "HP_Gene.h"
#include <sstream>

using namespace std;

HP_Record::HP_Record() {
	seqid = string();
	source = string();
	type = string();
	start = 0;
	end = 0;
	score = 0;
    strand = '.';
    phase = -1;
	ID = string();
    name = string();
    parent = string();
}

string HP_Record::toString() const {
	stringstream sstm;
	sstm << seqid << "\t" << source << "\t" << type << "\t" << start << "\t" << end << "\t";
	if (score < 0) {
		sstm << ".";
	} else {
		sstm << score;
	}
	sstm << "\t" << strand << "\t";
	if (phase < 0) {
		sstm << ".";
	} else {
		sstm << phase;
	}
	sstm << "\t";
	bool hasAttribute = false;
	if (!ID.empty()) {
		sstm  << "ID=" << ID;
		hasAttribute = true;
	}
	if (!name.empty()) {
		if (hasAttribute) {
			sstm << ";";
		}
		sstm << "Name=" << name;
		hasAttribute = true;
	}
	if (!parent.empty()) {
		if (hasAttribute) {
			sstm << ";";
		}
		sstm << "Parent=" << parent;
		hasAttribute = true;
	}
	return sstm.str();
}

string HP_Exon::toString() const {
	stringstream sstm;
	sstm << HP_Record::toString() << endl;
	return sstm.str();
}

string HP_MRNA::toString() const{
	stringstream sstm;
	sstm << HP_Record::toString() << endl;
	for (vector<HP_Exon>::const_iterator ii = exons.begin(); ii != exons.end(); ii++) {
		sstm << (*ii).toString();
	}
	return sstm.str();
}

int HP_MRNA::getLength() const {
	int len = 0;
	for (vector<HP_Exon>::const_iterator ii = exons.begin(); ii != exons.end(); ii++) {
		len += ii->end-ii->start+1;
	}
	return len;
}

string HP_Gene::toString() const{
	stringstream sstm;
	sstm << HP_Record::toString() << endl;
	for (vector<HP_MRNA>::const_iterator ii = mRNAs.begin(); ii != mRNAs.end(); ii++) {
		sstm << (*ii).toString();
	}
	return sstm.str();
}

// return the possible bounds for this gene
// beg and end are 1 based
void HP_Gene::getBounds(int &beg, int &end) const {
	beg = -1;
	end = -1;
	for (vector<HP_MRNA>::const_iterator ii = mRNAs.begin();
		 ii != mRNAs.end(); ii++) {
		if (beg == -1 || ii->start < beg) {
			beg = ii->start;
		}
		if (end == -1 || ii->end > end) {
			end = ii->end;
		}
	}
}
