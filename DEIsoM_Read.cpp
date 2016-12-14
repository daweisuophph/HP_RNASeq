/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#include <sstream>
#include "DEIsoM_Read.h"

using namespace std;

DEIsoM_Read::DEIsoM_Read() {
	pos = -1;
	name = string();
	len = 0;
	cigar = vector<int32_t>();
}

string DEIsoM_Read::toString() const{
	stringstream sstm;
	sstm << name << " " << pos << " " << len  << " ";
	for (int i = 0; i < cigar.size(); i++) {
		int op = cigar[i] & 0xF;
		int count = cigar[i] >> 4;
		const char *opSam = "MIDNSHP=X";
		if (op < 9) {
			sstm << opSam[op] << count;
		}
	}
	return sstm.str();
}

bool DEIsoM_Read::doesAlignTo(const DEIsoM_MRNA &mRNA) const {
	int leftLen = len;
	int curr = pos;
	vector<int32_t> isoCigar;
	for (list<DEIsoM_Exon>::const_iterator ii = mRNA.exons.begin();
		 ii != mRNA.exons.end() && leftLen > 0; ii++) {
		if (curr < ii->start) {
			int count = ii->start - curr;
			int32_t op = (count << 4) + 3; //'N'
			isoCigar.push_back(op);
			curr = ii->start;
		}
		if (curr >= ii->start && curr <= ii->end) {
			int count = 0;
			if (ii->end+1-curr > leftLen) {
				count = leftLen;
			} else {
				count = ii->end+1-curr;
			}
			int32_t op = (count << 4) + 0; //'M'
			isoCigar.push_back(op);
			leftLen -= count;
			curr += count;
		}
	}
	
	 //cout << toString() << endl;
	 //cout << mRNA.toString() << endl;
	 
	if (isoCigar.size() != cigar.size()) {
		return false;
	} else {
		for (int i = 0; i < isoCigar.size(); i++) {
			if (isoCigar[i] != cigar[i]) {
				return false;
			}
		}
	}
	return true;
}

bool DEIsoM_Read::isOverhangOK(int overhangLen) const {
	for (int i = 0; i < cigar.size(); i++) {
		int op = cigar[i] & 0xF;
		if (op == '0') {
			if ((cigar[i] >> 4) < overhangLen) {
				return false;
			}
		}
	}
	return true;
}

int DEIsoM_Read::getRelativePosOn(const DEIsoM_MRNA &mRNA) const{
	int relPos = 1; //1-based
	for (list<DEIsoM_Exon>::const_iterator ii = mRNA.exons.begin();
		 ii != mRNA.exons.end(); ii++) {
		if (pos >= ii->start && pos <= ii->end) {
			relPos += pos - ii->start;
			return relPos;
		} else if (pos < ii->start) {
			return 0; // before the isoform
		} else {
			relPos += ii->end-ii->start+1;
		}
	}
	return relPos;
}

int DEIsoM_GetInsertedLength(const DEIsoM_Read &read1, const DEIsoM_Read &read2, const DEIsoM_MRNA &mRNA) {
	int relPos1 = read1.getRelativePosOn(mRNA);
	int relPos2 = read2.getRelativePosOn(mRNA);
	return abs(relPos1-relPos2)+read1.len;
}