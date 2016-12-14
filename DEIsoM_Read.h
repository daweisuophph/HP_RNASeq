/*
 Author: Hao Peng (pengh@purdue.edu)
 Date: May 8, 2013
 Version: 1.0v
 */
#ifndef _DEIsoM_READ
#define _DEIsoM_READ

#include <string>
#include <list>
#include <vector>
#include <boost/cstdint.hpp>
#include "DEIsoM_Gene.h"

using namespace std;

class DEIsoM_Read {
public:
	int pos; // 1 based start position of read
	string name; // read name
	int len; // read length
	vector<int32_t> cigar;//cigar sequence (high 28 bits is count, low 4 bits is the operation);
	
	DEIsoM_Read();
	string toString() const;
	bool doesAlignTo(const DEIsoM_MRNA &mRNA) const;
	bool isOverhangOK(int overhangLen) const;
	int getRelativePosOn(const DEIsoM_MRNA &mRNA) const;
};

int DEIsoM_GetInsertedLength(const DEIsoM_Read &read1, const DEIsoM_Read &read2, const DEIsoM_MRNA &mRNA);

#endif