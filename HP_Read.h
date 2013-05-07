#ifndef _HP_READ
#define _HP_READ

#include <string>
#include <list>
#include <vector>
#include <boost/cstdint.hpp>
#include "HP_Gene.h"

using namespace std;

class HP_Read {
public:
	int pos; // 1 based start position of read
	string name; // read name
	int len; // read length
	vector<int32_t> cigar;//cigar sequence (high 28 bits is count, low 4 bits is the operation);
	
	HP_Read();
	string toString() const;
	bool doesAlignTo(const HP_MRNA &mRNA) const;
	bool isOverhangOK(int overhangLen) const;
	int getRelativePosOn(const HP_MRNA &mRNA) const;
};

int HP_GetInsertedLength(const HP_Read &read1, const HP_Read &read2, const HP_MRNA &mRNA);

#endif