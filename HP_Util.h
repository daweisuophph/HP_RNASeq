#ifndef _HP_UTIL
#define _HP_UTIL
#include <fstream>

using namespace std;
// test whether a file exists

bool inline fileExist(char *filename) {
	ifstream ifs(filename);
	if (ifs) {
		ifs.close();
		return true;
	} else {
		return false;
	}
}


#endif