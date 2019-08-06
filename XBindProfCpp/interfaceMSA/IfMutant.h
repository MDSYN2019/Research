
#ifndef _IfMutant_
#define _IfMutant_

#include <vector>
#include <string>
#include "IfSingleMutant.h"
#include "../utils/StringTools.h"


using namespace std;

class IfMutant
{
public:
	vector<IfSingleMutant> mutList;
	int mutNum;

public:
	IfMutant(const string& line);
	~IfMutant();
};

#endif
