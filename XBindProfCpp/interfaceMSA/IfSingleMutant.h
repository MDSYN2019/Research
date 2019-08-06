#ifndef _IfSingleMutant_
#define _IfSingleMutant_
#include <stdlib.h>
#include <stdio.h>
#include <string>
using namespace std;

class IfSingleMutant
{
public:
	int resID;
	char chainID;
	char natType;
	char mutType;

public:
	IfSingleMutant(string s);
	~IfSingleMutant();
	IfSingleMutant& operator=(const IfSingleMutant& other);
	string toString();


};

#endif
