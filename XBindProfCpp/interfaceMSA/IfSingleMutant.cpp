#include "IfSingleMutant.h"


#include <iostream>


IfSingleMutant::IfSingleMutant(string s)
{
	int len = s.length();
	this->natType = s.at(0);
	this->mutType = s.at(len - 1);
	this->chainID = s.at(1);
	string subs = s.substr(2, len - 3);
	this->resID = atoi(subs.c_str());

}

IfSingleMutant::~IfSingleMutant()
{
}

IfSingleMutant & IfSingleMutant::operator=(const IfSingleMutant & other)
{
	if (this == &other)
		return *this;
	natType = other.natType;
	mutType = other.mutType;
	chainID = other.chainID;
	resID = other.resID;
	return *this;
}

string IfSingleMutant::toString()
{
	char p[32] = { 0, };
	sprintf(p, "%c-%c-%d-%c", natType, chainID, resID, mutType);
	string s = p;
	
	return s;
}
