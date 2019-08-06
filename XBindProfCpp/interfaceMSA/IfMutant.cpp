#include "IfMutant.h"



IfMutant::IfMutant(const string& line)
{
	string s = line.substr(0, line.length()-1);
	vector<string> spt; 
	string x = ",";
	StringTools::splitString(s, x, &spt);
	this->mutNum = spt.size();
	for (int i = 0; i < mutNum; i++)
	{
		this->mutList.push_back(IfSingleMutant(spt.at(i)));
	}
}




IfMutant::~IfMutant()
{
}
