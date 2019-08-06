#ifndef _ResName_
#define _ResName_

#include <stdlib.h>
#include <string>
#include <map>
using namespace std;

class ResName
{
public:
	string aaSeq;
	map<char, int> sinToIntMap;
public:
	ResName();
	ResName& operator=(const ResName& other);
	char intToSin(int i);
	int sinToInt(char c);
	~ResName();


};

#endif
