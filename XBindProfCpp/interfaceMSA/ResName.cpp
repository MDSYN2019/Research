#include "ResName.h"



ResName::ResName()
{
        this->aaSeq = "ACDEFGHIKLMNPQRSTVWY";
	for (int i = 0; i < 20; i++)
	{
		char c = this->aaSeq[i];
		this->sinToIntMap[c] = i;
	}
}

ResName & ResName::operator=(const ResName & other)
{
	if (this == &other)
	{
		return *this;
	}

	sinToIntMap = other.sinToIntMap;
	return *this;
	// TODO: 在此处插入 return 语句
}

char ResName::intToSin(int i)
{
	if (i > 19 || i < 0)
		return 'X';
	return this->aaSeq[i];
}

int ResName::sinToInt(char c)
{
	map<char, int>::iterator p = this->sinToIntMap.find(c);
	if (p != this->sinToIntMap.end())
	{
		return sinToIntMap[p->first];
	}
	return -1;
}


ResName::~ResName()
{
}
