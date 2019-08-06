
#include "StringTools.h"



string StringTools::trimString(string& s)
{
	if (s.empty()) {
		return s;
	}
	int begin=0, end=s.length()-1;
	for (int i = 0; i < s.length(); i++)
	{
		if (s.at(i) == ' ')
			begin++;
		else
			break;
	}
	for (int i = s.length() - 1; i >= 0; i--)
	{
		if (s.at(i) == ' ')
			end--;
		else
			break;
	}

	return s.substr(begin, end - begin + 1);
}

void StringTools::splitString(std::string& s, std::string delim, std::vector< std::string >* ret)
{
	size_t last = 0;
	size_t index = s.find_first_of(delim, last);
	while (index != std::string::npos)
	{
		if(index > last)
			ret->push_back(s.substr(last, index - last));
		last = index + 1;
		index = s.find_first_of(delim, last);
	}
	if (index > last)
	{
		ret->push_back(s.substr(last, index - last));
	}
}



