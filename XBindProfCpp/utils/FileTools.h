#ifndef _FileTools_
#define _FileTools_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

namespace FileTools
{
	vector<string> FileToStringList(string fileName);
	vector<double> FileToDoubleList(string fileName);
}


#endif
