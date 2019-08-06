#include "FileTools.h"

vector<string> FileTools::FileToStringList(string fileName)
{
	vector<string> lines;
	fstream f(fileName.c_str());
	if (!f)
	{
		cerr << "can't open file: " << fileName << endl;
	}
	string line;
	while (!f.eof())
	{
		getline(f, line);
		lines.push_back(line);
	}
	f.close();
	return lines;
}

vector<double> FileTools::FileToDoubleList(string fileName)
{
	vector<double> dataList;
	fstream f(fileName.c_str(),ios::in);
	if (!f)
	{
		cerr << "can't open file: " << fileName << endl;
	}
	double d;
	while (!f.eof())
	{
		f >> d;
		dataList.push_back(d);
	}
	f.close();
	return dataList;
}
