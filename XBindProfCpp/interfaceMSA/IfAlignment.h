#ifndef _IfAlignment_
#define _IfAlignment_

#include <map>
#include <vector>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include "../utils/FileTools.h"
#include "../utils/StringTools.h"

using namespace std;

class IfAlignment
{
private:
	map<int, int> resIDToIfIDA;
	map<int, int> resIDToIfIDB;
	vector<string> seqsA;
	vector<string> seqsB;
	char chainA;
	char chainB;
	vector<string> pdbIDs;
	vector<double> iscores;

public:
	IfAlignment(string file, double cutoff);

	string getSiteAASeq(int site, char chainID);
	string getSiteAASeq(int site, char chainID, double scoreCutoff);
	~IfAlignment();
};

#endif
