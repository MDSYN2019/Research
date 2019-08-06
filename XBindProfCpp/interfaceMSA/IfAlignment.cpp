#include "IfAlignment.h"


IfAlignment::IfAlignment(string file, double cutoff)
{
	vector<string> lines;

	fstream f(file.c_str());
	if (!f)
	{
		cerr << "can't open file" << file  << endl;
		exit(0);
	}
	string line;
	while (!f.eof())
	{
		getline(f, line);
		lines.push_back(line);
	}
	f.close();

	if (lines.size() < 8)
	{
		cerr << "invalid file: " << file << endl;
		exit(0);
	}
	this->chainA = lines.at(1).at(12);
	this->chainB = lines.at(3).at(12);

	stringstream ss;
	vector<string> spt;
	StringTools::splitString(lines.at(2), " ", &spt);

	unsigned int i;
	for (i = 1; i < spt.size(); i++)
	{
		int resID = atoi(spt[i].c_str());
		int pos = i - 1;
		this->resIDToIfIDA[resID] = pos;
	}
	spt.clear();

	StringTools::splitString(lines.at(4), " ", &spt);
	for (i = 1; i < spt.size(); i++)
	{
		int resID = atoi(spt[i].c_str());
		int pos = i - 1;
		this->resIDToIfIDB[resID] = pos;
	}
	spt.clear();

	for (i = 7; i < lines.size(); i++)
	{
		StringTools::splitString(lines.at(i), " ", &spt);
		if (spt.size() != 8)
			continue;
                double score = atof(spt.at(5).c_str());
                if(score < cutoff)
                {
                   spt.clear();
                   continue;
                }
		this->pdbIDs.push_back(spt.at(0));
		this->iscores.push_back(atof(spt.at(5).c_str()));
		this->seqsA.push_back(spt.at(6));
		this->seqsB.push_back(spt.at(7));
		spt.clear();
	}

//        cout << this->seqsA.size() << endl;
//	cout << "resA number: " << this->resIDToIfIDA.size() << endl;
//	cout << "resB number: " << this->resIDToIfIDB.size() << endl;
//	cout << "seq Num: " << this->seqsA.size() << endl;
//	for (i = 0; i < this->seqsA.size(); i++)
//	{
//		cout << seqsA.at(i) << " " << seqsB.at(i) << endl;
//	}

}

string IfAlignment::getSiteAASeq(int site, char chainID)
{
	map<int, int>::iterator it;
	int ifPos;
	string s("");
	if (chainID == this->chainA)
	{
		it = this->resIDToIfIDA.find(site);
		if (it == this->resIDToIfIDA.end())
			return s;
		ifPos = it->second;
		if (ifPos >= this->seqsA.at(0).length())
			return s;
		for (int i = 0; i < this->seqsA.size(); i++)
		{
			s.append(1, seqsA.at(i).at(ifPos));
		}
		return s;
	}
	else if (chainID == this->chainB)
	{
		it = this->resIDToIfIDB.find(site);
		if (it == this->resIDToIfIDB.end())
			return s;
		ifPos = it->second;
		if (ifPos >= this->seqsB.at(0).length())
			return s;
		for (int i = 0; i < this->seqsB.size(); i++)
		{
			s.append(1, seqsB.at(i).at(ifPos));
		}
		return s;
	}


	return string();
}

string IfAlignment::getSiteAASeq(int site, char chainID, double scoreCutoff)
{
	return string();
}

IfAlignment::~IfAlignment()
{
}

