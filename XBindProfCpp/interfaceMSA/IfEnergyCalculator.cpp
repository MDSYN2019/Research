#include "IfEnergyCalculator.h"

double IfEnergyCalculator::getEnergy(string posAASeq, char aaType, double w1, double w2 ,double w3)
{
	/*
	w1: fix pseudo count number
	w2: evo pseudo count number
	w3: blank pseudo count number
	*/
	int blankNum = 0;
	int aaNum = 0;
	int count[20] = { 0 };
	ResName rn;
	int intType = rn.sinToInt(aaType);
	for (int i = 0; i < posAASeq.length(); i++)
	{
		

		int type = rn.sinToInt(posAASeq.at(i));
		if (type >= 0 && type < 20)
		{
			count[type] ++;
			aaNum++;
		}
		else
			blankNum++;
	}

	AAProbabilityArray pa(count);

	SubstitutionMatrix sub;
	AAProbabilityArray pb = pa.applySubstitution(sub);

/*
	cout << "count Num: " << count[intType] << endl;
	cout << "evoNum: " << w2*pb.getValue(intType) << endl;
	cout << "blankNum: " << w3*blankNum << endl;
*/
	double e = -15.0*log(count[intType] + w2*pb.getValue(intType) + w1 + blankNum*w3);
	return e;
}

double IfEnergyCalculator::getEnergy(string posAASeq, char aaType)
{
	return getEnergy(posAASeq, aaType, 25.0, 5.0, 15.0);
}

double IfEnergyCalculator::getDDG(IfAlignment align, IfSingleMutant mut)
{
	string aaSeq = align.getSiteAASeq(mut.resID, mut.chainID);
	char natType = mut.natType;
	char mutType = mut.mutType;
	return getEnergy(aaSeq, mutType) - getEnergy(aaSeq, natType);
}

double IfEnergyCalculator::getDDG(IfAlignment align, IfSingleMutant mut, double w1, double w2, double w3)
{
	string aaSeq = align.getSiteAASeq(mut.resID, mut.chainID);
	char natType = mut.natType;
	char mutType = mut.mutType;
	return getEnergy(aaSeq, mutType, w1, w2, w3) - getEnergy(aaSeq, natType, w1, w2, w3);
}


