#include "AAProbabilityArray.h"
#include <iostream>

AAProbabilityArray::AAProbabilityArray(double pa[])
{
	for (int i = 0; i < 20; i++)
		this->pa[i] = pa[i];
        this->sampleNum = -1;
	normalize();
}

AAProbabilityArray::AAProbabilityArray(int count[])
{
	this->sampleNum = 0;
	for (int i = 0; i < 20; i++)
	{
		this->pa[i] = count[i];
		this->sampleNum += count[i];
	}
	normalize();
}

AAProbabilityArray::AAProbabilityArray(double pa[], double sampleNum)
{
	for (int i = 0; i < 20; i++)
		this->pa[i] = pa[i];
	this->sampleNum = sampleNum;
	normalize();
}



AAProbabilityArray::~AAProbabilityArray()
{
}

void AAProbabilityArray::normalize()
{
	double tot = 0;
	for (int i = 0; i < 20; i++)
	{
		tot += pa[i];
	}
	if (tot == 0)
		return;
	for (int i = 0; i < 20; i++)
	{
		pa[i] /= tot;
	}
}

void AAProbabilityArray::addPseudoCount(AAProbabilityArray & pb, double psNum)
{
	if (this->sampleNum < 0)
	{
		cerr << "sample number unknow, can't add pseudo count" << endl;
		return ;
	}
	int totNum = this->sampleNum + psNum;
	for (int i = 0; i < 20; i++)
	{
		this->pa[i] = (this->pa[i] * this->sampleNum + pb.pa[i] * psNum)/totNum;
	}
	this->sampleNum = -1;
	normalize();
}

AAProbabilityArray AAProbabilityArray::applySubstitution(SubstitutionMatrix& sub)
{
	double p[20];
	for (int i = 0; i < 20; i++)
	{
		p[i] = 0;
		for (int j = 0; j < 20; j++)
		{
			p[i] += this->pa[j] * sub.getValue(i, j);
		}
	}
	return AAProbabilityArray(p);
}

double AAProbabilityArray::getValue(int i)
{
	if (i < 0 || i > 19)
	{
		cerr << "invalid index : " << i << endl;
		exit(0);
	}
	return this->pa[i];
}
