#ifndef _AAProbabilityArray_
#define _AAProbabilityArray_

#include <math.h>
#include "../interfaceMSA/ResName.h"
#include "SubstitutionMatrix.h"

using namespace std;

class AAProbabilityArray
{
protected: 
	double pa[20];
	double sampleNum;


public:
	AAProbabilityArray(double pa[]);
	AAProbabilityArray(int count[]);
	AAProbabilityArray(double pa[], double sampleNum);

	~AAProbabilityArray();


	void normalize();
	void addPseudoCount(AAProbabilityArray& pb, double num);

	AAProbabilityArray applySubstitution(SubstitutionMatrix& sub);
	

	double getValue(int i);

};

#endif
