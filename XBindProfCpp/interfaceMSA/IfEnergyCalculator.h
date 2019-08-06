#ifndef _IfEnergyCalculator_
#define _IfEnergyCalculator_

#include "IfAlignment.h"
#include "../math/AAProbabilityArray.h"
#include "IfSingleMutant.h"
#include "ResName.h"
#include <math.h>

namespace IfEnergyCalculator
{
	double getEnergy(string posAASeq, char aaType, double w1, double w2 , double w3);
	double getEnergy(string posAASeq, char aaType);
	double getDDG(IfAlignment align, IfSingleMutant mut);
	double getDDG(IfAlignment align, IfSingleMutant mut, double w1, double w2, double w3);

}



#endif

