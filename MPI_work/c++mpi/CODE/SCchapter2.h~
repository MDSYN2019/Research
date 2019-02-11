/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter2.h   (header file)
* Definitions:  SCchapter2.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER2_
#define _SCCHAPTER2_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"


float  FloatMachineEps();
double DoubleMachineEps();

SCstatus GramSchmidt(SCVector * x, SCVector * q);
SCstatus GramSchmidt(SCVector * x, SCVector * q, SCMatrix &r);
SCstatus ModifiedGramSchmidt(SCVector * x, SCVector * q, SCMatrix &r);
SCstatus QRDecomposition(SCMatrix A, SCMatrix &Q, SCMatrix &R);
#endif

