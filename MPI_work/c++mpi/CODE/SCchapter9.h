/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter9.h   (header file)
* Definitions:  SCchapter9.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER9_
#define _SCCHAPTER9_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"


void GaussElimination(SCMatrix &A, SCVector &b, int pivotflag);
double HouseholderTrans(SCVector &x, SCVector &w);
void HouseholderQR(SCMatrix &A, SCVector &v);
void Hessenberg(SCMatrix &A);
void ModifiedArnoldi(int m, const SCMatrix &A, SCMatrix &H, SCMatrix &V);
void ModifiedArnoldi(int m, const SCVector &x, const SCMatrix &A, SCMatrix &H, SCMatrix &V);
void GMRES(int m, const SCMatrix &A, const SCVector &b, SCVector &x);
#endif
 
