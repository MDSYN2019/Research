/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter10.h   (header file)
* Definitions:  SCchapter10.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER10_
#define _SCCHAPTER10_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"

double PowerMethod(SCMatrix &A, SCVector &x);
double PowerMethod(SCMatrix &A, SCVector &x, double tolerance);
double ISPowerMethod(SCMatrix &A, SCVector &x);
void TDQREigensolver(int N, double *a, double *b, double *lambda, double **Q);
void SolveSecularEq(double bm, int N, double *d, double *xi, double * lambda);
double SecularEq(double bm, int N, double *d, double *xi, double x);
#endif

