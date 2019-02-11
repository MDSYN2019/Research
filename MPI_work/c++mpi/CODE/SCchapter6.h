/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter6.h   (header file)
* Definitions:  SCchapter6.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER6_
#define _SCCHAPTER6_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"

/* Thomas Algorithm for Solving the Tri-diag system [b a c]*x=q */
void ThomasAlgorithm(int N, double b, double a, double c, double *x, double *q);
void ThomasAlgorithm(int N, double *b, double *a, double *c, double *x, double *q);

void ThomasAlgorithmLU(int N, double *b, double *a, double *c, double *l, double *u, double *d);
void ThomasAlgorithmSolve(int N, double *l, double *u, double *d, double *x, double *q);

/* Periodic Thomas Algorithm which uses Matrix Condensation */
void ThomasAlgorithm_per(int N, double b, double a, double c, double *x, double *q);
void ThomasAlgorithm_per(int N, double *b, double *a, double *c, double *x, double *q);


void ThomasAlgorithm_P(int mynode, int numnodes, 
		       int N, double *a, double *b, double *c, double *x, double *q);


#endif
