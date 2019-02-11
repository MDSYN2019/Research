/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter5.h   (header file)
* Definitions:  SCchapter5.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER5_
#define _SCCHAPTER5_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"

/***********************************************************************/
/* First Derivative Operators */

/* Second Order First Derivative in 1D */
void SO_FirstDeriv_1D    (int npts, double dx, double *u, double *u_x);

/* Second Order First Derivative in 1D, periodic boundaries */
void SO_FirstDeriv_1Dper (int npts, double dx, double *u, double *u_x);

/* Second Order First Derivative in 1D -- Parallel*/
void SO_FirstDeriv_1DP    (int npts, double dx, double *u, double *u_x, int mynode, int totalnodes);

/* Second Order First Derivative in 1D, periodic boundaries -- Parallel */
void SO_FirstDeriv_1DperP (int npts, double dx, double *u, double *u_x, int mynode, int totalnodes);
/***********************************************************************/


/***********************************************************************/
/* Second Derivative Operators */

/* Second Order Second Derivative in 1D */
void SO_SecondDeriv_1D    (int npts, double dx, double *u, double *u_x);

/* Second Order Second Derivative in 1D, periodic boundaries */
void SO_SecondDeriv_1Dper (int npts, double dx, double *u, double *u_x);

/* Second Order Second Derivative in 1D -- Parallel*/
void SO_SecondDeriv_1DP    (int npts, double dx, double *u, double *u_x, int mynode, int totalnodes);

/* Second Order Second Derivative in 1D, periodic boundaries -- Parallel */
void SO_SecondDeriv_1DperP (int npts, double dx, double *u, double *u_x, int mynode, int totalnodes);
/***********************************************************************/


void FornbergWeights(double xi, double *x, int m, int n, double ***C);

double AdamsBashforth(int order, double u_old, double dt, double * RHS);
void   AdamsBashforth(int order, int N, double *u_old, double *u_new, double dt, double **RHS);

double RungeKutta4(double uold, double time, double dt, 
		   double (*rkfunc)(double,double));
double RungeKutta(int order, double uold, double (*rkfunc)(double));

#endif

