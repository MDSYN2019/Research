/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter4.h   (header file)
* Definitions:  SCchapter4.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER4_
#define _SCCHAPTER4_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"

double SquareRoot(double value, double guess, int iterations);

double NewtonRaphson(double x0, double (*func)(double), 
		     double (*func_der)(double), int max_iter, 
		     int multiplicity=1);

double NewtonRaphson(double x0, double (*func)(double), 
		     double (*func_der)(double), 
		     double (*func_secondder)(double), 
		     int max_iter, int multiplicity=1);


SCVector ConjugateGradient(SCMatrix A, SCVector b, SCVector x0);

double MidpointRule(int level, double xleft, double xright, 
		    double (*func)(double));
double TrapezoidRule(int level, double xleft, double xright, 
		     double (*func)(double));
double Romberg(int m, int k, double xleft, double xright, 
	       double (*func)(double));
double SimpsonsRule(int m, double xleft, double xright, 
		    double (*func)(double));


/* Jacobi Polynomials */
double JacobiPoly(int degree, double x, double alpha, double beta);
double JacobiPolyDerivative(int degree, double x, double alpha, double beta);
void   JacobiZeros(int degree, double *z, double alpha, double beta);
void   JacobiZW(int degree, double * z, double *w, double alpha, double beta);

/* Hermite Polynomials */
double HermitePoly(int degree, double x);
double HermitePolyDerivative(int degree, double x);
void   HermiteZeros(int degree, double *z);
void   HermiteZW(int degree, double *z, double *w);

/* Laguerre Polynomials */
double LaguerrePoly(int degree, double x);
double LaguerrePolyDerivative(int degree, double x);
void   LaguerreZeros(int degree, double *z);
void   LaguerreZW(int degree, double * z, double *w);

#endif
