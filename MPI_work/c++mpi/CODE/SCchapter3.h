/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter3.h   (header file)
* Definitions:  SCchapter3.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER3_
#define _SCCHAPTER3_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"


double NewtonInterpolant(double x, int npts, double * xpts, double * coeffs);
void   NewtonDiffTable(int npts, double *xpts, double *funcvals, 
			 double * coeff);
double NewtonDiffFunction(int start_index, int ending_index, 
			  double * xpts, double * funcvals);


double LagrangePoly(double x, int pt, int npts, double * xpts);
double LagrangeInterpolant(double x, int npts, double *xpts, 
			   double * funcvals);

double ChebyshevPoly(int degree, double x);
void   ChebyshevPoints(int npts, double * x);


void  CreateGrid_EvenlySpaced(int npts, double *x, double a, double b);
void  CreateGrid_ChebyshevPts(int npts, double *x, double a, double b);

void LS_ComputeCoeffs(int npts, double *xpts, double *funcval, 
		      int deg, double *alpha, double *beta, double *lscoeffs);

double LS_OrthoPoly(int j, double x, double *alpha, double *beta);

double LSApproximatingPoly(int ndeg, double x, double *alpha, 
			   double *beta, double *lscoeffs);


double Square_2dInterpolant(SCPoint x, int npts,  double *funcvals);


class LSPoly{
 private:
  int ndeg;
  double *alpha, *beta, *lscoeffs;

  double LSPolyOrtho(int j, double x);
  
 public:
  LSPoly();
  ~LSPoly();

  void PrintCoeffs();
  int Initialize(int npts, int in_ndeg, double * xpts, double * funcvals);
  double Evaluate(double x);
};

#endif
