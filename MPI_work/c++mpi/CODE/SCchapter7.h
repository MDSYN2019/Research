/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter7.h   (header file)
* Definitions:  SCchapter7.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER7_
#define _SCCHAPTER7_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"

void Diffusion_EF_CentralDifference(int N, double DN, double *uold, double *unew);
void Diffusion_EB_CentralDifference(int N, double DN, double *uold, double *unew);

int Diffusion_Jacobi(int N, double dx, double dt, 
		     double **A, double **q, double abstol);
int Diffusion_GaussSeidel(int N, double dx, double dt, 
			  double **A, double **q, double abstol);

/* Jacobi Algorithm for solving Ax=b */
int Jacobi(int N, double **A, double *x, double *b, double abstol);
void Jacobi(int N, double **A, double *x, double *b, int iterations);

/* Parallel Implementations */

int Jacobi_P(int mynode, int numnodes, int N, double **A, double *x, 
	     double *b, double abstol);
void Jacobi_P(int mynode, int numnodes, int N, double **A, double *x, 
	      double *b, int iterations);


/* Gauss-Seidel Algorithm for solving Ax=b */
int GaussSeidel(int N, double **A, double *x, double *b, double abstol);
void GaussSeidel(int N, double **A, double *x, double *b, int iterations);


/* SOR (Successive Over Relaxaion  Algorithm for solving Ax=b */
int SOR(double omega, int N, double **A, double *x, double *b, double abstol);
void SOR(double omega, int N, double **A, double *x, double *b, int iterations);

#endif
