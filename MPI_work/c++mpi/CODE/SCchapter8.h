/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter8.h   (header file)
* Definitions:  SCchapter8.cpp (definition file)
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#ifndef _SCCHAPTER8_
#define _SCCHAPTER8_

#include <iostream.h>
#include <math.h>
#include "SCmathlib.h"

/* Euler Forward Central Difference Scheme - Periodic */
void EF_CentralDifference(int N, double CFL, double *uold, double *unew);

/* Euler Forward Central Difference Scheme - Periodic */
void EF_FirstOrderUpwind(int N, double CFL, double *uold, double *unew);

/* Lax-Friedrichs Scheme - Periodic */
void LaxFriedrichs(int N, double CFL, double *uold, double *unew);

/* Lax-Friedrichs-Tadmor Scheme - Periodic */
void LaxFriedrichsTadmor(int N, double CFL, double *uold, double *unew);
void LaxFriedrichsTadmor(int N, double CFL, double *uold, double *unew, double alpha);

/* Lax-Wendroff Scheme - Periodic */
void LaxWendroff(int N, double CFL, double *uold, double *unew);

/* Euler Foward Second Order Upwind Scheme - Periodic */
void EF_SecondOrderUpwind(int N, double CFL, double *uold, double *unew);

/* Leap-Frog Central Difference Scheme - Periodic */
void LeapFrog_CentralDifference(int N, double CFL, double **uold, double *unew);

/* Second Order Adams-Bashforth Central Difference Scheme - Periodic */
void AB2_CentralDifference(int N, double CFL, double **uold, double *unew);

/* Second Order Adams-Bashforth Central Difference Scheme - Periodic */
void AB3_CentralDifference(int N, double CFL, double **uold, double *unew);

/* Crank-Nicolson Central Difference Scheme - Periodic */
void CrankNicolson_CentralDifference(int N, double CFL, double *uold, double *unew);

/* Second Order Adams-Bashforth Central Difference Scheme - Non-Periodic */
void AB2_CentralDifferenceNP(int N, double CFL, double **uold, double *unew);

/* Leap-Frog Central Difference Scheme - Non-Periodic */
void LeapFrog_CentralDifferenceNP(int N, double CFL, double **uold, double *unew);

/* Crank-Nicolson Central Difference Scheme - Advection/Diffusion */
void CrankNicolson_CentralDifference_AdvectionDiffusion(int N, double CFL,
		  double DN, double *uold, double *unew);



/* MinMod functions used in Lax-Friedrich-Tadmor schemes */
double MinMod(double x, double y);
double MinMod(double x, double y, double z);

#endif
