#ifndef __MD_function__
#define __MD_function__

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>

double LJ (double dist_ab, double eps_ab, double sigma_ab);  // What does this calculate?
double distance(double ax,double ay,double az, double bx, double by, double bz); 
double minDist(double ax,double ay,double az, double bx, double by, double bz, double bxx, double byy, double bzz);
void minimumImage(double *a,double *b,double *c,double bxx,double byy,double bzz); // implementing the minimum image conventon
double centreOfMass(int Oa, int H1a, int H2a,int Ob, int H1b, int H2b,int Oc, int H1c, int H2c);

#endif
