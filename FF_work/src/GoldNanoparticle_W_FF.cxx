/*
Algorithm for obtaining a Coarse-grained force-field under the SDK implementation, where each 
water molecule is coarse-grained and described in a 3-molecule blob

 |-----------|
 |H20     H2O|
 |    H2O    | <--------> Au      <---- Illustration of the water bead interaction with a single gold atom
 |-----------|

*/

// This work is incomplete - the FF produced from this work is difficult to fit to interpolate with a
// any form of standard LJ-like function. I'v tried to fit it to a Morse function but that is even more
// tricky. 

// #include <iostream>
// #include <cmath>

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <complex>
#include <vector>

#include "FFtools.hpp"
      
double FFtools::LJ (double dist_ab, double eps_ab, double sigma_ab) {   
  double LJ;
 
  LJ = 4*eps_ab*(pow((sigma_ab/dist_ab),12)-pow((sigma_ab/dist_ab),6));
  
  if (LJ < 0.0 && dist_ab < 2.0) {
    std::cout << "WARNING! \n" std::endl;
    exit(0);
  }  
  return LJ;
}

  
double FFtools::distance(double ax,double ay,double az, double bx, double by, double bz) {
    double dist;
    dist = sqrt((bx-ax)*(bx-ax) + (by-ay)*(by-ay) + (bz-az)*(bz-az));
    return dist;
}

double FFtools::minDist(double ax,double ay,double az, double bx, double by, double bz,double bxx, double byy, double bzz) {
    double dist; 
    double dx; // Delta X
    double dy; // Delta Y 
    double dz; // Delta Z

    dx = bx - ax;
    dy = by - ay;
    dz = bz - az;
    
    minimum_image(&dx, &dy, &dz, bxx, byy, bzz);
    dist = sqrt(dx*dx + dy*dy + dz*dz);
    return dist;
}

  
void FFtools::minimumImage(double *a, double *b, double *c, double bxx, double byy, double bzz) {
  /* This function applies the minimum image convention to a vector */
  double half_bxx = bxx/2.0;
  double half_byy = byy/2.0;
  double half_bzz = bzz/2.0;
  double newa; 
  double newb;
  double newc;
  
  /*changing the distance between points to account for the minimum image criterion*/
     
  if (*a >= half_bxx) {
    *a = *a - bxx;
  }
  else if (*a < -half_bxx) {
      *a = *a + bxx;
  }

   if (*b >= half_byy) {
    *b = *b - byy;
  }
  else if (*b < -half_byy) {
      *b = *b + byy;
  }

  if (*c >= half_bzz) {
    *c = *c - bzz;
  }
  else if (*c < -half_bzz) {
      *c = *c + bzz;
  }
}

void FFtools::continuity(double bx,double by,double bz,double bxx,double byy,double bzz) {

  /* This function applies the minimum image convention to a vector */

  double half_bxx = bxx/2.0;
  double half_byy = byy/2.0;
  double half_bzz = bzz/2.0;
  double newa; 
  double newb;
  double newc;

 /* changing the distance between points to account for the minimum image criterion */
    
  if (bx < -half_bxx) {
    bx = bx + bxx;
  }
  else if (bx >= half_bxx) {
    bx = bx - bxx;
  }   
  if (by < -half_byy) {
    by = by + byy;
  }
  else if (by >= half_byy) {
    by = by - byy;
  }
  if (bz < -half_bzz) {
    bz = bz + bzz;
 }
  else if (bz >= half_bzz) {
    bz = bz - bzz;
  }
}

double FFtools::centreofmass(int Oa, int H1a,int H2a,int Ob, int H1b, int H2b,int Oc, int H1c, int H2c){

  double Oadum = Oa; 
  double Obdum = Ob;
  double Ocdum = Oc;
  double H1adum = H1a; 
  double H2adum = H2a;
  double H1bdum = H1b;
  double H2bdum = H2b; 
  double H1cdum = H1c;
  double H2cdum = H2c;
  double COM; 

  COM = (Oadum*16 + H1adum*1.0 + H2adum*1.0 + Obdum*16 + H1bdum*1.0 + H2bdum*1.0 + Ocdum*16 + H1cdum*1.0 + H2cdum*1.0)/(16*3 + 6); 

  /* With this COM definition we now know the COM in each cartesian coordinate.*/ 
  /* We now have to calculate the difference in distance between each COM point */
  
  return COM; 
}  
