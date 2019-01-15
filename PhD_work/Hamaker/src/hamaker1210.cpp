#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>
#include "HM1210.hpp"

//double LeJ(double sigma, double epsilon, double r, double pi, double goldden,double R);  
//double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R);

  /* Script to calculate the Hamaker constant */ 
using namespace HM1210;

void printPotential() {
  
  double r; 
  double sigma; 
  double epsilon;
  double LJ;
  double pi = 4.0*atan(1.0);
  double goldden = 0.059;
  double R = 20.0;
  double inter;
  double gradient; 
  int j = 0;
  double i;
  
  for (i = 20.10; i <= 35.0; i+= 0.01) { 
    r = i; 
    LJ = LeJ(3.5601,0.06205,r,pi,goldden,R); 
    inter = interaction(3.5601,0.06205,r,pi,goldden,R);
    gradient = -(LeJ(3.5601,0.06205,r+0.02,pi,goldden,R) - LeJ(3.5601,0.06205,r,pi,goldden,R))/0.02;
    j++;
    printf("%d %lf %lf %lf %lf \n",j,r,LJ,inter,gradient); 
  }
  
}
    
