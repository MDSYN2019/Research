#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>
#include <string.h>
#include <iostream>
/*

  Hydrophobic CM (epsilon: 0.4200 (kcal/mol), sigma: 4.5060 (\AA))

  Hydrophilic W (epsilon: 0.8950 (kcal/mol), sigma: 4.3710 (\AA))

 */


  int r; 
  double sigma; 
  double epsilon;
  double LJ;
  double i;
  double pi = 4.0*atan(1.0);
  double goldden = 0.059;
  double R = 20.0;
  double inter; 
  double gradient; 
  int j = 0;
  int a = 0;

  double P(double d, double r) {

    double p_particle;
    
    p_particle = 525*pow(d,12) + 10500*r*pow(d,9) + 93660*pow(r,2)*pow(d,8) + 490560*pow(r,3)*pow(d,7) + 1674456*pow(r,4)*pow(d,6) + 3911712*pow(r,5)*pow(d,5) + 6376832*pow(r,6)*pow(d,4) + 7200256*pow(r,7)*pow(d,3) + 5392384*pow(r,8)*pow(d,2) + 2416640*pow(r,9)*d + 491520*pow(r,10);

      return p_particle;
   }

double Upp(double d, double r, double sigma, double epsilon, double P) {

    double pp;

    double p1 = 0.059;
    double p2 = 0.059;
    
    pp = (64 * pow(pi,2) * p1 * p2 * epsilon * pow(sigma,12) * pow(r,6) * P)/(4725 * pow(d,7) * pow((2*r + d),8) * pow((4*r +d),7))
      + ((2 * pow(pi,2) * p1 * p2 * epsilon * pow(sigma,6))/3)*((r/2) * ((1/(4*r +d)) - (1/d) - (4*r/pow((2*r + d),2))) + log10(pow((2*r +d),2)/(d*(4*r+d))));
      return pp;
      }

double Pe;

int main () {
  
  for (i=10.0;i<=25.0;i+=0.01) { 
    
    Pe = 0;
    r = i;

    printf("%i %f %lf %.15lf \n",a, i, Upp(i-10.0, 10.0 , 4.5060, 0.42, P(i,10.0))+18.855768, (Upp(i-10.0, 10.0, 4.5060 ,0.42,P(i+0.02,10.0)) - Upp(i-10.0,10.0,4.5060,0.42,P(i,10.0)))/0.02);
    //LJ = LeJ(3.2529,0.666381,r,pi,goldden,R);
    //inter = interaction(3.2529,0.666381,r,pi,goldden,R);  
    a = a + 1;
  }
  return 0;
}
