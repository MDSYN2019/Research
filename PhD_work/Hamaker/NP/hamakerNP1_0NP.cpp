#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>
#include <string.h>


namespace HM128 {

double LeJ(double sigma, double epsilon, double r, double pi, double goldden,double R);
double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R);

  double r; 
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

  epsilon=atof(argv[1]);
  sigma=atof(argv[2]);
  

  for (i=20.10;i<=35.0;i+=0.01) { 
    
    r = i;
    //LJ = LeJ(3.2529,0.666381,r,pi,goldden,R);
    //inter = interaction(3.2529,0.666381,r,pi,goldden,R);
    LJ=LeJ(sigma,epsilon,r,pi,goldden,R);
    inter=interaction(sigma,epsilon,r,pi,goldden,R);
    gradient = -(LeJ(sigma,epsilon,r+0.02,pi,goldden,R) - LeJ(sigma,epsilon,r,pi,goldden,R))/0.02;
    j++; 

    printf("%d %lf %lf %lf %lf \n",j,r,LJ,inter,gradient);
  
  }


  /*  
     Ned to make sure that these functions for NP-NP partices work correctly 
   */
  double P(double d, double r) {

    double p_particle;
    
    p = 525.0*pow(d,12) + 10500*r*pow(d,9) + 93660*pow(r,2)*pow(d,8) + 490560*pow(r,3)*pow(d,7) + 1674456*pow(r,4)*pow(d,6) + 391171*pow(r,5)*pow(d,5) + 6376832*pow(r,6)*pow(d,4) + 7200256*pow(r,7)*pow(d,3) + 5392384*pow(r,8)*pow(d,2) + 2416640*pow(r,9)*d + 491520*pow(r,10)

      return p_particle;
   }

  double Upp(double d, double r, double sigma, double epsilon, double P) {

    double upp;

    double p1 = 0.113;
    double p2 = 0.113;
    
    upp = (64 * pow(pi,2) * p1 * p2 * epsilon * power(sigma,12) * pow(r,6) * P)/(4725 * pow(d,7) * pow((2*r + d),8) * pow((4*r +d),7))
      + (2 * pow(pi,2) * p1 * p2 * epsilon * pow(sigma,6)/3)*[(r/2) * (1/(4*r +d) - (1/d) - (4*r/pow((2*r + d),2))) + log10(pow((2*r +d),2)/d*(4*r+d))]
      return upp
  }
  
double LeJ(double sigma, double epsilon, double r,double pi, double goldden, double R)

  {
    double HM; 
    
    
    HM = (3*pi*goldden*epsilon*pow(sigma,12)*pow(R,3)/5)*(5*pow(R,6) + 45*pow(R,4)*pow(r,2) + 63*pow(R,2)*pow(r,4) + 15*pow(r,6))/pow(pow(r,2)-pow(R,2),9)
      - (9*pi*goldden*epsilon*pow(sigma,8)/5)*((pow(R,3)*(3*pow(R,2) + 5*pow(r,2))))/pow(pow(r,2) - pow(R,2),5);

    return HM;

  }


double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R)

{ 

  double interaction; 

  interaction = (108*pi*epsilon*goldden*pow(sigma,12)*pow(R,3)/5)*(5*pow(R,6) + 27*pow(R,4)*pow(r,2) + 27*pow(R,2)*pow(r,4) + 5*pow(r,6))*r/pow(pow(r,2) - pow(R,2),10)
              - 72*pi*goldden*epsilon*pow(sigma,8)*pow(R,3)*(pow(r,2) + pow(R,2))*r/pow(pow(r,2) - pow(R,2),6); 


  return interaction;

}
  
}
