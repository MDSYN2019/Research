/*Edit:The interaction potentials has been checked and it is correct!*/

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <complex.h>
#include <string.h>

double LeJ(double sigma, double epsilon, double r, double pi, double goldden,double R);
double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R);

int main(int argc, char *argv[]) 
{

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
  

 /*defining the logarithmic part*/ 
 
 
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

    /*defining the logarithmic part*/ 


    printf("%d %lf %lf %lf %lf \n",j,r,LJ,inter,gradient); 
  
  }
  
  return(0);
}



double LeJ(double sigma, double epsilon, double r,double pi, double goldden, double R)

  {
    double HM; 
 
    
    HM = (9*pi*goldden*epsilon*pow(sigma,9)/35)*((pow(R,3)*(3*pow(R,4) + 42*pow(R,2)*pow(r,2) + 35*pow(r,4)))/(r*pow(pow(r,2)-pow(R,2),6))) - (9*pi*goldden*epsilon*pow(sigma,6)*pow(R,3))/(pow(pow(r,2)-pow(R,2),3)); 

    return HM;

  }


double interaction(double sigma, double epsilon, double r, double pi, double goldden,double R)

{ 

  /*This still needs editing*/
  double interaction;
  

  interaction = ((27*pi*epsilon*goldden*pow(sigma,9))/35)*((pow(R,3)*(-pow(R,6) + 27*pow(R,4)*pow(r,2) + 189*pow(R,2)*pow(r,4) + 105*pow(r,6)))/(pow(r,2)*pow(pow(r,2) - pow(R,2),7))) - (54*pi*goldden*epsilon*pow(sigma,6))*((pow(R,3)*r)/pow(pow(r,2)-pow(R,2),4)); 


  return interaction;

}  
