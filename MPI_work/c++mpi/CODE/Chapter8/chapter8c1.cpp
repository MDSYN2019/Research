
#include <iostream.h>
#include <iomanip.h>
#include "SCchapter8.h"

double func(double x);
double truesoln(double x, double t);

int main(int argc, char * argv[]){
  int i,j;
  int npts = 1000; 
  int maxit = 10000;
  double CFL = .5;
  
  double dx = 1.0/npts;
  double * uold1 = new double[npts];
  double * unew1 = new double[npts];
  double * uold2 = new double[npts];
  double * unew2 = new double[npts];
  double * uold3 = new double[npts];
  double * unew3 = new double[npts];
  double * uold4 = new double[npts];
  double * unew4 = new double[npts];

  for(i=0;i<npts;i++){
    uold1[i] = func(i*dx);
    uold2[i] = func(i*dx);
    uold3[i] = func(i*dx);
    uold4[i] = func(i*dx);
  }

  for(i=0;i<maxit;i++){
    LaxFriedrichs(npts,CFL,uold1,unew1);
    LaxFriedrichsTadmor(npts,CFL,uold2,unew2,1.0);
    LaxFriedrichsTadmor(npts,CFL,uold3,unew3,2.0);
    LaxFriedrichsTadmor(npts,CFL,uold4,unew4,3.0);
    
    for(j=0;j<npts;j++){
      uold1[j] = unew1[j];
      uold2[j] = unew2[j];      
      uold3[j] = unew3[j];      
      uold4[j] = unew4[j];      
    }
  }
  
  cout << setiosflags(ios::scientific) << setprecision(5);
  for(i=0;i<npts;i++){
    cout << i*dx << "\t" << truesoln(i*dx,dx*CFL*maxit);
    cout << "\t" << unew1[i] << "\t" << unew2[i];
    cout << "\t" << unew3[i] << "\t" << unew4[i] << endl;
  }

  delete[] uold1;
  delete[] uold2;
  delete[] uold3;
  delete[] uold4;

  delete[] unew1;
  delete[] unew2;
  delete[] unew3;
  delete[] unew4;
}


double func(double x){
  double value = 0.0;

  while(x>1.0)
    x = x - 1.0;
  while(x<0.0)
    x = x + 1.0;

  if( (x >= 0.4) && (x <= 0.6))
    value = 100.0*(x-0.4)*(0.6-x);
  
  return(value);
}

double truesoln(double x, double time){
  return(func(x-time));
}


