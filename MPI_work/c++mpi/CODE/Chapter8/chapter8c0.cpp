
#include <iostream.h>
#include <iomanip.h>
#include "SCchapter8.h"

double func(double x);
double truesoln(double x, double t);

int main(int argc, char * argv[]){
  int i,j;
  int npts = 20; 
  int maxit = 150;
  double CFL = .1;

  double dx = 1.0/npts;
  double * uold = new double[npts];
  double * unew = new double[npts];

  for(i=0;i<npts;i++)
    uold[i] = func(i*dx);

  for(i=0;i<maxit;i++){
    EF_FirstOrderUpwind(npts,CFL,uold,unew);

    for(j=0;j<npts;j++)
      uold[j] = unew[j];
  }
  
  cout << setiosflags(ios::scientific) << setprecision(5);
  for(i=0;i<npts;i++)
    cout << i*dx << "\t" << truesoln(i*dx,dx*CFL*maxit) << "\t" << unew[i] << endl;

  delete[] unew;
  delete[] uold;

}


double func(double x){
  double value;
  value = sin(2.0*M_PI*x);
  return(value);
}

double truesoln(double x, double time){
  return(func(x-time));
}
