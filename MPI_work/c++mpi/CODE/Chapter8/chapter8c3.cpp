
#include <iostream.h>
#include <iomanip.h>
#include "SCchapter8.h"

double func(double x);
double truesoln(double x, double t);

int main(int argc, char * argv[]){
  int i,j;
  int npts = 46;
  int maxit = 25;
  double CFL = 1.0;
  double dx = 9.0/(npts-1);
  double DN = 0.01*CFL/dx;

  double *uold,*unew; 

  uold = new double[npts];
  unew = new double[npts];

  for(i=0;i<npts;i++)
    uold[i] = func(-3+i*dx);

  for(i=0;i<maxit;i++){
    CrankNicolson_CentralDifference_AdvectionDiffusion(npts,CFL,DN,uold,unew);
    unew[0] = 0.0;

    for(j=0;j<npts;j++){
      uold[j] = unew[j];
    }
  }
  
  cout << setiosflags(ios::scientific) << setprecision(5);
  for(i=0;i<npts;i++)
    cout << -3.0+i*dx << "\t" << truesoln(-3.0+i*dx,dx*CFL*maxit) << "\t" << unew[i] << endl;
  
  delete[] unew;
  delete[] uold;

}


double func(double x){
  double value;
  value = exp(-2.0*x*x);
  return(value);
}

double truesoln(double x, double time){
  return(func(x-time));
}
