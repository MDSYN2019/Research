
#include <iostream.h>
#include <iomanip.h>
#include "SCchapter8.h"

double func(double x);
double truesoln(double x, double t);

int main(int argc, char * argv[]){
  int i,j;
  int npts = 300;
  int maxit = 30000;
  double CFL = 0.1;

  double dx = 1.0/npts;
  double *uold[3],*unew; 

  uold[0] = new double[npts];
  uold[1] = new double[npts];
  uold[2] = new double[npts];
  unew = new double[npts];

  for(i=0;i<npts;i++){
    uold[0][i] = func(i*dx);
    uold[1][i] = func(i*dx);
    uold[2][i] = func(i*dx);
  }

  for(i=0;i<maxit;i++){
    AB3_CentralDifference(npts, CFL, uold, unew);

    for(j=0;j<npts;j++){
      uold[2][j] = uold[1][j];
      uold[1][j] = uold[0][j];
      uold[0][j] = unew[j];
    }
  }
  
  cout << setiosflags(ios::scientific) << setprecision(5);
  for(i=0;i<npts;i++)
    cout << i*dx << "\t" << truesoln(i*dx,dx*CFL*maxit) << "\t" << unew[i] << endl;

  delete[] unew;
  delete[] uold[0];
  delete[] uold[1];

}


double func(double x){
  double value;
  value = sin(2.0*M_PI*x);
  return(value);
}

double truesoln(double x, double time){
  return(func(x-time));
}
