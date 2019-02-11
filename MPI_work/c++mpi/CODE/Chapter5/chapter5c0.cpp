
#include <iostream.h>
#include <iomanip.h>
#include "SCchapter5.h"

double func(double x);
double func_first_der(double x);
double func_second_der(double x);
 

int main(int argc, char * argv[]){
  const int levels = 10;
  const double a = 0.0;
  const double b = 1.0;
  int i,j,npts;
  double dx,dxp,ux_error,uxx_error;
  double *u,*u_x,*u_xx;

  cout << "npts\tError (First Deriv)\tError (Second Deriv)\n";
  for(i=2;i<levels+2;i++){
    npts = (int) pow(2.0,i);
    dx = 1.0/(npts-1);
    
    u = new double[npts];
    u_x = new double[npts];
    u_xx = new double[npts];
    
    for(j=0;j<npts;j++)
      u[j] = func(j*dx);
    
    SO_FirstDeriv_1D  (npts,dx,u,u_x);
    SO_SecondDeriv_1D (npts,dx,u,u_xx);

    ux_error=0.0;
    uxx_error=0.0;
    for(j=0;j<npts;j++){
      ux_error  += dx*pow((u_x[j]-func_first_der(j*dx)),2);
      uxx_error += dx*pow((u_xx[j]-func_second_der(j*dx)),2);
    }
    
    cout << setprecision(10) << setiosflags(ios::scientific);
    cout << npts << "\t" << sqrt(ux_error) << "\t" << sqrt(uxx_error) << endl;

    delete[] u;
    delete[] u_x;
    delete[] u_xx;
  }
  
}

double func(double x){
  return(x*x*x*x);
}

double func_first_der(double x){
  return(4*x*x*x);
}

double func_second_der(double x){
  return(12*x*x);
}

