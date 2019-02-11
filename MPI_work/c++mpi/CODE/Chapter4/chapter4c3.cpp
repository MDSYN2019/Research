
#include <iostream.h>
#include <iomanip.h>
#include <math.h>
#include "SCchapter4.h"


double func(double x);

int main(int * argc, char ** argv[]){
  const double a=0.0,b=4.0;
  const int level = 3;

  cout << setprecision(10);
  cout << "The value of x*exp(x) integrated on [0,4] using \n Trapezoidal Rule (" << pow(2,level)+1 << " equidistant points) is: " ;
  cout << TrapezoidRule(level,a,b,func) << "  " << endl;
}

double func(double x){
  return(x*exp(x));
}


