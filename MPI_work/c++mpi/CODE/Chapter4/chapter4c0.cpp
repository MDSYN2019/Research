
#include <iostream.h>
#include "SCchapter4.h"

const double alpha = 2.0;

double func(double x);
double func_der(double x);
double func_secondder(double x);

int main(int argc, char * argv[]){
  int max_iter = 100;
  int multiplicity = 1;
  double x0 = 1;
 

  double x =  NewtonRaphson(x0,func,func_der,func_secondder,
			    max_iter,multiplicity);
  
  cout << x << " " << fabs(func(x)) <<endl;;
}



double func(double x){
  double value = (x-alpha)*(x-alpha)*(x-alpha)*(x+alpha);
  return(value);
}

double func_der(double x){
  double value = 3*(x-alpha)*(x-alpha)*(x+alpha)+ 
    (x-alpha)*(x-alpha)*(x-alpha);
  
  return(value);
}

double func_secondder(double x){
  double value = 6*(x-alpha)*(x+alpha)+6*(x-alpha)*(x-alpha);
  return(value);
}

