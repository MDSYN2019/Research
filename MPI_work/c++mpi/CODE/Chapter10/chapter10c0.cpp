
#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"
#include "SCchapter10.h"

const int size = 3;

int main(int argc, char *argv[]){
  int i;
  SCMatrix A(size);
  SCVector b(size),x(size);
  
  A(0,0) = 1.0; A(0,1) = 0.5; A(0,2) = 1.0/3.0;
  A(1,0) = 0.5; A(1,1) = 1.0/3.0; A(1,2) = 0.25;
  A(2,0) = 1./3.; A(2,1) = 0.25; A(2,2) = 1.0/5.0;
  
  A.Print();

  cout << "----------------" << endl;

  x.Initialize(1.0);
 
  double ff = ISPowerMethod(A,x);

  cout << "eigenvalue = " << ff << endl;

  x.Print();

}






