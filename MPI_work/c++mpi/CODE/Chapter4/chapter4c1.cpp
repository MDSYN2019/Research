

#include <iostream.h>
#include "SCmathlib.h"
#include "SCchapter4.h"


int main(int  argc, char * argv[]){
  int dim = 4;
  SCVector x(dim),b(dim),x0(dim);
  SCMatrix A(dim);
  
  // Set our initial guess
  x0(0) = x0(1) = x0(2) = x0(3) = 1.0;
  
  for(int i=0;i<dim;i++){
    for(int j=0;j<dim;j++){
      A(i,j) = 1+(i+1)*(j+1);
      
      /* We do this to make sure that the symmetric matrix that
	 we create has a determinant which is non-zero */
      
      if(i==3 && j == 2)
	A(i,j) = 12;
      if(i==2 && j == 3)
	A(i,j) = 12;
    }
  }

  cout << "The Matrix A that we are using: " << endl;
  A.Print();
  cout << endl;

  SCVector y(dim);
  y(0) = 2.;
  y(1) = -3.;
  y(2) = 5.43;
  y(3) = -22.56;

  cout << "The exact solution is: " << endl;
  y.Print();
  cout << endl;

  b = A*y;

  cout << "The right hand side, b, of the expression Ax=b: " << endl;
  b.Print();
  cout << endl;

  x = ConjugateGradient(A,b,x0);
  
  cout << "The approximate solution using Conjugate Gradient is: " << endl;
  x.Print();
  cout << endl;
 
}

