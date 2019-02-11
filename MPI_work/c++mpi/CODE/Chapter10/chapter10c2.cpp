

#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"
#include "SCchapter10.h"

const int N = 20;

int main(int argc, char *argv[]){
  int i,j;

  double * a = new double[N];
  double * b = new double[N];
  double * lambda = new double[N];
  double ** Q = CreateMatrix(N,N);


  for(i=0;i<N;i++){
    a[i] = 1.5 + i;
    b[i] = 0.3;
    lambda[i] = 0.0;
  }

  TDQREigensolver(N,a,b,lambda,Q);

  cout << "The eigenvalues are: " << endl;
  for(i=0;i<N;i++)
    cout << lambda[i] << endl;

  delete[] a;
  delete[] b;
  delete[] lambda;
  DestroyMatrix(Q,N,N);
}


