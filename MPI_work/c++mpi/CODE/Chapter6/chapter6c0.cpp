
#include <iostream.h>
#include <iomanip.h>
#include "SCchapter6.h"


int main(int argc, char * argv[]){
  int N = 8,i;
  double *a,*b,*c,*q,*x;

  a = new double[N];
  b = new double[N];
  c = new double[N];
  q = new double[N];
  x = new double[N];

  for(i=0;i<N;i++){
    a[i] = -2.0;
    b[i] = c[i] = 1.0;
    q[i] = (double)(i+1);
  }
  
  ThomasAlgorithm(N,b,a,c,x,q);

  for(i=0;i<N;i++)
    cout << x[i] << endl;

  delete[] a;
  delete[] b;
  delete[] c;
  delete[] q;
  delete[] x;
}

