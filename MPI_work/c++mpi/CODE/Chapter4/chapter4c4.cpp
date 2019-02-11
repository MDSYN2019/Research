

#include<iostream.h>
#include<iomanip.h>
#include <math.h>
#include "SCmathlib.h"
#include "SCchapter4.h"

int main(int  argc, char * argv[]){
  int i;
  int degree = 2;
  double alpha = 0.0,beta = 0.0;

  cout << " Enter degree: ";
  cin >> degree;
  
  double * z = new double[degree];
  double * w = new double[degree];

  cout << setprecision(10);

  cout << "---------------------------------------" << endl;
  cout << "Legendre Quadrature Zeros and Weights "  << endl;

  JacobiZW(degree,z,w,alpha,beta);
  for(i=0;i<degree;i++)
    cout << z[i]  << "\t" << w[i] << endl;
  
  cout << "---------------------------------------" << endl;  
  cout << "Hermite Quadrature Zeros and Weights "  << endl;

  HermiteZW(degree,z,w);
  for(i=0;i<degree;i++)
    cout << z[i]  << "\t" << w[i] << endl;
  
  cout << "---------------------------------------" << endl;
  cout << "Laguerre Quadrature Zeros and Weights "  << endl;
  LaguerreZW(degree,z,w);
  for(i=0;i<degree;i++)
    cout << z[i]  << "\t" << w[i] << endl;
  cout << "---------------------------------------" << endl;

  delete[] z;
  delete[] w;
}


