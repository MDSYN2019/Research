#include <iostream.h>

int main(int argc, char * argv[]){
  int i,npts;
  double *x;
  double dx;

  cout << "Enter the number of points into which you want to discretize [0,1]: ";
  cin >> npts;

  x = new double[npts];
  dx = 1.0/(npts-1);

  for(i=0;i<npts;i++)
    x[i] = i*dx;

  for(i=0;i<npts;i++)
    cout << "x[" << i << "] = " << x[i] << endl;

  delete[] x;

}
