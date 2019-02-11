
#include <iostream.h>
#include "SCchapter3.h" //contains declaration of
                        //CreateGrid_EvenlySpaced

int main(int argc, char * argv[]){
  int i,npts;
  double *x;
  double dx;

  cout << "Enter the number of points in [0,1]: ";
  cin >> npts;

  x = new double[npts];
  
  CreateGrid_EvenlySpaced(npts, x, 0.0, 1.0);

  for(i=0;i<npts;i++)
    cout << "x[" << i << "] = " << x[i] << endl;

  delete[] x;

}
