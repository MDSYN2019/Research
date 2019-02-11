#include <iostream.h>
#include "SCchapter3.h" //contains declaration of
                        //CreateGrid_EvenlySpaced

double NewtonDiffFunction(int start_index, int ending_index, double * xpts, double * funcvals) {
  double val;

  int diff = ending_index - start_index;

  if (diff == 0) {
    val = funcvals[start_index]; 
  }

  else {
    val = (NewtonDiffFunction(start_index, ending_index-1, xpts, funcvals) - NewtonDiffFunction(start_index + 1, ending_index, xpts, funcvals)) / (xpts[start_index] - xpts[ending_index]) 
  }
  return val;
}

int main(int argc, char * argv[]){
  int i,npts;
  double *x;
  double dx;

  std::cout << "Enter the number of points in [0,1]: ";
  std::cin >> npts;

  x = new double[npts];
  
  CreateGrid_EvenlySpaced(npts, x, 0.0, 1.0);

  for(i=0;i<npts;i++)
    cout << "x[" << i << "] = " << x[i] << endl;

  delete[] x;

}
