#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"


const int size = 15;

int main(int argc, char *argv[]){
  int i,j,k;
  int index1,index2,offset;
  double alpha,gamma;

  /* Part 1 */

  double * x = new double[size];
  for(i=0;i<size;i++)
    x[i] = 0.0;
  double * F = new double[size];
  double ** A = new double*[size];


  for(i=0;i<size;i++){
    A[i] = new double[size];
    for(j=0;j<size;j++)
      A[i][j] = 0.;
    F[i] = (double)i;
  }

  A[0][0] = -2.0; A[0][1] = 1.0;
  A[size-1][size-2] = 1.0; A[size-1][size-1] = -2.0;

  for(i=1;i<size-1;i++){
    A[i][i] = -2.0;
    A[i][i-1] = 1.0;
    A[i][i+1] = 1.0;
  }


  /* Part 2 */
  
  for(i=0;i<log2(size+1)-1;i++){
    for(j=pow(2.0,i+1)-1;j<size;j=j+pow(2.0,i+1)){
      offset = pow(2.0,i);
      index1 = j - offset;
      index2 = j + offset;

      alpha = A[j][index1]/A[index1][index1];
      gamma = A[j][index2]/A[index2][index2];
      
      for(k=0;k<size;k++){
	A[j][k] -= (alpha*A[index1][k] + gamma*A[index2][k]);
      }
      F[j] -= (alpha*F[index1] + gamma*F[index2]);
    }
  }


  /* Part 3 */

  int index = (size-1)/2;
  x[index] = F[index]/A[index][index];

  for(i=log2(size+1)-2;i>=0;i--){
    for(j=pow(2.0,i+1)-1;j<size;j=j+pow(2.0,i+1)){
      offset = pow(2.0,i);
      index1 = j - offset;
      index2 = j + offset;

      x[index1] = F[index1];
      x[index2] = F[index2];
      for(k=0;k<size;k++){
	if(k!= index1)
	  x[index1] -= A[index1][k]*x[k];
	if(k!= index2)
	  x[index2] -= A[index2][k]*x[k];
      }

      x[index1] = x[index1]/A[index1][index1];
      x[index2] = x[index2]/A[index2][index2];
    }
  }
  
  for(i=0;i<size;i++){
    cout << x[i] << endl;
  }
  
  delete[] x;
  delete[] F;
  for(i=0;i<size;i++)
    delete[] A[i];
  delete[] A;
}


