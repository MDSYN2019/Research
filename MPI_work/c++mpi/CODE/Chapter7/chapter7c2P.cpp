
#include <iostream.h>
#include <iomanip.h>
#include <mpi.h>
#include "SCmathlib.h"
#include "SCchapter7.h"

int main(int argc, char * argv[]){
  int i,j, N = 20;
  double **A,*x,*q;
  int totalnodes,mynode;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  if(mynode==0){
    A = CreateMatrix(N,N);
    x = new double[N];
    q = new double[N];

    for(i=0;i<N;i++){
      q[i] = i+1;
      A[i][i] = -2.0;
      if(i<N-1){
	A[i][i+1] = 1.0;
	A[i+1][i] = 1.0;
      }
    }
  }

  Jacobi_P(mynode,totalnodes,N,A,x,q,1.0e-14);
  
  if(mynode==0){
    for(i=0;i<N;i++)
      cout << x[i] << endl;
    DestroyMatrix(A,N,N);
    delete[] x;
    delete[] q;
  }

  MPI_Finalize();

}

