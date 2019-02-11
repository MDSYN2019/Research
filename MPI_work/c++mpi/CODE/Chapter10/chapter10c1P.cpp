
#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"
#include "SCchapter3.h"
#include<mpi.h>

void ChebyVandermonde(int npts, double *A, int row);

// Global variable to set size of the system
const int size = 10;


int main(int argc, char *argv[]){
  int i,j,k,index,cnt;
  int mynode, totalnodes;
  double tmps1,tmps2,sum,total;
  const int maxit = 1000; //maximum number of iterations for power method

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  int numrows = size/totalnodes;
  double **A_local   = new double*[numrows];
  double **A_local_t = new double*[numrows];
  double * xold = new double[size];
  double * xnew = new double[numrows];
  double * tmparray = new double[size];

  for(i=0;i<numrows;i++){
    A_local[i] = new double[size];
    A_local_t[i] = new double[size];
    index = mynode*numrows + i;
    ChebyVandermonde(size,A_local[i],index);
  }
  
  /********************************************************/
  /* Use the power method to obtain the first eigenvector */
  /********************************************************/

  // Initialize starting vector to 1.0
  for(i=0;i<size;i++)
    xold[i] = 1.0;  

  for(int it=0;it<maxit;it++){
    // Matrix-vector multiplication
    for(i=0;i<numrows;i++)
      xnew[i] = dot(size,A_local[i],xold);

    // Compute Euclidian norm of new vector
    sum = 0.0;
    for(i=0;i<numrows;i++)
      sum += xnew[i]*xnew[i];
    MPI_Allreduce(&sum,&total,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    total = sqrt(total);

    // Scale vector by its norm
    for(i=0;i<numrows;i++)
      xnew[i] = xnew[i]/total;
    
    // Gather (Allgather) new vector to all processors
    MPI_Allgather(xnew,numrows,MPI_DOUBLE,tmparray,numrows,
		   MPI_DOUBLE,MPI_COMM_WORLD);
    
    // Compute difference between old and new vector
    sum = 0.0;
    for(i=0;i<size;i++){
      sum += (xold[i]-tmparray[i])*(xold[i]-tmparray[i]);
      xold[i] = tmparray[i]; //replace old with new
    }
    
    if(sqrt(sum) < 1.0e-7) //termination condition
      break;
  }

  // Modify eigenvector per Householder transformation
  xold[0] = xold[0] + Sign(xold[0]);

  /**************************************************************/
  /* Compute S = A*H (S is a temporary state stored in A_local) */
  /**************************************************************/
  
  tmps1 = dot(size,xold,xold);
  for(i=0;i<numrows;i++){
    tmps2 = dot(size,A_local[i],xold);
    for(j=0;j<size;j++)
      A_local[i][j] = A_local[i][j] - 2.0*tmps2*xold[j]/tmps1;
  }


  /******************************/
  /* Tranpose temporary state S */
  /******************************/
  
  for(i=0;i<numrows;i++){
    MPI_Alltoall(A_local[i],numrows,MPI_DOUBLE,tmparray,numrows,
		 MPI_DOUBLE,MPI_COMM_WORLD);
    cnt = 0;
    for(k=0;k<totalnodes;k++)
      for(j=0;j<numrows;j++)
	A_local_t[j][i+k*numrows] = tmparray[cnt++];
  }


  /***************************/
  /* Compute G = H*S = H*A*H */
  /***************************/


  tmps1 = dot(size,xold,xold);
  for(i=0;i<numrows;i++){
    tmps2 = dot(size,A_local_t[i],xold);
    for(j=0;j<size;j++)
      A_local_t[i][j] = A_local_t[i][j] - 2.0*tmps2*xold[j]/tmps1;
  }


  /********************************************************/
  /* Tranpose G so that it is stored by rows acrossed the */ 
  /* processors, just as the original A was stored        */
  /********************************************************/
  
  for(i=0;i<numrows;i++){
    MPI_Alltoall(A_local_t[i],numrows,MPI_DOUBLE,tmparray,numrows,
		 MPI_DOUBLE,MPI_COMM_WORLD);
    cnt = 0;
    for(k=0;k<totalnodes;k++)
      for(j=0;j<numrows;j++)
	A_local[j][i+k*numrows] = tmparray[cnt++];
  }



  if(mynode == 0){
    for(i=0;i<numrows;i++){
      for(j=0;j<size;j++)
	cout << A_local[i][j] << " ";
      cout << endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(mynode == 1){
    for(i=0;i<numrows;i++){
      for(j=0;j<size;j++)
	cout << A_local[i][j] << " ";
      cout << endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);

  for(i=0;i<numrows;i++){
    delete[] A_local[i];
    delete[] A_local_t[i];
  }
  delete[] A_local;
  delete[] A_local_t;
  delete[] tmparray;
  delete[] xold;
  delete[] xnew;

  MPI_Finalize();
}



void ChebyVandermonde(int npts, double *A, int row){
  int i,j;
  double * x = new double[npts];
  
  ChebyshevPoints(npts,x);
  
  for(j=0;j<npts;j++)
    A[j] = pow(x[row],j);
  
  delete[] x;
}




