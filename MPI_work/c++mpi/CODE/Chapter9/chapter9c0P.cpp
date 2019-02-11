/* PART 1*/
#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"
#include "SCchapter3.h"
#include<mpi.h>

void ChebyVandermonde(int npts, double *A, int row);

// Global variable to set size of the system
const int size = 10;

int main(int argc, char *argv[]){
  int i,j,k,index;
  int mynode, totalnodes;
  double scaling;
  MPI_Status status;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  int numrows = size/totalnodes;
  double **A_local = new double*[numrows];
  int * myrows = new int[numrows];



  /* PART 2 */
  double * xpts = new double[size];
  ChebyshevPoints(size,xpts);
  
  for(i=0;i<numrows;i++){
    A_local[i] = new double[size+1];
    index = mynode + totalnodes*i;
    myrows[i] = index;
    ChebyVandermonde(size,A_local[i],index);

    // Set-up right-hand-side as the Runge function
    A_local[i][size] = 1.0/(1.0+25.0*xpts[index]*xpts[index]);
  }
  delete[] xpts;


  double * tmp = new double[size+1];
  double * x = new double[size];

  /* PART 3 */
  /* Gaussian Elimination of the augmented matrix */
  int cnt = 0;

  for(i=0;i<size-1;i++){
    if(i == myrows[cnt]){
      MPI_Bcast(A_local[cnt],size+1,MPI_DOUBLE,mynode,MPI_COMM_WORLD);
      for(j=0;j<size+1;j++)
	tmp[j] = A_local[cnt][j];
      cnt++;
    }
    else{
      MPI_Bcast(tmp,size+1,MPI_DOUBLE,i%totalnodes,MPI_COMM_WORLD);
    }
    for(j=cnt;j<numrows;j++){
      scaling = A_local[j][i]/tmp[i];
      for(k=i;k<size+1;k++)
	A_local[j][k] = A_local[j][k] - scaling*tmp[k]; 
    }
  }

  /* PART 4 */

  /* On each processor, initialize the value of x as equal to
     the modified (by Gaussian elimination) right-hand-side if
     that information is on the processor, otherwise initialize
     to zero. 
  */

  cnt = 0;
  for(i=0;i<size;i++){
    if(i==myrows[cnt]){
      x[i] = A_local[cnt][size];
      cnt++;
    }
    else
      x[i] = 0;
  }
  

  /* PART 5 */
  /* Backsolve to find the solution x */
  cnt = numrows-1;
  for(i=size-1;i>0;i--){
    if(cnt>=0){
      if(i == myrows[cnt]){
	x[i] = x[i]/A_local[cnt][i];
	MPI_Bcast(x+i,1,MPI_DOUBLE,mynode,MPI_COMM_WORLD);
	cnt--;
      }
      else
	MPI_Bcast(x+i,1,MPI_DOUBLE,i%totalnodes,MPI_COMM_WORLD);
    }
    else
      MPI_Bcast(x+i,1,MPI_DOUBLE,i%totalnodes,MPI_COMM_WORLD);
    
    for(j=0;j<=cnt;j++)
      x[myrows[j]] = x[myrows[j]] - A_local[j][i]*x[i];
  }
  
  if(mynode==0){
    x[0] = x[0]/A_local[cnt][0];
    MPI_Bcast(x,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  }
  else
    MPI_Bcast(x,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
  
  
  /* PART 6 */

  if(mynode==0){
    for(i=0;i<size;i++)
      cout << x[i] << endl;
  }
  
  delete[] tmp;
  delete[] myrows;
  for(i=0;i<numrows;i++)
    delete[] A_local[i];
  delete[] A_local;

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





