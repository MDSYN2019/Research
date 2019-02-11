

#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"
#include "SCchapter10.h"
#include<mpi.h>


// Global variable to set size of the system
const int size = 50;

int main(int argc, char *argv[]){
  int i,j,k,ll,m,isum,ioffset,cnt,N1,N2;
  int mynode, totalnodes;
  MPI_Status status;
  double bn,*tmpd,**tmpdd;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  // Set up storage on each process
  double * a = new double[size];
  double * b = new double[size];
  double * lambda = new double[size];
  double ** Q =  CreateMatrix(size,size);
  double ** Q1 = CreateMatrix(size,size);
  double ** Q2 = CreateMatrix(size,size);
  double * xi = new double[size];
  double * d  = new double[size];

  int ** index = ICreateMatrix(totalnodes,totalnodes);
  double ** adjust = CreateMatrix(totalnodes,totalnodes);

  //Form the Matrix A
  for(i=0;i<size;i++){
    a[i] = i+1.5;
    b[i] = 0.3;
  }

  //Set up recursive partitioning of the matrix.
  //Each process will solve for some subset of the problem

  index[0][0] = size; 
  for(i=0;i<log2(totalnodes);i++){
    isum = 0;
    for(j=0;j<pow(2.0,i);j++){
      index[i+1][2*j] = index[i][j]/2;
      index[i+1][2*j+1] = index[i][j] - index[i+1][2*j];
      isum += index[i+1][2*j];
      adjust[i][j] = b[isum-1];
      a[isum-1] = a[isum-1] - b[isum-1];
      a[isum]   = a[isum]   - b[isum-1];
      isum += index[i+1][2*j+1];
    }
  }


  // Each process solves recursively for its subpart of
  // the problem.

  ioffset = (int) log2(totalnodes);
  isum = 0;
  for(k=0;k<mynode;k++)
    isum += index[ioffset][k];
  
  TDQREigensolver(index[ioffset][mynode],&a[isum],&b[isum],d,Q1);


  // Fan-in algorithm to finish solving system

  for(i=0;i<log2(totalnodes);i++){
    isum = 0; cnt = 0;
    for(j=0;j<totalnodes;j+=(int)pow(2.0,i)){
      if(mynode == j){
	if(isum%2==0){
	  MPI_Recv(d+index[ioffset][isum],index[ioffset][isum+1],MPI_DOUBLE,
		   j+(int)pow(2.0,i),1,MPI_COMM_WORLD, &status);
	  for(k=0;k<index[ioffset][isum+1];k++)
	    MPI_Recv(Q2[k],index[ioffset][isum+1],MPI_DOUBLE,
		     j+(int)pow(2.0,i),1,MPI_COMM_WORLD, &status);
	  N1 = index[ioffset][isum];
	  N2 = index[ioffset][isum+1];
	  bn = adjust[ioffset-1][cnt++];
	}
	else{
	  MPI_Send(d,index[ioffset][isum],MPI_DOUBLE,j-(int)pow(2.0,i),
		   1,MPI_COMM_WORLD);
	  for(k=0;k<index[ioffset][isum];k++)
	    MPI_Send(Q1[k],index[ioffset][isum],MPI_DOUBLE,j-(int)pow(2.0,i),
		     1,MPI_COMM_WORLD);
	}
      }
      isum++;
    }
    
    for(j=0;j<totalnodes;j+=(int)pow(2.0,i+1)){
      if(mynode == j){
	cnt = 0;
	for(k=0;k<N1;k++)
	  xi[cnt++] = Q1[N1-1][k];
	for(k=0;k<N2;k++)
	  xi[cnt++] = Q2[0][k];

	// Solve for the secular equation to 
	// obtain eigenvalues
	SolveSecularEq(bn,N1+N2,d,xi,lambda);
    
	// Form the Q matrix from Q1 and Q2
	for(k=0;k<N1;k++){
	  for(ll=0;ll<N1+N2;ll++){
	    Q[k][ll] = 0.0;
	    for(m=0;m<N1;m++)
	      Q[k][ll] += Q1[k][m]*xi[m]/(d[m]-lambda[ll]);
	  }
	}
	for(k=0;k<N2;k++){
	  for(ll=0;ll<N1+N2;ll++){
	    Q[N1+k][ll] = 0.0;
	    for(m=0;m<N2;m++)
	      Q[k+N1][ll] += Q2[k][m]*xi[N1+m]/(d[N1+m]-lambda[ll]);
	  }
	}


	// Normalize the Q matrix so that each eigenvector
	// has length one
	double sum;
	for(k=0;k<N1+N2;k++){
	  sum = 0.0;
	  for(ll=0;ll<N1+N2;ll++)
	    sum+= Q[ll][k]*Q[ll][k];
	  sum = sqrt(sum);
	  for(ll=0;ll<N1+N2;ll++)
	    Q[ll][k] = Q[ll][k]/sum;
	}
	
	// Swap d and lambda arrays for use in the
	// next part of the fan-in algorithm
	tmpd = d;
	d = lambda;
	lambda = tmpd;

	// Swap Q and Q1 for use in the
	// next part of the fan-in algorithm	
	tmpdd = Q1;
	Q1 = Q;
	Q = tmpdd;
      }          
    }
    ioffset = ioffset - 1;
  }


  if(mynode==0){
    cout << "The eigenvalues are: " << endl;
    for(k=0;k<size;k++)
      cout << d[k] << endl;
  }


  MPI_Barrier(MPI_COMM_WORLD);

  delete[] a;
  delete[] b;
  delete[] lambda;
  delete[] xi;
  delete[] d;

  DestroyMatrix(Q,size,size);
  DestroyMatrix(Q1,size,size);
  DestroyMatrix(Q2,size,size);
  DestroyMatrix(adjust,totalnodes,totalnodes);
  IDestroyMatrix(index,totalnodes,totalnodes);

  MPI_Finalize();
}


