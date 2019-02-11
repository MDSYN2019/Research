/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter7.cpp  (definition file)
* Header file:  SCchapter7.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/

#include "SCmathlib.h"
#include "SCchapter6.h"
#include "SCchapter7.h"
#ifdef MPISRC
#include <mpi.h>
#endif


void Diffusion_EF_CentralDifference(int N, double DN, double *uold, 
				    double *unew){
  
  for(int i=1;i<N-1;i++)
    unew[i] = uold[i] + DN*(uold[i+1] - 2.0*uold[i] + uold[i-1]);
  
  return;
}


void Diffusion_EB_CentralDifference(int N, double DN, double *uold, 
				    double *unew){
  
  ThomasAlgorithm(N-2,-DN,1.0+2.0*DN,-DN,&unew[1],&uold[1]);
  
  return;
}


int Diffusion_Jacobi(int N, double dx, double dt, 
		      double **A, double **q, double abstol){
  int i,j,k;
  int maxit = 100000;
  double sum;
  double ** Aold = CreateMatrix(N,N);
  
  double D = dt/(dx*dx);
  
  for(i=1; i<N-1; i++)
    for(j=1;j<N-1;j++)
      Aold[i][j] = 1.0;
  
  /* Boundary Conditions -- all zeros */
  for(i=0;i<N;i++){
    A[0][i] = 0.0;
    A[N-1][i] = 0.0;
    A[i][0] = 0.0;
    A[i][N-1] = 0.0;
  }

  for(k=0; k<maxit; k++){
    for(i = 1; i<N-1; i++){
      for(j=1; j<N-1; j++){
	A[i][j] = dt*q[i][j] + Aold[i][j] +
	  D*(Aold[i+1][j] + Aold[i][j+1] - 4.0*Aold[i][j] + 
	     Aold[i-1][j] + Aold[i][j-1]);
      }
    }
    
    sum = 0.0;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
	sum += (Aold[i][j]-A[i][j])*(Aold[i][j]-A[i][j]);
	Aold[i][j] = A[i][j];
      }
    }

    if(sqrt(sum)<abstol){
      DestroyMatrix(Aold,N,N);
      return k;
    }

  }
      
  cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
  DestroyMatrix(Aold,N,N);
  
  return maxit;
}



int Diffusion_GaussSeidel(int N, double dx, double dt, 
		      double **A, double **q, double abstol){
  int i,j,k;
  int maxit = 100000;
  double sum;
  double ** Aold = CreateMatrix(N,N);
  
  double D = dt/(dx*dx);
  
  for(i=1; i<N-1; i++)
    for(j=1;j<N-1;j++)
      Aold[i][j] = 1.0;
  
  /* Boundary Conditions -- all zeros */
  for(i=0;i<N;i++){
    A[0][i] = 0.0;
    A[N-1][i] = 0.0;
    A[i][0] = 0.0;
    A[i][N-1] = 0.0;
  } 

  for(k=0; k<maxit; k++){
    for(i = 1; i<N-1; i++){
      for(j=1; j<N-1; j++){
	A[i][j] = dt*q[i][j] + Aold[i][j] +
	  D*(Aold[i+1][j] + Aold[i][j+1] - 4.0*Aold[i][j] + 
	     A[i-1][j] + A[i][j-1]);
      }
    }
    
    sum = 0.0;
    for(i=0;i<N;i++){
      for(j=0;j<N;j++){
	sum += (Aold[i][j]-A[i][j])*(Aold[i][j]-A[i][j]);
	Aold[i][j] = A[i][j];
      }
    }

    if(sqrt(sum)<abstol){
      DestroyMatrix(Aold,N,N);
      return k;
    }

  }
      
  cerr << "Gauss Seidel: Maximum Number of Interations Reached Without Convergence\n";
  DestroyMatrix(Aold,N,N);
  
  return maxit;
}


int Jacobi(int N, double **A, double *x, double *b, double abstol){
  int i,j,k;
  int maxit = 100000;
  double sum1,sum2;
  double * xold = new double[N];
  
  //set initial guess to all 1.0
  for(i=0; i<N; i++){
    xold[i] = 1.0;
  }
  
  for(k=0; k<maxit; k++){
    for(i = 0; i<N; i++){
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i; j++)
	sum1 = sum1 + A[i][j]*xold[j];
      for(j=i+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (-sum1 - sum2 + b[i])/A[i][i];
    }

    if(sqrt(dot(N,x,xold))<abstol){
      delete[] xold;
      return k;
    }
    
    for(i=0; i<N; i++)
      xold[i] = x[i];
  }

  cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
  delete[] xold;
  
  return maxit;
}

#ifdef MPISRC
int Jacobi_P(int mynode, int numnodes, int N, double **A, double *x, double *b, double abstol){
  int i,j,k,i_global;
  int maxit = 100000;
  int rows_local,local_offset,last_rows_local,*count,*displacements;
  double sum1,sum2,*xold;
  double error_sum_local, error_sum_global;
  MPI_Status status;

  rows_local = (int) floor((double)N/numnodes);
  local_offset = mynode*rows_local;
  if(mynode == (numnodes-1)) 
    rows_local = N - rows_local*(numnodes-1);

  /*Distribute the Matrix and R.H.S. among the processors */
  if(mynode == 0){
    for(i=1;i<numnodes-1;i++){
      for(j=0;j<rows_local;j++)
	MPI_Send(A[i*rows_local+j],N,MPI_DOUBLE,i,j,MPI_COMM_WORLD);
      MPI_Send(b+i*rows_local,rows_local,MPI_DOUBLE,i,rows_local,
	       MPI_COMM_WORLD);
    }
    last_rows_local = N-rows_local*(numnodes-1);
    for(j=0;j<last_rows_local;j++)
      MPI_Send(A[(numnodes-1)*rows_local+j],N,MPI_DOUBLE,numnodes-1,j,
	       MPI_COMM_WORLD);
    MPI_Send(b+(numnodes-1)*rows_local,last_rows_local,MPI_DOUBLE,numnodes-1,
	     last_rows_local,MPI_COMM_WORLD);
  }
  else{
    A = CreateMatrix(rows_local,N);
    x = new double[rows_local];    
    b = new double[rows_local];
    for(i=0;i<rows_local;i++)
      MPI_Recv(A[i],N,MPI_DOUBLE,0,i,MPI_COMM_WORLD,&status);
    MPI_Recv(b,rows_local,MPI_DOUBLE,0,rows_local,MPI_COMM_WORLD,&status);
  }


  xold = new double[N];
  count = new int[numnodes];
  displacements = new int[numnodes];


  //set initial guess to all 1.0
  for(i=0; i<N; i++){
    xold[i] = 1.0;
  }

  for(i=0;i<numnodes;i++){
    count[i] = (int) floor((double)N/numnodes);
    displacements[i] = i*count[i];
  }
  count[numnodes-1] = N - ((int)floor((double)N/numnodes))*(numnodes-1);
  
  for(k=0; k<maxit; k++){
    error_sum_local = 0.0;
    for(i = 0; i<rows_local; i++){
      i_global = local_offset+i;
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i_global; j++)
	sum1 = sum1 + A[i][j]*xold[j];
      for(j=i_global+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (-sum1 - sum2 + b[i])/A[i][i_global];
      error_sum_local += (x[i]-xold[i_global])*(x[i]-xold[i_global]);
    }
    
    MPI_Allreduce(&error_sum_local,&error_sum_global,1,MPI_DOUBLE,
		  MPI_SUM,MPI_COMM_WORLD);
    MPI_Allgatherv(x,rows_local,MPI_DOUBLE,xold,count,displacements,
		   MPI_DOUBLE,MPI_COMM_WORLD);
    
    if(sqrt(error_sum_global)<abstol){
      if(mynode == 0){
	for(i=0;i<N;i++)
	  x[i] = xold[i];
      }
      else{
	DestroyMatrix(A,rows_local,N);
	delete[] x;
        delete[] b;
      }
      delete[] xold;
      delete[] count;
      delete[] displacements;
      return k;
    }
  }

  cerr << "Jacobi: Maximum Number of Interations Reached Without Convergence\n";
  if(mynode == 0){
    for(i=0;i<N;i++)
      x[i] = xold[i];
  }
  else{
    DestroyMatrix(A,rows_local,N);
    delete[] x;
    delete[] b;
  }
  delete[] xold;
  delete[] count;
  delete[] displacements;
  
  return maxit;
}


void Jacobi_P(int mynode, int numnodes, int N, double **A, double *x, double *b, int iterations){
  int i,j,k,i_global;
  int rows_local,local_offset,last_rows_local,*count,*displacements;
  double sum1,sum2,*xold;
  double error_sum_local, error_sum_global;
  MPI_Status status;

  rows_local = (int) floor((double)N/numnodes);
  local_offset = mynode*rows_local;
  if(mynode == (numnodes-1)) 
    rows_local = N - rows_local*(numnodes-1);

  /*Distribute the Matrix and R.H.S. among the processors */
  if(mynode == 0){
    for(i=1;i<numnodes-1;i++){
      for(j=0;j<rows_local;j++)
	MPI_Send(A[i*rows_local+j],N,MPI_DOUBLE,i,j,MPI_COMM_WORLD);
      MPI_Send(b+i*rows_local,rows_local,MPI_DOUBLE,i,rows_local,
	       MPI_COMM_WORLD);
    }
    last_rows_local = N-rows_local*(numnodes-1);
    for(j=0;j<last_rows_local;j++)
      MPI_Send(A[(numnodes-1)*rows_local+j],N,MPI_DOUBLE,numnodes-1,j,
	       MPI_COMM_WORLD);
    MPI_Send(b+(numnodes-1)*rows_local,last_rows_local,MPI_DOUBLE,numnodes-1,
	     last_rows_local,MPI_COMM_WORLD);
  }
  else{
    A = CreateMatrix(rows_local,N);
    x = new double[rows_local];    
    b = new double[rows_local];
    for(i=0;i<rows_local;i++)
      MPI_Recv(A[i],N,MPI_DOUBLE,0,i,MPI_COMM_WORLD,&status);
    MPI_Recv(b,rows_local,MPI_DOUBLE,0,rows_local,MPI_COMM_WORLD,&status);
  }


  xold = new double[N];
  count = new int[numnodes];
  displacements = new int[numnodes];


  //set initial guess to all 1.0
  for(i=0; i<N; i++){
    xold[i] = 1.0;
  }

  for(i=0;i<numnodes;i++){
    count[i] = (int) floor((double)N/numnodes);
    displacements[i] = i*count[i];
  }
  count[numnodes-1] = N - ((int)floor((double)N/numnodes))*(numnodes-1);
  
  for(k=0; k<iterations; k++){
    error_sum_local = 0.0;
    for(i = 0; i<rows_local; i++){
      i_global = local_offset+i;
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i_global; j++)
	sum1 = sum1 + A[i][j]*xold[j];
      for(j=i_global+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (-sum1 - sum2 + b[i])/A[i][i_global];
      error_sum_local += (x[i]-xold[i_global])*(x[i]-xold[i_global]);
    }
    
    MPI_Allreduce(&error_sum_local,&error_sum_global,1,MPI_DOUBLE,
		  MPI_SUM,MPI_COMM_WORLD);
    MPI_Allgatherv(x,rows_local,MPI_DOUBLE,xold,count,displacements,
		   MPI_DOUBLE,MPI_COMM_WORLD);
  }

  if(mynode == 0){
    for(i=0;i<N;i++)
      x[i] = xold[i];
  }
  else{
    DestroyMatrix(A,rows_local,N);
    delete[] x;
    delete[] b;
  }
  delete[] xold;
  delete[] count;
  delete[] displacements;
  
  return;
}
#endif



void Jacobi(int N, double **A, double *x, double *b, int iterations){
  int i,j,k;
  double sum1,sum2;
  double * xold = new double[N];
  
  //set initial guess
  for(i=0; i<N; i++){
    xold[i] = 1.0;
  }
  
  for(k=0; k<iterations; k++){
    for(i = 0; i<N; i++){
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i; j++)
	sum1 = sum1 + A[i][j]*xold[j];
      for(j=i+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (-sum1 - sum2 + b[i])/A[i][i];
    }
    for(i=0; i<N; i++)
      xold[i] = x[i];
  }

  delete[] xold;
  return;
}



int GaussSeidel(int N, double **A, double *x, double *b, double abstol){
  int i,j,k;
  int maxit = 100000;
  double sum1,sum2;
  double * xold = new double[N];
  
  //set initial guess
  for(i=0; i<N; i++){
    x[i] =    1.0e0;
    xold[i] = 1.0e0;
  }
  
  for(k=0; k<maxit; k++){
    for(i = 0; i<N; i++){
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i; j++)
	sum1 = sum1 + A[i][j]*x[j];
      for(j=i+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (-sum1 - sum2 + b[i])/A[i][i];
    }
    
    if(sqrt(dot(N,x,xold))<abstol){
      delete[] xold;
      return k;
    }
    
    for(i=0; i<N; i++)
      xold[i] = x[i];
  }
  
  cerr << "GaussSeidel: Maximum Number of Interations Reached Without Convergence\n";
  delete[] xold;

  return maxit;
}



void GaussSeidel(int N, double **A, double *x, double *b, int iterations){
  int i,j,k;
  double sum1,sum2;
  double * xold = new double[N];
  
  //set initial guess
  for(i=0; i<N; i++){
    x[i] = 1.0e0;
    xold[i] = 1.0e0;
  }
  
  for(k=0; k<iterations; k++){
    for(i = 0; i<N; i++){
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i; j++)
	sum1 = sum1 + A[i][j]*x[j];
      for(j=i+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (-sum1 - sum2 + b[i])/A[i][i];
    }    
    for(i=0; i<N; i++)
      xold[i] = x[i];
  }
  
  delete[] xold;
  return;
}



int SOR(double omega, int N, double **A, double *x, double *b, double abstol){
  int i,j,k;
  int maxit = 100000;
  double sum1,sum2;
  double * xold = new double[N];
  
  //set initial guess
  for(i=0; i<N; i++){
    x[i] = 1.0e0;
    xold[i] = 1.0e0;
  }
  
  for(k=0; k<maxit; k++){
    for(i = 0; i<N; i++){
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i; j++)
	sum1 = sum1 + A[i][j]*x[j];
      for(j=i+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (1.0-omega)*xold[i] + omega*(-sum1 - sum2 + b[i])/A[i][i];
    }
    
    if(sqrt(dot(N,x,xold))<abstol){
      delete[] xold;
      return k;
    }
    
    for(i=0; i<N; i++)
      xold[i] = x[i];
  }
  
  cerr << "SOR: Maximum Number of Interations Reached Without Convergence\n";
  delete[] xold;

  return maxit;
}



void SOR(double omega, int N, double **A, double *x, 
	 double *b, int iterations){
  int i,j,k;
  double sum1,sum2;
  double * xold = new double[N];
  
  //set initial guess
  for(i=0; i<N; i++){
    x[i] = 1.0e0;
    xold[i] = 1.0e0;
  }
  
  for(k=0; k<iterations; k++){
    for(i = 0; i<N; i++){
      sum1 = 0.0; sum2 = 0.0;
      for(j=0; j < i; j++)
	sum1 = sum1 + A[i][j]*x[j];
      for(j=i+1; j < N; j++)
	sum2 = sum2 + A[i][j]*xold[j];
      
      x[i] = (1.0-omega)*xold[i] + omega*(-sum1 - sum2 + b[i])/A[i][i];
    }    
    for(i=0; i<N; i++)
      xold[i] = x[i];
  }
  
  delete[] xold;
  return;
}










