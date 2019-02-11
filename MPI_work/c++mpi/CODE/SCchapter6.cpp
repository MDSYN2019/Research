/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter6.cpp  (definition file)
* Header file:  SCchapter6.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/

#include "SCchapter6.h"

#ifdef MPISRC
#include <mpi.h>
#endif

void ThomasAlgorithm(int N, double b, double a, double c, double *x, double *q){
  int i;
  double *l,*u,*d,*y;

  l = new double[N];
  u = new double[N];
  d = new double[N];
  y = new double[N];

  /* LU Decomposition */
  d[0] = a; 
  u[0] = c;
  for(i=0;i<N-2;i++){
    l[i] = b/d[i];
    d[i+1] = a - l[i]*u[i];
    u[i+1] = c;
  }
  l[N-2] = b/d[N-2];
  d[N-1] = a - l[N-2]*u[N-2];

  /* Forward Substitution [L][y] = [q] */
  y[0] = q[0];
  for(i=1;i<N;i++)
    y[i] = q[i] - l[i-1]*y[i-1];

  /* Backward Substitution [U][x] = [y] */
  x[N-1] = y[N-1]/d[N-1];
  for(i=N-2;i>=0;i--)
    x[i] = (y[i] - u[i]*x[i+1])/d[i];

  delete[] l;
  delete[] u;
  delete[] d;
  delete[] y;
  return;
}


void ThomasAlgorithm(int N, double *b, double *a, double *c, double *x, double *q){
  double *l,*u,*d;

  l = new double[N];
  u = new double[N];
  d = new double[N];

  ThomasAlgorithmLU(N,b,a,c,l,u,d);
  ThomasAlgorithmSolve(N,l,u,d,x,q);

  delete[] l;
  delete[] u;
  delete[] d;
  return;
}


void ThomasAlgorithmLU(int N, double *b, double *a, double *c, double *l, double *u, double *d){
  int i;

  /* LU Decomposition */
  d[0] = a[0]; 
  u[0] = c[0];
  for(i=0;i<N-2;i++){
    l[i] = b[i]/d[i];
    d[i+1] = a[i+1] - l[i]*u[i];
    u[i+1] = c[i+1];
  }
  l[N-2] = b[N-2]/d[N-2];
  d[N-1] = a[N-1] - l[N-2]*u[N-2];

  return;
}

void ThomasAlgorithmSolve(int N, double *l, double *u, double *d, double *x, double *q){
  int i;
  double *y = new double[N];

  /* Forward Substitution [L][y] = [q] */
  y[0] = q[0];
  for(i=1;i<N;i++)
    y[i] = q[i] - l[i-1]*y[i-1];

  /* Backward Substitution [U][x] = [y] */
  x[N-1] = y[N-1]/d[N-1];
  for(i=N-2;i>=0;i--)
    x[i] = (y[i] - u[i]*x[i+1])/d[i];
  
  delete[] y;
  return;
}


void ThomasAlgorithm_per(int N, double b, double a, double c, double *x, double *q){
  int i;
  double *x1,*x2,*q2;

  x1 = new double[N-1];
  x2 = new double[N-1];
  q2 = new double[N-1];

  /* Prepare secondary q */
  for(i=0;i<N-1;i++)
    q2[i] = 0.0;
  q2[0] = -b;
  q2[N-2] = -c;
  
  ThomasAlgorithm(N-1,b,a,c,x1,q);
  ThomasAlgorithm(N-1,b,a,c,x2,q2);
  
  x[N-1] = (q[N-1] - c*x1[0] - b*x1[N-2])/
    (a + c*x2[0] + b*x2[N-2]);
  
  for(i=0;i<N-1;i++)
    x[i] = x1[i] + x2[i]*x[N-1];

  delete[] x1;
  delete[] x2;
  delete[] q2;
}

void ThomasAlgorithm_per(int N, double *b, double *a, double *c, double *x, double *q){
  int i;
  double *x1,*x2,*q2;

  x1 = new double[N-1];
  x2 = new double[N-1];
  q2 = new double[N-1];

  /* Prepare secondary q */
  for(i=0;i<N-1;i++)
    q2[i] = 0.0;
  q2[0] = -b[N-1];
  q2[N-2] = -c[N-2];
  
  ThomasAlgorithm(N-1,b,a,c,x1,q);
  ThomasAlgorithm(N-1,b,a,c,x2,q2);
  
  x[N-1] = (q[N-1] - c[N-1]*x1[0] - b[N-2]*x1[N-2])/
    (a[N-1] + c[N-1]*x2[0] + b[N-2]*x2[N-2]);
  
  for(i=0;i<N-1;i++)
    x[i] = x1[i] + x2[i]*x[N-1];

  delete[] x1;
  delete[] x2;
  delete[] q2;
}



#ifdef MPISRC
void ThomasAlgorithm_P(int mynode, int numnodes, 
		       int N, double *b, double *a, double *c, double *x, double *q){
  int i;
  int rows_local,local_offset;
  double S[2][2],T[2][2],s1tmp,s2tmp;
  double *l,*d,*y;
  MPI_Status status;

  l = new double[N];
  d = new double[N];
  y = new double[N];
  
  for(i=0;i<N;i++)
    l[i] = d[i] = y[i] = 0.0;

  S[0][0] = S[1][1] = 1.0;
  S[1][0] = S[0][1] = 0.0;

  rows_local = (int) floor((double)N/numnodes);
  local_offset = mynode*rows_local;
  
  if(mynode==0){
    s1tmp = a[local_offset]*S[0][0];
    S[1][0] = S[0][0];
    S[1][1] = S[0][1];
    S[0][1] = a[local_offset]*S[0][1];
    S[0][0] = s1tmp;
    for(i=1;i<rows_local;i++){
      s1tmp = a[i+local_offset]*S[0][0] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][0]; 
      s2tmp = a[i+local_offset]*S[0][1] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][1];
      S[1][0] = S[0][0];
      S[1][1] = S[0][1];
      S[0][0] = s1tmp;
      S[0][1] = s2tmp;
    }
  }
  else{
    for(i=0;i<rows_local;i++){
      s1tmp = a[i+local_offset]*S[0][0] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][0]; 
      s2tmp = a[i+local_offset]*S[0][1] - b[i+local_offset-1]*c[i+local_offset-1]*S[1][1];
      S[1][0] = S[0][0];
      S[1][1] = S[0][1];
      S[0][0] = s1tmp;
      S[0][1] = s2tmp;
    }
  }

  for(i=0; i<=log2(numnodes);i++){
    if(mynode+pow(2.0,i) < numnodes)
      MPI_Send(S,4,MPI_DOUBLE,int(mynode+pow(2.0,i)),0,MPI_COMM_WORLD);
    if(mynode-pow(2.0,i)>=0){
      MPI_Recv(T,4,MPI_DOUBLE,int(mynode-pow(2.0,i)),0,MPI_COMM_WORLD,&status);
      s1tmp = S[0][0]*T[0][0] + S[0][1]*T[1][0];
      S[0][1] = S[0][0]*T[0][1] + S[0][1]*T[1][1];
      S[0][0] = s1tmp;
      s1tmp = S[1][0]*T[0][0] + S[1][1]*T[1][0];
      S[1][1] = S[1][0]*T[0][1] + S[1][1]*T[1][1];
      S[1][0] = s1tmp;
    }
  }
  
  d[local_offset+rows_local-1] = (S[0][0] + S[0][1])/(S[1][0] + S[1][1]);
  if(mynode == 0){
    MPI_Send(&d[local_offset+rows_local-1],1,MPI_DOUBLE,1,0,MPI_COMM_WORLD);
  }
  else{
    MPI_Recv(&d[local_offset-1],1,MPI_DOUBLE,mynode-1,0,MPI_COMM_WORLD,&status);
    if(mynode != numnodes-1)
      MPI_Send(&d[local_offset+rows_local-1],1,MPI_DOUBLE,mynode+1,0,MPI_COMM_WORLD);
  }

  
  if(mynode == 0){
    l[0] = 0;
    d[0] = a[0];
    for(i=1;i<rows_local-1;i++){
      l[local_offset+i] = b[local_offset+i-1]/d[local_offset+i-1];
      d[local_offset+i] = a[local_offset+i] - l[local_offset+i]*c[local_offset+i-1];
    }
    l[local_offset+rows_local-1] = b[local_offset+rows_local-2]/d[local_offset+rows_local-2];
  }
  else{
    for(i=0;i<rows_local-1;i++){
      l[local_offset+i] = b[local_offset+i-1]/d[local_offset+i-1];
      d[local_offset+i] = a[local_offset+i] - l[local_offset+i]*c[local_offset+i-1];
    }
    l[local_offset+rows_local-1] = b[local_offset+rows_local-2]/d[local_offset+rows_local-2];
  }

  /***********************************************************************************/
  
  if(mynode>0)
    d[local_offset-1] = 0;
  
  double * tmp = new double[N];
  for(i=0;i<N;i++)
    tmp[i] = d[i];
  MPI_Allreduce(tmp,d,N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(i=0;i<N;i++)
    tmp[i] = l[i];
  MPI_Allreduce(tmp,l,N,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  delete[] tmp;

  if(mynode ==0){    
    /* Forward Substitution [L][y] = [q] */
    y[0] = q[0];
    for(i=1;i<N;i++)
      y[i] = q[i] - l[i]*y[i-1];
    
    /* Backward Substitution [U][x] = [y] */
    x[N-1] = y[N-1]/d[N-1];
    for(i=N-2;i>=0;i--)
      x[i] = (y[i] - c[i]*x[i+1])/d[i];
    
  }
  
  delete[] l;
  delete[] y;
  delete[] d;
  return;
}
#endif



