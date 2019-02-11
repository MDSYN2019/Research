
#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"
#include<mpi.h>

const int rows_per_proc = 40;
const double c1 = 20.0;
const double c2 = 1000.0;
const double tol = 1.0e-14;

int main(int argc, char *argv[]){
  int i,j,k;
  int mynode, totalnodes, totalsize, offset;
  MPI_Status status;
  double sum,local_sum,c,d,alpha,beta;
  double ** A, *q, *x, *grid;

  double *p,*z,*r,*mr;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);

  totalsize = totalnodes*rows_per_proc;

  p = new double[totalsize];
  z = new double[rows_per_proc];
  r = new double[rows_per_proc];
  mr = new double[rows_per_proc];

  x = new double[rows_per_proc];
  q = new double[rows_per_proc];
  grid = new double[rows_per_proc];

  double dx = 1.0/(totalsize+1);
  for(i=0;i<rows_per_proc;i++){
    grid[i] = dx*(1+rows_per_proc*mynode+i);
    q[i] = -dx*dx*sin(2.0*M_PI*grid[i])*
      (-4.0*M_PI*M_PI - c2*exp(c1*(grid[i]-0.5)*(grid[i]-0.5)));
    x[i] = 1.0;
  }

  /* Part 2 */

  A = new double*[rows_per_proc];
  for(i=0;i<rows_per_proc;i++){
    A[i] = new double[totalsize];
    for(j=0;j<totalsize;j++)
      A[i][j] = 0.0;
  }

  if(mynode==0){
    A[0][0] = 2.0 + dx*dx*c2*exp(c1*(grid[0]-0.5)*(grid[0]-0.5)); 
    A[0][1] = -1.0;
    for(i=1;i<rows_per_proc;i++){
      A[i][i] = 2.0 + dx*dx*c2*exp(c1*(grid[i]-0.5)*(grid[i]-0.5)); 
      A[i][i-1] = -1.0;
      A[i][i+1] = -1.0;
    }
  }
  else if(mynode == (totalnodes-1)){
    A[rows_per_proc-1][totalsize-1] = 2.0 + 
        dx*dx*c2*exp(c1*(grid[rows_per_proc-1]-0.5)*(grid[rows_per_proc-1]-0.5));  
    A[rows_per_proc-1][totalsize-2] = -1.0;
    for(i=0;i<rows_per_proc-1;i++){
      offset = rows_per_proc*mynode + i;
      A[i][offset] = 2.0 + dx*dx*c2*exp(c1*(grid[i]-0.5)*(grid[i]-0.5)); 
      A[i][offset-1] = -1.0;
      A[i][offset+1] = -1.0;
    }
  }
  else{
    for(i=0;i<rows_per_proc;i++){
      offset = rows_per_proc*mynode + i;
      A[i][offset] = 2.0 + dx*dx*c2*exp(c1*(grid[i]-0.5)*(grid[i]-0.5)); 
      A[i][offset-1] = -1.0;
      A[i][offset+1] = -1.0;
    }
  }

  /* Part 3 */

  offset = mynode*rows_per_proc;
 
  for(i=0;i<totalsize;i++)
    p[i] = 1.0;

  for(i=0;i<rows_per_proc;i++){
    r[i] = q[i] - dot(totalsize,A[i],p); //calculation of residual
    mr[i] = r[i]/A[i][offset+i];   //calculation of modified residual
  }
  
  local_sum = dot(rows_per_proc,mr,r);
  MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  c = sum;
  
  MPI_Allgather(mr,rows_per_proc,MPI_DOUBLE,p,rows_per_proc,MPI_DOUBLE,MPI_COMM_WORLD);


  /* Part 4 */

  for(k=0;k<totalsize;k++){

    for(i=0;i<rows_per_proc;i++)
      z[i] = dot(totalsize,A[i],p);
    
    local_sum = dot(rows_per_proc,z,p+offset);

    MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    alpha = c/sum;

    for(i=0;i<rows_per_proc;i++){
      x[i] = x[i] + alpha*p[offset+i];
      r[i] = r[i] - alpha*z[i];
    }

    /* Preconditioning Stage */
    for(i=0;i<rows_per_proc;i++)
      mr[i] = r[i]/A[i][offset+i];
    
    local_sum = dot(rows_per_proc,mr,r);

    MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    
    d = sum; //contains inner product of residual and modified residual

    local_sum = dot(rows_per_proc,r,r);

    MPI_Allreduce(&local_sum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    //sum now contains inner product of residual and residual

    if(mynode == 0)
      cout << k << "\t" << "dot(mr,r) = " << d << "\t" << "dot(r,r) = " << sum << endl;
 
    if(fabs(d) < tol) break;
    if(fabs(sum) < tol) break;
    
    beta = d/c;

    for(i=0;i<rows_per_proc;i++)
      z[i] = mr[i] + beta*p[i+offset];

    MPI_Allgather(z,rows_per_proc,MPI_DOUBLE,p,rows_per_proc,MPI_DOUBLE,MPI_COMM_WORLD);

    c = d;
    
  }

  delete[] p;
  delete[] z;
  delete[] r;
  delete[] mr;
  delete[] x;
  delete[] q;
  delete[] grid;
  for(i=0;i<rows_per_proc;i++)
    delete[] A[i];
  delete[] A;
  
  MPI_Finalize();
}






