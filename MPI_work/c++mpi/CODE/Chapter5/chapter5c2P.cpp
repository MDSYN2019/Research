
#include <iostream.h>
#include <iomanip.h>
#include <mpi.h>
#include "SCchapter5.h"

double func(double x);
double func_first_der(double x);
double func_second_der(double x);
 

int main(int argc, char * argv[]){
  const int levels = 10;
  const double global_a = 0.0;
  const double global_b = 1.0;
  int i,j,global_npts,local_npts;
  int mynode, totalnodes;
  double dx,dxp,ux_error,uxx_error;
  double *u,*u_x,*u_xx;
  double local_a;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes);
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode);
  

  if(mynode == 0)
    cout << "npts\tError (First Deriv)\tError (Second Deriv)\n";

  for(i=2;i<levels+2;i++){
    global_npts = (int) pow(2,i);    
    local_npts  = global_npts/totalnodes;
    global_npts = local_npts*totalnodes;
    dx = (global_b-global_a)/(global_npts-1);

    local_a = global_a + dx*local_npts*mynode;
    
    u = new double[local_npts];
    u_x = new double[local_npts];
    u_xx = new double[local_npts];
    
    for(j=0;j<local_npts;j++)
      u[j] = func(local_a + j*dx);
    
    SO_FirstDeriv_1DP  (local_npts,dx,u,u_x,mynode,totalnodes);
    SO_SecondDeriv_1DP (local_npts,dx,u,u_xx,mynode,totalnodes);

    ux_error=0.0;
    uxx_error=0.0;
    for(j=0;j<local_npts;j++){
      ux_error  += dx*pow((u_x[j]-func_first_der(local_a + j*dx)),2);
      uxx_error += dx*pow((u_xx[j]-func_second_der(local_a + j*dx)),2);
    }
    
    double answer1,answer2;
    MPI_Reduce(&ux_error,&answer1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&uxx_error,&answer2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    
     if(mynode==0){
       cout << setprecision(10) << setiosflags(ios::scientific);
       cout << global_npts << "\t" << sqrt(answer1) << "\t" << sqrt(answer2) << endl;
    }
    delete[] u;
    delete[] u_x;
    delete[] u_xx;
  }

  MPI_Finalize();

}

double func(double x){
  return(x*x*x*x);
}

double func_first_der(double x){
  return(4*x*x*x);
}

double func_second_der(double x){
  return(12*x*x);
}

