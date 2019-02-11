
#include <iostream.h>
#include <iomanip.h>
#include "SCmathlib.h"
#include "SCchapter7.h"


int main(int argc, char * argv[]){
  int i,j;
  int npts = 40,itcount;
  double **A,**q,dx,dt,truesoln,sum;

  dx = 1.0/(npts-1);
  dt = 0.25*dx*dx;

  A = CreateMatrix(npts,npts);
  q = CreateMatrix(npts,npts);

  for(i=0;i<npts;i++)
    for(j=0;j<npts;j++)
      q[i][j] = -2.0*M_PI*M_PI*sin(M_PI*dx*i)*sin(M_PI*dx*j);
    
  itcount = Diffusion_Jacobi(npts,dx,dt,A,q,1.0e-14);
 
  sum = 0.0;
  for(i=0;i<npts;i++){
    for(j=0;j<npts;j++){
      truesoln = -sin(M_PI*dx*i)*sin(M_PI*dx*j);
      sum += (A[i][j] - truesoln)*(A[i][j] - truesoln);
    }
  }
  
  cout << setprecision(5) << setiosflags(ios::scientific);
  cout << "Jacobi: L2 Error is " << sqrt(sum) << " in " << itcount << " iterations" << endl;
  
  DestroyMatrix(A,npts,npts);
  DestroyMatrix(q,npts,npts);
}

