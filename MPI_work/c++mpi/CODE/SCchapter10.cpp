/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter10.cpp  (definition file)
* Header file:  SCchapter10.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#include "SCchapter10.h"
#include "SCchapter9.h"


double PowerMethod(SCMatrix &A, SCVector &x){
  const double tolerance = 1.0e-7;
  return(PowerMethod(A,x,tolerance));
}


double PowerMethod(SCMatrix &A, SCVector &x, double tolerance){
  const int maxit = 1000;
  double eigval,eold=1.0;
  SCVector xnew(x.Dimension());

  x.Initialize(1.0);

  for(int i=0; i<maxit; i++){
    xnew = A*x;
    eigval = xnew.ElementofMaxMod();

    if(fabs(eigval-eold)<tolerance)
      break;

    x = xnew/eigval;
    eold = eigval;
  }

  return eigval;
}

double ISPowerMethod(SCMatrix &A, SCVector &x){
  const int maxit = 100;
  const double tol = 1.0e-7;
  double tmp,sigma1,sigma2;
  SCMatrix Amod(A.Rows());
  SCVector y(A.Rows());

  // Calculation of the Rayleigh Quotient
  y = A*x;
  sigma1 = dot(y,x)/dot(x,x);

  for(int k=0;k<maxit;k++){
    Amod = A;
     
    for(int i=0;i<Amod.Rows();i++)
      Amod(i,i) = Amod(i,i) - sigma1;
    
    // Use LU solver from Chapter 9 to obtain solution x 
    GaussElimination(Amod, x, 1);

    // Normalize with the element of maximum modulus 
    tmp = x.ElementofMaxMod();
    x = x/tmp;
    
    y = A*x;
    sigma2 = dot(y,x)/dot(x,x);
    
    if(fabs(sigma1-sigma2) < tol)
      return sigma2;
    
    sigma1 = sigma2;
  }
  
  cout << "ISPowerMethod - Maximum number of iterations reached\n";
  return sigma2;
}


void TDQREigensolver(int N, double *a, double *b, 
		     double *lambda, double **Q){

  if(N == 1){
    Q[0][0] = 1.0;
    lambda[0] = a[0];
  }
  else{
    int i,j,k,N1,N2,cnt;
    N1 = N/2;
    N2 = N-N1;
    double * d = new double[N]; 
    double ** Q1 = CreateMatrix(N1,N1);
    double ** Q2 = CreateMatrix(N2,N2);
    
    a[N1-1] = a[N1-1] - b[N1-1];
    a[N1]   = a[N1]   - b[N1-1];
    TDQREigensolver(N1,a,b,d,Q1);
    TDQREigensolver(N2,&a[N1],&b[N1],&d[N1],Q2);

    double * xi = new double[N];
    cnt = 0;
    for(i=0;i<N1;i++)
      xi[cnt++] = Q1[N1-1][i];
    for(i=0;i<N2;i++)
      xi[cnt++] = Q2[0][i];

    SolveSecularEq(b[N1-1],N,d,xi,lambda);
    
    for(i=0;i<N1;i++){
      for(j=0;j<N;j++){
	Q[i][j] = 0.0;
	for(k=0;k<N1;k++)
	  Q[i][j] += Q1[i][k]*xi[k]/(d[k]-lambda[j]);
      }
    }
    for(i=0;i<N2;i++){
      for(j=0;j<N;j++){
	Q[N1+i][j] = 0.0;
	for(k=0;k<N2;k++)
	  Q[i+N1][j] += Q2[i][k]*xi[N1+k]/(d[N1+k]-lambda[j]);
      }
    }

    double sum;
    for(i=0;i<N;i++){
      sum = 0.0;
      for(j=0;j<N;j++)
	sum+= Q[j][i]*Q[j][i];
      sum = sqrt(sum);
      for(j=0;j<N;j++)
	Q[j][i] = Q[j][i]/sum;
    }
    
    delete[] xi;
    delete[] d;
    DestroyMatrix(Q1,N1,N1);
    DestroyMatrix(Q2,N2,N2);
  }
  return;
}


void SolveSecularEq(double bm, int N, double *d, 
		    double *xi, double * lambda){
  
  double xl,xr,xm;
  double yl,yr,ym;
  const double offset = 1.0e-5;
  const double tol = 1.0e-6;


  for(int i=0;i<N-1;i++){
    xl = d[i] + offset;
    yl = SecularEq(bm,N,d,xi,xl);
    xr = d[i+1] - offset;
    yr = SecularEq(bm,N,d,xi,xr);
    xm = 0.5*(xl+xr);
    ym = SecularEq(bm,N,d,xi,xm);
   
    if(yl*yr > 0){
      lambda[i] = xl;
      continue;
    }
    
    while(fabs(ym)>tol){
      if(yl*ym<0)
	xr = xm;
      else
	xl = xm;
      xm = 0.5*(xl+xr);
      ym = SecularEq(bm,N,d,xi,xm);
    }
    lambda[i] = xm;
  }

  xl = d[N-1] + offset;
  yl = SecularEq(bm,N,d,xi,xl);
  xr = 2*d[N-1];
  yr = SecularEq(bm,N,d,xi,xr);
  xm = 0.5*(xl+xr);
  ym = SecularEq(bm,N,d,xi,xm);
  
  if(yl*yr > 0){
    lambda[N-1] = xl;
  }
  else{
    while(fabs(ym)>tol){
      if(yl*ym<0)
	xr = xm;
      else
	xl = xm;
      xm = 0.5*(xl+xr);
      ym = SecularEq(bm,N,d,xi,xm);
    }
    lambda[N-1] = xm;
  }
}


double SecularEq(double bm, int N, double *d, double *xi, double x){
  double sum = 0.0e0;
  for(int i=0;i<N;i++)
    sum += xi[i]*xi[i]/(d[i]-x);
  sum = bm*sum + 1.0e0;
  return sum;
}



