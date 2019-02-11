/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter4.cpp  (definition file)
* Header file:  SCchapter4.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#include "SCchapter4.h"

double SquareRoot(double value, double guess, int iterations){
  int i;
  double xn = guess;

  for(i=0; i<iterations;i++)
    xn = 0.5*(xn + value/xn);
  
  return xn;
}


double NewtonRaphson(double x0, double (*func)(double), 
		     double (*func_der)(double), 
		     int max_iter, int multiplicity){
  double x,p = (double) multiplicity;
  const double tolerance = 1.0e-14;

  for(int i=0;i<max_iter;i++){
    x = x0 - p*func(x0)/func_der(x0);
    if(fabs(func(x))<tolerance) break;
    x0 = x;
  }
  
  return x;
}


double NewtonRaphson(double x0, double (*func)(double), 
		     double (*func_der)(double), 
		     double (*func_secondder)(double), int max_iter, 
		     int multiplicity){
  double x,p = (double) multiplicity;
  const double tolerance = 1.0e-14;
  
  cout << "-----------------------------------" << endl;
  cout << "x \t f(x) \t g'(x)" << endl;
  
  for(int i=0;i<max_iter;i++){
    cout << x0 << "\t" << func(x0) << "\t" <<  (1-p) + 
      p*func(x0)*func_secondder(x0)/pow(func_der(x0),2) << endl;;
    x = x0 - p*func(x0)/func_der(x0);
    if(fabs(func(x))<tolerance) break;
    x0 = x;
  }
  
  return x;
}


double NewtonRaphson(double x0, double (*func)(double), 
		     double (*func_der)(double), int iterations){
  double x, p = 1; //here, p stands for the multiplicity; 
                   //we assume standard p=1
  const double tolerance = 1.0e-14;

  for(int i=0;i<iterations;i++){
    x = x0 - p*func(x0)/func_der(x0);
    if(fabs(func(x))<tolerance) break;
    x0 = x;
  }
  
  return x;
}



SCVector ConjugateGradient(SCMatrix A, SCVector b, SCVector x0){
  int dim = x0.Dimension();
  const double tolerance = 1.0e-14;
  SCVector x(dim),r(dim),v(dim),z(dim);
  double c,t,d;

  x = x0;
  r = b - A*x;
  v = r;
  c = dot(r,r);
 

  for(int i=0;i<dim;i++){
    if(sqrt(dot(v,v))<tolerance){
      cerr << "An error has occurred in ConjugateGradient: execution of function terminated" << endl;
      break;
    }
    z = A*v;
    t = c/dot(v,z);
    x = x + t*v;
    r = r - t*z;
    d = dot(r,r);
    if(sqrt(d) < tolerance)
      break;
    v = r + (d/c)*v;
    c = d;
  }
  
  return x;
} 


double MidpointRule(int level, double xleft, double xright, double (*func)(double)){
  int i, nsteps = (int) pow(2.0,level)-1;
  double h = (xright-xleft)/pow(2.0,level);
  double sum = 0.0;

  for(i=0;i<=nsteps;i++)
    sum += func(xleft + (i+0.5)*h);
  sum *= h;

  return sum;
}


double TrapezoidRule(int level, double xleft, double xright, 
		     double (*func)(double)){
  int i, nsteps = (int) pow(2.0,level)-1;
  double h = (xright-xleft)/pow(2.0,level);
  double sum = 0.0;

  for(i=1;i<=nsteps;i++)
    sum += func(xleft + i*h);
  sum *= 2;

  /* Add the first and the last point to the summation */
  sum += func(xleft) + func(xright);
  sum *= 0.5*h;

  return sum;
}


double Romberg(int m, int k, double xleft, double xright, 
	       double (*func)(double)){
  double RI,I1,I2;
  double coeff = pow(4.0,m);

  if(k < m){
    cerr << "ROMBERG::Value of k must be >= m; setting k=m\n";
    k = m;
  }

  if(m==0){
    RI = TrapezoidRule(k,xleft,xright,func);
  }
  else{
    I1 = Romberg(m-1,k,  xleft,xright,func);
    I2 = Romberg(m-1,k-1,xleft,xright,func);
    RI = (coeff*I1 - I2)/(coeff-1.0);
  }
  
  return RI;
}


double SimpsonsRule(int m, double xleft, double xright, 
		    double (*func)(double)){
  int i; 
  double h = (xright-xleft)/(2.0*m);
  double sum = 0.0, x;

  for(i=1;i<=m;i++){
    x = xleft + (2*i-1)*h;
    sum += func(x-h) + 4.0*func(x) + func(x+h);
  }
  sum *= (h/3.0);

  return sum;
}



double JacobiPoly(int degree, double x, double alpha, double beta){
  double value;
  double tmp,degm1;
  double a1=0.,a2=0.,a3=0.,a4=0.;

  switch(degree){
  case 0:
    value = 1.0;
    break;
  case 1:
    value = 0.5*(alpha-beta+(alpha+beta+2.0)*x);
    break;
  default:
    degm1 = degree-1.0; 
    tmp = 2.0*degm1+alpha+beta;
    a1= 2.0*(degm1+1)*(degm1+alpha+beta+1)*tmp;
    a2= (tmp+1)*(alpha*alpha-beta*beta);
    a3= tmp*(tmp+1.0)*(tmp+2.0);
    a4= 2.0*(degm1+alpha)*(degm1+beta)*(tmp+2.0);

    value = ((a2+a3*x)*JacobiPoly(degree-1,x,alpha,beta)-
	     a4*JacobiPoly(degree-2,x,alpha,beta))/a1;
  }
  
  return value;

}


double JacobiPolyDerivative(int degree, double x, double alpha, double beta){
  double value;
  double tmp;
  double b1,b2,b3;

  switch(degree){
  case 0:
    value = 0.0;
    break;
  default:
    tmp = 2.0*degree+alpha+beta;
    b1 = tmp*(1.0-x*x);
    b2 = degree*(alpha-beta-tmp*x);
    b3 = 2.0*(degree+alpha)*(degree+beta);

    value = (b2*JacobiPoly(degree,x,alpha,beta)+ b3*JacobiPoly(degree-1,x,alpha,beta))/b1;
  }
  
  return value;
   
}


void JacobiZeros(int degree, double *z, double alpha, double beta){  
  int i,j,k;
  const int maxit = 30;
  const double EPS = 1.0e-14;
  double   dth = M_PI/(2.0*degree);
  double   poly,pder,rlast=0.0;
  double   sum,delr,r;
  double one = 1.0, two = 2.0;

  // If the degree of the polynomial is zero (or less), then there are no roots
  if(degree<=0)
    return;

  
  for(k = 0; k < degree; k++){
    r = -cos((two*k + one) * dth);
    if(k) r = 0.5*(r + rlast);
    
    for(j = 1; j < maxit; ++j){
      poly = JacobiPoly(degree,r,alpha,beta);
      pder = JacobiPolyDerivative(degree,r,alpha,beta);
      
      sum = 0.0;
      for(i = 0; i < k; ++i) 
	sum += one/(r - z[i]);
      
      delr = -poly / (pder - sum * poly);
      r   += delr;
      if( fabs(delr) < EPS ) break;
    }
    z[k]  = r;
    rlast = r;
  }

  return;
}



void JacobiZW(int degree, double * z, double *w, double alpha, double beta){
  int i;
  double fac, one = 1.0, two = 2.0, apb = alpha + beta;

  JacobiZeros(degree, z, alpha, beta);  

  for(i=0;i<degree;i++)
    w[i] = JacobiPolyDerivative(degree,z[i],alpha,beta);
      
  fac  = pow(two,apb + one)*GammaF(alpha + degree + one)*GammaF(beta + degree + one);
  fac /= GammaF(degree + one)*GammaF(apb + degree + one);
  
  for(i = 0; i < degree; ++i) 
    w[i] = fac/(w[i]*w[i]*(one-z[i]*z[i]));
  
  return;
}


double HermitePoly(int degree, double x){
  double value;

  switch(degree){
  case 0:
    value = 1.0;
    break;
  case 1:
    value = 2.0*x;
    break;
  default:
    value = 2.0*(x*HermitePoly(degree-1,x)-(degree-1.0)*HermitePoly(degree-2,x));
  }
  return value;
}


double HermitePolyDerivative(int degree, double x){
  double value;

  switch(degree){
  case 0:
    value = 0.0;
    break;
  case 1:
    value = 2.0;
    break;
  default:
    value = 2.0*degree*HermitePoly(degree-1,x);
  }
  return value;
}

void HermiteZeros(int degree, double *z){  
  int i,j,k;
  const int maxit = 30;
  const double EPS = 1.0e-14;
  double   dth = M_PI/(2.0*degree);
  double   poly,pder,rlast=0.0;
  double   sum,delr,r;
  double one = 1.0, two = 2.0;

  // If the degree of the polynomial is zero (or less), then there are no roots
  if(degree<=0)
    return;
  
  for(k = 0; k < degree; k++){
    r = -cos((two*k + one) * dth);
    if(k) r = 0.5*(r + rlast);
    
    for(j = 1; j < maxit; ++j){
      poly = HermitePoly(degree,r);
      pder = HermitePolyDerivative(degree,r);
      
      sum = 0.0;
      for(i = 0; i < k; ++i) 
	sum += one/(r - z[i]);
      
      delr = -poly / (pder - sum * poly);
      r   += delr;
      if( fabs(delr) < EPS ) break;
    }

    if(j==maxit) 
      cerr << "Maximum iteration of Newton Iteration has been reached in function: HermiteZeros\n";

    z[k]  = r;
    rlast = r;
  }

  return;
}


void HermiteZW(int degree, double * z, double *w){
  int i;
  double tmp;
  double spi = sqrt(M_PI);
  double nfac = Factorial(degree);
  double twonp1 = pow(2.0,degree+1);
  
  HermiteZeros(degree, z);  
  
  for(i = 0; i < degree; ++i){ 
    tmp = HermitePolyDerivative(degree,z[i]);
    w[i] = spi*nfac*twonp1/(tmp*tmp);  
  }
  return;
}


double LaguerrePoly(int degree, double x){
  double value;
  double degm1 = degree-1.0;

  switch(degree){
  case 0: 
    value = 1.0;
    break;
  case 1:
    value = 1-x;
    break;
  default:
    value = ((2.0*degm1+1.0-x)*LaguerrePoly(degree-1,x) - degm1*LaguerrePoly(degree-2,x))/degree;
  }
  
  return value;
}


double LaguerrePolyDerivative(int degree, double x){
  double value;

  switch(degree){
  case 0:
    value = 0.0;
    break;
  case 1:
    value = -1.0;
    break;
  default:
    value = LaguerrePolyDerivative(degree-1,x) - LaguerrePoly(degree-1,x);
  }
  return value;
}

void LaguerreZeros(int degree, double *z){  
  int i,j,k;
  const int maxit = 30;
  const double EPS = 1.0e-14;
  double   dth = M_PI/(2.0*degree);
  double   poly,pder,rlast=0.0;
  double   sum,delr,r;
  double one = 1.0, two = 2.0;

  // If the degree of the polynomial is zero (or less), then there are no roots
  if(degree<=0)
    return;
  
  for(k = 0; k < degree; k++){
    r = -cos((two*k + one) * dth);
    if(k) r = 0.5*(r + rlast);
    
    for(j = 1; j < maxit; ++j){
      poly = LaguerrePoly(degree,r);
      pder = LaguerrePolyDerivative(degree,r);
      
      sum = 0.0;
      for(i = 0; i < k; ++i) 
	sum += one/(r - z[i]);
      
      delr = -poly / (pder - sum * poly);
      r   += delr;
      if( fabs(delr) < EPS ) break;
    }

    if(j==maxit) 
      cerr << "Maximum iteration of Newton Iteration has been reached in function: HermiteZeros\n";

    z[k]  = r;
    rlast = r;
  }

  return;
}


void LaguerreZW(int degree, double * z, double *w){
  int i;
  double prefactor = 1.0/degree;

  LaguerreZeros(degree, z);  
  
  for(i = 0; i < degree; ++i){
    w[i] =  -prefactor/(LaguerrePoly(degree-1,z[i])*LaguerrePolyDerivative(degree,z[i]));
  }  

  return;
}










