/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter3.cpp  (definition file)
* Header file:  SCchapter3.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#include "SCchapter3.h"


void NewtonDiffTable(int npts, double *xpts, double *funcvals, 
			 double * newton_coeffs){
  int i;  
  for(i=0;i<npts;i++)
    newton_coeffs[i] = NewtonDiffFunction(0,i, xpts, funcvals);

  return;
}


double NewtonDiffFunction(int start_index, int ending_index, 
			  double * xpts, double * funcvals){
  double val;
  
  int diff = ending_index-start_index;
  
  if(diff == 0){
    val = funcvals[start_index];
  }
  else{
    val = (NewtonDiffFunction(start_index,ending_index-1,xpts,funcvals) -
	   NewtonDiffFunction(start_index+1,ending_index,xpts,funcvals))/
      (xpts[start_index]-xpts[ending_index]);
  }
  
  return val;
}


double NewtonInterpolant(double x, int npts, double * xpts, 
			 double * newton_coeffs){
  int i,j;
  double sum = 0.0, xval;
  
  for(i=0;i<npts;i++){
    xval = 1.0;
    for(j=0;j<i;j++)
      xval = xval*(x-xpts[j]);
    sum = sum + newton_coeffs[i]*xval;
  }
  
  return sum;
}


double LagrangeInterpolant(double x, int npts, double *xpts, 
			   double * funcvals){
  int i;
  double sum = 0.0;

  for(i=0;i<npts;i++){
    sum = sum + funcvals[i]*LagrangePoly(x,i,npts,xpts);
  }
  return sum;
}


double LagrangePoly(double x, int pt, int npts, double * xpts){
  int i;
  double h=1.0;
  
  for(i=0;i<pt;i++)
      h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);

  for(i=pt+1;i<npts;i++)
    h = h * (x - xpts[i])/(xpts[pt]-xpts[i]);
  
  return h;
}


double ChebyshevPoly(int degree, double x){
  double value;

  switch(degree){
  case 0: 
    value = 1.0;
    break;
  case 1:
    value = x;
    break;
  default:
    value = 2.0*x*ChebyshevPoly(degree-1,x) - ChebyshevPoly(degree-2,x);
  }
  
  return value;
}


void ChebyshevPoints(int npts, double * x){
  int k;
  double c1 = 4.0*atan(1.0)/(npts*2.0);

  for(k=0;k<npts;k++)
    x[k] = cos(c1*(2.0*k+1.0));

  return;
}

void  CreateGrid_EvenlySpaced(int npts, double *x, double a, double b){
  double dx = (b-a)/(npts-1.0);
  
  for(int i=0;i<npts;i++)
    x[i] = a + i*dx;

  return;
}

void  CreateGrid_ChebyshevPts(int npts, double *x, double a, double b){
  double bma = b-a;
  double apb = a+b;

  /* First get the Chebyshev points on [-1,1] */
  ChebyshevPoints(npts, x); 

  /* Scale the points to the interval [a,b] */
  for(int i=0;i<npts;i++)
    x[i] = 0.5*bma*x[i] + apb;

  return;

}

void LS_ComputeCoeffs(int npts, double *xpts, double *funcvals, 
		      int ndeg, double *alpha, double *beta, double *lscoeffs){
  int i,j;
  double xi,tmpd;
  double * gamma = new double[ndeg+1];  

  //////////////////////////
  // Compute average first
  xi = 0.0;
  for(i=0;i<npts;i++){
    xi += xpts[i];
  }
  xi /= (double) npts;
  /////////////////////////

  gamma[0] = npts;
  alpha[0] = beta[0] = 0.0;
  alpha[1] = xi;
  
  for(j=1;j<=ndeg-1;j++){
    gamma[j]   = 0.0;
    alpha[j+1] = 0.0;
    for(i=0;i<npts;i++){
      tmpd = LS_OrthoPoly(j,xpts[i],alpha,beta);
      gamma[j] += tmpd*tmpd;
      alpha[j+1] += xpts[i]*tmpd*tmpd;
    }
    alpha[j+1] /= gamma[j];
    beta[j] = gamma[j]/gamma[j-1];
  }
  
  gamma[ndeg] = 0.0;
  
  for(i=0;i<npts;i++){
    tmpd = LS_OrthoPoly(ndeg,xpts[i],alpha,beta);
    gamma[ndeg] += tmpd*tmpd;
  }
  
  beta[ndeg] = gamma[ndeg]/gamma[ndeg-1];
  
  for(j=0;j<=ndeg;j++){
    lscoeffs[j] = 0.0;
    for(i=0;i<npts;i++)
      lscoeffs[j] = lscoeffs[j] + funcvals[i]*LS_OrthoPoly(j,xpts[i],alpha,beta);
    lscoeffs[j] /= gamma[j];
  }
  
  delete[] gamma;

  return;
}


double LS_OrthoPoly(int j, double x, double *alpha, double *beta){
  double value;
  
  switch(j){
  case 0:
    value = 1.0;
    break;
  case 1:
    value = x - alpha[j];
    break;
  default:
    value = (x-alpha[j])*LS_OrthoPoly(j-1,x,alpha,beta) - 
      beta[j-1]*LS_OrthoPoly(j-2,x,alpha,beta);
    break;
  }
  
  return value;
}

double LSApproximatingPoly(int ndeg, double x, double *alpha, 
			   double *beta, double *lscoeffs){
  double value = 0.0;

  for(int i=0;i<=ndeg;i++)
    value += lscoeffs[i]*LS_OrthoPoly(i, x, alpha, beta);


  return value;

}


///////////////////////////////////////////////////////////////////////////
LSPoly::LSPoly(){
  ndeg = 0;
  alpha = NULL;
  beta = NULL;
  lscoeffs = NULL;
}
    
LSPoly::~LSPoly(){
  delete[] alpha;
  delete[] beta;
  delete[] lscoeffs;
  
  ndeg = 0;
  alpha = NULL;
  beta = NULL;
  lscoeffs = NULL;
}

double LSPoly::LSPolyOrtho(int j,double x){
  double value;
  
  switch(j){
  case 0:
    value = 1.0;
    break;
  case 1:
    value = x - alpha[j];
    break;
  default:
    value = (x-alpha[j])*LSPolyOrtho(j-1,x) - 
      beta[j-1]*LSPolyOrtho(j-2,x);
    break;
  }
  
  return value;
}

int LSPoly::Initialize(int npts, int in_ndeg, double * xpts, double * funcvals){
  int i,j;
  double xi,tmpd;

  if(alpha!=NULL){
    cerr << "An Error Has Occurred:: LSPoly has already been Initialized" << endl;
    return 0;
  }
  
  ndeg = in_ndeg;
  
  /* Storage for the this object */
  lscoeffs = new double[ndeg+1];
  alpha = new double[ndeg+1];
  beta =  new double[ndeg+1];  
  

  /* Storage for the just this method */
  double * gamma = new double[ndeg+1];  

  //////////////////////////
  // Compute average first
  xi = 0.0;
  for(i=0;i<npts;i++){
    xi += xpts[i];
  }
  xi /= (double) npts;
  /////////////////////////

  gamma[0] = npts;
  alpha[0] = beta[0] = 0.0;
  alpha[1] = xi;
  
  for(j=1;j<=ndeg-1;j++){
    gamma[j]   = 0.0;
    alpha[j+1] = 0.0;
    for(i=0;i<npts;i++){
      tmpd = LS_OrthoPoly(j,xpts[i],alpha,beta);
      gamma[j] += tmpd*tmpd;
      alpha[j+1] += xpts[i]*tmpd*tmpd;
    }
    alpha[j+1] /= gamma[j];
    beta[j] = gamma[j]/gamma[j-1];
  }
  
  gamma[ndeg] = 0.0;
  
  for(i=0;i<npts;i++){
    tmpd = LSPolyOrtho(ndeg,xpts[i]);
    gamma[ndeg] += tmpd*tmpd;
  }
  
  beta[ndeg] = gamma[ndeg]/gamma[ndeg-1];
  
  for(j=0;j<=ndeg;j++){
    lscoeffs[j] = 0.0;
    for(i=0;i<npts;i++)
      lscoeffs[j] = lscoeffs[j] + funcvals[i]*LSPolyOrtho(j,xpts[i]);
    lscoeffs[j] /= gamma[j];
  }
  
  delete[] gamma;

  return 1;
}


void LSPoly::PrintCoeffs(){
  cout << endl;
  cout << "*********************************" << endl;
  cout << "i\talpha\tbeta\tlscoeffs" << endl;

  for(int j=0;j<=ndeg;j++){
    cout << j << "\t" << alpha[j] << "\t" << beta[j] << "\t" << lscoeffs[j] << endl;
  }
 
  cout << "*********************************" << endl << endl;
 
  return;
}



double LSPoly::Evaluate(double x){
  double value = 0.0;

  for(int i=0;i<=ndeg;i++)
    value += lscoeffs[i]*LSPolyOrtho(i,x);

  return value;
}

///////////////////////////////////////////////////////////////////////////


double Square_2dInterpolant(SCPoint x, int npts, double *funcvals){
  double value = 0.;
  double h[4];

  if(npts != 4){
    cerr << "An Error Has occurred in Square_2dInterpolant -- Invalid npts given" << endl;
    return value;
  }
  
  h[0] = 0.5*(1.0-x(0));
  h[1] = 0.5*(1.0+x(0));
  h[2] = 0.5*(1.0-x(1));
  h[3] = 0.5*(1.0+x(1));

  value = funcvals[0]*h[0]*h[2] + funcvals[1]*h[1]*h[2] + 
    funcvals[2]*h[1]*h[3] + funcvals[3]*h[0]*h[3];

  return value;
}








