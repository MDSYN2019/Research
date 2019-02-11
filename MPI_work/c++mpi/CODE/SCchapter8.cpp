/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter8.cpp  (definition file)
* Header file:  SCchapter8.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/

#include "SCchapter6.h"
#include "SCchapter8.h"

/* Euler Forward Central Difference Scheme - Periodic */
void EF_CentralDifference(int N, double CFL, double *uold, double *unew){

  for(int i=1; i<N-1; i++){
    unew[i] = uold[i] - 0.5*CFL*(uold[i+1]-uold[i-1]);
  }
  
  unew[0] = uold[0] - 0.5*CFL*(uold[1]-uold[N-1]);
  unew[N-1] = uold[N-1] - 0.5*CFL*(uold[0]-uold[N-2]);

  return;
}


/* Euler Forward Central Difference Scheme - Periodic */
void EF_FirstOrderUpwind(int N, double CFL, double *uold, double *unew){
  
  for(int i=1; i<N; i++)
    unew[i] = uold[i] - CFL*(uold[i]-uold[i-1]);
  
  unew[0] = uold[0] - CFL*(uold[0]-uold[N-1]);
}


/* Lax-Friedrichs Scheme - Periodic */
void LaxFriedrichs(int N, double CFL, double *uold, double *unew){
  
  for(int i=1; i<N-1; i++){
    unew[i] = 0.5*(uold[i+1] + uold[i-1] - CFL*(uold[i+1]-uold[i-1]));
  }
  
  unew[0] = 0.5*(uold[1] + uold[N-1] - 0.5*CFL*(uold[1]-uold[N-1]));
  unew[N-1] = 0.5*(uold[0] + uold[N-2] - 0.5*CFL*(uold[0]-uold[N-2]));

  return;
}


double MinMod(double x, double y){
  double value = fabs(x);
  value = (value < fabs(y))?value:fabs(y);
  value *= 0.5*(Sign(x)+Sign(y));
  return value;
}

double MinMod(double x, double y, double z){
  double value = fabs(x);
  value = (value < fabs(y))?value:fabs(y);
  value = (value < fabs(z))?value:fabs(z);
  value *= (Sign(x)+Sign(y)+Sign(z))/3.0;
  return value;
}

void LaxFriedrichsTadmor(int N, double CFL, double *uold, double *unew){
  
  for(int i=2;i<N-2;i++){
    unew[i] = 0.5*(uold[i+1]+uold[i-1]) + 
      0.25*(1.0-CFL*CFL)*(MinMod(uold[i]-uold[i-1],uold[i-1]-uold[i-2])-
			  MinMod(uold[i+2]-uold[i+1],uold[i+1]-uold[i])) -
      0.5*CFL*(uold[i+1]-uold[i-1]);
  }
  
  unew[0] = 0.5*(uold[1]+uold[N-1]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(uold[0]-uold[N-1],uold[N-1]-
			       uold[N-2])-MinMod(uold[2]-uold[1],uold[1]-uold[0])) -
    0.5*CFL*(uold[1]-uold[N-1]);
  unew[1] = 0.5*(uold[2]+uold[0]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(uold[1]-uold[0],uold[0]-uold[N-1])-
			MinMod(uold[3]-uold[2],uold[2]-uold[1])) -
    0.5*CFL*(uold[2]-uold[0]);
  unew[N-1] = 0.5*(uold[0]+uold[N-2]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(uold[N-1]-uold[N-2],
			       uold[N-2]-uold[N-3])-
			MinMod(uold[1]-uold[0],uold[0]-uold[N-1])) -
    0.5*CFL*(uold[0]-uold[N-2]);
  unew[N-2] = 0.5*(uold[N-1]+uold[N-3]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(uold[N-2]-uold[N-3],
			       uold[N-3]-uold[N-4])-
			MinMod(uold[0]-uold[N-1],
			       uold[N-1]-uold[N-2])) -
    0.5*CFL*(uold[N-1]-uold[N-3]);
  
}

void LaxFriedrichsTadmor(int N, double CFL, double *uold, double *unew, double alpha){
  for(int i=2;i<N-2;i++){
    unew[i] = 0.5*(uold[i+1]+uold[i-1]) + 
      0.25*(1.0-CFL*CFL)*(MinMod(alpha*(uold[i]-uold[i-1]),
				 alpha*(uold[i-1]-uold[i-2]),
				 0.5*(uold[i]-uold[i-2]))-
			  MinMod(alpha*(uold[i+2]-uold[i+1]),
				 alpha*(uold[i+1]-uold[i]),
				 0.5*(uold[i+2]-uold[i]))) -
      0.5*CFL*(uold[i+1]-uold[i-1]);
  }
  
  unew[0] = 0.5*(uold[1]+uold[N-1]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(alpha*(uold[0]-uold[N-1]),
			       alpha*(uold[N-1]-uold[N-2]),
			       0.5*(uold[0]-uold[N-2]))-
			MinMod(alpha*(uold[2]-uold[1]),
			       alpha*(uold[1]-uold[0]),
			       0.5*(uold[2]-uold[0]))) -
    0.5*CFL*(uold[1]-uold[N-1]);
  unew[1] = 0.5*(uold[2]+uold[0]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(alpha*(uold[1]-uold[0]),
			       alpha*(uold[0]-uold[N-1]),
			       0.5*(uold[1]-uold[N-1]))-
			MinMod(alpha*(uold[3]-uold[2]),
			       alpha*(uold[2]-uold[1]),
			       0.5*(uold[3]-uold[1]))) -
    0.5*CFL*(uold[2]-uold[0]);
  unew[N-1] = 0.5*(uold[0]+uold[N-2]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(alpha*(uold[N-1]-uold[N-2]),
			       alpha*(uold[N-2]-uold[N-3]),
			       0.5*(uold[N-1]-uold[N-3]))-
			MinMod(alpha*(uold[1]-uold[0]),
			       alpha*(uold[0]-uold[N-1]),
			       0.5*(uold[1]-uold[N-1]))) -
    0.5*CFL*(uold[0]-uold[N-2]);
  unew[N-2] = 0.5*(uold[N-1]+uold[N-3]) + 
    0.25*(1.0-CFL*CFL)*(MinMod(alpha*(uold[N-2]-uold[N-3]),
			       alpha*(uold[N-3]-uold[N-4]),
			       0.5*(uold[N-2]-uold[N-4]))-
			MinMod(alpha*(uold[0]-uold[N-1]),
			       alpha*(uold[N-1]-uold[N-2]),
			       0.5*(uold[N-1]-uold[N-2]))) -
    0.5*CFL*(uold[N-1]-uold[N-3]);
  
}





/* Euler Foward Second Order Upwind Scheme - Periodic */
void EF_SecondOrderUpwind(int N, double CFL, double *uold, double *unew){

  for(int i=2; i<N; i++){
    unew[i] = (1.0e0 - CFL*((3.0/2.0) - CFL/2.0))*uold[i] + CFL*(2.0-CFL)*uold[i-1] 
      + (CFL/2.0)*(CFL-1.0)*uold[i-2];
  }

  unew[0] = (1.0 - CFL*((3.0/2.0) - CFL/2.0))*uold[0] + CFL*(2.0-CFL)*uold[N-1] + 
    (CFL/2.0)*(CFL-1.0)*uold[N-2];

  unew[1] = (1.0 - CFL*((3.0/2.0) - CFL/2.0))*uold[1] + CFL*(2.0-CFL)*uold[0] + 
    (CFL/2.0)*(CFL-1.0)*uold[N-1];

  return;
}


/* Lax-Wendroff Scheme - Periodic */
void LaxWendroff(int N, double CFL, double *uold, double *unew){
  
  for(int i=1; i<N-1; i++){
    unew[i] = uold[i] - (CFL/2.0)*(uold[i+1]-uold[i-1])+(CFL*CFL/2.0)*
      (uold[i+1]-2.0*uold[i]+uold[i-1]);
  }

  unew[0] = uold[0] - (CFL/2.0)*(uold[1]-uold[N-1])+(CFL*CFL/2.0)*
    (uold[1]-2.0*uold[0]+uold[N-1]);

  unew[N-1] = uold[N-1] - (CFL/2.0)*(uold[0]-uold[N-2])+(CFL*CFL/2.0)*
    (uold[0]-2.0*uold[N-1]+uold[N-2]);
  
  return;
}

/* Leap-Frog Central Difference Scheme - Periodic */
void LeapFrog_CentralDifference(int N, double CFL, double **uold, double *unew){

  for(int i=1; i<N-1; i++)
    unew[i] = uold[1][i] - CFL*(uold[0][i+1]-uold[0][i-1]); 

  unew[0] = uold[1][0] - CFL*(uold[0][1]-uold[0][N-1]);
  unew[N-1] = uold[1][N-1] - CFL*(uold[0][0]-uold[0][N-2]);
  
  return;
}


/* Second Order Adams-Bashforth Central Difference Scheme - Periodic */
void AB2_CentralDifference(int N, double CFL, double **uold, double *unew){

  for(int i=1; i<N-1; i++){
    unew[i] = uold[0][i] - 0.75*CFL*(uold[0][i+1]-uold[0][i-1]) + 0.25*CFL*
      (uold[1][i+1]-uold[1][i-1]);
  }

  unew[0] = uold[0][0] - 0.75*CFL*(uold[0][1]-uold[0][N-1]) + 0.25*CFL*
    (uold[1][1]-uold[1][N-1]);
  
  unew[N-1] = uold[0][N-1] - 0.75*CFL*(uold[0][0]-uold[0][N-2]) + 0.25*CFL*
    (uold[1][0]-uold[1][N-2]);

  return;
}

/* Second Order Adams-Bashforth Central Difference Scheme - Periodic */
void AB3_CentralDifference(int N, double CFL, double **uold, double *unew){
  double c1 = 23.0/24.0,
    c2 = 16.0/24.0,
    c3 = 5.0/24.0;
  
  for(int i=1; i<N-1; i++){
    unew[i] = uold[0][i] - c1*CFL*(uold[0][i+1]-uold[0][i-1]) + 
      c2*CFL*(uold[1][i+1]-uold[1][i-1]) - c3*CFL*(uold[2][i+1]-uold[2][i-1]);
  }

  unew[0] = uold[0][0] - c1*CFL*(uold[0][1]-uold[0][N-1]) + 
    c2*CFL*(uold[1][1]-uold[1][N-1]) - c3*CFL*(uold[2][1]-uold[2][N-1]);
  
  unew[N-1] = uold[0][N-1] - c1*CFL*(uold[0][0]-uold[0][N-2]) + 
    c2*CFL*(uold[1][0]-uold[1][N-2]) - c3*CFL*(uold[2][0]-uold[2][N-2]);

  return;
}

/* Crank-Nicolson Central Difference Scheme - Periodic */
void CrankNicolson_CentralDifference(int N, double CFL, double *uold, double *unew){
  double c = 0.25*CFL;
  double *q = new double[N];
  
  for(int i=1; i<N-1; i++)
    q[i] = uold[i] - 0.25*CFL*(uold[i+1]-uold[i-1]);
  
  q[0] =   uold[0] - 0.25*CFL*(uold[1]-uold[N-1]);
  q[N-1] = uold[N-1] - 0.25*CFL*(uold[0]-uold[N-2]);

  ThomasAlgorithm_per(N,-c,1.0e0,c,unew,q);
  
  delete[] q;

  return;
}

/* Second Order Adams-Bashforth Central Difference Scheme - Non-Periodic */
void AB2_CentralDifferenceNP(int N, double CFL, double **uold, double *unew){

  for(int i=1; i<N-1; i++){
    unew[i] = uold[0][i] - 0.75*CFL*(uold[0][i+1]-uold[0][i-1]) + 0.25*CFL*
      (uold[1][i+1]-uold[1][i-1]);
  }
  
  unew[N-1] = uold[0][N-1] - 1.5*CFL*(uold[0][N-1]-uold[0][N-2]) + 0.5*CFL*
    (uold[1][N-1]-uold[1][N-2]);
  
  return;
}

/* Leap-Frog Central Difference Scheme - Non-Periodic */
void LeapFrog_CentralDifferenceNP(int N, double CFL, double **uold, double *unew){

  for(int i=1; i<N-1; i++)
    unew[i] = uold[1][i] - CFL*(uold[0][i+1]-uold[0][i-1]); 
  
  unew[N-1] = uold[1][N-1] - 2.0*CFL*(uold[0][N-1]-uold[0][N-2]);
  
  return;
}


/* Crank-Nicolson Central Difference Scheme - Periodic */
void CrankNicolson_CentralDifference_AdvectionDiffusion(int N, double CFL,
			     double DN, double *uold, double *unew){
  double c1 = 0.25*CFL - 0.5*DN,
    c2 = -0.25*CFL - 0.5*DN;
  double *q = new double[N-1];
  
  for(int i=1; i<N-1; i++)
    q[i-1] = uold[i] - 0.25*CFL*(uold[i+1]-uold[i-1]) 
      + 0.5*DN*(uold[i+1] - 2.0*uold[i] + uold[i-1]);
  
  ThomasAlgorithm_per(N-1,c2,1.0e0+DN,c1,&unew[1],q);
  
  unew[N-1] = unew[N-2];

  delete[] q;

  return;
}
