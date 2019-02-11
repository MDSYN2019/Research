/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter5.cpp  (definition file)
* Header file:  SCchapter5.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/

#include "SCchapter5.h"
#ifdef MPISRC
#include <mpi.h>
#endif


/* Second Order First Derivative in 1D */
void SO_FirstDeriv_1D    (int npts, double dx, double *u, double *u_x){
  double two_invdx = 1.0/(2.0*dx);
  
  for(int i=1;i<npts-1;i++)
    u_x[i] = (u[i+1]-u[i-1])*two_invdx;
    
  // Foward Differencing
  u_x[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx;
  
  // Backward Differencing
  u_x[npts-1] = (3.0*u[npts-1] - 4.0*u[npts-2] + u[npts-3])*two_invdx;
  
  return;
}


/* Second Order First Derivative in 1D, periodic boundaries */
void SO_FirstDeriv_1Dper (int npts, double dx, double *u, double *u_x){
  double two_invdx = 1.0/(2.0*dx);
  
  for(int i=1;i<npts-1;i++)
    u_x[i] = (u[i+1]-u[i-1])*two_invdx;
    
  // Left Endpoint
  u_x[0] = (u[1]-u[npts-1])*two_invdx;
  // Right Endpoint
  u_x[npts-1] = (u[0]-u[npts-2])*two_invdx;
  
  return;
}


#ifdef MPISRC
/* Second Order First Derivative in 1D -- Parallel*/
void SO_FirstDeriv_1DP    (int npts, double dx, double *u, double *u_x, 
			   int mynode, int totalnodes){
  double two_invdx = 1.0/(2.0*dx);
  double mpitemp;
  MPI_Status status;

  if(mynode == 0)
    u_x[0] = (-3.0*u[0] + 4.0*u[1] - u[2])*two_invdx;
  
  if(mynode == (totalnodes-1))
    u_x[npts-1] = (3.0*u[npts-1] - 4.0*u[npts-2] + u[npts-3])*two_invdx;
    
  for(int i=1;i<npts-1;i++)
    u_x[i] = (u[i+1]-u[i-1])*two_invdx;
  
  if(mynode == 0){
#ifdef MPI_SENDRECV
    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD);

    MPI_Recv(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD, &status); 
    u_x[npts-1] = (mpitemp - u[npts-2])*two_invdx;  
#else
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,1,1,1,1, MPI_COMM_WORLD, &status);
    u_x[npts-1] = (mpitemp - u[npts-2])*two_invdx;  
    
#endif
  }
  else if(mynode == (totalnodes-1)){
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_x[0] = (u[1]-mpitemp)*two_invdx;

    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,
			 mynode-1,1, MPI_COMM_WORLD, &status);
    u_x[0] = (u[1]-mpitemp)*two_invdx;
#endif
  }
  else{
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_x[0] = (u[1]-mpitemp)*two_invdx;
    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);

    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD);
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD, &status);     
    u_x[npts-1] = (mpitemp-u[npts-2])*two_invdx;
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,
			 mynode-1,1, MPI_COMM_WORLD, &status);
    u_x[0] = (u[1]-mpitemp)*two_invdx;
    
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode+1,1,
			 mynode+1,1, MPI_COMM_WORLD, &status);
    u_x[npts-1] = (mpitemp-u[npts-2])*two_invdx;
#endif
  }
  
  return;
}

/* Second Order First Derivative in 1D, periodic boundaries -- Parallel */
void SO_FirstDeriv_1DperP (int npts, double dx, double *u, double *u_x, 
			   int mynode, int totalnodes){
  double two_invdx = 1.0/(2.0*dx);
  double mpitemp;
  MPI_Status status;
  
  for(int i=1;i<npts-1;i++)
    u_x[i] = (u[i+1]-u[i-1])*two_invdx;
  
  if(mynode == 0){
#ifdef MPI_SENDRECV
    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD);
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD, &status); 
    u_x[npts-1] = (mpitemp - u[npts-2])*two_invdx;  

    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,totalnodes-1,1,MPI_COMM_WORLD);
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,totalnodes-1,1,MPI_COMM_WORLD, &status);     
    u_x[0] = (u[1]-mpitemp)*two_invdx;
#else
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,1,1,1,1, MPI_COMM_WORLD, &status);
    u_x[npts-1] = (mpitemp - u[npts-2])*two_invdx;  

    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,totalnodes-1,1,
			 totalnodes-1,1, MPI_COMM_WORLD, &status);
    u_x[0] = (u[1]-mpitemp)*two_invdx;  
#endif
  }
  else if(mynode == (totalnodes-1)){
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_x[0] = (u[1]-mpitemp)*two_invdx;

    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);

    MPI_Recv(&mpitemp,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD, &status); 
    u_x[npts-1] = (mpitemp-u[npts-2])*two_invdx;
    
    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,
			 mynode-1,1, MPI_COMM_WORLD, &status);
    u_x[0] = (u[1]-mpitemp)*two_invdx;

    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,0,1,0,1, MPI_COMM_WORLD, &status);
    u_x[npts-1] = (mpitemp - u[npts-2])*two_invdx;
#endif
  }
  else{
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_x[0] = (u[1]-mpitemp)*two_invdx;
    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);

    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD);
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD, &status);     
    u_x[npts-1] = (mpitemp-u[npts-2])*two_invdx;
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,
			 mynode-1,1, MPI_COMM_WORLD, &status);
    u_x[0] = (u[1]-mpitemp)*two_invdx;
    
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode+1,1,
			 mynode+1,1, MPI_COMM_WORLD, &status);
    u_x[npts-1] = (mpitemp-u[npts-2])*two_invdx;
#endif
  }
  
  return;
}
#endif



/* Second Order Second Derivative in 1D */
void SO_SecondDeriv_1D    (int npts, double dx, double *u, double *u_xx){
  int i;
  double inv_dx2 = 1.0/(dx*dx);
  
  // Forward differencing
  u_xx[0] = (2.0*u[0]-5.0*u[1]+4.0*u[2]-u[3])*inv_dx2;
  
  // Central differencing
  for(i=1;i<npts-1;i++)
    u_xx[i] = (u[i+1]-2.0*u[i]+u[i-1])*inv_dx2;
  
  // Backward differencing
  u_xx[npts-1] = (2.0*u[npts-1]-5.0*u[npts-2]+4.0*u[npts-3]-u[npts-4])*inv_dx2;
  
  return;
}



/* Second Order Second Derivative in 1D, periodic boundaries */
void SO_SecondDeriv_1Dper (int npts, double dx, double *u, double *u_xx){
  int i;
  double inv_dx2 = 1.0/(dx*dx);

  u_xx[0] = (u[1]-2.0*u[0]+u[npts-1])*inv_dx2;
  
  // Central differencing
  for(i=1;i<npts-1;i++)
    u_xx[i] = (u[i+1]-2.0*u[i]+u[i-1])*inv_dx2;
  
  u_xx[npts-1] = (u[0]-2.0*u[npts-1]+u[npts-2])*inv_dx2;
  
  return;
}


#ifdef MPISRC
/* Second Order Second Derivative in 1D -- Parallel*/
void SO_SecondDeriv_1DP    (int npts, double dx, double *u, double *u_xx, 
			    int mynode, int totalnodes){
  int i;
  double inv_dx2 = 1.0/(dx*dx);
  double mpitemp;
  MPI_Status status;

  if(mynode == 0)
    u_xx[0] = (2.0*u[0]-5.0*u[1]+4.0*u[2]-u[3])*inv_dx2;
    
  if(mynode == (totalnodes-1))
    u_xx[npts-1] = (2.0*u[npts-1]-5.0*u[npts-2]+4.0*u[npts-3]-u[npts-4])*inv_dx2;
  
  for(i=1;i<npts-1;i++)
    u_xx[i] = (u[i+1]-2.0*u[i]+u[i-1])*inv_dx2;
  
  if(mynode == 0){
#ifdef MPI_SENDRECV
    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD);

    MPI_Recv(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD, &status); 
    u_xx[npts-1] = (mpitemp - 2.0*u[npts-1] + u[npts-2])*inv_dx2;  
#else
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,1,1,1,1, MPI_COMM_WORLD, &status);
    u_xx[npts-1] = (mpitemp - 2.0*u[npts-1] + u[npts-2])*inv_dx2;  
#endif
  }
  else if(mynode == (totalnodes-1)){
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_xx[0] = (u[1] - 2.0*u[0] + mpitemp)*inv_dx2;

    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,mynode-1,1, 
			 MPI_COMM_WORLD, &status);
    u_xx[0] = (u[1] - 2.0*u[0] + mpitemp)*inv_dx2;
#endif
  }
  else{
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_xx[0] = (u[1] -2.0*u[0] + mpitemp)*inv_dx2;
    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);

    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD);
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD, &status);     
    u_xx[npts-1] = (mpitemp -2.0*u[npts-1] + u[npts-2])*inv_dx2;
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,
			 mynode-1,1, MPI_COMM_WORLD, &status);
    u_xx[0] = (u[1] -2.0*u[0] + mpitemp)*inv_dx2;
    
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode+1,1,
			 mynode+1,1, MPI_COMM_WORLD, &status);
    u_xx[npts-1] = (mpitemp -2.0*u[npts-1] + u[npts-2])*inv_dx2;
#endif
  }
  
  return;
}

/* Second Order Second Derivative in 1D, periodic boundaries -- Parallel */
void SO_SecondDeriv_1DperP (int npts, double dx, double *u, double *u_xx, 
			    int mynode, int totalnodes){
  int i;
  double inv_dx2 = 1.0/(dx*dx);
  double mpitemp;
  MPI_Status status;
  
  for(i=1;i<npts-1;i++)
    u_xx[i] = (u[i+1]-2.0*u[i]+u[i-1])*inv_dx2;
  
  if(mynode == 0){
#ifdef MPI_SENDRECV
    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD);

    MPI_Recv(&mpitemp,1,MPI_DOUBLE,1,1,MPI_COMM_WORLD, &status); 
    u_xx[npts-1] = (mpitemp - 2.0*u[npts-1] + u[npts-2])*inv_dx2;  

    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,totalnodes-1,1,MPI_COMM_WORLD);
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,totalnodes-1,1,MPI_COMM_WORLD, &status);     
    u_xx[0] = (u[1]-2.0*u[0]+mpitemp)*inv_dx2;

#else
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,1,1,1,1, MPI_COMM_WORLD, &status);
    u_xx[npts-1] = (mpitemp - 2.0*u[npts-1] + u[npts-2])*inv_dx2;  

    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,totalnodes-1,1,
			 totalnodes-1,1, MPI_COMM_WORLD, &status);
    u_xx[0] = (u[1] - 2.0*u[0] + mpitemp)*inv_dx2;  
#endif
  }
  else if(mynode == (totalnodes-1)){
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_xx[0] = (u[1] - 2.0*u[0] + mpitemp)*inv_dx2;

    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);

    MPI_Recv(&mpitemp,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD, &status); 
    u_xx[npts-1] = (mpitemp-2.0*u[npts-1]+u[npts-2])*inv_dx2;
    
    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,mynode-1,1, 
			 MPI_COMM_WORLD, &status);
    u_xx[0] = (u[1] - 2.0*u[0] + mpitemp)*inv_dx2;

    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,0,1,0,1, 
			 MPI_COMM_WORLD, &status);
    u_xx[npts-1] = (mpitemp -  2.0*u[npts-1] + u[npts-2])*inv_dx2;
#endif
  }
  else{
#ifdef MPI_SENDRECV
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD, &status); 
    u_xx[0] = (u[1] -2.0*u[0] + mpitemp)*inv_dx2;
    mpitemp = u[0];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode-1,1,MPI_COMM_WORLD);

    mpitemp = u[npts-1];
    MPI_Send(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD);
    MPI_Recv(&mpitemp,1,MPI_DOUBLE,mynode+1,1,MPI_COMM_WORLD, &status);     
    u_xx[npts-1] = (mpitemp -2.0*u[npts-1] + u[npts-2])*inv_dx2;
#else
    mpitemp = u[0];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode-1,1,
			 mynode-1,1, MPI_COMM_WORLD, &status);
    u_xx[0] = (u[1] -2.0*u[0] + mpitemp)*inv_dx2;
    
    mpitemp = u[npts-1];
    MPI_Sendrecv_replace(&mpitemp,1,MPI_DOUBLE,mynode+1,1,
			 mynode+1,1, MPI_COMM_WORLD, &status);
    u_xx[npts-1] = (mpitemp -2.0*u[npts-1] + u[npts-2])*inv_dx2;
#endif
  }
  
  return;
}
#endif



/**************************************************************************************/
void CD_SecondDeriv(int npts, double dx, double dy, double **u, double **u_xx_yy){
  double inv_dx2 = 1.0/(dx*dx);
  double inv_dy2 = 1.0/(dy*dy);
  
  for(int i=1;i<npts-1;i++)
    for(int j=1;j<npts-1;j++){
      u_xx_yy[i][j] = (u[i-1][j]-2.0*u[i][j] +u[i+1][j])*inv_dx2 +
	(u[i][j-1]-2.0*u[i][j] +u[i][j+1])*inv_dy2;
    }
  
  return;
}

void CrossDerivative(int npts, double dx, double dy, double **u, double **u_xy){
  double inv_dxdy = 1.0/(4.0*dx*dy);

  for(int i=1;i<npts-1;i++)
    for(int j=1;j<npts-1;j++){
      u_xy[i][j] = inv_dxdy*(u[i+1][j+1]-u[i+1][j-1]-u[i-1][j+1]+u[i-1][j-1]);
    }

  return;
}

void FornbergWeights(double xi, double *x, int m, int n, double ***C){
  int i,j,k,mn;
  double C1,C2,C3;
  
  C[0][0][0] = 1.0;
  C1 = 1.0;
  
  for(j=1;j<=n;j++){
    if(j<m)
      mn = j;
    else
      mn = m;
    C2 = 1.0;
    
    for(k=0;k<=(j-1);k++){
      C3 = x[j]-x[k];
      C2 = C2*C3;
      if(j<=m) C[j][j-1][k]=0.;
      C[0][j][k] = (x[j]-xi)*C[0][j-1][k]/C3;
      for(i=1;i<=mn;i++)
	C[i][j][k] = ((x[j]-xi)*C[i][j-1][k]-i*C[i-1][j-1][k])/C3;
    }

    C[0][j][j] = -C1*(x[j-1]-xi)*C[0][j-1][j-1]/C2;
    
    for(i=1;i<=mn;i++)
      C[i][j][j] = C1*(i*C[i-1][j-1][j-1]-(x[j-1]-xi)*C[i][j-1][j-1])/C2;
    C1 = C2;
  }
  return;
}


double AdamsBashforth(int order, double u_old, double dt, double * RHS){
  double answer;
  
  switch(order){
  case 1: /* 1st Order Adams-Bashforth -- Euler Forward */
    answer = u_old + dt*RHS[0];
    break;
  case 2: /* 2nd Order Adams-Bashforth */
    answer = u_old + dt*(1.5*RHS[0] - 0.5*RHS[1]);
    break;
  case 3: /* 3rd Order Adams-Bashforth */
    answer = u_old + dt*( (23./12.)*RHS[0] - (4./3.)*RHS[1] + (5./12.)*RHS[2]);
    break;
  default: /* default is Euler Forward */
    answer = u_old + dt*RHS[0];
  }
  
  return answer;
  
}


void AdamsBashforth(int order, int N, double *u_old, double *u_new, double dt, double ** RHS){
  int i;
  
  switch(order){
  case 1: /* 1st Order Adams-Bashforth -- Euler Forward */
    for(i=0;i<N;i++)
      u_new[i] = u_old[i] + dt*RHS[i][0];
    break;
  case 2: /* 2nd Order Adams-Bashforth */
    for(i=0;i<N;i++)
      u_new[i] = u_old[i] + dt*(1.5*RHS[i][0] - 0.5*RHS[i][1]);
    break;
  case 3: /* 3rd Order Adams-Bashforth */
    for(i=0;i<N;i++)
      u_new[i] = u_old[i] + dt*( (23./12.)*RHS[i][0] - (4./3.)*RHS[i][1] + (5./12.)*RHS[i][2]);
    break;
  default: /* default is Euler Forward */
    for(i=0;i<N;i++)
      u_new[i] = u_old[i] + dt*RHS[i][0];
  }
  
  return;
}


double RungeKutta4(double uold, double time, double dt, double (*rkfunc)(double,double)){
  double unew, hold[4];

  hold[0] = rkfunc(uold,time);
  hold[1] = rkfunc(uold+0.5*hold[0],time+0.5*dt);
  hold[2] = rkfunc(uold+0.5*hold[1],time+0.5*dt);
  hold[3] = rkfunc(uold+hold[2],time+dt);

  unew = uold + (1.0/6.0)*(hold[0]+2.0*(hold[1]+hold[2])*hold[3]);

  return unew;
}


double RungeKutta(int order, double dt, double uold, double (*rkfunc)(double)){
  double unew = uold;

  for(int k=order;k>=1;k--)
    unew = uold + (dt/k)*rkfunc(unew);

  return unew;
}







