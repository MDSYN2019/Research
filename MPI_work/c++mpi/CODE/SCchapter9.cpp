/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter9.cpp  (definition file)
* Header file:  SCchapter9.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#include "SCchapter9.h"


void GaussElimination(SCMatrix &A, SCVector &b, int pivotflag){
  int k,pivot;
  int N = A.Rows();

  /* NOTE: The values contained in both the matrix A and
     the vector b are modified in this routine.  Upon 
     returning, A contains the upper triangular matrix
     obtained from its LU decomposition, and b contains
     the solution of the system Ax=b*/

  // Steps (1) and (2) (decomposition and solution of Ly = b)
  switch(pivotflag){
  case 1: // Case in which pivoting is employed

    for(k=0; k < N-1; k++){
      pivot = A.MaxModInColumnindex(k,k);
      A.RowSwap(pivot,k);
      Swap(b(pivot),b(k)); 
      for(int i=k+1;i<N;i++){
	double l_ik = A(i,k)/A(k,k);
	for(int j=k;j<N;j++)
	  A(i,j) = A(i,j) - l_ik*A(k,j);
	b(i) = b(i) - l_ik*b(k); 
      }
    }
    break;

  case 0:  // Case 0/default in which no pivoting is used
  default:

    for(k=0; k < N-1; k++){
      for(int i=k+1;i<N;i++){
	double l_ik = A(i,k)/A(k,k);
	for(int j=k;j<N;j++)
	  A(i,j) = A(i,j) - l_ik*A(k,j);
	b(i) = b(i) - l_ik*b(k); 
      }
    }
  }

  // Step (3) (backsolveing to solve Ux=y)  
  b(N-1) = b(N-1)/A(N-1,N-1);
  for(k=N-2;k>=0;k--){
    for(int j=k+1;j<N;j++)
      b(k) -= A(k,j)*b(j);
    b(k) = b(k)/A(k,k);
  }
}


double HouseholderTrans(SCVector &x, SCVector &w){
  double alpha;
  double xm = x.MaxMod();
  
  for(int i=0;i<x.Dimension();i++)
    w(i) = x(i)/xm;

  alpha = Sign(w(0))*w.Norm_l2();

  w(0) = w(0) + alpha;

  alpha = -alpha*xm;

  return alpha;
}


void HouseholderQR(SCMatrix &A, SCVector &v){
  int i,j,k,q;
  double beta,gamma;
  int N = A.Rows();
  SCVector *x = new SCVector(N),
    *w = new SCVector(N);

  for(k=0;k<N-1;k++){
    A.GetColumn(k,*x,k);
    A(k,k) = HouseholderTrans(*x,*w);
    beta = 2.0/dot(N-k,*w,*w);
    v(k) = (*w)(0);
    for(i=k+1;i<N;i++)
      A(i,k) = (*w)(i-k);
    for(j=k+1;j<N;j++){
      gamma = 0.0;
      for(q=k;q<N;q++)
	gamma += beta*(*w)(q-k)*A(q,j);
      for(i=k;i<N;i++){
	A(i,j) = A(i,j) - gamma*(*w)(i-k);
      }
    }
  }
  
  delete x;
  delete w;
}


void Hessenberg(SCMatrix &A){
  int i,j,k,q;
  double beta,gamma;
  int N = A.Rows();
  SCVector *x = new SCVector(N),
    *w = new SCVector(N);

  for(k=0;k<N-2;k++){
    A.GetColumn(k,*x,k+1);
    A(k+1,k) = HouseholderTrans(*x,*w);
    beta = 2.0/dot(N-k-1,*w,*w);
   
    for(i=k+2;i<N;i++)
      A(i,k) = (*w)(i-k);
    for(j=k+1;j<N;j++){
      gamma = 0.0;
      for(q=k+1;q<N;q++)
	gamma += beta*(*w)(q-k-1)*A(q,j);
      for(i=k+1;i<N;i++){
	A(i,j) = A(i,j) - gamma*(*w)(i-k-1);
      }
    }
    for(i=0;i<N;i++){
      gamma = 0.0;
      for(q=k+1;q<N;q++)
	gamma += beta*(*w)(q-k-1)*A(i,q); 
      for(j=k+1;j<N;j++){
	A(i,j) = A(i,j) - gamma*(*w)(j-k-1);
      }
    }
  }
  
  delete x;
  delete w;
}


void ModifiedArnoldi(int m, const SCMatrix &A, SCMatrix &H, 
		     SCMatrix &V){
  SCVector v(A.Rows()),w(A.Rows());
  v.Initialize(0.0);
  v(0) = 1.0;

  V.PutColumn(0,v);

  for(int j=0;j<m;j++){
    w = A*v;
    for(int i=0;i<=j;i++){
      V.GetColumn(i,v);
      H(i,j) = dot(w,v);
      w = w - H(i,j)*v;
    }
    H(j+1,j) = w.Norm_l2();
    v = w/H(j+1,j);
    V.PutColumn(j+1,v);
  }
}


void ModifiedArnoldi(int m, const SCVector &x, const SCMatrix &A, SCMatrix &H, 
		     SCMatrix &V){
  SCVector v(A.Rows()),w(A.Rows());
  v = x;

  V.PutColumn(0,v);
  
  for(int j=0;j<m;j++){
    w = A*v;
    for(int i=0;i<=j;i++){
      V.GetColumn(i,v);
      H(i,j) = dot(w,v);
      w = w - H(i,j)*v;
    }

    H(j+1,j) = w.Norm_l2();
    v = w/H(j+1,j);
    V.PutColumn(j+1,v);
  }

}


void GMRES(int m, const SCMatrix &A, const SCVector &b, SCVector &x){
  int i,j,k,ll,nr;
  int N = A.Rows();
  SCMatrix H(m+1,m),V(N,m+1);
  SCVector w(N),r(N),y(m+1),z(N);
  double * c = new double[m+1];
  double * s = new double[m+1];
  const int maxit = 1000;
  const double tol = 1.0e-7;
  double delta,rho,tmp;

  x.Initialize(0.0);

  r = b - A*x;

  for(j=0;j<maxit;j++){
    y.Initialize(0.0);
    y(0) = r.Norm_l2();
    r.Normalize();

    ModifiedArnoldi(m,r,A,H,V);

    /* Givens Rotation to accomplish QR factorization */
    for(i=0;i<m;i++){
      for(k=1;k<=i;k++){
	tmp = H(k-1,i);
	H(k-1,i) = c[k-1]*H(k-1,i) + s[k-1]*H(k,i);
	H(k,i) = -s[k-1]*tmp + c[k-1]*H(k,i);
      }

      delta = sqrt(H(i,i)*H(i,i)+H(i+1,i)*H(i+1,i));     
      c[i] = H(i,i)/delta;
      s[i] = H(i+1,i)/delta;

      H(i,i) = c[i]*H(i,i) + s[i]*H(i+1,i);
      
      for(k=i+1;k<m+1;k++)
	H(k,i) = 0.0;

      y(i+1) = -s[i]*y(i);
      y(i)   =  c[i]*y(i);
      rho = fabs(y(i+1));
      if(rho < tol){
	nr = i;
	break;
      }
    }

    /* Backsolve to obtain coefficients */
    z.Initialize(0.0);
    if(i>=(m-1)){ 
      nr = m; 
      z(nr-1) = y(nr-1)/H(nr-1,nr-1);
    }

    for(k=nr-2;k>=0;k--){
      z(k) = y(k);
      for(ll=k+1;ll<nr;ll++)
	z(k) -= H(k,ll)*z(ll);
      z(k) = z(k)/H(k,k);
    }
    
    /* Linear combination of basis vectors of the Krylov space */
    for(i=0;i<nr;i++){
      V.GetColumn(i,r);
      x = x + z(i)*r;
    }

    if(rho<tol)
      break;

    r = b - A*x;
  }
  delete[] c;
  delete[] s;
}



 
