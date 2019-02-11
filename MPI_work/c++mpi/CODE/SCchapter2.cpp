/**************************************************************************
*
* Name:         Scientific Computing Mathematics Library (SCmathlib)
*
* File name:    SCchapter2.cpp  (definition file)
* Header file:  SCchapter2.h    (header file)      
*
* Developers:   R.M. Kirby and G.E. Karniadakis
*
*************************************************************************/


#include "SCchapter2.h"


float FloatMachineEps(){
  float  fmachine_e, ftest;
  fmachine_e = 1.0;
  
  ftest = 1.0 + fmachine_e;
  while(1.0 != ftest){
    fmachine_e = fmachine_e/2.0;
    ftest = 1.0 + fmachine_e;
  }

  return fmachine_e;
}


double DoubleMachineEps(){
  double dmachine_e, dtest;
  dmachine_e = 1.0;
  
  dtest = 1.0 + dmachine_e;
  while(1.0 != dtest){
    dmachine_e = dmachine_e/2.0;
    dtest = 1.0 + dmachine_e;
  }

  return dmachine_e;
}



SCstatus GramSchmidt(SCVector * x, SCVector * q){
  int i,j;
  int dim = x[0].Dimension();
  SCVector y(dim);
  SCMatrix r(dim);

  r(0,0) = x[0].Norm_l2();

  if(r(0,0)==0.0)
    return(FAIL);
  else
    q[0] = x[0]/r(0,0);


  for(j=1;j<dim;j++){
    for(i=0;i<=j-1;i++)
      r(i,j) = dot(q[i],x[j]);
    
    y = x[j];
    for(i=0;i<=j-1;i++)
      y = y - r(i,j)*q[i];

    r(j,j) = y.Norm_l2();
    
    if(r(j,j) == 0.0)
      return(FAIL);
    else
      q[j] = y/r(j,j);
  }

  return(SUCCESS);
}

SCstatus GramSchmidt(SCVector * x, SCVector * q, SCMatrix &r){
  int i,j;
  int dim = x[0].Dimension();
  SCVector y(dim);

  r(0,0) = x[0].Norm_l2();

  if(r(0,0)==0.0)
    return(FAIL);
  else
    q[0] = x[0]/r(0,0);


  for(j=1;j<dim;j++){
    for(i=0;i<=j-1;i++)
      r(i,j) = dot(q[i],x[j]);
    
    y = x[j];
    for(i=0;i<=j-1;i++)
      y = y - r(i,j)*q[i];

    r(j,j) = y.Norm_l2();
    
    if(r(j,j) == 0.0)
      return(FAIL);
    else
      q[j] = y/r(j,j);
  }

  return(SUCCESS);
}


SCstatus ModifiedGramSchmidt(SCVector * x, SCVector * q, SCMatrix &r){
  int i,j;
  int dim = x[0].Dimension();
  SCVector y(dim);

  r(0,0) = x[0].Norm_l2();

  if(r(0,0)==0.0)
    return(FAIL);
  else
    q[0] = x[0]/r(0,0);


  for(j=1;j<dim;j++){    
    y = x[j];
    for(i=0;i<=j-1;i++){
      r(i,j) = dot(q[i],y);
      y = y - r(i,j)*q[i];
    }

    r(j,j) = y.Norm_l2();
    
    if(r(j,j) == 0.0)
      return(FAIL);
    else
      q[j] = y/r(j,j);
  }

  return(SUCCESS);
}


SCstatus QRDecomposition(SCMatrix A, SCMatrix &Q, SCMatrix &R){
  int i,j;
  int num_vecs = A.Rows();
  int dim = A.Columns();
  SCstatus scflag;

  SCVector * q = new SCVector[num_vecs], 
    * v = new SCVector[num_vecs];
  
  for(i=0;i<num_vecs;i++){
    q[i].Initialize(dim);
    v[i].Initialize(dim);
  }
  

  for(i=0;i<num_vecs;i++){
    for(j=0;j<dim;j++)
      v[i](j) = A(j,i);
  }

  scflag = ModifiedGramSchmidt(v,q,R);

  for(i=0;i<num_vecs;i++)
    for(j=0;j<dim;j++)
      Q(j,i) = q[i](j); 

  return scflag;
}







