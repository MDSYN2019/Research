#ifndef __GSLMAT__
#define __GSLMAT__

#include <iostream>
#include <vector>
#include "gsl.hpp"

class initiateMatrix {
public:
  initiateMatrix(); // constructor
  initiateMatrix(gsl_matrix * m, int i, int j); // constructor
  gsl_matrix * M;
  int ii;
  int jj;
  void set_matrix();
  ~initiateMatrix(); // destructor 
};

void cross_product(const gsl_vector *u, const gsl_vector *v, gsl_vector *product);
double triple_scalar_product(const gsl_vector *u, const gsl_vector *v, const gsl_vector *w);

//class cross_triple {
//
//public:
//  cross_triple(const gsl_vector *a, const gsl_vector *b, const gsl_vector *c, int N);
//  ~cross_triple();
//  void cross_product(const gsl_vector *u, const gsl_vector *v, const gsl_vector *product);
//  void triple_scalar_product(const gsl_vector *u, const gsl_vector *v, const gsl_vector *w); 
//};

#endif
