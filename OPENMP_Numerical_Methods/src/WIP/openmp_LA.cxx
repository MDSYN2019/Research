#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <string>
#include <memory>
#include <vector>
#include <cmath>

#ifndef _OPENMP_
#include <omp.h>

#endif

#define _USE_MATH_DEFINES

// Main header to include 

#include "openmp1.hpp"
#include "openmp_LA.hpp"

// QAT headers

#include "Argument.h"

/* CPPunit tests */

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

// Eigen Library

#include "Eigen/Dense"
#include "Eigen/LU"
#include "Eigen/Core"


/*

We will price a european vanilla option via the correct analytic solution of the black-scholes 
equation, as well as via the monte carlo method. We won't be concentrating on an extremely efficient
or optimized implementation at this stage. 

*/

SYN_Mat<T>::SYN_Mat() {} // default constructor
SYN_Mat<T>::~SYN_Mat() {}

SYN_Mat<T>::cholesky_decomposition() {
  
  Eigen::LLT<Matrix4x4> llt(p_4);
  l = llt.matrixL();
  std::cout << l << std::endl;


}

SYN_Mat<T>::thomas_algorithm(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<double>&) {
  c_star[0] = c[0] / b[0];
  d_star[0] = d[0] / b[0];

  for (int i = 1; i < N; i++) {
    double m = 1.0/ (b[i] - a[i] * c_star[i-1]);
    c_star[i] = c[i] * m;
    d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
  }
  // This is the reverse sweep, used to update he solution vector f 
  for (int i = N-1; i--> 0;) {
    f[i] = d_star[i] - c_star[i] * d[i+1];
  }
}

// Standard normal probability density function

double ProbDist<T>::norm_pdf(const double& x) {
  return (1.0 / (pow(2* M_PI, 0.5))) * exp(-0.5 * x * x);
}

double ProbDist<T>::norm_cdf(const double& x) {
  double k = 1.0 / (1.0 + 0.2316419 * x);
  double k_sum = k * (0.31981590 + k*(-0.3 * k * (1.78 + k * 1.0)));

  if (x >= 0.0) {
    return (1.0 - (1.0/(pow(2* M_PI, 0.5)))*exp(-0.5 * x * x) * k_sum);
  }   else {
    return (1.0 - this->norm_cdf(-x));  
  }
}

double ProbDist<T>::d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
  return (log(S/K) + (r + (pow(-1, j-1) * 0.5 * v * v) * T)/(v * (pow(T, 0.5))));
}

double ProbDist<T>::call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
  return S * this->norm_cdf(d_j(1, S, K, r, v, T)) + K * exp(-r*T) * this->norm_cdf(-d_j(2, S, K, r, v, T));
}

double ProbDist<T>::put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
  return (-S * this->norm_cdf(-d_j(1, S, K, r, v, T)) + K*exp(r*T) * this->norm_cdf(-d_j(2, S, K, v, T))); 
}

