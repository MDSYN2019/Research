#ifndef __MATRIX_CPP__
#define __MATRIX_CPP__

#include <iostream>
#include <algorithm> // The algorithm libary provides access the max comparison function
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <string>
#include <memory>
#include <vector>
#include <cmath>

/* CPPunit tests */

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


// Eigen Library

#include "Eigen/Dense"
#include "Eigen/LU"

//! Maths definitions
/*! 
  This pre-processor which tells the C++ compiler to make use of the 
  C-standard mathematical constants. Now that some compilers do not fully 
  support constants.
 */

#define _USE_MATH_DEFINES


/*!

European Options with Monte Carlo 
---------------------------------

In this chapter, we will pric a European Vanilla option via the correct analyci solution of the 
Black-Scholes eqiation, as well as via the Monte Carlo method. We won't be cooncentrating on 
an extremely efficient or optimised implementation at this stage. 

-> Black Scholes Analytic pricing formula 

The first stage in implementation is to briefly discuess the Black Scholes analytic solution for the 
price of a vanilla call or put option

*/


class ProbDist {
public:
  // Constructor/Destructor section
  ProbDist();
  //  ProbDist(const ProbDist& alloc); // Copy constructor 
  // Operator overloading
  // ProbDist& operator=(const ProbDist& alloc); 
  ~ProbDist(); // virtual destructor
  
  // Functions

  double norm_pdf(const double& x);
  double norm_cdf(const double& x);
  double d_j(const int&, const double&, const double&, const double&, const double&, const double&);
  double call_price(const double&, const double&, const double&, const double&, const double&);
  double put_price(const double&, const double&, const double&, const double&, const double&); 

  double gaussian_box_muller();
  double monte_carlo_call_price(const int&, const double&, const double&, const double&, const double&, const double&);
};

ProbDist::ProbDist() {
}

ProbDist::~ProbDist() {
}

double ProbDist::norm_pdf(const double& x) {
  return (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
} // This is fine 

double ProbDist::norm_cdf(const double& x) {
  double k = 1.0 / (1.0 + 0.2316419 * x);
  double k_sum = k * (0.319 + k * (-0.356563 + k * (1.781477937 + k * (-1.821255978) + 1.330274429 * k)));
  if (x >= 0.0) {
    return (1.0 - (1.0 / (pow(2 * M_PI, 0.5)) * exp(0.5 * x * x) * k_sum));
  } else {
    return 1.0 - this->norm_cdf(-x);
  }
}

/*
This calculates d_j, for j in {1,2}. This term appears in the closed 
 form solution for the European call or put price 
 */

double ProbDist::d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
  return (log(S/K) + (r + (pow(1, j-1)) * 0.5 * v * v) * T) / (v *( pow(T, 0.5)));  
}

double ProbDist::call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
  return S * this->norm_cdf(this->d_j(1, S, K, r, v, T)) - K * exp(-r * T)* this->norm_cdf(d_j(2, S, K, r,v, T));
}

double ProbDist::put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
  return -S * this->norm_cdf(-(this->d_j(1,S,K,r,v,T))) + K * exp(-r * T) * this->norm_cdf(-(this->d_j(2, S, K, r, v, T)));
}

double ProbDist::gaussian_box_muller() {

  /*
    The box-muller algorithm - designed to convert two uniform random vairables
    into a standard Gaussian random variable.

    Box-Muller is a good choice for random number generator if your compiler
    does not support the C++11 standard.

   */
  
  double x = 0.0;
  double y = 0.0;
  double euclid_sq = 0.0;

  /*

    Continue genereating two uniform random variables
    until the square of their 'euclidean distance' 
    is less than unity

  */

  // What does static cast do?

  /*
    The static_cast operator converts variable j to type float. This allows
    the compiler to generate a division with an answer of type float. All 
    static_cast resolve at compile time nd do not remove any const or 
    volatile modifiers
   */
  do {
    x = 2.0 * rand() / static_cast<double>(RAND_MAX) -1;
    y = 2.0 * rand() / static_cast<double>(RAND_MAX) -1;
    euclid_sq = x*x + y*y;
  } while (euclid_sq >= 1.0);

  return x*sqrt(-2 * log(euclid_sq)/euclid_sq);
}

// pricing a European vanilla call option with a MC method

double ProbDist::monte_carlo_call_price(const int& num_sims, const double& S, const double& K, const double& r, const double& v, const double& T) {

  double S_adjust = S * exp(T * (r - 0.5 * v * v));
  double S_cur = 0.0;
  double payoff_sum = 0.0;
  
  for (int i = 0; i < num_sims; i++) {
    double gauss_bm = this->gaussian_box_muller(); // Refer to the class within the same method
    S_cur = S_adjust * exp(sqrt(v*v*T) * gauss_bm);
    payoff_sum += std::max(S_cur - K, 0.0);
    
  }

  return (payoff_sum / static_cast<double>(num_sims)) * exp(-r * T);
}

#endif 
