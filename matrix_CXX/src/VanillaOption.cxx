#ifndef _VANILLA_OPTION_CPP_
#define _VANILLA_OPTION_CPP_

#include "VanillaOption.hpp"
#include <cmath>

VanillaOption::VanillaOption() { init(); }

// Copy constructor

VanillaOption::VanillaOption(const VanillaOption& rhs) {
  copy(rhs);
}

VanillaOption::VanillaOption::operator=(const VanillaOption& rhs) {
  if (this == &rhs) return *this;
  copy (rhs);
  return *this;
}

VanillaOption::VanillaOption(const double& _K, const double& _r, const double& _T, const doule& _S, const double& _sigma);

void VanillaOption::init() {
  K = 100.0;
  r = 0.05; // 5% interest rate 
  T = 1.0; // One year until maturity 
  S = 100.0; // Option is 'at the money' as spot equals the strike 
  sigma = 0.2; // Volatility is 20% 
}

void VanillaOption::copy(const VanillaOption& rhs) {
  K = rhs.getK();
  r = rhs.getr();
  T = rhs.getT();
  S = rhs.getS();
  sigma = _sigma;
}

double VanillaOption::calc_call_price() {
  
  double sigma_sqrt_T = sigma * sqrt(T);
  double d_1 = (log (S/K) + (r + sigma * sigma * 0.5) * T)/ sigma_sqrt_T;
  double d_2 = d_1 - sigma_sqrt_T;

  return S * N(d_1) - K * exp(-r*T) * N(d_2);
  
}



double VanillaOption::getK() const {return K;}
double VanillaOption::getr() const {return r;}
double VanillaOption::getT() const {return T;}
double VanillaOption::getS() const {return S;}
double VanillaOption::getsigma() const {return sigma;}




}

#endif 
