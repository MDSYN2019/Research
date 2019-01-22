#ifndef _VANILLA_OPTION_CPP_
#define _VANILLA_OPTION_CPP_

#include "VanillaOption.hpp"
#include <cmath>

VanillaOption::VanillaOption() { init(); }
VanillaOption::VanillaOption(const double& _K, const double& _r, const double& _T, const doule& _S, const double& _sigma);

void VanillaOption::init() {
  K = 100.0;
  r = 0.05; // 5% interest rate 
  T = 1.0; // One year until maturity 
  S = 100.0; // Option is 'at the money' as spot equals the strike 
  sigma = 0.2; // Volatility is 20% 
}

void VanillaOption::copy(const VanillaOption& rhs) {
}


#endif 
