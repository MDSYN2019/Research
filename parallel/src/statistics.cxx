#ifndef __STAT_CXX__
#define __STAT_CXX__

#include "statistics.h"
#include <iostream>
#include <cmath>

#define _USE_MATH_DEFINES

/* Define _USE_MATH_DEFINES before including math.h to expose these macro
 * definitions for common math constants.  These are placed under an #ifdef
 * since these commonly-defined names are not part of the C/C++ standards.
 */

/* Definitions of useful mathematical constants
 * M_E        - e
 * M_LOG2E    - log2(e)
 * M_LOG10E   - log10(e)
 * M_LN2      - ln(2)
 * M_LN10     - ln(10)
 * M_PI       - pi
 * M_PI_2     - pi/2
 * M_PI_4     - pi/4
 * M_1_PI     - 1/pi
 * M_2_PI     - 2/pi
 * M_2_SQRTPI - 2/sqrt(pi)
 * M_SQRT2    - sqrt(2)
 * M_SQRT1_2  - 1/sqrt(2)
 */

//! Statistical Distribution
/*!
  A more elaborate class definition
 */

StatisticalDistribution::StatisticalDistribution() {} //! Constructor
StatisticalDistribution::~StatisticalDistribution() {}


// Constructor/Destructor for the inherited class

StandardNormalDistribution::StandardNormalDistribution() {}
StandardNormalDistribution::~StandardNormalDistribution() {}

// Probability density function

double StandardNormalDistribution::pdf(const double& x) const {
  return (1.0 / sqrt(2.0 * M_PI)) * exp(-0.5 * x * x);
}

// Cumulative density function

double StandardNormalDistribution::cdf(const double& x) const {
  double k = 1.0 / (1.0 + 0.2316419 * x);
  double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255 + 1.330274429 * k))));

  if (x >= 0.0) {
    return (1.0 - (1.0 / (pow(2 * M_PI,0.5))) * exp(-0.5 * x * x) * k_sum);
  } else {
    return 1.0 - cdf(-x); 
  }
}

// Inverse cumulative distribution function (aka the probit function)
double StandardNormalDistribution::inv_cdf(const double& quantile) const {
  // This is the Beasley-Moro algorithm which can
  // be found in Glasserman. We won't go into the
  // details here

  static double a[4] = {2.5066, -18.615, 41.3911, -25.4410};
  static double b[4] = {-8.473, 23.083, -21.06, 3.1308};
  static double c[9] = {0.3347, 0.976, 0.160, 0.0276, 0.0038, 0.00039, 0.000321, 0.000002, 0.000003};

  if (quantile >= 0.5 && quantile <= 0.92) {

    double num = 0.0;
    double denom = 1.0;

    for (int i = 0; i < 4; i++) {

      num += a[i] * pow((quantile - 0.5), 2*i + 1);
      denom += b[i] * pow((quantile - 0.5), 2 * i); 
    }

    return num/denom;
  } else if (quantile > 0.92 && quantile < 1) {

    double num = 0.0;

    for (int i = 0; i < 9; i++) {
      num += c[i] * pow((log(-log(1-quantile))), i);
    }

    return num;
  }
}


// Expectation/mean
double StandardNormalDistribution::mean() const {return 0.0;}
// Variance
double StandardNormalDistribution::var() const {return 1.0;}
// Standard Deviation
double StandardNormalDistribution::stddev() const {return 1.0;}


// Obtain a seqeuence of random draws from this distribution
void StandardNormalDistribution::random_draws(const std::vector<double>& uniform_draws, std::vector<double>& dist_draws) {

  /*
    The simplest method is to calculate this is with the Box-Muller method, 
    which has been used procedurally in many other chapters

    Check that the uniform draws and dist_draws are the same size 
    and have an even number of elements 
   */

  if (uniform_draws.size() != dist_draws.size()) {
    std::cout << "" << std::endl;
    return;

  }

  // Slow, but easy to implement
  for (int i = 0; i< uniform_draws.size()  / 2; i++) {
    dist_draws[2*i] = sqrt(-2.0 * log(uniform_draws[2*i])) * sin(2 * M_PI * uniform_draws[2*i+1]);
    dist_draws[2 * i+1] = sqrt(-2.0 * log(uniform_draws[2 * i])) * cos (2*M_PI*uniform_draws[2*i+1]);
  }
  return;
}


#endif 
