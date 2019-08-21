#ifndef __STAT_CXX__
#define __STAT_CXX__

#include "statistics.h"
#include <iostream>

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
  return (1.0 / sqrt(2.0 * M_PI);
}

#endif 
