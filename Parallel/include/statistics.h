#ifndef __STAT__
#define __STAT__

#include <cmath>
#include <vector>


/*!

One of the most common examples of concepts in quantitiative finance is that of a statistical distribtion.
Random variables play a huge part in quantitive financial modelling. Derivatives, pricing, cash-flow 
forceasting and quantitive trading all make use of statitiscal methods in some fashion


Many of the chapters within this book have made use of random number generators
in order to carry out pricing tasks 


 */

/**
 *
 */
class StatisticalDistribution {

 public:
  StatisticalDistribution();
  virtual ~StatisticalDistribution();  

  // Distribution functions
  virtual double pdf(const double& x) const = 0;
  virtual double cdf(const double& x) const = 0;

 private:
};

#endif
