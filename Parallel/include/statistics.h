#ifndef __STAT__
#define __STAT__

#include <cmath>
#include <vector>


/*!

One of the most common examples of concepts in quantitiative finance is that of a statistical distribtion.
Random variables play a huge part in quantitive financial modelling. Derivatives, pricing, cash-flow 
forceasting and quantitive trading all make use of statitiscal methods in some fashion

Many of the chapters within this book have made use of random number generators
in order to carry out pricing tasks.

In a nutshell, we are splitting the generation of (uniform integer) random numbers from 
draws of specific statistical distribution,s such taht we can use the statics classes
elsewhere withut bringing along the heavy random number generation functions.

Equally useful is the fact taht we will be able to "swap out" different random number 
generators for out statisics classes for reliability, extensibility and efficiency

 */

/**
 *
 */
//! A test class - find and replace from this template for future class definitions
/*!
  A more elaborate class definition
 */

class QTstyle_Test {
};
//! Statistical Distribution Class

/*! 
  We've specified pure virtual methods for the probability density function (pdf), cumulative
  density function (cdf), inverse cdf (inv_cdf), as well as descriptive statistics functions such as 
  as mean, var (variance) and stdev. 

  Finally, we have a method that takes in a vector of uniform random variables on the open interval (0,1), then fills 
  a vector of identical length with draws from the distribution

 */

class StatisticalDistribution {

 public:
  //! A constructor
  /*! 
    Statistical Distribution constructor 
   */
  StatisticalDistribution();

  //! Virtual destructor
  /*!
    A more elaborate explanation here 
   */
  virtual ~StatisticalDistribution(); // Virtual function - useful for polymorphism when inherting with other classes

  // Distribution functions
  virtual double pdf(const double& x) const = 0;
  virtual double cdf(const double& x) const = 0;

  // Inverse cumulative distribution functions
  virtual double inv_cdf(const double& quantile) const = 0;

  // Descriptive stats
  virtual double mean() const = 0; //! Variable 1 
  virtual double var() const = 0; //! Varable 2 
  virtual double stdev() const = 0; //! Variable 3

  // obtain a sequence of random draws from this distribution
  virtual void random_draws(const std::vector<double>& uniform_draws, std::vector<double>& dist_draws) = 0;

};

  //! Standard Normal Distribution Implementation
  /*!
    A more elaborate explanation here 
   */
class StandardNormalDistribution : public StatisticalDistribution {
 public:
  StandardNormalDistribution();
  virtual ~StandardNormalDistribution();

  // Distribution functions
  virtual double pdf(const double& x) const;
  virtual double cdf(const double& x) const;

  // Inverse cumulative distribution function (aka probit functions)
  virtual double inv_cdf(const double& quantile) const;


  // Descriptive stats
  virtual double mean() const; //! Equal to 0 
  virtual double var() const; //! Equal to 1
  virtual double stddev() const; //! Variable 1

  // Obtain a sequence of random draws from the standard normal distribution
  virtual void random_draws(const std::vector<double>& uniform_draws, std::vector<double>& dist_draws);
  
  

};

#endif
