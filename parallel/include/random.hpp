#ifndef __RANDOM__
#define __RANDOM__

#include <vector>

#define M_E        2.71828182845904523536
#define M_LOG2E    1.44269504088896340736
#define M_LOG10E   0.434294481903251827651
#define M_LN2      0.693147180559945309417
#define M_LN10     2.30258509299404568402
#define M_PI       3.14159265358979323846
#define M_PI_2     1.57079632679489661923
#define M_PI_4     0.785398163397448309616
#define M_1_PI     0.318309886183790671538
#define M_2_PI     0.636619772367581343076
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2    1.41421356237309504880
#define M_SQRT1_2  0.707106781186547524401

/*
------------------------------------------------------
Random Number Generators and Statistical Distributions
------------------------------------------------------


In this chapter, we are going to construct classes to help us encapsulate the generation of random
numbers. Random number generators (RNGs) are an essential tool in quantitiative finance 
as they are necessary for Monte Carlo simulations that power option pricing techniques

Overview
--------

Why write our own random numnber generators?

- It allows us to make use of psuedo-random numbers. These are sequences of numbers that posess the correct 
  statistical properties to "emulate" random numbers in order to improve the convergence of Monte Carlo simulations. 
  The interferace for random numbers and pseudo-random numners is identical and we can hide away the details in the 
  specific classes. In particular, we can implement low discrepancy numbers and anti-theric sampling in this manner

- Replying on the rand function 


 */


class RandomNumberGenerator {
  
protected:
  /*
    Since we're creaint an abstract base class, is it a good idea to be using a protected class?

    This is actually a contentious issue - Sometimes protected variables are frowned upon. Instead, it is arguemed taht 
    all data should be private and that accessor methods should be used.


   */

  
  unsigned long init_seed; // Initial random seed value
  unsigned long cur_seed; // Current random seed value
  unsigned long num_draws; // Dimensionality of the RNG
public:
  RandomNumberGenerator(unsigned long _num_draws, unsigned long _init_seed) : num_draws(_num_draws), init_seed(_init_seed), cur_seed(_init_seed) {};
  virtual ~RandomNumberGenerator() {};

  /*

    We then have four separate access and reset methods (all virtual), which get, set, and reset the random 
    seed and another which resets the number of random draws. They are all directly implemented in the header file, once again 
    stopping us from need to create a random.cpp source file

   */
  virtual unsigned long get_random_seed() const {return cur_seed;}
  virtual void set_random_seed(unsigned long _seed) {cur_seed = _seed;}

  virtual void set_num_draws(unsigned long _num_draws) {
    num_draws = _num_draws;
  }

  // Obtain a random integer (needed for creating random uniforms)
  virtual unsigned long get_random_integer() = 0;
};

#endif

