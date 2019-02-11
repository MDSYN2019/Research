#ifndef __MD_Algorithm__
#define __MD_Algorithm__

#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cstdio>

// Object flips
class flips {

public:

  flips(int); // constructor flips
  ~flips(); // destructor flips
  void printFlips();  
  void printDiff();
  // private variables
  
private:
  int TT;
  const gsl_rng_type * T; // random number generators
  gsl_rng * r;
  int i, n;
  double mu;
  
  int heads;
  int tails;
  int difference; 
};

#endif 
