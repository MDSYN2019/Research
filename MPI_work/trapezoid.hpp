/*  
  Implementation of the trapezium rule  - A form of numerical integration 
*/
#ifndef __TRAP__
#define __TRAP__

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <array>
#include <cassert>
#include "mpi.h"

// A simple string class

class Str {
public:
  typedef Vec<char>::size_type size_type;
  
};

template <class T> class {
};

/*
1. First we called a serial algorithm for solving our problem - we studied the trapezoidal rule
   for estimating a definite integral 

2. In order to parallelize the serial algorithm, we simply partitioned the data 
 */
class Trap {
public:
  void computeTrapezium();
  void read();
private:
  float integral; // Store result in integral
  float a, b; // 
  float h;
  float x;
  int i;
  std::stringstream ss; 
};

class MPITrap : private Trap {
public:
  MPITrap();
  ~MPITrap();
  void Trap();
  void getData();
  
private:
  int my_rank; // used in constructor 
  int p;
  float local_a;
  float local_b;
  float x;
  int i;
  int local_n;
  float total;
  int source;
  int dest = 0;
  int tag = 0;
  MPI_Status status;
};

#endif 

