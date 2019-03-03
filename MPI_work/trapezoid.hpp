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
  
private:
  int my_rank;
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

