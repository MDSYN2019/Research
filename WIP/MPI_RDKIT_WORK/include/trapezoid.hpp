/*  
  Implementation of the trapezium rule  - A form of numerical integration 
*/


#ifndef __TRAP_h_
#define __TRAP_h_

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
  Trap();
  ~Trap();
  void read();
  void computeTrapezium();
private:
  float integral; // Store result in integral
  float a, b; //
  int n;
  float h;
  float x;
  int i;
  std::stringstream ss; 
};



#endif 

