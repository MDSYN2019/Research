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



#endif 

