#ifndef __placeholder__
#define __placeholder__

#include <iostream>
#include <iomanip>
#include <string>
#include <cassert>
#include <cmath>
#include <array>
#include <cassert> 
#include <deque>

#include "mpi.h"

template <class ItemType>
class InitiateVectorMethod {
public:
  InitiateVectorMethod() {};
  virtual ~InitiateVectorMethod() {};
  // Void methods
  void setup(int* ); 
  void traits();

private:
  int my_rank, comm_sz;
  std::deque<ItemType> A;
};

#endif

