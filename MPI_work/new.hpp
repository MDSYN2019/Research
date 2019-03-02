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
  InitiateVectorMethod(int, int) {};
  virtual ~InitiateVectorMethod() {};

  // void public  methods
  void setup(int* ); 
  void traits();
  void SendVector();

private:
  int my_rank, comm_sz;
  std::deque<ItemType> A;
  int var1, var2;
};

#endif

