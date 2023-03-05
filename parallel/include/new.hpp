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

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

template <class ItemType>
class InitiateVectorMethod {
public:
  InitiateVectorMethod(int, int) {};
  virtual ~InitiateVectorMethod() {};

  // void public  methods

  void setup(int*); 
  void traits();
  void SendVector();
  void GetData();
  
private:
  MPI_Group  group_world;  
  MPI_Group  first_row_group;                                        
  MPI_Comm   first_row_comm;
  MPI_Status status;
  int my_rank, comm_sz;
  std::deque<ItemType> A;
  int var1, var2;  
};

class VectorMethodTest : public CppUnit::TestFixture {
private: // Need to check if this is indeed private 
  CPPUNIT_TEST_SUITE();
  CPPUNIT_TEST();
  CPPUNIT_TEST_SUITE_END();  
public:
  void setUp();
  void tearDown();
  void testConstructor();
};

#endif
