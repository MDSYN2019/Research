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

// TODO : Testing 
MPI_Datatype getType(std::string* inputString) { 
  switch (inputString) {
  case 'i':
    return MPI_INT;
    break;
  case 'f':
    return MPI_FLOAT;
    break;
  case 'j':
    return MPI_UNSIGNED;
    break;
  case 'd':
    return MPI_DOUBLE;
    break;
  case 'c':
    return MPI_CHAR;
    break;
  case 's':
    return MPI_SHORT;
    break;
  case 'l':
    return MPI_LONG;
    break;
  case 'm':
    return MPI_UNSIGNED_LONG;
    break;
  case 'b':
    return MPI_BYTE;
    break;
  }
}

template <class ItemType>
class InitiateVectorMethod {
public:
  InitiateVectorMethod(int, int) {};
  virtual ~InitiateVectorMethod() {};

  // void public  methods

  void setup(int*); 
  void traits();
  void SendVector();
private:
  MPI_Group  group_world;                                                                                                                           
  MPI_Group  first_row_group;                                                                                                                       
  MPI_Comm   first_row_comm;       
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

/*

Pack/Unpack

An alternative approach to grouping data is provided by the MPI functions MPI_Pack 
and MPI_unpack functions. MPI_Pack allows one to explicitly store noncintiguous data in 
contiguous memory locations 

 */
