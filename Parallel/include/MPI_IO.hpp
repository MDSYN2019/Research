/*
MPI Types

MPI_CHAR - signed char
MPI_SHORT - signed short int
MPI_INT - signed int
MPI_LONG - signed long int
MPI_UNSIGNED_CHAR - unsigned char
MPI_UNSIGNED_SHORT - unsigned short int
MPI_UNSIGNED - unsigned int
MPI_UNSIGNED_LONG - unsigned long int
MPI_FLOAT - float
MPI_DOUBLE - double
MPI_LONG_DOUBLE - long double

MPI_BYTE
MPI_PACKED


On most parallel systems, the processes involved in the execution of a parallel program are identified 
by a sequence of nonnegative integers. If there are p processes executing a program, hey will have ranks

0, 1, ... p-1. 

So one possibility here is for each process other than 0 to send a message to process 0. 



 */

#ifndef __MPI_IO__
#define __MPI_IO__

#include <iostream>
#include <iomanip>
#include <vector>
#include <mpi.h>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


class MPIInput : public CppUnit::TestCase { // Inherit from Cppunittest

public:
  //  MPIInput(); /*!< Default constructor - at the moment disabled*/  
  MPIInput(int, int); /*!< int int constructor */
  void getData(int*, int*, int*);
  void bubbleSort(int*, int);
  void getDataPack(float*, float*, int*);
  ~MPIInput(); 
  // Test functions
  void test1();
  
private:
  int source = 0; /*!< process sending integral */ 
  int dest = 0; /*!< All messages go to 0 */
  int tag = 0; /*!< The tag int is used for marking MPI messages */
  int my_rank;
  int p;
  
  int* start_inp;
  int* end_inp;
  int* n_ptr;

  // MPI pack pointers
  float* a_ptr;
  float* b_ptr;
  int B = 23;
  int C = 24;
  MPI_Status status;

};

#endif 
