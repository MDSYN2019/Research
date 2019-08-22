/*
  One obvious problem with our program is its lack of generality. The function, 
  f(x), and the input data, a, b, n, are hardwired. So if we want to change any 
  of these, we must edit and recompule the program. 
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

/
class MPIInput {

public:
  MPIInput(); /*!< Default constructor */ 
  MPIInput(int, int); /*!< int int constructor */
  void MPIStart();
  void getData();
  void bubbleSort();
  void oddEvenSort();
  void I_send();

  virtual ~MPIInput();
  
private:
  int* n_ptr;
  int source = 0;
  int dest;
  int tag;
  int my_rank;
  int p;
  int stage;
  //  
  float* a_ptr;
  float* b_ptr;
  MPI_Status status;
};

#endif 
