/*
  Tree structured-program looking at MPI_Broadcast 
  A communiation pattern that involves all the processes in a communiator is a 
  collective communication. As a consequence, a collective commmunication involves more 
  than two processes. 

  A broadcast is a collective communcation in which a single process
  sends the same data to every proceess in the communicator.
*/

#ifndef __MPI_BC__
#define __MPI_BC__

#include <iostream>
#include <vector>
#include <algorithm>
#include "mpi.h"
#include <numeric>

// cppunit tests

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

/*
  For the constructors that take a size, we will allocate 
*/

template <typename T, typename Q, typename R>  class MPI_BC_Generic {
public:
  MPI_BC_Generic();
  MPI_BC_Generic(std::size_t n);
private:
  T* data;
  T* limit;
};

/*
  Inheritance in C++

  1. Public mode: If we derive a sub class from a public base class, then the public members 
     of the base class will become public in the derived class and protected members of the 
     base class will become protected in the derived class (public -> public) 

  2. Protected mode: If we derive a sub class from a protected base class, then both public 
     member and protected members of the base class will become protected in derived 
     class

  3. Private mode: If we derive a sub class from a Private base class. Then both public 
  member and protected members of the base class 

     
 */

class MPI_BC : public CppUnit::TestCase {
public:
  MPI_BC(); // Default constructor 
  //  MPI_BC(int) : CppUnit::TestCase(name) {} // Test constructor  
  virtual ~MPI_BC();   
  void packData(); // Using MPI_Pack/MPI_unpack
  void time_ellapsed(); // Total time for the MPI program to execute
  void broadcast_input(); // broadcasting values
  void broadcast_vector(); // broadcasting vector
  void buildMpiType(double*, double*, int*, MPI_Datatype*); // Allocates the data type into a struct and broadcasts it 
  void Send(float, float, int, int); // Standard send/receive pair 
  void SendVector(); // Standard send/receive pair 
  void Receive(float*, float*, int*, int); // Standard send/receive pair 
  // void GetData(float*, float*, int*, int, int); // WIP
  void parallelAllocateVec(double*, double*, int, std::vector<int>*, MPI_Datatype*);
 
  // Functions yet to be implemented 
  // void add_vector();
  // void vectorTest(); // unit test function  
private:
  float* a_p;
  float* b_p;
  int* n_p;
  int my_rank, comm_sz;
  int lenOfVec;
  int* pointerToArray;
  double start, finish; // integers for measuring the start and finish of the function called from this class
  MPI_Aint aint; // What does MPI_Aint mean?
  MPI_Status status;
  MPI_Datatype AA;

  // Vectors 
  std::vector<int> vectorOfBlockLengths;
  std::vector<int> MPItype;
  std::vector<int> v;
  std::vector<MPI_Datatype> MPIDatatype;
  std::vector<MPI_Aint> MPIdisplacements;
  std::vector<std::vector<int> > TdVector;
  // Testing 
};

#endif
