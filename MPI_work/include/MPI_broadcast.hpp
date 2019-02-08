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
#include <mpi.h>
#include <numeric>


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


class MPI_BC {
public:
  MPI_BC(); // Default constructor 
  explicit MPI_BC(std::size_t n,   ); 
  // MPI_BC(int); // Custom constructor 
  virtual ~MPI_BC();   
  void packData(); // Using MPI_Pack/MPI_unpack
  void time_ellapsed(); // Total time for the MPI program to execute

  void buildMpiType(double*, double*, int*, MPI_Datatype*); // Allocates the data type into a struct and broadcasts it 
  void Send(float, float, int, int); // Standard send/receive pair 
  void Receive(float*, float*, int*, int); // Standard send/receive pair 
  // void GetData(float*, float*, int*, int, int); // WIP
  void parallelAllocateVec(double*, double*, int, std::vector<int>*, MPI_Datatype*);
private:
  int my_rank, comm_sz;
  int lenOfVec;
  int* pointerToArray;
  double start, finish; // integers for measuring the start and finish of the function called from this class
  MPI_Aint aint; // What does MPI_Aint mean?
  std::vector<int> vectorOfBlockLengths;
  std::vector<int> MPItype;
  std::vector<MPI_Datatype> MPIDatatype;
  std::vector<MPI_Aint> MPIdisplacements;
};

#endif
