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

class MPI_BC {
public:
  MPI_BC();
  virtual ~MPI_BC();
  
  void Get_Input(int, int, double*, double*, int*); // WIP
  void buildMpiType(double*, double*, int*, MPI_Datatype*); // WIP

  void Send(float, float, int, int); //
  void Receive(float*, float*, int*, int);

  void GetData(float*, float*, int*, int, int); // WIP
private:
  int my_rank, comm_sz;
  MPI_Aint aint; // What does MPI_Aint mean?
  std::vector<int> vectorOfBlockLengths;
};

#endif
