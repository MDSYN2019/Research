//! Comments on Code

/*!

Performance evaluation of MPI programs

We're usually not interested in the time taken from the sstart of program execution 
to the end of program execution. 

We're only interested in the time it takes to do the actual multiplicaton, so we need to modify our source code by adding in calls to a function that will tell us the amount of time that elapses from the beginning to the end of the actual matrix

*/


#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include "mpi.h"
#include "MPI_broadcast.hpp"

// cppunit tests

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

// Trying to get this vector understood

std::map<std::string, std::string> typeConvDict; // TODO


//! bcast definition
/*! my_bcast function does the following..

 */

void my_bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator) {

  int world_rank; /*!< int world rank */
  int world_size; /*!< int world size */
  
  MPI_Comm_rank(communicator, &world_rank);
  MPI_Comm_size(communicator, &world_size);

  if (world_rank == root) {
    
    // If we are the root process, send our data to everyone
   
    for (int i = 0; i < world_size; i++) {
      if (i != world_rank) {
	MPI_Send(data, count, datatype, i, 0, communicator);
      }
    }
  } else {
    // If we are a receiver process, receive the data from the root
    MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
  }
}

//! MPI_BC class template filling 
/*! 
  Placeholder
*/

MPI_BC::MPI_BC() {

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

}

void MPI_BC::parallelAllocateVec(double* aa, double* bb, int lenOfVec, std::vector<int>* vecpart, MPI_Datatype* input_mpi_t_p) {

  std::iota(MPItype.begin(), MPItype.end(), 1); // Vector allocation of types
  std::iota(MPIDatatype.begin(), MPIDatatype.end(), MPI_INT); // Vector allocation of MPI_INt
  std::iota(MPIdisplacements.begin(), MPIdisplacements.end(), sizeof(int));  // vector allocation of the size of the vector 

  pointerToArray = &MPItype[0];  
  MPI_Get_address(&vecpart[0], &aint);
  //  MPI_Type_create_struct(lenOfVec, pointerToArray, MPIdisplacements, MPItype, input_mpi_t_p);
  MPI_Type_commit(input_mpi_t_p);
  finish = MPI_Wtime();
}


void MPI_BC::packData() {
}


void MPI_BC::buildMpiType(double* a_p, double* b_p, int* n_p, MPI_Datatype* input_mpi_t_p) {
  /*  
    A derived datatype can be used to represent any collection 
    of data items by storing both the types of items and their 
    relative locations in memory. 

    If a function that sends data knows the types and the relative 
    locations in memory of a collection of data items,   
  */
  
  start = MPI_Wtime();
  int array_of_blocklengths[3] = {1,1,1};

  MPI_Datatype array_of_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint a_addr, b_addr, n_addr;
  MPI_Aint array_of_displacements[3] = {0};  

  MPI_Get_address(a_p, &a_addr);
  MPI_Get_address(b_p, &b_addr);
  MPI_Get_address(n_p, &n_addr);

  array_of_displacements[1] = b_addr - a_addr;
  array_of_displacements[2] = n_addr - a_addr;
  //MPI_Type_create_struct();
  //MPI_Type_commmit(input_mpi_t_p);
}

// Build MPI type

/*
void MPI_BC::Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) { // input, input, input, output, output
  if (my_rank == 0) {
    std::cout << "Enter a, b and n \n";
    scanf("%lf %lf %d", a_p, b_p, n_p); 
  }
  MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Get Input
}

void MPI_BC::Get_input2(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) { // input, input, input, output, output

  MPI_Datatype input_mpi_t;
  buildMpiType(a_p, b_p, n_p, &input_mpi_t); // TODO
  // Get Input
  if (my_rank == 0) {

  }
}

*/
void MPI_BC::Send(float a, float b, int n, int dest) {
   start = MPI_Wtime();
   MPI_Send(&a, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
   MPI_Send(&b, 1, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
   MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
} /* Send */


void MPI_BC::SendVector() {
  
  for (unsigned int i = 0; i < 10; i++) {
    v.push_back(i);
  }

  //MPI_Datatype column_mpi_t;
  //MPI_Type_vector(10, 1, 10, MPI_INT, &column_mpi_t);
  //MPI_Type_commit(&column_mpi_t);

  //if (my_rank == 0) {
  //  MPI_Send(&v[0][1], 1, column_mpi_t, 1, 0, MPI_COMM_WORLD);
  //} else {
  //  MPI_Recv(&v[0][1], 1, column_mpi_t, 0, 0, MPI_COMM_WORLD, &status);
  // }
  // MPI_Send(&a, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
  // MPI_Send(&b, 1, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
  // MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
} /* Send */

void MPI_BC::Receive(float* a_ptr, float* b_ptr, int* n_ptr, int source) {  
  MPI_Recv(a_ptr, 1, MPI_FLOAT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(b_ptr, 1, MPI_FLOAT, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(n_ptr, 1, MPI_INT, source, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  finish = MPI_Wtime();
} /* Receive */

MPI_BC::~MPI_BC() {
  MPI_Finalize();
} // destructor 

// Test methods

//MPI_BC::vectorTest() {
//  CPPUNIT_ASSERT(v.size() == v.size());
//}
