#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <mpi.h>
#include "MPI_broadcast.hpp"

std::map<std::string, std::string> typeConvDict; // TODO

enum MPI_TYPE {MPI_CHAR, MPI_SHORT, MPI_INT, MPI_LONG, MPI_LONG_LONG, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_SHORT, MPI_UNSIGNED, MPI_UNSIGNED_LONG, MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE, MPI_BYTE, MPI_PACKED};

MPI_BC::MPI_BC() {
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
} // constructor 

MPI_BC::~MPI_BC() {
  MPI_Finalize();
} // destructor 

void MPI_BC::parallelAllocateVec(double* aa, double* bb, int lenOfVec, std::vector<int>& vecpart) {
  int array_of_blocklengths =
  // Vector related operations
    
  MPI_Datatype type = {MPI_DOUBLE};
  MPI_Get_address(&vecpart, ..);
  vectorofBlocklengths.push_back(lenOfVec); // Add in length of vector  
  MPI_Get_address(&vectorOfBlockLengths[0], aint);
  MPI_Datatype type = {MPI_INT};
  MPI_Type_create_struct();
  
}

void MPI_BC::buildMpiType(double* a_p, double* b_p, int* n_p, MPI_Datatype input_mpi_t_p) {
  /*
    
    A derived datatype can bbe used to represent any collection 
    of data items by storing both the types of items and their 
    relative locations in memory. 

    If a function that sends data knows the types and the relative 
    locations in memory of a collection of data items, 
    
  */
  
  int array_of_blocklengths[3] = {1,1,1};

  MPI_Datatype array_of_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint a_addr, b_addr, n_addr;
  MPI_Aint array_of_displacements[3] = {0};  

  MPI_Get_address(a_p, &a_addr);
  MPI_Get_address(b_p, &b_addr);
  MPI_Get_address(n_p, &n_addr);

  array_of_displacements[1] = b_addr - a_addr;
  array_of_displacements[2] = n_addr - a_addr;
  MPI_Type_create_struct();
  MPI_Type_commmit(input_mpi_t_p);
}
// Build MPI type

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

void MPI_BC::Send(float a, float b, int n, int dest) {
    MPI_Send(&a, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
    MPI_Send(&b, 1, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
    MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
} /* Send */

void MPI_BC::Receive(float* a_ptr, float* b_ptr, int* n_ptr, int     source) {
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(n_ptr, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
} /* Receive */

void MPI_BC::GetData() {
  if (my_rank == 0) {
    std::cout << "Enter a, b, and n \n";
    scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
  }
  
  
}
