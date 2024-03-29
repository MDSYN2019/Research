

#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <string>
#include "mpi.h"
#include "MPI_broadcast.hpp"


/*
  An alternative approach to grouping data is provided by the MPI functions
  MPI_Pack and MPI_Unpack, MPI_Pack allows one to explicitly store noncontiguous 
  data in contiguous memory locations, and MPI_unpack can be used to copy ..
 */

/*
void Pack_unpackGet(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank) {
  std::string buffer; // Keep data in the buffer 
  int position; 

  if (my_rank == 0) {
    std::cout << "Enter a, b, and n \n";
    position = 0;
    MPI_Pack(a_ptr, 1, MPI_FLOAT, buffer, 100, &position, MPI_COMM_WORLD); 
    MPI_Pack(b_ptr, 1, MPI_FLOAT, buffer, 100, &position, MPI_COMM_WORLD);
    
  }

}
*/

std::map<std::string, std::string> typeConvDict; // TODO
void my_bcast(void* data, int count, MPI_Datatype datatype, int root, MPI_Comm communicator) {
  int world_rank;
  MPI_Comm_rank(communicator, &world_rank);
  int world_size;
  MPI_Comm_size(communicator, &world_size);

  if (world_rank == root) {
    // If we are the root process, send our data to everyone
    int i;
    for (i = 0; i < world_size; i++) {
      if (i != world_rank) {
        MPI_Send(data, count, datatype, i, 0, communicator);
      }
    }
  } else {
    // If we are a receiver process, receive the data from the root
    MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
  }
}

MPI_BC::MPI_BC() {
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
} // constructor 


MPI_BC::~MPI_BC() {
  MPI_Finalize();
} // destructor 


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

void MPI_BC::buildMpiType(double* a_p, double* b_p, int* n_p, MPI_Datatype* input_mpi_t_p) {
  /*  
    A derived datatype can bbe used to represent any collection 
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
} 

void MPI_BC::broadcast_input() { // input, input, input, output, output
  if (my_rank == 0) {
    std::cout << "Enter a, b and n \n";
    scanf("%lf %lf %d", a_p, b_p, n_p); 
  }
  MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);  
}

void MPI_BC::add_vector() {  
  /*
    [ 1 2 3 4 5 6 7 8 9 10 ] 
    [ 1 2 3 4 5 6 7 8 9 10 ] 
    [ 1 2 3 4 5 6 7 8 9 10 ] 

    ... for a 10 X 10 matrix
   */

  for (int i = 0; i < 10; i++) {
    TdVector[i].push_back(i);
  }
  MPI_Type_vector(10, 1, 10, MPI_INT, &AA); // count, block_length, stride, element_type, new_mpi_t 
  MPI_Type_commit(&AA);  

  if (my_rank == 0) {
    MPI_Send(&(TdVector[0][1]), 1, AA, 1, 0, MPI_COMM_WORLD);  
  } else {
    MPI_Recv(&(TdVector[0][1]), 1, AA, 0, 0, MPI_COMM_WORLD, &status);
  }
}


void MPI_BC::broadcast_vector() { // input, input, input, output, output

  int* vecData = v.data();
  int count = v.size();
  int tag;
  // std::cout << v[0] << std::endl;
  if (my_rank == 0) {
    // Broadcast vector data - recv does not need to be called as 
    MPI_Bcast(v.data(), v.size() * sizeof(decltype(v)::value_type), MPI_BYTE, 0, MPI_COMM_WORLD);
  }
  else {
    //  MPI_Recv(v.data(), v.size() * sizeof(decltype(v)::value_type), MPI_BYTE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); 
    std::cout << " " << v[0] << " "  << " " << my_rank << " " << std::endl; 
  }
}

void MPI_BC::Receive(float* a_ptr, float* b_ptr, int* n_ptr, int   source) {  
  MPI_Recv(a_ptr, 1, MPI_FLOAT, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(b_ptr, 1, MPI_FLOAT, source, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  MPI_Recv(n_ptr, 1, MPI_INT, source, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  finish = MPI_Wtime();
} /* Receive */


void MPI_BC::vectorTest() {
//  CPPUNIT_ASSERT(v.size() == v.size());
}
