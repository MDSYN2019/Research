#include <iostream>
#include <iomanip>
#include <vector>
#include "mpi.h"

#include "MPI_broadcast.hpp"

MPI_BC::MPI_BC() {
} // constructor 

MPI_BC::~MPI_BC() {
} // destructor 

void MPI_BC::build_mpi_type(double* a_p, double* b_p, int* n_p, MPI_Datatype input_mpi_t_p) {

  int array_of_blocklengths[3] = {1,1,1};
  MPI_Datatype array_of_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint a_addr, b_addr, n_addr;
  MPI_Aint array_of_displacements[3] = {0};
  
  
}

void MPI_BC::Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) { // input, input, input, output, output
  if (my_rank == 0) {
    std::cout << "Enter a, b and n \n";
    scanf("%lf %lf %d", a_p, b_p, n_p);
    MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(n_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }
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


void MPI_BC::Get_data2() {
  if (my_rank == 0) {
    std::cout << "Enter a, b, and n \n";
    scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
  }  
}
