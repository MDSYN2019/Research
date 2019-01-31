#include <iostream>
#include <iomanip>
#include <vector>
#include <mpi.h>

// Custom headers
#include "MPI_IO.hpp"

MPI_input::MPI_input() {
} // constructor 

MPI_input::MPI_input(int mr, int pe) {
  my_rank = mr;
  p = pe;
}

void MPI_input::MPI_start() { 
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);  
}

void MPI_input::Get_data() {
  if (my_rank == 0) {
    std::cout << "Enter a, b and n \n";
    scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
    for (int dest = 1; dest < p; dest++) {
      tag = 0;
      MPI_Send(a_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 1;
      MPI_Send(b_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 2;
      MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
  } else {
    tag = 0;
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 1;
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 2;
    MPI_Recv(n_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  }
}

MPI_input::~MPI_input() {
  MPI_Finalize();
} // destructor 

