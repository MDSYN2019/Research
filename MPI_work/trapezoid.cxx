/*  
    Implementation of the trapezium rule  - A form of numerical integration 

*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <cmath>
#include <array>
#include <cassert>
#include "mpi.h"

#include "trapezoid.hpp"

Trap::Trap() {}
Trap::~Trap() {}

void Trap::read() {
  std::cout << "Please input Values:" << std::endl;
  std::cin >> a >> b >> n;
  std::cout << a << b << n << std::endl;
}

void Trap::computeTrapezium() {

  h = (b - a)/n;
  integral = ((float)(a) + (float)(b))/2.0;
  x = a;
  
  for (unsigned int i = 1; i <= n -1; i++) {
    x = x + h;
    integral = integral + (float)(x);
  }  

  integral = integral * h;
}
// S
MPITrap::MPITrap() {
  MPI_Init(&argc, &argv); // 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // allocate rank to address of my_rank 
  MPI_Comm_size(MPI_COMM_WORLD, &p); // allocate size to p    
}

MPITrap::~MPITrap() {
  MPI_Finalize();
}

void MPITrap::Trap() {
  integral = ((float)(local_a) + (float)(local_b)) / 2.0;
  x = local_a;

  for (unsigned int i = 1; i <= local_n-1; i++) {
    x = x + h;
    integral = integral + (float)(x);
  }
  integral = integral * h;
}

void MPITrap::getData() {

}

void MPITrap::AddIntegral() {
  if (my_rank == 0) {
    total = integral;
    for (int source = 1; source < p; source++) {
      tag = 0;
      MPI_Send(a_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 1;
      MPI_Send(b_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 2;
      MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 1;
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 2;
    MPI_Recv(n_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  } // Get data 
}


