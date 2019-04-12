#include <iostream>
#include <cassert>
#include <cstdio>
#include <cmath>
#include "mpi.h"
#include "MPI_ping_pong.hpp"

MPI_PP::MPI_PP() {
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank); // What does this represent?
  MPI_Comm_size(MPI_COMM_WORLD, &world_size); 
  /*! If the number of processors is less than 2, then terminate */
  if (world_size != 2) {
    fprintf(stderr, "World size must be two for %s\n", argv[0]);
    MPI_Abort(MPI_COMM_WORLD, 1);
    return 1;
  }
}

MPI_PP::~MPI_PP() {
  // Nothing to do here - not deallocating anything 
}

MPI_PP::serial_dot(float*, float*, int ) {
  for (int i = 0; i < n; i++) {
    sum = sum + x[i] * y[i];
  }
  return sum;
}
