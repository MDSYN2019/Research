#include <Eigen/Dense>
#include <cstdio>
#include <string>
#include <iostream>
#include <mpi.h>
#include <map> // most commonly used associative container
#include <vector> // most commonly used incremental container
#include "MPI_functions.hpp"

/*

1. Partition the program solution into tasks

2. Identify the communcation channels between the tasks

3. Aggregate the tasks into composite tasks

4. Map the composite tasks to cores

 */


const int MAX_STRING = 100;

int main (void) {

  int n = 1024; 
  char greeting[MAX_STRING];
  int comm_sz; // number of processes
  intt my_rank; // my process rank

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  if (my_rank != 0) { // If the CPU->Memory rank is not 0, we have to send it to the 0th rank

    sprintf(greeting, "Greetings from process %d of %d!", my_rank, comm_sz);

    MPI_Send(greeting, strlen(greeting)+1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
  }

  else {
    std::cout << "Greetings from process %d of %d\n" <<  my_rank << comm_sz;
    for (int q = 1; q < comm_sz; q++) {
      MPI_Recv(greeting, MAX_STRING, MPI_CHAR, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("%s\n", greeting); 
    }
  }
  
  MPI_Finalize();
  return 0;
}
