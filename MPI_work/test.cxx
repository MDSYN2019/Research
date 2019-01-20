#include <iostream>
#include <iomanip>
#include <mpi.h>

void Get_input(int my_rank, int comm_sz, double* a_p, double b_p, int* n_p) {
}

int main(void) {

  int my_rank, comm_sz;

  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  std::cout << "Proc %d of %d > Does anyone have a toothpick \n" << my_rank << comm_sz;

  MPI_Finalize();
  return 0;
}
