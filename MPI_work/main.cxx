#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cassert>
#include <array>

#include "mpi.h"

double total_d, local_d;

void Build_mpi_type (double* a_p, double* b_p, int* n_p, MPI_Datatype* input_mpi_t_p) {
  int array_of_blocklengths[3] = {1,1,1};
  MPI_Datatype array_of_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint array_of_displacements[3] = {0};
  MPI_Aint a_addr, b_addr, n_addr;

  MPI_Get_address(a_p, &a_addr);
  MPI_Get_address(b_p, &b_addr);
  MPI_Get_address(n_p, &n_addr);
  
  array_of_displacements[1] = b_addr - a_addr;
  array_of_displacements[2] = n_addr - a_addr;
  MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements, array_of_types, input_mpi_t_p);
  MPI_Type_commit(input_mpi_t_p);

}

void Get_Input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) {
  MPI_Datatype input_mpi_t;  
  Build_mpi_type(a_p, b_p, n_p, &input_mpi_t);
  if (my_rank == 0) {
    std::cout << "Enter a, b and n \n";
    scanf("%lf %lf %d", a_p, b_p, n_p);
    MPI_Bcast(a_p, 1, input_mpi_t, 0, MPI_COMM_WORLD);
  } else {    
    total_d = *a_p;
    for (int source = 1; source < comm_sz; source++) {
      MPI_Recv(&a_p, 1, input_mpi_t, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      std::cout << total_d << std::endl;
      total_d += total_d;
    }
  }

  MPI_Type_free(&input_mpi_t);


}


// MPI headers

#include "MPI_IO.hpp"
//#include "MPI_broadcast.hpp"
//#include "MPI_functions.hpp"

std::vector<int> bb;

int main(int argc, char** argv) {
  double a = 3.0;
  double b = 4.0;
  int c = 5;

  int my_rank, comm_sz;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  Get_Input(0, comm_sz, &a, &b, &c);

  for (int q = 1; q < comm_sz; q++) {
    std::cout << q << a << b << c << std::endl;
  }

  std::cout << total_d << std::endl;
  
  MPI_Finalize();
  //MPI_input a(3, 5);
  // a.Get_data();
  
  /*
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0) {
    int value = 17;
    int result = MPI_Send(&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

    if (result == MPI_SUCCESS)
      std::cout << "Rank 0 OK!" << std::endl;
  } else if (rank == 1) {
    int value;
    int result = MPI_Recv(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (result == MPI_SUCCESS && value == 17)
      std::cout << "Rank 1 OK!" << std::endl;
  }
  
  MPI_Finalize();
  */
 

  return 0;
}
