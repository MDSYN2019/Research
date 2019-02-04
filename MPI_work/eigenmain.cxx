#include <iostream>
#include <Eigen/Dense>
#include <mpi.h>
using std::cin;
using std::cout;
using std::endl;
using namespace Eigen;

int main(int argc, char ** argv) {

  MatrixXd a = MatrixXd::Ones(3, 4);
  int myrank;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Status status;
  if (0 == myrank)
    {
      MPI_Send(a.data(), 12, MPI_DOUBLE, 1, 99, MPI_COMM_WORLD);
    }
  else if (1 == myrank)
    {
      MPI_Recv(a.data(), 12, MPI_DOUBLE, 0, 99, MPI_COMM_WORLD, &status);
      cout << "RANK " << myrank << endl;
      cout << a << endl;
    }
  MPI_Finalize();
  return 0;
}
