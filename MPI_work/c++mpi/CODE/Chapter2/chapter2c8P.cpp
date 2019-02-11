

#include <iostream.h>
#include <math.h>
#include <mpi.h>


void main(int argc, char **argv){

  MPI_Init(&argc,&argv);

  cout << "HELLO!!" << endl;

  MPI_Finalize();
}


