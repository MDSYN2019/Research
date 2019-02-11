#include <iostream>
#include <ctime>
#include <iomanip>
#include <mpi.h>

/*


*/

void timestamp();

int main (int argc, char **argv) {

  int mynode, totalnodes;
  int datasize; // number of data units to be sent/received
  int sender; // process number of the sending process
  int receiver; // process number of the receiving process 
  int tag; // integer message tag

  // --
  int i;
  double x[100];
  double dx = 1.0/99.0;

  for (int i = 0; i <= 100; i++) {
    std::cout << 
    
  }
  
  // --
  MPI_Status status;
  
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &totalnodes); // Passes in the number of processes the parallel program is running with 
  MPI_Comm_rank(MPI_COMM_WORLD, &mynode); // identifies the rank 

  // Determine datasize 

  double * databuffer = new double[datasize];

  // Fill in sender, receiver and tag on sender/receiver processes
  // and fill in databuffer on the sender process

  if (mynode == sender) {
    MPI_Send(databuffer, datasize, MPI_DOUBLE, receiver, tag, MPI_COMM_WORLD);
  }

  if (mynode == receiver) {

    MPI_Recv(databuffer, datasize, MPI_DOUBLE, sender, tag, MPI_COMM_WORLD, &status);
    
  }
  
  // Senc Recv complete
  
  std::cout << "Hello world from process " << mynode;
  std::cout << " of " << mynode;
  MPI_Finalize();
  return 0;  
}
