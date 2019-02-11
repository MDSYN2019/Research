#include <iostream>
#include <sstream>
#include <string>
#include <mpi.h>
#include <cstring>

const char *MY_NAME = "Zed A. Shaw";


void scope_demo(int count) {


}


const int MAX_STRING = 100;


int main(void) {
  std::string greeting;
  int comm_sz; // number of processes
  int my_rank; // my process rank

  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_size(MPI_COMM_WORLD, &my_rank);

  if (my_rank != 0) {

    sprintf(greeting, "Greetings from process %d of %d", my_rank, comm_sz);

    MPI_Send(greeting, strlen(greeting) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);

  } else {

    printf("Greetings from process %d of %d \n", my_rank, comm_sz);
    for (int q = 1; q < comm_sz; q++) {

      MPI_Recv(greeting, MAX_STRING, MPI_CHAR, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      printf("%s\n", greeting);
      
    }
  }

  MPI_Finalize();

  return 0;
  
}
