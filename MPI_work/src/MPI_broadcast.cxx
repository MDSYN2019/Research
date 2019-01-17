#include <iostream>
#include <iomanip>
#include <vector>
#include "MPI_broadcast.hpp"

MPI_BC::MPI_BC() {
} // constructor 

MPI_BC::~MPI_BC() {
} // destructor 

int MPI_BC::Ceiling_log2(int x) {
    /* Use unsigned so that right shift will fill                                                                                                 
     * leftmost bit with 0 
     */
    unsigned temp = (unsigned) x - 1;
    int result = 0;

    while (temp != 0) {
         temp = temp >> 1;
         result = result + 1 ;
    }
    return result;
} 

void MPI_BC::Send(float a, float b, int n, int dest) {
    MPI_Send(&a, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
    MPI_Send(&b, 1, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
    MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
} /* Send */


void MPI_BC::Receive(float* a_ptr, float* b_ptr, int* n_ptr, int     source) {
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(n_ptr, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
} /* Receive */


void MPI_BC::Get_data2() {
  if (my_rank == 0) {
    std::cout << "Enter a, b, and n \n";
    scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
  }
  
}
