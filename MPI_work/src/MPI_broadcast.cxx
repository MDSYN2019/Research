#include <iostream>
#include <iomanip>
#include <vector>
#include "MPI_broadcast.hpp"

MPI_BC::MPI_BC() {
} // constructor 

MPI_BC::~MPI_BC() {
} // destructor 

int MPI_BC::Ceiling_log2(int  x  /* in */) {
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

void MPI_BC::Send(float  a     /* in */,
        float  b     /* in */,
	int    n     /* in */,
        int    dest  /* in */) {

    MPI_Send(&a, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
    MPI_Send(&b, 1, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
    MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
} /* Send */


void MPI_BC::Receive(
        float*  a_ptr  /* out */,
	float*  b_ptr  /* out */,
        int*    n_ptr  /* out */,
        int     source /* in  */) {

    MPI_Status status;

    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, 0,
	MPI_COMM_WORLD, &status);
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, 1,
        MPI_COMM_WORLD, &status);
    MPI_Recv(n_ptr, 1, MPI_INT, source, 2,
        MPI_COMM_WORLD, &status);
} /* Receive */
