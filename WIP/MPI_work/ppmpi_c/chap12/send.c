/* send.c -- wrapper for MPI_Send that takes timing information.
 *     Should be compiled and linked with ./bcast.c
 *
 * Input: none
 * Output: Print_send_time prints total elapsed time spent in MPI_Send 
 *
 * See Chap 12, pp. 271 & ff in PPMPI.
 */
#include <stdio.h>
#include "mpi.h"

static double send_time = 0.0;

int MPI_Send(void* buffer, int count,
        MPI_Datatype datatype, int dest, int tag, MPI_Comm comm) {
    double start_time;
    double finish_time;
    int return_val;

    start_time = MPI_Wtime();
    return_val = PMPI_Send(buffer, count, datatype,
                      dest, tag, comm);
    finish_time = MPI_Wtime();
    send_time = send_time + finish_time - start_time;
    return return_val;
}

void Print_send_time(void) {
    int my_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    printf("Process %d > Total time in MPI_Send = %e\n", my_rank, send_time);
}
