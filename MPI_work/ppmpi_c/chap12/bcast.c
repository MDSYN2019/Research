/* bcast.c -- linear loop broadcast, illustrates use of MPI's 
 *     profiling interface.
 * 
 * Input: none
 * Output: 
 *     Total time spent by each process in MPI_Send
 *     Result of broadcast -- should be "1" printed by each
 *
 * Note:  Should be linked with ./send.c and the MPI profiling
 *     library.
 *
 * See Chap 12, pp. 271 & ff in PPMPI.
 */
#include <stdio.h>
#include "mpi.h"

extern void Print_send_time(void);

main(int argc, char* argv[]) {
    int           p;
    int           my_rank;
    int           root = 0;
    MPI_Comm      comm;
    int           proc;
    MPI_Status    status;
    int           size = 1;
    int           tag = 0;
    MPI_Datatype  datatype = MPI_INT;
    int*          x;
    int           y;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    if (my_rank == root)
        y = 1;
    else
        y = 0;
    x = &y;

    if (my_rank == root) {
        for (proc = 0; proc < p; proc++)
            MPI_Send(x, size, datatype, proc, tag, comm);
    }
    MPI_Recv(x, size, datatype, root, tag, comm, &status);

    Print_send_time();
    printf("Process %d > *x = %d\n", my_rank, *x);
    MPI_Finalize();
}  /* main */
