/* arg_test.c -- test whether an MPI implementation distributes command
 *     line arguments among the processes
 *
 * Input:  command line arguments
 * Output:  command line arguments accessible to each process
 *
 * Notes: Link with cio.o
 *
 * See Chap 8, pp. 154 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"

/* Header file for our I/O library */
#include "cio.h"

main(int argc, char* argv[]) {
    MPI_Comm  io_comm;
    int       i;

    MPI_Init(&argc, &argv);

    /* Set up communicator for I/O */
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);

    for (i = 0; i < argc; i++)
        Cprintf(io_comm,"","argv[%d] = %s", i, argv[i]);

    MPI_Finalize();
}
