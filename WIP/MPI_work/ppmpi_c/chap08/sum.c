/* sum.c -- add two vectors using cyclic distribution of arrays.  Program
 *     to illustrate use of cyclic_io functions.
 *
 * Input: 
 *     n:  order of vectors
 *     x, y:  the vectors being added
 *
 * Output:
 *     z: the sum vector
 *
 * Notes:  Compile with Makefile.sum
 *
 * See Chap 8, pp. 170 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"

/* Header file for the basic I/O functions */
#include "cio.h"

/* Header file for the cyclic array I/O functions */
#include "cyclic_io.h"

main(int argc, char* argv[]) {
    CYCLIC_ARRAY_STRUCT  x;
    CYCLIC_ARRAY_STRUCT  y;
    CYCLIC_ARRAY_STRUCT  z;
    int                  n;
    MPI_Comm             io_comm;
    int                  i;

    MPI_Init(&argc, &argv);

    /* Build communicator for I/O */
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    if (Cache_io_rank(MPI_COMM_WORLD, io_comm) ==
            NO_IO_ATTR)
        MPI_Abort(MPI_COMM_WORLD, -1);

    /* Get n */
    Cscanf(io_comm, "Enter the array order", "%d", &n);

    /* Initialize scalar members.  Calls  */
    /* function for building derived type */
    Initialize_params(&io_comm, n, &x);
    Initialize_params(&io_comm, n, &y);
    Initialize_params(&io_comm, n, &z);

    /* Get vector elements */
    Read_entries("Enter elements of x", &x);
    Read_entries("Enter elements of y", &y);

    /* Add local entries */
    for (i = 0; i < Local_size(&x); i++)
        Local_entry(&z,i) =
            Local_entry(&x,i) + Local_entry(&y,i);

    /* Print z */
    Print_entries("x + y =", &z);

    MPI_Finalize();
}
