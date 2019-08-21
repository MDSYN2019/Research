/* cache_test.c -- cache and retrieve a process rank with a communicator
 *
 * Input: none
 * Output: Message from the process whose rank was cached
 *
 * See Chap 8, pp. 139 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"
#include <stdlib.h>

main(int argc, char* argv[]) {
    int p;
    int my_rank;
    MPI_Comm  io_comm;      /* Communicator for I/O      */
    int       IO_KEY;       /* I/O process attribute key */
    int*      io_rank_ptr;  /* Attributes are pointers   */
    void*     extra_arg;    /* Unused                    */
    int       flag;
    int*      io_rank_attr;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Get a separate communicator for I/O functions by */
    /*     duplicating MPI_COMM_WORLD                   */
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);

    /* Create the attribute key */
    MPI_Keyval_create(MPI_DUP_FN, MPI_NULL_DELETE_FN,
        &IO_KEY, extra_arg);

    /* Allocate storage for the attribute content */
    io_rank_ptr = (int*) malloc(sizeof(int));

    /* Set the attribute content */
    *io_rank_ptr = 0;

    /* Cache the attribute with io_comm */
    MPI_Attr_put(io_comm, IO_KEY, io_rank_ptr);

    /* Retrieve the I/O process rank  */
    MPI_Attr_get(io_comm, IO_KEY, &io_rank_attr, &flag);

    /* If flag == 0, something went wrong: */
    /*     there's no attribute cached.    */
    if ((flag != 0) && (my_rank == *io_rank_attr)) {
        printf("Greetings from the I/O Process!\n");
        printf("My rank is %d\n", *io_rank_attr);
    }

    MPI_Finalize();
}  /* main */
