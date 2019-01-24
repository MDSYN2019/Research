/* comm_create.c -- builds a communicator from the first q processes
 *     in a communicator containing p = q^2 processes.
 *
 * Input: none
 * Output: q -- program tests correct creation of new communicator
 *     by broadcasting the value 1 to its members -- all other 
 *     processes have the value 0 -- global sum computed across
 *     all the processes.
 *
 * Note:  Assumes that MPI_COMM_WORLD contains p = q^2 processes
 *
 * See Chap 7, pp. 117 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"
#include <math.h>
#include <stdlib.h>

main(int argc, char* argv[]) {
    int        p;
    int        q; /* = sqrt(p) */
    int        my_rank;
    MPI_Group  group_world;
    MPI_Group  first_row_group;
    MPI_Comm   first_row_comm;
    int*       process_ranks;
    int        proc;
    int        test = 0;
    int        sum;
    int        my_rank_in_first_row;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    q = (int) sqrt((double) p);

    /* Make a list of the processes in the new
     * communicator */
    process_ranks = (int*) malloc(q*sizeof(int));
    for (proc = 0; proc < q; proc++)
        process_ranks[proc] = proc;

    /* Get the group underlying MPI_COMM_WORLD */
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);

    /* Create the new group */
    MPI_Group_incl(group_world, q, process_ranks,
        &first_row_group);

    /* Create the new communicator */
    MPI_Comm_create(MPI_COMM_WORLD, first_row_group,
        &first_row_comm);

    /* Now check whether we can do collective ops in first_row_comm */
    if (my_rank < q) {
        MPI_Comm_rank(first_row_comm, &my_rank_in_first_row);
        if (my_rank_in_first_row == 0) test = 1;
        MPI_Bcast(&test, 1, MPI_INT, 0, first_row_comm);
    }
    MPI_Reduce(&test, &sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("q = %d, sum = %d\n", q, sum);
    }

    MPI_Finalize();
}  /* main */
