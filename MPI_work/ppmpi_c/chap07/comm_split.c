/* comm_split.c -- build a collection of q communicators using MPI_Comm_split
 *
 * Input: none
 * Output:  Results of doing a broadcast across each of the q communicators.
 *     
 * Note:  Assumes the number of processes, p = q^2
 *
 * See Chap. 7, pp. 120 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"
#include <math.h>

main(int argc, char* argv[]) {
    int       p;
    int       my_rank;
    MPI_Comm  my_row_comm;
    int       my_row;
    int       q;
    int       test;
    int       my_rank_in_row;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    q = (int) sqrt((double) p);

    /* my_rank is rank in MPI_COMM_WORLD.
     * q*q = p */
    my_row = my_rank/q;
    MPI_Comm_split(MPI_COMM_WORLD, my_row, my_rank,
        &my_row_comm);

    /* Test the new communicators */
    MPI_Comm_rank(my_row_comm, &my_rank_in_row);
    if (my_rank_in_row == 0)
        test = my_row;
    else
        test = 0;

    MPI_Bcast(&test, 1, MPI_INT, 0, my_row_comm);

    printf("Process %d > my_row = %d, my_rank_in_row = %d, test = %d\n", 
        my_rank, my_row, my_rank_in_row, test);

    MPI_Finalize();
}  /* main */
