/* top_fcns.c -- test basic topology functions
 *
 * Input: none
 * Output: results of calls to various functions testing topology
 *     creation
 *
 * Algorithm:
 *     1.  Build a 2-dimensional Cartesian communicator from
 *         MPI_Comm_world
 *     2.  Print topology information for each process
 *     3.  Use MPI_Cart_sub to build a communicator for each
 *         row of the Cartesian communicator
 *     4.  Carry out a broadcast across each row communicator
 *     5.  Print results of broadcast
 *     6.  Use MPI_Cart_sub to build a communicator for each
 *         column of the Cartesian communicator
 *     7.  Carry out a broadcast across each column communicator
 *     8.  Print results of broadcast
 *
 * Note: Assumes the number of processes, p, is a perfect square
 *
 * See Chap 7, pp. 121 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"
#include <math.h>

main(int argc, char* argv[]) {
    int       p;
    int       my_rank;
    int       q;
    MPI_Comm  grid_comm;
    int       dim_sizes[2];
    int       wrap_around[2];
    int       reorder = 1;
    int       coordinates[2];
    int       my_grid_rank;
    int       grid_rank;
    int       free_coords[2];
    MPI_Comm  row_comm;
    MPI_Comm  col_comm;
    int       row_test;
    int       col_test;


    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    q = (int) sqrt((double) p);

    dim_sizes[0] = dim_sizes[1] = q;
    wrap_around[0] = wrap_around[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim_sizes,
        wrap_around, reorder, &grid_comm);

    MPI_Comm_rank(grid_comm, &my_grid_rank);
    MPI_Cart_coords(grid_comm, my_grid_rank, 2,
        coordinates);

    MPI_Cart_rank(grid_comm, coordinates, &grid_rank);

    printf("Process %d > my_grid_rank = %d, coords = (%d,%d), grid_rank = %d\n",
        my_rank, my_grid_rank, coordinates[0], coordinates[1], grid_rank);

    free_coords[0] = 0;
    free_coords[1] = 1;
    MPI_Cart_sub(grid_comm, free_coords, &row_comm);
    if (coordinates[1] == 0)
        row_test = coordinates[0];
    else
        row_test = -1;
    MPI_Bcast(&row_test, 1, MPI_INT, 0, row_comm);
    printf("Process %d > coords = (%d,%d), row_test = %d\n",
         my_rank, coordinates[0], coordinates[1], row_test);

    free_coords[0] = 1;
    free_coords[1] = 0;
    MPI_Cart_sub(grid_comm, free_coords, &col_comm);
    if (coordinates[0] == 0)
        col_test = coordinates[1];
    else
        col_test = -1;
    MPI_Bcast(&col_test, 1, MPI_INT, 0, col_comm);
    printf("Process %d > coords = (%d,%d), col_test = %d\n",
         my_rank, coordinates[0], coordinates[1], col_test);

    MPI_Finalize();
}  /* main */
