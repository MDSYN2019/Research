/* ub.c -- build a derived type that uses MPI_UB
 *
 * Input: none
 * Output: result of MPI_Scatter using derived type
 *
 * Note:  Scattered array has dimension 12.  So the number of processes
 *     should divide 12.
 *
 * See Chap 8, pp. 164 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"

#define n 12

main(int argc, char* argv[]) {
    int p;
    int my_rank;
    int           block_lengths[2];
    MPI_Aint      displacements[2];
    MPI_Datatype  types[2];
    MPI_Datatype  vector_mpi_t;
    MPI_Datatype  ub_mpi_t;
    int count;
    int padded_size;
    int stride;
    float input_data[n]; 
    float local_data[n];
    int i;
    MPI_Comm io_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    padded_size = n;

    /* First create vector_mpi_t                    */
    count = padded_size/p;
    stride = p;
    MPI_Type_vector(count, 1, stride, MPI_FLOAT,
        &vector_mpi_t);

    /* The first type is vector_mpi_t               */
    types[0] = vector_mpi_t;

    /* The second type is MPI_UB                    */
    types[1] = MPI_UB;

    /* vector_mpi_t starts at displacement 0        */
    displacements[0] = 0;

    /* MPI_UB starts at displacement sizeof(float)  */
    displacements[1] = sizeof(float);

    /* The derived type will have 1 element of type */
    /*     vector_mpi_t and 1 element of type       */
    /*     MPI_UB                                   */
    block_lengths[0] = block_lengths[1] = 1;

    /* Now create the full type */
    MPI_Type_struct(2, block_lengths, displacements,
        types, &ub_mpi_t);
    MPI_Type_commit(&ub_mpi_t);

    if (my_rank == 0)
        for (i = 0; i < padded_size; i++)
            input_data[i] = (float) i;

    MPI_Scatter(input_data, 1, ub_mpi_t,
        local_data, padded_size/p, MPI_FLOAT,
        0, io_comm);

    printf("Process %d > ", my_rank);
    for (i = 0; i < padded_size/p; i++)
        printf("%4.1f ", local_data[i]);
    printf("\n");

    MPI_Finalize();
}  /* main */
