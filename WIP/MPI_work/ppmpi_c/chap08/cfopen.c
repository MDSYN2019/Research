/* cfopen.c -- open a file and write ints received from each process to it
 *
 * Input: none
 * Output: process ranks in file "testfile"
 *
 * Note: Link with cio.o
 *
 * See Chap 8, p. 157 in PPMPI
 */
#include <stdio.h>
#include "mpi.h"
#include "cio.h"

#define MAX 32


FILE* Cfopen(
        char*     filename  /* in */,
        char*     mode      /* in */,
        MPI_Comm  io_comm   /* in */);

main(int argc, char* argv[]) {
    int       p;
    int       my_rank;
    int       io_rank;
    int       list[MAX];
    MPI_Comm  io_comm;
    int       i;
    FILE*     fp;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);

    fp = Cfopen("testfile", "w", io_comm);

    Get_io_rank(io_comm, &io_rank);
    MPI_Gather(&my_rank, 1, MPI_INT, list, 1, MPI_INT,
        io_rank, io_comm);
    if (my_rank == io_rank) {
        for (i = 0; i < p; i++)
            fprintf(fp,"%d ", list[i]);
        fprintf(fp,"\n");
        fclose(fp);
    }

    MPI_Finalize();
}  /* main */


/************************************************************/
FILE* Cfopen(
        char*     filename  /* in */,
        char*     mode      /* in */,
        MPI_Comm  io_comm   /* in */) {

    int    root;
    int    my_io_rank;
    FILE*  fp;

    Get_io_rank(io_comm, &root);

    MPI_Comm_rank(io_comm, &my_io_rank);

    if (my_io_rank == root) {
        fp = fopen(filename, mode);
        return fp;
    } else {
        return NULL;
    }
}  /* Cfopen */
