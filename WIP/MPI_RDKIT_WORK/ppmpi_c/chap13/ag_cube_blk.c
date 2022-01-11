/* ag_cube_blk.c -- hypercube allgather using blocking sends and receives
 *
 * Input: series of blocksizes for allgather, 0 to stop.
 * Output: Contents of gathered array on each process -- list of
 *     process ranks, each rank appearing in a block of size blocksize.
 *
 * Note:  array sizes are hardwired in MAX and LOCAL_MAX.
 *
 * See Chap 13, pp. 280 & ff, in PPMPI.
 */

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "cio.h"

#define MAX 128
#define LOCAL_MAX 128

void  Allgather_cube(
         float    x[]        /* in  */, 
         int      blocksize  /* in  */, 
         float    y[]        /* out */, 
         MPI_Comm comm       /* in  */);

void Print_arrays(MPI_Comm io_comm, char* title, 
         float y[], int blocksize);


/********************************************************************/
main(int argc, char* argv[]) {
    int       p;
    int       my_rank;
    float     x[LOCAL_MAX];
    float     y[MAX];
    int       blocksize;
    MPI_Comm  io_comm;
    int       i;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);

    Cscanf(io_comm,"Enter the local array size","%d", &blocksize);

    while(blocksize > 0) {
        for (i = 0; i < blocksize; i++)
            x[i] = (float) my_rank;
        Allgather_cube(x, blocksize, y, MPI_COMM_WORLD);
        Print_arrays(io_comm, "Gathered_arrays", y, blocksize);
        /* Enter 0 to stop. */
        Cscanf(io_comm,"Enter the local array size",
            "%d", &blocksize);
    }

    MPI_Finalize();
}  /* main */


/********************************************************************/
void Print_arrays(
         MPI_Comm  io_comm    /* in */,
         char*     title      /* in */, 
         float     y[]        /* in */, 
         int       blocksize  /* in */) {

    char item[16];
    char list[4*MAX];
    int  i;
    int  p;

    MPI_Comm_size(io_comm, &p);

    list[0] = '\0';
    for (i = 0; i < blocksize*p; i++) {
        sprintf(item, "%3.1f ", y[i]);
        strcat(list, item);
    }
    Cprintf(io_comm, title, "%s", list);
}  /* Print_arrays */


/********************************************************************/
int log_base2(int p) {
/* Just counts number of bits to right of most significant
 * bit.  So for p not a power of 2, it returns the floor
 * of log_2(p).
 */
    int return_val = 0;
    unsigned q;

    q = (unsigned) p;
    while(q != 1) {
        q = q >> 1;
        return_val++;
    }
    return return_val;
}  /* log_base2 */


/********************************************************************/
void  Allgather_cube(
         float    x[]        /* in  */, 
         int      blocksize  /* in  */, 
         float    y[]        /* out */, 
         MPI_Comm comm       /* in  */) {

    int           i, d, p, my_rank;
    unsigned      eor_bit;
    unsigned      and_bits;
    int           stage, partner;
    MPI_Datatype  hole_type;
    int           send_offset, recv_offset;
    MPI_Status    status;

    int log_base2(int p);

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    /* Copy x into correct location in y */
    for (i = 0; i < blocksize; i++)
        y[i + my_rank*blocksize] = x[i];

    /* Set up */
    d = log_base2(p);
    eor_bit = 1 << (d-1);
    and_bits = (1 << d) - 1;

    for (stage = 0; stage < d; stage++) {
        partner = my_rank ^ eor_bit;
        send_offset = (my_rank & and_bits)*blocksize;
        recv_offset = (partner & and_bits)*blocksize;

        MPI_Type_vector(1 << stage, blocksize, 
            (1 << (d-stage))*blocksize, MPI_FLOAT,
            &hole_type);
        MPI_Type_commit(&hole_type);

        MPI_Send(y + send_offset, 1, hole_type,
            partner, 0, comm);
        MPI_Recv(y + recv_offset, 1, hole_type,
            partner, 0, comm, &status);

        MPI_Type_free(&hole_type);  /* Free type so we  */
                 /* can build new type during next pass */
        eor_bit = eor_bit >> 1; 
        and_bits = and_bits >> 1;
    }
} /* Allgather_cube */
