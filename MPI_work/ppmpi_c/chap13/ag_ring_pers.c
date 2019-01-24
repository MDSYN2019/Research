/* ag_ring_pers.c -- ring allgather using persistent communication requests
 *
 * Input: series of blocksizes for allgather, 0 to stop.
 * Output: Contents of gathered array on each process -- list of
 *     process ranks, each rank appearing in a block of size blocksize.
 *
 * Note:  array sizes are hardwired in MAX, LOCAL_MAX, and MAX_BYTES.
 *
 * See Chap 13, pp. 301 & ff, in PPMPI.
 */

#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "cio.h"

#define MAX 128
#define LOCAL_MAX 128
#define MAX_BYTES MAX*sizeof(float)

void Allgather_ring(float x[], int blocksize,
         float y[], MPI_Comm ring_comm);

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
        Allgather_ring(x, blocksize, y, MPI_COMM_WORLD);
        Print_arrays(io_comm, "Gathered_arrays", y, blocksize);
        /* Enter 0 to stop. */
        Cscanf(io_comm,"Enter the local array size",
            "%d", &blocksize);
    }

    MPI_Finalize();
}  /* main */


/********************************************************************/
void Allgather_ring(
         float     x[]        /* in  */, 
         int       blocksize  /* in  */, 
         float     y[]        /* out */, 
         MPI_Comm  ring_comm  /* in  */) {

    int         i, p, my_rank;
    int         successor, predecessor;
    int         send_offset, recv_offset;
    MPI_Status  status;
    MPI_Request send_request;
    MPI_Request recv_request;
    float       send_buf[MAX];
    float       recv_buf[MAX];
    int         position;
    
    MPI_Comm_size(ring_comm, &p);
    MPI_Comm_rank(ring_comm, &my_rank);

    /* Copy x into correct location in y */
    for (i = 0; i < blocksize; i++)
        y[i + my_rank*blocksize] = x[i];

    successor = (my_rank + 1) % p;
    predecessor = (my_rank - 1 + p) % p;

    MPI_Send_init(send_buf, blocksize*sizeof(float), 
        MPI_PACKED, successor, 0, ring_comm, 
        &send_request);
    MPI_Recv_init(recv_buf, blocksize*sizeof(float), 
        MPI_PACKED, predecessor, 0, ring_comm, 
        &recv_request );

    send_offset = my_rank*blocksize;
    for (i = 0; i < p - 1; i++) {
        position = 0;
        MPI_Pack(y+send_offset, blocksize, MPI_FLOAT, 
            send_buf, MAX_BYTES, &position, ring_comm);
        MPI_Start(&send_request);
        MPI_Start(&recv_request);
        recv_offset = send_offset = 
           ((my_rank - i - 1 + p) % p)*blocksize;
        position = 0;
        MPI_Wait(&send_request, &status);
        MPI_Wait(&recv_request, &status);
        MPI_Unpack(recv_buf, MAX_BYTES, &position, 
            y+recv_offset, blocksize, MPI_FLOAT, ring_comm);
    }
    MPI_Request_free(&send_request);
    MPI_Request_free(&recv_request);
} /* Allgather_ring */


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
