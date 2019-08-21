/* ping_pong.c -- two-process ping-pong -- send from 0 to 1 and send back
 *     from 1 to 0
 *
 * Input: none
 * Output: time elapsed for each ping-pong
 *
 * Notes:
 *     1.  Size of message is MAX_ORDER floats.
 *     2.  Number of ping-pongs is MAX.
 *
 * See Chap 12, pp. 267 & ff. in PPMPI.
 */

#include <stdio.h>
#include "mpi.h"

#define MAX_ORDER 100
#define MAX 2

main(int argc, char* argv[]) {
    int    p;
    int    my_rank;
    int    test;
    int    min_size = 0;
    int    max_size = 16;
    int    incr = 8;
    float  x[MAX_ORDER];
    int    size;
    int    pass;
    MPI_Status  status;
    int         i;
    double      wtime_overhead;
    double      start, finish;
    double      raw_time;
    MPI_Comm    comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &comm);

    wtime_overhead = 0.0;
    for (i = 0; i < 100; i++) {
        start = MPI_Wtime();
        finish = MPI_Wtime();
        wtime_overhead = wtime_overhead + (start - finish);
    }
    wtime_overhead = wtime_overhead/100.0;

    if (my_rank == 0) {
        for (test = 0, size = min_size;
                size <= max_size; size = size + incr, test++) {
            for (pass = 0; pass < MAX; pass++) {
                MPI_Barrier(comm);
                start = MPI_Wtime();
                MPI_Send(x, size, MPI_FLOAT, 1, 0, comm);
                MPI_Recv(x, size, MPI_FLOAT, 1, 0, comm,
                    &status);
                finish = MPI_Wtime();
                raw_time = finish - start - wtime_overhead;
                printf("%d %f\n", size, raw_time);
            }
        }
    } else { /* my_rank == 1 */
        for (test = 0, size = min_size; size <= max_size;
                 size = size + incr, test++) {
            for (pass = 0; pass < MAX; pass++) {
		MPI_Barrier(comm); 
                MPI_Recv(x, size, MPI_FLOAT, 0, 0, comm,
		    &status); 
                MPI_Send(x, size, MPI_FLOAT, 0, 0, comm);
	    } 
        } 
    }


    MPI_Finalize();
}  /* main */
