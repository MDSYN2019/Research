/* comm_time.c  
 * Version 5:  Process 0 starts off the ring pass.
 *     Fixed erroneous calc of time_array_order in Initialize.
 *     Changed number of message sizes to 2.
 *     Removed debug output from Print_results.
 *     Corrected calculation of source process rank.
 *     Try to check sends and receives using Cprintfs in loop.
 * Time communication around a ring of processes.
 *     Guaranteed to have bugs.
 * 
 * Input: None (see notes).
 *
 * Output:  Average, minimum, and maximum time for messages 
 *    of varying sizes to be forwarded around a ring of 
 *    processes.
 *
 * Algorithm:
 *    1.  Allocate and initialize storage for messages 
 *        and communication times
 *    2.  Compute ranks of neighbors in ring.
 *    3.  Foreach message size
 *    3b.     Foreach test
 *    3a.         Start clock
 *    3c.         Send message around loop
 *    3d.         Add elapsed time to running sum
 *    3e.         Update max/min elapsed time
 *    4.  Print times.
 *
 * Functions:
 *    Initialize:  Allocate and initialize arrays
 *    Print_results:  Send results to I/O process
 *        and print.
 *
 * Notes:  
 *    1. Due to difficulties some MPI implementations 
 *       have with input, the number of tests, the max 
 *       message size, the min message size, and the size 
 *       increment are hardwired.
 *    2. We assume that the size increment evenly divides
 *       the difference max_size - min_size
 *
 * See Chap 9, pp. 192 & ff and pp. 208 & ff in PPMPI.
 */
#include <stdio.h>
#include "mpi.h"
#include "cio.h"

void Initialize(int max_size, int min_size, int size_incr, 
    int my_rank, float** x_ptr, double** times_ptr, 
    double** max_times_ptr, double** min_times_ptr, 
    int* order_ptr);

void Print_results(MPI_Comm io_comm, int my_rank, 
    int min_size, int max_size, int size_incr, 
    int time_array_order, int test_count, 
    double* times, double* max_times, double* min_times);

main(int argc, char* argv[]) {
    int         test_count = 2;    /* Number of tests  */
    int         max_size = 1000;   /* Max msg. length  */
    int         min_size = 0;      /* Min msg. length  */
    int         size_incr = 1000;  /* Increment for    */
                                   /*     msg. sizes   */
    float*      x;                 /* Message buffer   */
    double*     times;             /* Elapsed times    */
    double*     max_times;         /* Max times        */
    double*     min_times;         /* Min times        */
    int         time_array_order;  /* Size of timing   */
                                   /*     arrays.      */
    double      start;             /* Start time       */
    double      elapsed;           /* Elapsed time     */
    int         i, test, size;     /* Loop variables   */
    int         p, my_rank, source, dest; 
    MPI_Comm    io_comm;
    MPI_Status  status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);

    Cprintf(io_comm,"","Before Initialize, p = %d, my_rank = %d",
        p, my_rank);
    Initialize(max_size, min_size, size_incr, my_rank,
        &x, &times, &max_times, &min_times, 
        &time_array_order);

    source = (my_rank -  1 + p) % p;
    dest = (my_rank + 1) % p;

    /* For each message size, find average circuit time */
    /*     Loop var size = message size                 */
    /*     Loop var i = index into arrays for timings   */
    for (size = min_size, i = 0; size <= max_size;
            size = size + size_incr, i++) {
        Cprintf(io_comm,"",
            "Before if, my_rank = %d, source = %d, dest = %d",
            my_rank, source, dest);
        if (my_rank == 0) {
            times[i] =0.0;
            max_times[i] = 0.0;
            min_times[i] = 1000000.0;
            for (test = 0; test < test_count; test++) {
                start = MPI_Wtime();
        Cprintf(io_comm,"",
            "Before send x[0] = %f, size = %d, dest = %d",
            my_rank, x[0], size, dest);
                MPI_Send(x, size, MPI_FLOAT, dest, 0,
                    MPI_COMM_WORLD);
                MPI_Recv(x, size, MPI_FLOAT, source, 0,
                    MPI_COMM_WORLD, &status);
        Cprintf(io_comm,"",
            "After recv x[0] = %f, size = %d, source = %d",
            my_rank, x[0], size, source);
                elapsed = MPI_Wtime() - start;
                times[i] = times[i] + elapsed;
                if (elapsed > max_times[i])
                    max_times[i] = elapsed;
                if (elapsed < min_times[i])
                    min_times[i] = elapsed;
            }
        } else { /* my_rank != 0 */
            for (test = 0; test < test_count; test++) {
                MPI_Recv(x, size, MPI_FLOAT, source, 0,
                     MPI_COMM_WORLD, &status);
        Cprintf(io_comm,"",
            "After recv x[0] = %f, size = %d, source = %d",
            my_rank, x[0], size, source);
        Cprintf(io_comm,"",
            "Before send x[0] = %f, size = %d, dest = %d",
            my_rank, x[0], size, dest);

                MPI_Send(x, size, MPI_FLOAT, dest, 0,
                     MPI_COMM_WORLD);
            }
        }
    } /* for size . . . */

    Print_results(io_comm, my_rank, min_size, max_size, 
        size_incr, time_array_order, test_count, times, 
        max_times, min_times);

    MPI_Finalize();
} /* main */

/********************************************************/
void Initialize(int max_size, int min_size, int size_incr, 
    int my_rank, float** x_ptr, double** times_ptr, 
    double** max_times_ptr, double** min_times_ptr, 
    int* order_ptr) {
    int i;

    *x_ptr = (float *) malloc(max_size*sizeof(float));

    *order_ptr = (max_size - min_size)/size_incr + 1;
    *times_ptr = 
       (double *) malloc((*order_ptr)*sizeof(double));
    *max_times_ptr = 
       (double *) malloc((*order_ptr)*sizeof(double));
    *min_times_ptr = 
       (double *) malloc((*order_ptr)*sizeof(double));

    /* Initialize buffer -- why this? */
    for (i = 0; i < max_size; i++) 
        (*x_ptr)[i] = (float) my_rank;
} /* Initialize */


/********************************************************/
/* Send results from process 0 in MPI_COMM_WORLD to     */
/* I/O process in io_comm, which prints the results.    */
void Print_results(MPI_Comm io_comm, int my_rank, 
    int min_size, int max_size, int size_incr, 
    int time_array_order, int test_count, double* times, 
    double* max_times, double* min_times) {
    int i;
    int size;
    MPI_Status status;
    int io_process;
    int io_rank;

    Get_io_rank(io_comm, &io_process);
    MPI_Comm_rank(io_comm, &io_rank);

    if (my_rank == 0) {
        MPI_Send(times, time_array_order, MPI_DOUBLE, 
            io_rank, 0, io_comm);
        MPI_Send(max_times, time_array_order, MPI_DOUBLE, 
            io_process, 0, io_comm);
        MPI_Send(min_times, time_array_order, MPI_DOUBLE, 
            io_process, 0, io_comm);
    }
    if (io_rank == io_process) {
        MPI_Recv(times, time_array_order, MPI_DOUBLE,
            MPI_ANY_SOURCE, 0, io_comm, &status);
        MPI_Recv(max_times, time_array_order, MPI_DOUBLE, 
            MPI_ANY_SOURCE, 0, io_comm, &status);
        MPI_Recv(min_times, time_array_order, MPI_DOUBLE, 
            MPI_ANY_SOURCE, 0, io_comm, &status);

        printf("Message size (floats):  ");
        for (size = min_size; 
             size <= max_size; size += size_incr)
            printf("%10d ", size);
        printf("\n");

        printf("Avg circuit time (ms):  ");
        for (i = 0; i < time_array_order; i++)
            printf("%10f ",1000.0*times[i]/test_count);
        printf("\n");
        printf("Max circuit time (ms):  ");
        for (i = 0; i < time_array_order; i++)
            printf("%10f ",1000.0*max_times[i]);
        printf("\n");
        printf("Min circuit time (ms):  ");
        for (i = 0; i < time_array_order; i++)
            printf("%10f ",1000.0*min_times[i]);
        printf("\n\n");
        fflush(stdout);
    }
} /* Print_results */
