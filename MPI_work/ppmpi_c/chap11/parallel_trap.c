/* parallel_trap.c -- parallel trapezoidal rule with timing functions
 *
 * Input:
 *     a: left endpoint
 *     b: right endpoint
 *     n: number of trapezoids
 * Output:
 *     Integral
 *     Elapsed time in seconds (excluding I/O).
 *
 * Note:  f(x) is hardwired.
 *
 * See Chap 11, pp. 254 & ff in PPMPI.
 */
#include <stdio.h>
#include "mpi.h"
#include "cio.h"

main(int argc, char* argv[]) {
    int       p;
    int       my_rank;
    float     a;
    float     b;
    int       n;
    float     h;
    float     integral;
    float     total = 0.0;
    int       local_n;
    float     local_a;
    float     local_b;
    MPI_Comm  io_comm;
    int       i;
    double    overhead;
    double    start, finish;

    float Trap(float local_a, float local_b, int local_n,
              float h);    /* Calculate local integral  */

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);
    
    Cscanf(io_comm,"Enter a, b, and n","%f %f %d", &a, &b, &n);

    /* Estimate overhead */
    overhead = 0.0;
    for (i = 0; i < 100; i++) {
        MPI_Barrier(MPI_COMM_WORLD); 
        start = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        finish = MPI_Wtime();
        overhead = overhead + (finish - start);
    }
    overhead = overhead/100.0;

    MPI_Barrier(MPI_COMM_WORLD);
    start = MPI_Wtime();

    h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoidals */

    /* Length of each process' interval of
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;

    /* Call the serial trapezoidal function */
    integral = Trap(local_a, local_b, local_n, h);

    /* Add up the integrals calculated by each process */
    MPI_Reduce(&integral, &total, 1, MPI_FLOAT,
        MPI_SUM, 0, MPI_COMM_WORLD);

    MPI_Barrier(MPI_COMM_WORLD); 
    finish = MPI_Wtime(); 
    Cprintf(io_comm,"Our estimate is","%f",total);
    Cprintf(io_comm,"Elapsed time in seconds","%e",
        (finish - start) - overhead);

    MPI_Finalize();
}  /* main */


/********************************************************************/
float Trap(
          float  local_a   /* in */, 
          float  local_b   /* in */, 
          int    local_n   /* in */, 
          float  h         /* in */) { 

    float integral;   /* Store result in integral  */ 
    float x; 
    int i; 

    float f(float x); /* function we're integrating */

    integral = (f(local_a) + f(local_b))/2.0; 
    x = local_a; 
    for (i = 1; i <= local_n-1; i++) { 
        x = x + h; 
        integral = integral + f(x); 
    } 
    integral = integral*h; 
    return integral;
} /*  Trap  */


/********************************************************************/
float f(float x) { 
    float return_val; 
    /* Calculate f(x). */
    /* Store calculation in return_val. */ 
    return_val = x*x;
    return return_val; 
} /* f */


