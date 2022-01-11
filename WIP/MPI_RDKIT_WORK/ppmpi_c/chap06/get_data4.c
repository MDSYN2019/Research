/* get_data4.c -- Parallel Trapezoidal Rule.  Uses MPI_Pack/Unpack in
 *     distribution of input data.
 *
 * Input: 
 *    a, b: limits of integration.
 *    n: number of trapezoids.
 * Output:  Estimate of the integral from a to b of f(x) 
 *    using the trapezoidal rule and n trapezoids.
 *
 * Notes:  
 *    1.  f(x) is hardwired.
 *    2.  the number of processes (p) should evenly divide
 *        the number of trapezoids (n).
 *
 * See Chap 6., pp. 100 & ff in PPMPI
 */
#include <stdio.h>

/* We'll be using MPI routines, definitions, etc. */
#include "mpi.h"

main(int argc, char** argv) {
    int         my_rank;   /* My process rank           */
    int         p;         /* The number of processes   */
    float       a;         /* Left endpoint             */
    float       b;         /* Right endpoint            */
    int         n;         /* Number of trapezoids      */
    float       h;         /* Trapezoid base length     */
    float       local_a;   /* Left endpoint my process  */
    float       local_b;   /* Right endpoint my process */
    int         local_n;   /* Number of trapezoids for  */
                           /* my calculation            */
    float       integral;  /* Integral over my interval */
    float       total;     /* Total integral            */
    int         source;    /* Process sending integral  */
    int         dest = 0;  /* All messages go to 0      */
    int         tag = 0;
    MPI_Status  status;

    void Get_data4(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank);
    float Trap(float local_a, float local_b, int local_n,
              float h);    /* Calculate local integral  */

    /* Let the system do what it needs to start up MPI */
    MPI_Init(&argc, &argv);

    /* Get my process rank */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    /* Find out how many processes are being used */
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    Get_data4(&a, &b, &n, my_rank);

    h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoids */

    /* Length of each process' interval of 
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    integral = Trap(local_a, local_b, local_n, h);

    /* Add up the integrals calculated by each process */
    MPI_Reduce(&integral, &total, 1, MPI_FLOAT,
        MPI_SUM, 0, MPI_COMM_WORLD);

    /* Print the result */
    if (my_rank == 0) {
        printf("With n = %d trapezoids, our estimate\n", 
            n);
        printf("of the integral from %f to %f = %f\n", 
            a, b, total); 
    }

    /* Shut down MPI */
    MPI_Finalize();
} /*  main  */


/********************************************************************/
void Get_data4(
         float*  a_ptr    /* out */, 
         float*  b_ptr    /* out */, 
         int*    n_ptr    /* out */,
         int     my_rank  /* in  */) {

    char  buffer[100];  /* Store data in buffer        */
    int   position;     /* Keep track of where data is */    
                        /*     in the buffer           */

    if (my_rank == 0){
        printf("Enter a, b, and n\n");
        scanf("%f %f %d", a_ptr, b_ptr, n_ptr);

        /* Now pack the data into buffer.  Position = 0 */
        /* says start at beginning of buffer.           */
        position = 0;  

        /* Position is in/out */
        MPI_Pack(a_ptr, 1, MPI_FLOAT, buffer, 100,
            &position, MPI_COMM_WORLD);
        /* Position has been incremented: it now refer- */
        /* ences the first free location in buffer.     */

        MPI_Pack(b_ptr, 1, MPI_FLOAT, buffer, 100,
            &position, MPI_COMM_WORLD);
        /* Position has been incremented again. */

        MPI_Pack(n_ptr, 1, MPI_INT, buffer, 100,
            &position, MPI_COMM_WORLD);
        /* Position has been incremented again. */

        /* Now broadcast contents of buffer */
        MPI_Bcast(buffer, 100, MPI_PACKED, 0,
            MPI_COMM_WORLD);
    } else {
        MPI_Bcast(buffer, 100, MPI_PACKED, 0,
            MPI_COMM_WORLD);

        /* Now unpack the contents of buffer */
        position = 0;
        MPI_Unpack(buffer, 100, &position, a_ptr, 1,
            MPI_FLOAT, MPI_COMM_WORLD);
        /* Once again position has been incremented: */
        /* it now references the beginning of b.     */

        MPI_Unpack(buffer, 100, &position, b_ptr, 1,
            MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(buffer, 100, &position, n_ptr, 1,
            MPI_INT, MPI_COMM_WORLD);
    }
} /* Get_data4 */


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


