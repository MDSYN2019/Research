/* parallel_dot1.c -- Computes a parallel dot product.  Uses MPI_Allreduce.
 *
 * Input: 
 *     n: order of vectors
 *     x, y:  the vectors
 *
 * Output:
 *     the dot product of x and y as computed by each process.
 *
 * Note:  Arrays containing vectors are statically allocated.  Assumes that
 *     n, the global order of the vectors, is evenly divisible by p, the
 *     number of processes.
 *
 * See Chap 5, pp. 76 & ff in PPMPI.
 */
#include <stdio.h>
#include "mpi.h"

#define MAX_LOCAL_ORDER 100

main(int argc, char* argv[]) {
    float  local_x[MAX_LOCAL_ORDER];
    float  local_y[MAX_LOCAL_ORDER];
    int    n;
    int    n_bar;  /* = n/p */
    float  dot;
    int    p;
    int    my_rank;

    void Read_vector(char* prompt, float local_v[], int n_bar, int p,
             int my_rank);
    float Parallel_dot(float local_x[], float local_y[], int n_bar);
    void Print_results(float dot, int my_rank, int p);
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        printf("Enter the order of the vectors\n");
        scanf("%d", &n);
    }
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    n_bar = n/p;

    Read_vector("the first vector", local_x, n_bar, p, my_rank);
    Read_vector("the second vector", local_y, n_bar, p, my_rank);

    dot = Parallel_dot(local_x, local_y, n_bar);

    Print_results(dot, my_rank, p);

    MPI_Finalize();
}  /* main */
   

/*****************************************************************/
void Read_vector(
         char*  prompt     /* in  */,
         float  local_v[]  /* out */,
         int    n_bar      /* in  */,
         int    p          /* in  */,
         int    my_rank    /* in  */) {
    int i, q;
    float temp[MAX_LOCAL_ORDER];
    MPI_Status status;

    if (my_rank == 0) {
        printf("Enter %s\n", prompt);
        for (i = 0; i < n_bar; i++)
            scanf("%f", &local_v[i]);
        for (q = 1; q < p; q++) {
            for (i = 0; i < n_bar; i++)
                scanf("%f", &temp[i]);
            MPI_Send(temp, n_bar, MPI_FLOAT, q, 0, MPI_COMM_WORLD);
        }
    } else {
        MPI_Recv(local_v, n_bar, MPI_FLOAT, 0, 0, MPI_COMM_WORLD,
            &status);
    }
}  /* Read_vector */


/*****************************************************************/
float Serial_dot(
          float  x[]  /* in */, 
          float  y[]  /* in */, 
          int    n    /* in */) {

    int    i; 
    float  sum = 0.0;

    for (i = 0; i < n; i++)
        sum = sum + x[i]*y[i];
    return sum;
} /* Serial_dot */


/*****************************************************************/
float Parallel_dot(
          float  local_x[]  /* in */,
          float  local_y[]  /* in */,
          int  n_bar        /* in */) {

    float  local_dot;
    float  dot = 0.0;
    float  Serial_dot(float x[], float y[], int m);

    local_dot = Serial_dot(local_x, local_y, n_bar);
    MPI_Allreduce(&local_dot, &dot, 1, MPI_FLOAT,
        MPI_SUM, MPI_COMM_WORLD);
    return dot;
} /* Parallel_dot */


/*****************************************************************/
void Print_results(
         float  dot      /* in */, 
         int    my_rank  /* in */, 
         int    p        /* in */) {
    int         q;
    float       temp;
    MPI_Status  status;

    if (my_rank == 0) {
        printf("dot = \n");
        printf("Process 0 > %f\n", dot);
        for (q = 1; q < p; q++) {
            MPI_Recv(&temp, 1, MPI_FLOAT, q, 0, MPI_COMM_WORLD,
                &status);
            printf("Process %d > %f\n", q, temp);
        }
    } else {
        MPI_Send(&dot, 1, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
    }

}  /* Print_results */
