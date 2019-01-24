/* parallel_mat_vect.c -- computes a parallel matrix-vector product.  Matrix
 *     is distributed by block rows.  Vectors are distributed by blocks.
 *
 * Input:
 *     m, n: order of matrix
 *     A, x: the matrix and the vector to be multiplied
 *
 * Output:
 *     y: the product vector
 *
 * Notes:  
 *     1.  Local storage for A, x, and y is statically allocated.
 *     2.  Number of processes (p) should evenly divide both m and n.
 *
 * See Chap 5, p. 78 & ff in PPMPI.
 */

#include <stdio.h>
#include "mpi.h"

#define MAX_ORDER 100

typedef float LOCAL_MATRIX_T[MAX_ORDER][MAX_ORDER];

main(int argc, char* argv[]) {
    int             my_rank;
    int             p;
    LOCAL_MATRIX_T  local_A; 
    float           global_x[MAX_ORDER];
    float           local_x[MAX_ORDER];
    float           local_y[MAX_ORDER];
    int             m, n;
    int             local_m, local_n;

    void Read_matrix(char* prompt, LOCAL_MATRIX_T local_A, int local_m, int n,
             int my_rank, int p);
    void Read_vector(char* prompt, float local_x[], int local_n, int my_rank,
             int p);
    void Parallel_matrix_vector_prod( LOCAL_MATRIX_T local_A, int m, 
             int n, float local_x[], float global_x[], float local_y[],
             int local_m, int local_n);
    void Print_matrix(char* title, LOCAL_MATRIX_T local_A, int local_m,
             int n, int my_rank, int p);
    void Print_vector(char* title, float local_y[], int local_m, int my_rank,
             int p);

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    if (my_rank == 0) {
        printf("Enter the order of the matrix (m x n)\n");
        scanf("%d %d", &m, &n);
    }
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    local_m = m/p;
    local_n = n/p;

    Read_matrix("Enter the matrix", local_A, local_m, n, my_rank, p);
    Print_matrix("We read", local_A, local_m, n, my_rank, p);

    Read_vector("Enter the vector", local_x, local_n, my_rank, p);
    Print_vector("We read", local_x, local_n, my_rank, p);

    Parallel_matrix_vector_prod(local_A, m, n, local_x, global_x, 
        local_y, local_m, local_n);
    Print_vector("The product is", local_y, local_m, my_rank, p);

    MPI_Finalize();

}  /* main */


/**********************************************************************/
void Read_matrix(
         char*           prompt   /* in  */, 
         LOCAL_MATRIX_T  local_A  /* out */, 
         int             local_m  /* in  */, 
         int             n        /* in  */,
         int             my_rank  /* in  */, 
         int             p        /* in  */) {

    int             i, j;
    LOCAL_MATRIX_T  temp;

    /* Fill dummy entries in temp with zeroes */
    for (i = 0; i < p*local_m; i++)
        for (j = n; j < MAX_ORDER; j++)
            temp[i][j] = 0.0;

    if (my_rank == 0) {
        printf("%s\n", prompt);
        for (i = 0; i < p*local_m; i++) 
            for (j = 0; j < n; j++)
                scanf("%f",&temp[i][j]);
    }
    MPI_Scatter(temp, local_m*MAX_ORDER, MPI_FLOAT, local_A,
        local_m*MAX_ORDER, MPI_FLOAT, 0, MPI_COMM_WORLD);

}  /* Read_matrix */


/**********************************************************************/
void Read_vector(
         char*  prompt     /* in  */,
         float  local_x[]  /* out */, 
         int    local_n    /* in  */, 
         int    my_rank    /* in  */,
         int    p          /* in  */) {

    int   i;
    float temp[MAX_ORDER];

    if (my_rank == 0) {
        printf("%s\n", prompt);
        for (i = 0; i < p*local_n; i++) 
            scanf("%f", &temp[i]);
    }
    MPI_Scatter(temp, local_n, MPI_FLOAT, local_x, local_n, MPI_FLOAT,
        0, MPI_COMM_WORLD);

}  /* Read_vector */


/**********************************************************************/
/* All arrays are allocated in calling program */
/* Note that argument m is unused              */
void Parallel_matrix_vector_prod(
         LOCAL_MATRIX_T  local_A     /* in  */,
         int             m           /* in  */,
         int             n           /* in  */,
         float           local_x[]   /* in  */,
         float           global_x[]  /* in  */,
         float           local_y[]   /* out */,
         int             local_m     /* in  */,
         int             local_n     /* in  */) {

    /* local_m = m/p, local_n = n/p */

    int i, j;

    MPI_Allgather(local_x, local_n, MPI_FLOAT,
                   global_x, local_n, MPI_FLOAT,
                   MPI_COMM_WORLD);
    for (i = 0; i < local_m; i++) {
        local_y[i] = 0.0;
        for (j = 0; j < n; j++)
            local_y[i] = local_y[i] +
                         local_A[i][j]*global_x[j];
    }
}  /* Parallel_matrix_vector_prod */


/**********************************************************************/
void Print_matrix(
         char*           title      /* in */, 
         LOCAL_MATRIX_T  local_A    /* in */, 
         int             local_m    /* in */, 
         int             n          /* in */,
         int             my_rank    /* in */,
         int             p          /* in */) {

    int   i, j;
    float temp[MAX_ORDER][MAX_ORDER];

    MPI_Gather(local_A, local_m*MAX_ORDER, MPI_FLOAT, temp, 
         local_m*MAX_ORDER, MPI_FLOAT, 0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("%s\n", title);
        for (i = 0; i < p*local_m; i++) {
            for (j = 0; j < n; j++)
                printf("%4.1f ", temp[i][j]);
            printf("\n");
        }
    } 
}  /* Print_matrix */


/**********************************************************************/
void Print_vector(
         char*  title      /* in */, 
         float  local_y[]  /* in */, 
         int    local_m    /* in */, 
         int    my_rank    /* in */,
         int    p          /* in */) {

    int   i;
    float temp[MAX_ORDER];

    MPI_Gather(local_y, local_m, MPI_FLOAT, temp, local_m, MPI_FLOAT,
        0, MPI_COMM_WORLD);

    if (my_rank == 0) {
        printf("%s\n", title);
        for (i = 0; i < p*local_m; i++)
            printf("%4.1f ", temp[i]);
        printf("\n");
    } 
}  /* Print_vector */
