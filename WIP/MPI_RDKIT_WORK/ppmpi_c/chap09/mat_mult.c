/* mat_mult.c
 * Multiply a sequence of 2x2 matrices -- 1 factor from 
 *     each process.  Erroneous.
 *
 * Input: none
 *
 * Output: product of a sequence of 2x2 matrices
 *
 * Algorithm
 *    1. Generate local matrix
 *    2. Send local matrix to process 0
 *    3  if (my_rank == 0)
 *    3b.    for each process, receive matrix and 
 *              multiply by product
 *    3c.    print product
 *
 * Notes: 
 *    1. The matrices are stored as linear arrays.  The 
 *       correspondence is row major: Matrix[i][j] <-> 
 *       Array[2*i + j]
 *    2. Local matrices have the form
 *            [my_rank    my_rank+1]
 *            [my_rank+2  my_rank  ]
 *
 * See Chap 9, pp. 188 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"

#define MATRIX_ORDER 2
#define ARRAY_ORDER 4

#define Entry(mat,i,j) (mat[MATRIX_ORDER*(i) + (j)])

void Initialize(float my_matrix[], int my_rank);
void Mult(float product[], float factor[]);
void Print_matrix(char* title, float matrix[]);

main(int argc, char** argv) {
    float        my_matrix[ARRAY_ORDER];
    float        temp[ARRAY_ORDER];
    float        product[ARRAY_ORDER] = {1, 0, 0, 1};  
                     /* product is the identity matrix */
    int          p;
    int          my_rank;
    MPI_Status   status;
    int          i;
char title[100];

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    Initialize(my_matrix, my_rank);

    MPI_Send(my_matrix, ARRAY_ORDER, MPI_FLOAT, 0, 0,
        MPI_COMM_WORLD);

    if (my_rank == 0) {
        for (i = 0; i < p; i++) {
            MPI_Recv(temp, ARRAY_ORDER, MPI_FLOAT, 
                MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, 
                &status); 
            Mult(product, temp);
        }
        Print_matrix("The product is", product); 
    }

    MPI_Finalize();
} /* main */


/**********************************************************************/
void Initialize(
         float  my_matrix[]  /* out */,
         int    my_rank      /* in  */) {

    my_matrix[0] = my_matrix[3] = (float) my_rank;
    my_matrix[1] = (float) (my_rank + 1);
    my_matrix[2] = (float) (my_rank + 2);
}  /* Initialize */


/**********************************************************************/
void Mult(
         float  product[]  /* in/out */, 
         float  factor[]   /* in     */) {
    int    i, j, k;
    float  temp[ARRAY_ORDER];

    for (i = 0; i < MATRIX_ORDER; i++)
        for (j = 0; j < MATRIX_ORDER; j++) {
            Entry(temp,i,j) = 0.0;
            for (k = 0; k < MATRIX_ORDER; k++)
                Entry(temp,i,j) = Entry(temp,i,j) +
                    Entry(product,i,k)*Entry(factor,k,j);
        }

    for (i = 0; i < ARRAY_ORDER; i++)
        product[i] = temp[i];

}  /* Mult */


/**********************************************************************/
void Print_matrix(
         char*  title     /* in */, 
         float  matrix[]  /* in */) {

    int i, j;

    printf("%s\n", title);
    for (i = 0; i < MATRIX_ORDER; i++) {
        for (j = 0; j < MATRIX_ORDER; j++)
            printf("%f ", Entry(matrix,i,j));
        printf("\n");
    }
}  /* Print_matrix */
