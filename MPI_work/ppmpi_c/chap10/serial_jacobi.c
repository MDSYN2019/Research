/* serial_jacobi.c -- serial version of Jacobi's method for solving
 *     the linear system Ax = b.
 *
 * Input:
 *     n:  order of system
 *     tol:  convergence tolerance
 *     max_iter:  maximum number of iterations
 *     A:  coefficient matrix
 *     b:  right-hand side of system
 *
 * Output:
 *     x:  the solution if the method converges
 *     max_iter:  if the method fails to converge
 *
 * Notes:  
 *     1.  A should be strongly diagonally dominant in
 *         order to insure convergence.
 *     2.  A, x, and b are statically allocated.
 *
 * See Chap 10, pp. 218 & ff in PPMPI.
 */
#include <stdio.h>
#include <math.h>

#define MAX_DIM 12

typedef float MATRIX_T[MAX_DIM][MAX_DIM];

int Jacobi(MATRIX_T A, float x[], float b[], int n,
        float tol, int max_iter);
void Read_matrix(char* prompt, MATRIX_T A, int n);
void Read_vector(char* prompt, float x[], int n);
void Print_vector(char* title, float x[], int n);
void Print_matrix(char* title, MATRIX_T	A, int n);

main(int argc, char* argv) {
    MATRIX_T  A;
    float     x[MAX_DIM];
    float     b[MAX_DIM];
    int       n;
    float     tol;
    int       max_iter;
    int       converged;

    printf("Enter n, tolerance, and max number of iterations\n");
    scanf("%d %f %d", &n, &tol, &max_iter);
    Read_matrix("Enter the matrix", A, n);
    Read_vector("Enter the right-hand side", b, n);

    converged = Jacobi(A, x, b, n, tol, max_iter);

    if (converged)
        Print_vector("The solution is", x, n);
    else
        printf("Failed to converge in %d iterations\n", max_iter);
}  /* main */


/*********************************************************************/
/* Return 1 if iteration converged, 0 otherwise */
/* MATRIX_T is just a 2-dimensional array       */
int Jacobi(
        MATRIX_T  A         /* in  */, 
        float     x[]       /* out */, 
        float     b[]       /* in  */, 
        int       n         /* in  */, 
        float     tol       /* in  */, 
        int       max_iter  /* in  */) {
    int    i, j;
    int    iter_num;
    float  x_old[MAX_DIM];

    float Distance(float x[], float y[], int n);
    
    /* Initialize x */
    for (i = 0; i < n; i++)
        x[i] = b[i];

    iter_num = 0;
    do {
        iter_num++;
        
        for (i = 0; i < n; i++)
            x_old[i] = x[i];

        for (i = 0; i < n; i++){
            x[i] = b[i];
            for (j = 0; j < i; j++)
                x[i] = x[i] - A[i][j]*x_old[j];
            for (j = i+1; j < n; j++)
                x[i] = x[i] - A[i][j]*x_old[j];
            x[i] = x[i]/A[i][i];
        }
    } while ((iter_num < max_iter) && 
             (Distance(x,x_old,n) >= tol));

    if (Distance(x,x_old,n) < tol)
        return 1;
    else
        return 0;
} /* Jacobi */


/*********************************************************************/
float Distance(float x[], float y[], int n) {
    int i;
    float sum = 0.0;

    for (i = 0; i < n; i++) {
        sum = sum + (x[i] - y[i])*(x[i] - y[i]);
    }
    return sqrt(sum);
} /* Distance */


/*********************************************************************/
void Read_matrix(
         char*     prompt  /* in  */,
         MATRIX_T  A       /* out */,
         int       n       /* in  */) {
    int i, j;

    printf("%s\n", prompt);
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            scanf("%f", &A[i][j]);
}  /* Read_matrix */


/*********************************************************************/
void Read_vector(
         char*  prompt  /* in  */,
         float  x[]     /* out */,
         int    n       /* in  */) {
    int i;

    printf("%s\n", prompt);
    for (i = 0; i < n; i++)
        scanf("%f", &x[i]);
}  /* Read_vector */


/*********************************************************************/
void Print_vector(
         char*  title  /* in */,
         float  x[]    /* in */,
         int    n      /* in */) {
    int i;

    printf("%s\n", title);
    for (i = 0; i < n; i++)
        printf("%4.1f ", x[i]);
    printf("\n");
}  /* Print_vector */

/*********************************************************************/
void Print_matrix(
         char*     title  /* in */, 
         MATRIX_T  A      /* in */, 
         int       n      /* in */) {
    int i,j;

    printf("%s\n", title);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%f", A[i][j]);
        printf("\n");
    }
}  /* Print_matrix */
