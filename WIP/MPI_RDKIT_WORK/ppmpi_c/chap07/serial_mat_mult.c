/* serial_mat_mult.c -- multiply two square matrices on a single processor
 *
 * Input: 
 *     n: order of the matrices
 *     A,B: factor matrices
 *
 * Output:
 *     C: product matrix
 *
 * See Chap 7, pp. 111 & ff in PPMPI
 */
#include <stdio.h>

#define MAX_ORDER 10

typedef float MATRIX_T[MAX_ORDER][MAX_ORDER];

main() {
    int       n;
    MATRIX_T  A;
    MATRIX_T  B;
    MATRIX_T  C;

    void Read_matrix(char* prompt, MATRIX_T A, int n);
    void Serial_matrix_mult(MATRIX_T A, MATRIX_T B, MATRIX_T C, int n);
    void Print_matrix(char* title, MATRIX_T C, int n);

    printf("What's the order of the matrices?\n");
    scanf("%d", &n);

    Read_matrix("Enter A", A, n);
    Print_matrix("A = ", A, n);
    Read_matrix("Enter B", B, n);
    Print_matrix("B = ", B, n);
    Serial_matrix_mult(A, B, C, n);
    Print_matrix("Their product is", C, n);

}  /* main */


/*****************************************************************/
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


/*****************************************************************/
/* MATRIX_T is a two-dimensional array of floats */
void Serial_matrix_mult(
        MATRIX_T   A   /* in  */,
        MATRIX_T   B   /* in  */,
        MATRIX_T   C   /* out */,
        int        n   /* in  */) {

    int i, j, k;

    void Print_matrix(char* title, MATRIX_T C, int n);

    Print_matrix("In Serial_matrix_mult A = ", A, n);
    Print_matrix("In Serial_matrix_mult B = ", B, n);

    for (i = 0; i < n; i++)    
        for (j = 0; j < n; j++) {
            C[i][j] = 0.0;
            for (k = 0; k < n; k++)
                C[i][j] = C[i][j] + A[i][k]*B[k][j];
            printf("i = %d, j = %d, c_ij = %f\n", i, j, C[i][j]);
        }
}  /* Serial_matrix_mult */


/*****************************************************************/
void Print_matrix(
         char*     title  /* in  */,
         MATRIX_T  C       /* out */,
         int       n       /* in  */) {
    int i, j;

    printf("%s\n", title);
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            printf("%4.1f ", C[i][j]);
        printf("\n");
    }
}  /* Read_matrix */
