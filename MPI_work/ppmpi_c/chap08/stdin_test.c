/* stdin_test.c -- test whether an MPI implementation allows input from
 *     stdin
 *
 * Input: an int
 * Output: prompt for input and int read (if stdin is OK)
 *
 * See Chap 8, p. 154 in PPMPI.
 */
#include <stdio.h>
#include "mpi.h"

main(int argc, char* argv[]) {
    int x;

    MPI_Init(&argc, &argv);

    printf("Enter an integer\n");
    fflush(stdout);
    scanf("%d", &x);
    printf("We read x = %d\n", x);
    fflush(stdout);

    MPI_Finalize();
}
