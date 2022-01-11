/* multi_files.c -- each process opens and writes to a different file.
 *
 * Input: none
 * Output: To file file.X, "Greetings from Process X", X = 0, 1, . . ., p-1
 *
 * See Chap 8, pp. 157 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"

main(int argc, char* argv[]) {
    FILE*  my_fp;
    int    my_rank;
    char   filename[100];

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    sprintf(filename,"file.%d", my_rank);
    my_fp = fopen(filename, "w");

    fprintf(my_fp, "Greetings from Process %d!\n", 
        my_rank);
    
    fclose(my_fp);  
    MPI_Finalize();
}
