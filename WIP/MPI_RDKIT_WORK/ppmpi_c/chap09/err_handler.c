/* err_handler.c -- change default error handler in MPI to MPI_ERRORS_RETURN
 *
 * Input: none.
 * Output: Error message from each process.
 *
 * See Chap 9, pp. 210 & ff in PPMPI.
 */

#include <stdio.h>
#include "mpi.h"

main(int argc, char* argv[]) {
    int p;
    int my_rank;
    char error_message[MPI_MAX_ERROR_STRING];
    int message_length;
    int error_code;
    int x[2];
    int count = 2;
    MPI_Comm comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    error_code = MPI_Errhandler_set(MPI_COMM_WORLD,
                  MPI_ERRORS_RETURN);

    if (my_rank == 0) {
        x[0] = 1;
        x[1] = 2;
    }
    error_code = MPI_Bcast(x, count, MPI_INT, 0, comm);
    if (error_code != MPI_SUCCESS) {
        MPI_Error_string(error_code, error_message,
            &message_length);
        fprintf(stderr, "Error in call to MPI_Bcast = %s\n",
            error_message);
        fprintf(stderr, "Exiting from function XXX\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }


    MPI_Finalize();
}  /* main */
