/* queue.c -- functions for handling pending messages -- for use in
 *     parallel tree search
 *
 * See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
 *     search.
 */

#include "queue.h"
#include "mpi.h"
#include "node_stack.h"

extern int my_rank;

/*********************************************************************/
/* The only pending messages should be requests for work, rejects or */
/*     solutions.  Use stack as scratch area for receives.           */
/*     WARNING:  WILL OVERWRITE CONTENTS OF STACK!                   */
void Clean_up_queues(
         MPI_Comm  comm  /* in */) {
    int         done = FALSE;
    int         message_pending;
    STACK_T     stack;
    MPI_Status  status;

    Get_local_stack(&stack);
    while (!done) {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &message_pending,
            &status);
        if (message_pending) {
            MPI_Recv(Stack_list(stack), Max_size(stack), MPI_INT,
                MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &status);
        } else {
            done = TRUE;
        }
    }  /* while */
}  /* Clean_up_queues */


/*********************************************************************/
int Work_requests_pending(
        MPI_Comm  comm  /* in */) {
    MPI_Status  status;
    int         message_in_queue;

    MPI_Iprobe(MPI_ANY_SOURCE, REQUEST_TAG, comm, &message_in_queue,
        &status);

    if (message_in_queue) {
        return TRUE;
    } else {
        return FALSE;
    }
}  /* Work_requests_pending */


/*********************************************************************/
/* Only called when a message is pending */
int Get_dest(
        MPI_Comm  comm  /* in */) {
    int         node_request;
    MPI_Status  status;

    MPI_Recv(&node_request, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_TAG,
        comm, &status);
    return status.MPI_SOURCE;
}  /* Get_dest */
