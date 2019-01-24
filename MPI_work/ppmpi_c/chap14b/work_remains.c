/* work_remains.c -- controls exit from main loop in parallel tree
 *     search program.  Checks whether local stack is empty.  If not,
 *     returns TRUE.  If it is, executes the following algorithm.
 *
 *         Send "energy" to process 0 (see terminate.c, below).
 *         while(TRUE) {
 *             Send reject messages to all processes requesting work.
 *             if (the search is complete) {
 *                 Cancel outstanding work requests.
 *                 return FALSE.
 *             } else if (there's no outstanding request) {
 *                 Generate a process rank to which to send a request.
 *                 Send a request for work.
 *             } else if (a reply has been received) {
 *                 if (work has been received) return TRUE.
 *             }
 *         }
 *
 * See Chap 14, pp. 332 & ff, in PPMPI.
 */
#include <stdlib.h>
#include "work_remains.h"
#include "node_stack.h"
#include "terminate.h"
#include "service_requests.h"
#include "queue.h"
#ifdef STATS
#include "stats.h"
#endif

extern int my_rank;
extern int p;

MPI_Request posted_recv;

/*********************************************************************/
int Work_remains(
        STACK_T    local_stack  /* in/out */,
        MPI_Comm   comm         /* in     */ ) {

    int      work_available;
    int      work_request_process;
    int      request_sent;

#ifdef STATS
    Start_time(work_rem_time);
#endif 
   
#ifdef DEBUG
    printf("Process %d > In Work_remains \n", my_rank);
    fflush(stdout);
#endif

    if (!Empty(local_stack)) {
#ifdef STATS
        Finish_time(work_rem_time);
#endif 

#ifdef DEBUG
        printf("Process %d > In Work_remains, stack not empty\n", my_rank);
        fflush(stdout);
#endif
        return TRUE;
    } else {
        Return_energy(comm);
#ifdef DEBUG
        printf("Process %d > In Work_remains, stack empty, energy returned\n", 
           my_rank);
        fflush(stdout);
#endif
        request_sent = FALSE;
        while (TRUE) {
            Send_all_rejects(comm);
            if (Search_complete(comm)) {
#ifdef DEBUG
                printf("Process %d > In Work_remains, search complete\n", 
                    my_rank);
                fflush(stdout);
#endif

#ifdef STATS
                Finish_time(work_rem_time);
#endif 
                if (request_sent) Cancel_request();
                return FALSE;
            } else if (!request_sent) {
                work_request_process = New_request(comm);
                Send_request(work_request_process, comm);
                request_sent = TRUE;
#ifdef DEBUG
                printf("Process %d > In Work_remains, request sent to %d\n", 
                    my_rank, work_request_process);
                fflush(stdout);
#endif
            } else if (Reply_received(work_request_process,
                
                    &work_available, local_stack, comm)) {
#ifdef DEBUG
                printf("Process %d > In Work_remains, reply received from %d\n", 
                        my_rank, work_request_process);
                fflush(stdout);
#endif
                if (work_available) {
#ifdef DEBUG
                    printf("Process %d > In Work_remains, received work from %d\n", 
                        my_rank, work_request_process);
                    fflush(stdout);
#endif

#ifdef STATS
                    Finish_time(work_rem_time);
#endif 
                    return TRUE;
                } else {
#ifdef DEBUG
                    printf("Process %d > In Work_remains, no work at %d\n", 
                        my_rank, work_request_process);
                    fflush(stdout);
#endif
                    request_sent = FALSE;
                }
            }

        } /* while (TRUE) */

    } 
} /* Work_remains */


/*********************************************************************/
/* Generates a process rank to which a new request should be sent --
 *     simply rand() % p
 */
int New_request(
        MPI_Comm  comm  /* in */) {
    static int seeded = FALSE;
    int        rank;
    

    if (!seeded) {
        srand(my_rank);
        seeded = TRUE;
    }
 
    while ( (rank = rand() %p) == my_rank);

    return rank;
}  /* New_request */


/*********************************************************************/
void Send_request(
         int       work_request_process  /* in  */, 
         MPI_Comm  comm                  /* out */) {

    int     x = 0;
    STACK_T local_stack;

    MPI_Send(&x, 1, MPI_INT, work_request_process, REQUEST_TAG,
        comm); 

    /* Post nonblocking receive */
    Get_local_stack(&local_stack);
    MPI_Irecv(Stack_list(local_stack), Max_size(local_stack),
        MPI_INT, work_request_process, WORK_TAG, comm, 
        &posted_recv);
#ifdef STATS
    Incr_stat(requests_sent);
#endif
}  /* Send_request */


/*********************************************************************/
int Reply_received(
        int       work_request_process  /* in  */,
        int*      work_available        /* out */,
        STACK_T   local_stack           /* out */, 
        MPI_Comm  comm                  /* in  */) {
    int         reply_received;
    MPI_Status  status;
    int         count;

    MPI_Test(&posted_recv, &reply_received, &status);

    if (reply_received) {
        if (*Stack_list(local_stack) == -1) {
#ifdef STATS
            Incr_stat(rejects_recd);
#endif
            *work_available = FALSE;
        } else {
            MPI_Get_count(&status, MPI_INT, &count);
            Set_stack_scalars(count, local_stack);
#ifdef DEBUG
            Print_stack_list(count, local_stack, comm);
#endif
            Recv_half_energy(work_request_process, comm);
#ifdef STATS
            Incr_stat(work_recd);
#endif
            *work_available = TRUE;
        }
        return TRUE;
    } else {
        return FALSE;
    }
}  /* Reply_received */


/*********************************************************************/
/* Tree search has completed, but there is still an outstanding work */
/*    request.  Try to cancel it.                                    */
void Cancel_request(void) {
    MPI_Status status;

/*  Cancel not implemented in mpich prior to 1.0.13 */
    MPI_Cancel(&posted_recv);
    MPI_Wait(&posted_recv, &status);

}  /* Cancel_request */
