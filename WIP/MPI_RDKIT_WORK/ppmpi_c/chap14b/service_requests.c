/* service_requests.c -- functions for handling requests for work
 *    from other processes.
 *
 * See Chap 14, p. 332, in PPMPI. 
 */
#include <stdlib.h>
#include "service_requests.h"
#include "queue.h"
#include "terminate.h"
#include "node_stack.h"
#ifdef STATS
#include "stats.h"
#endif

extern int my_rank;
extern int cutoff_depth;
extern MPI_Comm io_comm;

MPI_Datatype send_stack_mpi_t;

static int   allocated_size;
static int*  block_lengths = (int*) NULL;
static int*  displacements = (int*) NULL; 
static NODE_T* node_list;
static int   node_count;


/*********************************************************************/
/* Allocate arrays for use in constructing MPI derived datatypes. 
 *     Return negative if malloc fails, 0 otherwise 
 */
int Allocate_type_arrays(
        int  cutoff_depth  /* in */,
        int  max_children  /* in */) {

    allocated_size = cutoff_depth*(max_children - 1)/2;

    block_lengths = (int*) malloc(allocated_size*sizeof(int));
    displacements = (int*) malloc(allocated_size*sizeof(int));
    node_list     = (NODE_T*) malloc(allocated_size*sizeof(NODE_T));

    if ( (block_lengths == (int*) NULL) ||
         (displacements == (int*) NULL) ||
         (node_list     == (NODE_T*) NULL) )
        return -1;
    else
        return 0;
}  /* Allocate_type_arrays */


/*********************************************************************/
void Free_type_arrays(void) {
    free(block_lengths);
    free(displacements);
    free(node_list);
}  /* Free_type_arrays */
    

/*********************************************************************/
/* Use MPI_Iprobe (via Work_requests_pending) to see whether other
 *     processes have requested work.  If so, check local stack to 
 *     see whether there are any nodes above cutoff_depth.  If so, 
 *     send half of them to the requesting process.  If not, send
 *     a "reject" message.  The splitting of the stack is done using
 *     MPI_Type_indexed.  After the data has been sent.  The stack
 *     is "compressed".
 */
void Service_requests(
         STACK_T   local_stack  /* in/out */,
         MPI_Comm  comm         /* in     */) {
    STACK_T send_stack;
    int     destination;

#ifdef DEBUG
    printf("Process %d > In Service_requests\n", my_rank);
    fflush(stdout);
#endif

    while (Work_requests_pending(comm)) {
#ifdef DEBUG
        printf("Process %d > In Service_requests, queue not empty\n", my_rank);
        fflush(stdout);
#endif
        destination = Get_dest(comm);
#ifdef DEBUG
        printf("Process %d > In Service_requests, dest = %d\n", 
               my_rank, destination);
        fflush(stdout);
#endif
        if (Nodes_available(local_stack)) {
#ifdef DEBUG
            printf("Process %d > In Service_requests, work available\n", 
               my_rank);
            fflush(stdout);
#endif
            Split(local_stack, &send_stack);
            Send_work(destination, send_stack, comm);
#ifdef DEBUG
            printf("Process %d > In Service_requests, work sent to %d\n", 
               my_rank, destination);
            Print_stack("After compress", local_stack, comm);
            fflush(stdout);
#endif
        } else {
            Send_reject(destination, comm);
#ifdef DEBUG
            printf("Process %d > In Service_requests, reject sent to %d\n", 
               my_rank, destination);
            fflush(stdout);
#endif
        }
    }

} /* Service_requests */


/*********************************************************************/
/* Count until at least two nodes above cutoff depth are found */
/*     Should probably be combined with Split                  */
int Nodes_available(
        STACK_T  local_stack  /* in */) {
    NODE_T node;
    int    count = 0;

    node = Top(local_stack);
    while ((node != NODE_NULL) && (count < 2)) {
        if (Depth(node) <= cutoff_depth)
            count++;
        node = Pred(local_stack, node);
    }
    
    if (count >= 2)
        return TRUE;
    else
        return FALSE;
}  /* Nodes_available */


/*********************************************************************/
/* Split builds a derived datatype that picks out the nodes to be */
/*     sent */
void Split(
         STACK_T   local_stack    /* in  */, 
         STACK_T*  send_stack_ptr /* out */) {
    int     index = 0;
    NODE_T  node;
    int     odd = 0;

    node = Stack_list(local_stack);
    node_count = 0;
    while ( (node != NODE_NULL) && (node <= Top(local_stack)) ) {
        if (Depth(node) <= cutoff_depth)
            if (odd) {
                block_lengths[node_count] = Size(node);
                displacements[node_count] = index;
                node_list[node_count] = node;
                node_count++;
                odd = 0;   /* skip the next one above cutoff */
            } else {
                odd = 1;   /* take the next one above cutoff */
            }
        index = index + Size(node);
        node = Next(local_stack, node);
    }

    MPI_Type_indexed(node_count, block_lengths, displacements, MPI_INT,
        &send_stack_mpi_t);
    MPI_Type_commit(&send_stack_mpi_t);

    Get_local_stack(send_stack_ptr);

}  /* Split */


/*********************************************************************/
void Send_work(
         int       destination  /* in */,
         STACK_T   send_stack   /* in */, 
         MPI_Comm  comm         /* in */) {

    MPI_Send(Stack_list(send_stack), 1, send_stack_mpi_t, destination,
        WORK_TAG, comm);
    MPI_Type_free(&send_stack_mpi_t);
    Send_half_energy(destination, comm);
#ifdef STATS
    Incr_stat(work_sent);
#endif
    Compress(send_stack);
}  /* Send_work */


/*********************************************************************/
void Compress(
         STACK_T  stack  /* in/out */) {
    NODE_T compress_point;  /* In general, not a node, just a point */
                            /*     in the stack                     */
    NODE_T delete_node;
    NODE_T save_node;
    int    new_in_use;
    int    i;

#ifdef DEBUG
    printf("Process %d > In Compress, raw stack_list = \n",my_rank);
    Print_stack_list(In_use(stack), stack, io_comm);
#endif

    compress_point = node_list[0];
    new_in_use = displacements[0];
    
    for (i = 0; i < node_count; i++) {
        delete_node = node_list[i];
        if (Next(stack, delete_node) <= Top(stack)) {
            save_node = Next(stack, delete_node);
            Copy_node(save_node, compress_point);
            new_in_use = new_in_use + Size(compress_point);
            compress_point = Next(stack, compress_point);
        }
    }

    Top(stack) = Pred(stack, compress_point);
    In_use(stack) = new_in_use;

#ifdef DEBUG
    printf("Process %d > After Compress, raw stack_list = \n",my_rank);
    Print_stack_list(In_use(stack), stack, io_comm);
    printf("Process %d > After Compress, Top(stack) = \n", my_rank);
    Print_node(Top(stack), stack, io_comm);
#endif
}  /* Compress */


/*********************************************************************/
void Send_reject(
         int       destination  /* in */,
         MPI_Comm  comm         /* in */){
    int x = -1;

    MPI_Send(&x, 1, MPI_INT, destination, WORK_TAG, comm);
#ifdef STATS
    Incr_stat(rejects_sent);
#endif
}  /* Send_reject */


/*********************************************************************/
void Send_all_rejects(
         MPI_Comm  comm) {
    int destination;

    while (Work_requests_pending(comm)) {
        destination = Get_dest(comm);
        Send_reject(destination, comm);
    }
}  /* Send_all_rejects */
