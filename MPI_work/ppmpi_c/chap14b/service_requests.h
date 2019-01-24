/* service_requests.h
 *
 * Definitions and declarations used by the service_requests functions 
 */
#ifndef SERVICE_REQUESTS_H
#define SERVICE_REQUESTS_H

#include "node_stack.h"
#include "mpi.h"

/* Called in Service_requests */
int   Allocate_type_arrays(int cutoff_depth, int max_children);
void  Free_type_arrays(void);
void  Service_requests(STACK_T local_stack, MPI_Comm comm);
int   Nodes_available(STACK_T local_stack);
void  Split(STACK_T local_stack, STACK_T* send_stack_ptr);
void  Send_work(int destination, STACK_T send_stack, MPI_Comm comm);
void  Compress(STACK_T stack);
void  Send_reject(int destination, MPI_Comm comm);
void  Send_all_rejects(MPI_Comm comm);

#endif
