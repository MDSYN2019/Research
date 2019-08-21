/* queue.h
 *
 * Definitions and declarations for accessing the message "queue".
 */
#ifndef QUEUE_H
#define QUEUE_H

#define REQUEST_TAG 10
#define WORK_TAG    20

#include "mpi.h"

void Clean_up_queues(MPI_Comm comm);
int  Work_requests_pending(MPI_Comm comm);
int  Get_dest(MPI_Comm comm);

#endif
