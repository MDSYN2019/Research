/* solution.h
 *
 * Declarations and definitions for manipulating the solutions
 */
#ifndef SOLUTION_H
#define SOLUTION_H

#include "mpi.h"
#include "node_stack.h"

typedef int COST_T;

typedef struct {
    COST_T soln;
    int    rank;
} UPDATE_STRUCT_T;

#define INFINITY 1000000
#define update_mpi_t MPI_2INT
#define cost_mpi_t MPI_INT

#define SOLUTION_TAG 100

COST_T Best_solution(MPI_Comm comm);
COST_T Local_best_solution(void);
void   Local_solution_update(COST_T cost, NODE_T node);
void   Bcast_solution(MPI_Comm comm);
int    Initialize_soln(int max_depth);
void   Update_solution(MPI_Comm comm);
void   Print_local_solution(MPI_Comm comm);
void   Print_solution(MPI_Comm io_comm);
void   Free_soln(void);

#endif
