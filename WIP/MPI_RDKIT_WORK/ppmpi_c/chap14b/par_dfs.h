/* par_dfs.h
 *
 * Declarations and definitions used in the par_dfs functions
 */
#ifndef PAR_DFS_H
#define PAR_DFS_H

#include "node_stack.h"
#include "solution.h"
#include "mpi.h"

#define MAX(a,b) ((a) >= (b) ? (a) : (b))

void    Par_dfs(STACK_T local_stack, MPI_Comm comm);
void    Bcast_solution(MPI_Comm comm);
int     Solution(NODE_T node);
COST_T  Evaluate(NODE_T node);
int     Feasible(NODE_T node, MPI_Comm comm);
void    Expand(NODE_T parent, STACK_T local_stack);

#endif
