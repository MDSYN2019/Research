/* par_tree_search.h
 *
 * Definitions and declarations used in par_tree_search.c
 */
#ifndef PAR_TREE_SEARCH_H
#define PAR_TREE_SEARCH_H

#include "mpi.h"
#include "node_stack.h"
#include "solution.h"

/* Local headers */
#include "node_stack.h"

void   Par_tree_search(NODE_T root, MPI_Comm comm);
void   Generate(NODE_T root, NODE_T** node_list, int p, MPI_Comm comm);
void   Scatter(NODE_T* node_list, NODE_T* node_ptr, MPI_Comm comm);
COST_T Local_evaluate(NODE_T node);
void   PTS_expand(NODE_T parent, STACK_T stack, int* stack_size_ptr);
void   Initialize(NODE_T node, STACK_T* stack_ptr);

#endif
