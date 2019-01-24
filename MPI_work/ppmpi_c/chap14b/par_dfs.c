/* par_dfs.c -- local depth-first search function.  Performs iterative
 *     depth-first search on tree until local_stack is exhausted or
 *     it has expanded max_work nodes.
 *
 * See Chap 14, pp. 330 & ff., in PPMPI.
 */
#include "par_dfs.h"
#include "node_stack.h"
#include "solution.h"
#ifdef STATS
#include "stats.h"
#endif

extern int max_work;
extern int max_depth;
extern int max_children;
extern int my_rank;

/*********************************************************************/
void Par_dfs(
         STACK_T   local_stack  /* in/out */, 
         MPI_Comm  comm         /* in     */) {
    int     count;
    NODE_T  node;
    float   temp_solution;


    /* Search local subtree for a while */
    count = 0;
    while (!Empty(local_stack) && (count < max_work)) {
        node = Pop(local_stack);
#ifdef STATS
        Incr_stat(nodes_expanded);
#endif
#ifdef DEBUG
        Print_node(node, local_stack, comm);
#endif
        if (Solution(node)) {
            temp_solution = Evaluate(node);
            if (temp_solution < Best_solution(comm)) {
                Local_solution_update(temp_solution, node);
                Bcast_solution(comm);
            }
        } else if (Feasible(node, comm)) {
            Expand(node, local_stack);
        }
        count++;
    }  /* while */

} /* Par_dfs */


/********************************************************************/
int Solution(
        NODE_T  node  /* in */) {
    if (Depth(node) == max_depth)
        return TRUE;
    else
        return FALSE;
}  /* Solution */


/********************************************************************/
COST_T Evaluate(
           NODE_T    node  /* in  */) {
    
    return (COST_T) Seed(node);
}  /* Evaluate */ 


/********************************************************************/
int Feasible(NODE_T node, MPI_Comm comm) {
    return TRUE;
} 

/********************************************************************/
void Expand(
         NODE_T  parent       /* in     */, 
         STACK_T local_stack  /* in/out */) {

    int    num_children;
    int    i;
    NODE_T temp_node;
    NODE_T temp_parent;
    int    local_seed;
    static int divisor = 0;
    int    quotient;
    int    min_children;
    int    error;

    if (divisor == 0) 
        divisor = MAX((max_depth + 1)/max_children, 1);

    local_seed = Seed(parent);
    num_children = My_rand(local_seed) % (max_children + 1);
    quotient = Depth(parent)/divisor;
    min_children = max_children - quotient - 1;
    if (num_children < min_children)
        num_children = min_children;
    
    temp_parent = Copy_to_scratch(parent);
    local_seed = My_rand(local_seed + My_rand(Depth(parent)));   
    for (i = num_children - 1; i >= 0; i--) {
        error = Allocate_node(local_stack, &temp_node);
        if (error < 0) {
            fprintf(stderr, "Process %d > Stack overflow, In_use(stack) = %d, Max_size(stack) = %d\n", 
                my_rank, In_use(local_stack), Max_size(local_stack));
            fprintf(stderr, "Quitting!\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        local_seed = My_rand(local_seed);
        Initialize_node(temp_node, temp_parent, i, local_seed);
        Push(temp_node, local_stack);
    }

}  /* Expand */
