/* par_tree_search.c -- controlling program for parallel tree search.
 *     Process 0 generates initial tree and distributes it among the
 *     processes.  Then each process runs through the basic loop
 *
 *         do {
 *             Local depth-first search for a while;
 *             Handle requests for work;
 *         } while (There's any work left on all processes);
 *
 * After exiting the loop, the global best solution is updated and
 * printed.
 * 
 * See Chap 14, pp. 328 & ff., in PPMPI.
 */
#include <stdio.h>
#include <malloc.h>
#include "mpi.h"
#include "cio.h"
#ifdef STATS
#include "stats.h"
#endif

/* Local headers */
#include "node_stack.h"
#include "solution.h"
#include "par_tree_search.h"
#include "par_dfs.h"
#include "service_requests.h"
#include "work_remains.h"
#include "terminate.h"

/* Global variables */
extern int       max_children;
extern int       max_depth;
extern MPI_Comm  io_comm;
extern int       p;
extern int       my_rank;



/*********************************************************************/
void Par_tree_search(
        NODE_T       root      /* in */, 
        MPI_Comm     comm      /* in */) {

    NODE_T*  node_list;
    NODE_T   node;
    STACK_T  local_stack;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    /* Generate initial set of nodes, 1 per process */
    if (my_rank == 0) {
        Generate(root, &node_list, p, comm);
    }
    Scatter(node_list, &node, comm);

    Initialize(node, &local_stack);

    do {
        /* Search for a while */
#ifdef STATS
        Start_time(par_dfs_time);
#endif
        Par_dfs(local_stack, comm);
#ifdef STATS
        Finish_time(par_dfs_time);
#endif

        /* Service requests for work. */
#ifdef STATS
        Start_time(svc_req_time);
#endif
        Service_requests(local_stack, comm);
#ifdef STATS
        Finish_time(svc_req_time);
#endif

        /* If local_stack isn't empty, return.          */
        /* If local_stack is empty, send                */
        /* requests until either we get work, or we     */
        /* receive a message terminating program.       */       
    } while(Work_remains(local_stack, comm));

    /* Get global best solution */
    Update_solution(comm);
    Print_solution(io_comm);
} /* Par_tree_search */


/*********************************************************************/
/* Called in Par_tree_search by process 0 to generate initial tree --
 *     basically a local depth-first search. */
void Generate(
         NODE_T    root       /* in  */, 
         NODE_T**  node_list  /* out */,
         int       size       /* in  */,
         MPI_Comm  comm       /* in  */) {
    int      stack_size = 1;
    NODE_T   node;
    COST_T   temp_solution;
    STACK_T  local_stack;

    Get_local_stack(&local_stack);

    while (!Empty(local_stack) && (stack_size < size)) {
        node = Pop(local_stack);
#ifdef STATS
        Incr_stat(nodes_expanded);
#endif

#ifdef DEBUG
        Print_node(node, local_stack, comm);
#endif
        stack_size = stack_size - 1;
        if (Solution(node)) {
            temp_solution = Local_evaluate(node);
            if (temp_solution < Local_best_solution()) {
                Local_solution_update(temp_solution, node);
            }
        } else if (Feasible(node, comm)) {
            PTS_expand(node, local_stack, &stack_size);
        }
    }  /* while */

#ifdef DEBUG
    printf("Process 0 > stack_size = %d\n", stack_size);
    Print_stack("Stack after Generate",local_stack, comm);
    fflush(stdout);
#endif

    if (stack_size < size) {
        fprintf(stderr, "Too few nodes in tree!\n");
        fprintf(stderr, "Aborting.\n");
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

} /* Generate */


/*********************************************************************/
COST_T Local_evaluate(
           NODE_T  node  /* in */) {

    return ((COST_T) My_rand(Seed(node)));
}  /* Local_evaluate */


/*********************************************************************/
/* Called by process 0 for expanding nodes in generating initial tree */
void PTS_expand(
         NODE_T   parent          /* in     */, 
         STACK_T  local_stack     /* in/out */, 
         int*     stack_size_ptr  /* in/out */) {

    int    num_children;
    int    i;
    NODE_T temp_node;
    NODE_T temp_parent;
    int    local_seed;
    int    divisor;
    int    quotient;
    int    min_children;
    int    error;

    /* MAX is defined in par_dfs.h */
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
        *stack_size_ptr = *stack_size_ptr + 1;
    }

} /*  PTS_expand  */



/*********************************************************************/
/* Process 0 scatters initial tree among processes */
void Scatter(
         NODE_T*   node_list  /* in  */, 
         NODE_T*   node_ptr   /* out */, 
         MPI_Comm  comm       /* in  */) {

    int        p, q;
    int        my_rank;
    STACK_T    local_stack;
    NODE_T     node;
    MPI_Status status;
    int        error;

    MPI_Comm_size(comm, &p);
    MPI_Comm_rank(comm, &my_rank);

    Get_local_stack(&local_stack);
    if (my_rank == 0) {
        /* Don't use MPI_Scatter because we will need an extra */
        /*     buffer to avoid overwriting local_stack         */
        for (q = p-1; q >= 1; q--) {
            node = Pop(local_stack);
            MPI_Send(node, Calc_size(node), node_mpi_t, q,
                0, comm);
        }
    } else {
        error = Allocate_node(local_stack, node_ptr);
        if (error < 0) {
            fprintf(stderr,"Process %d > Stack overflow, In_use(stack) = %d, Max_size(stack) = %d\n",
                my_rank, In_use(local_stack), Max_size(local_stack));
            fprintf(stderr, "Quitting!\n");
            MPI_Abort(MPI_COMM_WORLD, -1);
        }
        MPI_Recv(*node_ptr, MAX_NODE_SIZE, node_mpi_t, 0, 0,
            comm, &status);
    }
        
}  /* Scatter */


/*********************************************************************/
void Initialize(
         NODE_T    node       /* in  */, 
         STACK_T*  stack_ptr  /* out */) {

    Get_local_stack(stack_ptr);
    if (my_rank != 0) {
        Push(node, *stack_ptr);
    }
}  /* Initialize */
