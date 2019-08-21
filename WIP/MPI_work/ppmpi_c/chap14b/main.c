/* main.c 
 *
 * Main program that calls Par_tree_search.
 *
 * Input:
 *     1. max_depth: the maximum depth of nodes in the tree.
 *     2. cutoff_depth:  the maximum depth of nodes that can be sent
 *        to other processes in dynamic load balancing.
 *     3. max_work:  the maximum number of nodes expanded in a
 *        single call to Par_dfs.
 *     4. max_children:  the maximum number of children that a
 *        tree node can have.
 *     
 * Output:
 *     1. If STATS has been defined, statistics on the execution
 *        of the program.  See stats.c
 *
 * Algorithm:
 *     1. Start up MPI and get input.
 *     2. Call various setup functions.
 *     3. Process 0:  initialize root of tree
 *     4. Call Par_tree_search
 *     5. Print stats
 *     6. Clean up message queues
 *     7. Free storage
 *
 * Global Variables:
 *     max_work: maximum number of nodes expanded in a single call
 *         Par_dfs
 *     cutoff_depth: maximum depth of nodes redistributed in load
 *         balancing
 *     max_children: maximum number of children of a node
 *     max_depth: maximum depth of a node
 *     io_comm: communicator for I/O -- duplicate of MPI_COMM_WORLD
 *     p: number of processes in MPI_COMM_WORLD
 *     my_rank: process rank in MPI_COMM_WORLD
 *
 * See Chap 14, pp. 328 & ff., in PPMPI.
 */
#include <stdio.h>
#include <malloc.h>
#include "mpi.h"
#include "cio.h"

/* Local headers */
#include "node_stack.h"
#include "solution.h"
#include "par_tree_search.h"
#include "service_requests.h"
#include "terminate.h"
#include "queue.h"

#ifdef STATS
#include "stats.h"
#endif

/* Global variables */
int       max_work;
int       cutoff_depth;
int       max_children;
int       max_depth;

MPI_Comm  io_comm;
int       p;
int       my_rank;

void Get_root(NODE_T* root_ptr);

/*********************************************************************/
int
main(int argc, char* argv[]) {
    int       error;
    NODE_T    root;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);

    Cscanf(io_comm,
        "Enter max depth, cutoff depth, max work, and max children",
        "%d %d %d %d", &max_depth, &cutoff_depth, &max_work,
        &max_children);

    Setup_term_detect();

    error = Initialize_soln(max_depth);
    Cerror_test(io_comm, "Initialize_soln", error);

    error = Allocate_lists(max_depth, max_children);
    Cerror_test(io_comm, "Allocate_lists", error);

    error = Allocate_type_arrays(max_depth, max_children);
    Cerror_test(io_comm, "Allocate_type_arrays", error);

    if (my_rank == 0) {
        Get_root(&root);
    } else {
        root = NODE_NULL;
    }
   
    Par_tree_search(root, MPI_COMM_WORLD);

#ifdef STATS
    Print_stats(io_comm);
#endif

    Clean_up_queues(MPI_COMM_WORLD);

    Free_lists();
    Free_type_arrays();
    Free_soln();

    MPI_Finalize();

    return 0;
}  /* main */


/*********************************************************************/
void Get_root(
        NODE_T* root_ptr) {
    STACK_T stack;

    Allocate_root(root_ptr);
    Depth(*root_ptr) = 0;
    Sibling_rank(*root_ptr) = 0;
    Seed(*root_ptr) = 1;
    Size(*root_ptr) = Calc_size(*root_ptr);
    Get_local_stack(&stack);
    Push(*root_ptr, stack);
    
}  /* Get_root */
