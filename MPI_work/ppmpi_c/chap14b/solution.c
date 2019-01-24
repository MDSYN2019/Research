/* solution.c -- functions for keeping track of the best solution found
 *     so far -- for use in parallel tree search.
 *
 * See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
 *     search.
 */
#include <stdio.h>
#include "mpi.h"
#include "cio.h"
#include "node_stack.h"
#include "solution.h"

static int* best_solution;
static int* temp_solution;
static int* solution1;  
    /* first entry = cost, remaining entries, sibling */
    /* ranks on path from root to solution node       */
    /* Thus solution_size = 1 + max_depth + 1         */
static int* solution2;
static int  solution_size;

extern int max_depth;

/*********************************************************************/
/* Check queue for new solutions.  Return best solution found */
COST_T Best_solution(
           MPI_Comm  comm  /* in */) {
    MPI_Status  probe_status;
    MPI_Status  recv_status;
    int         done = FALSE;
    int         message_pending;
    int*        temp;

    while (!done) {
        MPI_Iprobe(MPI_ANY_SOURCE, SOLUTION_TAG, comm,
            &message_pending, &probe_status);
        if (message_pending) {
            MPI_Recv(temp_solution, solution_size, MPI_INT,
                probe_status.MPI_SOURCE, SOLUTION_TAG, comm,
                &recv_status);
            if (*temp_solution < *best_solution) {
                temp = temp_solution;
                temp_solution = best_solution;
                best_solution = temp;
            }
        } else {
            done = TRUE;
        }
    }  /* while */

    return Local_best_solution();
}  /* Best_solution */


/*********************************************************************/
COST_T Local_best_solution(void) { 

    return *best_solution;
} /* Local_best_solution */

/********************************************************************/ 
void Local_solution_update(
         COST_T  cost  /* in */,
         NODE_T  node  /* in */) {
    int   i;
    int*  bs_ptr;
    int*  a_ptr;

    *best_solution = cost;
    for (i = 0, bs_ptr = best_solution+1, a_ptr = Ancestors(node); 
             i < max_depth; i++, bs_ptr++, a_ptr++)
        *bs_ptr = *a_ptr;
    *bs_ptr = Sibling_rank(node);

}  /* Local_solution_update */


/********************************************************************/ 
/* Assumes plenty of system buffering.  So beware! */
void Bcast_solution(MPI_Comm comm){
    int my_rank;
    int p;
    int q;

    MPI_Comm_rank(comm, &my_rank);
    MPI_Comm_size(comm, &p);
    for (q = 0; q < my_rank; q++)
        MPI_Send(best_solution, solution_size, MPI_INT, q,
            SOLUTION_TAG, comm);
    for (q = my_rank+1; q < p; q++)
        MPI_Send(best_solution, solution_size, MPI_INT, q,
            SOLUTION_TAG, comm);
} /* Bcast_solution */

/*********************************************************************/
/* Return 0 if OK, -1 otherwise */
int Initialize_soln(
        int  max_depth  /* in */) {
    int i;

    solution_size = max_depth+2;
    solution1 = (int*) malloc(solution_size*sizeof(int));
    solution2 = (int*) malloc(solution_size*sizeof(int));
    if ((solution1 == (int*) NULL) || (solution2 == (int*) NULL))
        return -1;
 
    best_solution = solution1;
    temp_solution = solution2;
    best_solution[0] = INFINITY;  

    best_solution[1] = 0;        /* root is on best path */
    for (i = 1; i < solution_size; i++)
        best_solution[i] = -1;

    return 0;
}  /* Initialize_soln */


/*********************************************************************/
void Free_soln(void) {
    free(best_solution);
    free(temp_solution);
}  /* Free_soln */


/*********************************************************************/
/*  Assumes each process has access to stdout */
void Print_local_soln(
         MPI_Comm  comm  /* in */) {
    int i;
    int my_rank;

    MPI_Comm_rank(comm, &my_rank);

    printf("*************************************************************\n");
    printf("Process %d > Minimum cost = %d\n", my_rank, 
        Local_best_solution());
    printf("Min cost path = ");
    for (i = 1; i < solution_size; i++)
        printf("%d ", best_solution[i]);
    printf("\n");

} /* Print_local_soln */


/*********************************************************************/
void Print_solution(
         MPI_Comm  io_comm  /* in */) {
    int io_rank;
    int my_rank;
    int i;

    MPI_Comm_rank(io_comm, &my_rank);
    Get_io_rank(io_comm, &io_rank);
    if (my_rank == io_rank) {
        printf("Minimum cost = %d\n", Local_best_solution());
        printf("Path to minimum cost =\n");
        printf("    ");
        for (i = 1; i < solution_size; i++)
            printf("%d ", best_solution[i]);
        printf("\n");
    }
    
}  /* Print_soln */

/*********************************************************************/
/* Find owner of global best solution and broadcast */
void Update_solution(MPI_Comm comm){
    UPDATE_STRUCT_T local_data;
    UPDATE_STRUCT_T global_data;

    local_data.soln = Local_best_solution();
    MPI_Comm_rank(comm, &(local_data.rank));
    MPI_Allreduce(&local_data, &global_data, 1, update_mpi_t,
        MPI_MINLOC, comm);

    MPI_Bcast(best_solution, solution_size, MPI_INT,
        global_data.rank, comm);
}  /* Update_solution */

