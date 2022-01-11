/* node_stack.c -- functions for manipulating tree nodes and the stack.  For
 *     use in parallel tree search.
 *
 * See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
 *     search.
 */
#include <stdio.h>
#include <stdlib.h>  /* Get malloc and free */
#include "mpi.h"
#include "node_stack.h"

static STACK_STRUCT_T local_stack_struct;
STACK_T local_stack = &local_stack_struct;
static NODE_T scratch_node;

extern int max_depth;
extern MPI_Comm io_comm;
extern int my_rank;

/*********************************************************************/
void Set_stack_scalars(
         int      count  /* in  */,
         STACK_T  stack  /* out */) {
    int  subtract;

    In_use(stack) = count;
    subtract = Stack_int(stack, count-1);
    Top(stack) = &Stack_int(stack, count - subtract);
}  /* Set_stack_scalars */


/*********************************************************************/
/* Return index of stack element that precedes node.
NODE_T Pred(
           STACK_T  stack  /* in */,
           NODE_T   node   /* in */) {
    if (node > Stack_list(stack))
        return node - (*(node-1));
    else
        return NODE_NULL;
}  /* Pred */

/*********************************************************************/
void Copy_node(
         NODE_T  node1  /* in  */,
         NODE_T  node2  /* out */) {
    int   i;
    int*  n2_ptr;
    int*  n1_ptr;
    int   size;

    size = Size(node1);
    for (i = 0, n2_ptr = node2, n1_ptr = node1; 
         i < size; i++, n2_ptr++, n1_ptr++)
        *n2_ptr = *n1_ptr;

}  /* Copy_node */


/*********************************************************************/
NODE_T Copy_to_scratch(
         NODE_T  node  /* in */) {
    int   i;
    int*  s_ptr;
    int*  n_ptr;

    for (i = 0, s_ptr = scratch_node, n_ptr = node; 
         i < Size(node); i++, s_ptr++, n_ptr++)
        *s_ptr = *n_ptr;

    return scratch_node;
} /* Copy_to_scratch */


/*********************************************************************/
/* Return 0 if successful, negative otherwise */
int Allocate_lists(
        int  max_depth     /* in */,
        int  max_children  /* in */) {
    int max_size;

    max_size = (max_depth+1)*(NODE_MEMBERS-1);
    max_size = max_size + max_depth*(max_depth+1)/2;
    max_size = (max_children-1)*max_size;
    max_size = max_size + (NODE_MEMBERS - 1) + max_depth;
    Stack_list(local_stack) = (int*) malloc(max_size*sizeof(int));
    scratch_node = (NODE_T) malloc(MAX_NODE_SIZE*sizeof(int));

    if ( (Stack_list(local_stack) == ((int*) NULL)) ||
         (scratch_node == (NODE_T) NULL) ) 
        return -1;
    else
        Max_size(local_stack) = max_size;

    Top(local_stack) = NODE_NULL;
    In_use(local_stack) = 0;

#ifdef DEBUG
{
    int my_rank;
    MPI_Comm_rank(io_comm, &my_rank);
    printf("Process %d > In Allocate_lists, Max_size(local_stack) = %d\n", 
        my_rank, Max_size(local_stack));
    fflush(stdout);
}
#endif
    return 0;
        
}  /* Allocate_lists */


/*********************************************************************/
void Free_lists(void) {
    free(Stack_list(local_stack));
    free(scratch_node);
}  /* Free_lists */


/*********************************************************************/
int Allocate_root(
        NODE_T* root_ptr /* out */) {
    
    /* Top(comm_stack) == NODE_NULL */
    Top(local_stack) = Stack_list(local_stack);
    
    *root_ptr = Top(local_stack);
    return 0;

}  /* Allocate_root */


/*********************************************************************/
void Get_local_stack(STACK_T* local_stack_ptr) {
    *local_stack_ptr = local_stack;
}  /* Get_local_stack */


/*********************************************************************/
/* Returns 0 if successful, negative otherwise             */
/* Doesn't change Top of stack or In_use -- Push does that */
int Allocate_node(
        STACK_T stack    /* in/out */,
        NODE_T* node_ptr /* out    */) {
    
    if (Top(stack) == NODE_NULL) {
        *node_ptr = Stack_list(stack);
    } else if (In_use(stack) >= (Max_size(stack) - MAX_NODE_SIZE)) {
        return -1;
    } else {
        *node_ptr = Next(stack,Top(stack));
    }
    
    return 0;

}  /* Allocate_node */


/*********************************************************************/
void Initialize_node(
        NODE_T  node          /* out */, 
        NODE_T  parent        /* in  */, 
        int     sibling_rank  /* in  */,
        int     seed          /* in  */) {
    int  i;

    Depth(node) = Depth(parent) + 1;
    Sibling_rank(node) = sibling_rank;
    Seed(node) = seed;
    for (i = 0; i < Depth(parent); i++)
        Ancestor(node,i) = Ancestor(parent,i);
    Ancestor(node,i) = Sibling_rank(parent);
    Size(node) = Calc_size(node);
    
}  /* Initialize_node */


/*********************************************************************/
void Initialize_stack(NODE_T node, STACK_T* stack_ptr) {
}  /* Initialize_stack */ 


/*********************************************************************/
int Empty(STACK_T stack){ 
    if (Top(stack) == NODE_NULL)
        return TRUE; 
    else
        return FALSE;
}  /* Empty */


/*********************************************************************/
/* Assumes Empty has already been called */
NODE_T Pop(STACK_T stack) { 
    NODE_T temp;

    temp = Top(stack);
    Top(stack) = Pred(stack, temp);
    In_use(stack) = In_use(stack) - Size(temp); 
#ifdef DEBUG
{
/*
    int my_rank;
    MPI_Comm_rank(io_comm, &my_rank);
    printf("Process %d > In Pop, Max_size(stack) = %d, In_use(stack) = %d\n", 
        my_rank, Max_size(local_stack), In_use(stack));
    fflush(stdout);
*/
}
#endif

    return temp; 
}  /* Pop */


/*********************************************************************/
void Push(
         NODE_T   node   /* in     */, 
         STACK_T  stack  /* in/out */) {

    Top(stack) = node;
    In_use(stack) = In_use(stack) + Size(node);
#ifdef DEBUG
{
/*
    int my_rank;
    MPI_Comm_rank(io_comm, &my_rank);
    printf("Process %d > In Push, Max_size(stack) = %d, In_use(stack) = %d\n", 
        my_rank, Max_size(local_stack), In_use(stack));
    fflush(stdout);
*/
}
#endif
}  /* Push */


/*********************************************************************/
/* Assumes all processes have access to stdout -- for use in debugging */
void Print_node(
         NODE_T    node  /* in */, 
         STACK_T   stack /* in */,
         MPI_Comm  comm  /* in */) {
    int i;
    int my_rank;

    MPI_Comm_rank(comm, &my_rank);

    printf("--------------------------------------------------------------\n");
    printf("Process %d > Depth = %d, Sibling rank = %d, Seed = %d, Size = %d\n",
        my_rank, Depth(node), Sibling_rank(node), Seed(node), Size(node));
    printf("Process %d > Ancestors = ", my_rank);
    for (i = 0; i < Depth(node); i++)
        printf("%d ", Ancestor(node,i));
    printf("\n");
    fflush(stdout);
}  /* Print_node */
    

/*********************************************************************/
void Print_stack(
         char*     title,
         STACK_T   stack,
         MPI_Comm  comm) {
    NODE_T node;
    int    my_rank;
    
    MPI_Comm_rank(comm, &my_rank);
 
    printf("**************************************************************\n");
    printf("Process %d > %s.  Stack = \n", my_rank, title);
    node = Top(stack);
    while (node != NODE_NULL) {
        Print_node(node, stack, comm);
        node = Pred(stack, node);
    }
    printf("**************************************************************\n");
    fflush(stdout);
}  /* Print_stack */


/*********************************************************************/
void Print_stack_list(
         int       count  /* in */,
         STACK_T   stack  /* in */,
         MPI_Comm  comm   /* in */) {
    int  i;
    int* s_ptr;
    int  my_rank;

    MPI_Comm_rank(comm, &my_rank);
    printf("Process %d > count = %d, Stack_list = ", my_rank, count);
    for (i = 0, s_ptr = Stack_list(stack); i < count; i++, s_ptr++)
        printf("%d ", *s_ptr);
    printf("\n");
    fflush(stdout);

}  /* Print_stack_list */
