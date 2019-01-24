/* node_stack.h -- header file for node_stack.c
 *
 * Basic definitions and declarations for stack and tree node manipulation.
 */
#ifndef NODE_STACK_H
#define NODE_STACK_H
#include <stdio.h>  /* Get definition of NULL */
#include "mpi.h"

#define TRUE 1
#define FALSE 0

typedef int* NODE_T;
#define NODE_NULL  ((NODE_T) NULL)

#define NODE_MEMBERS       5   /* includes depth, sibling_rank, cost, */
                               /*    ancestors (1 member), size       */
#define Depth(node)        (*(node))
#define Sibling_rank(node) (*((node)+1))
#define Seed(node)         (*((node)+2))
#define Ancestors(node)    ((node)+3)
/* Number of ancestors = depth */
#define Ancestor(node,i)   (*((node)+3+(i)))
#define Size(node)         (*((node)+ NODE_MEMBERS + Depth(node) - 2))
#define Calc_size(node)    (NODE_MEMBERS + Depth(node) - 1)

extern int max_depth;
#define MAX_NODE_SIZE (NODE_MEMBERS + max_depth - 1)

#define node_mpi_t MPI_INT

/*
#define MODULUS 65536
#define MULTIPLIER 25173
#define INCREMENT 13849
*/
#define MODULUS 65537
#define MULTIPLIER 16383
#define INCREMENT 1001
#define My_rand(seed)  ((MULTIPLIER*(seed) + INCREMENT) % MODULUS)

typedef struct {
    int     max_size;    /* size of stack_list    */
    int     in_use;      /* number of ints in use */
    NODE_T  stack_top;   
    int*    stack_list;
} STACK_STRUCT_T;

typedef STACK_STRUCT_T* STACK_T;
/* Top references an actual record unless the stack is empty -- */
/*    in which case it's NULL                                   */
#define Max_size(stack)    ((stack)->max_size)
#define In_use(stack)      ((stack)->in_use)
#define Top(stack)         ((stack)->stack_top)
#define Stack_list(stack)  ((stack)->stack_list)
#define Stack_int(stack,i) (*(((stack)->stack_list) + (i)))

#define Next(stack,node)   ((node) + Size(node))
/*
#define Pred(stack,node)   ((node) > Stack_list(stack) ? \
                            ((node) - (*(node-1))) : \
                            NODE_NULL)
*/
void    Set_stack_scalars(int count, STACK_T stack);
NODE_T  Pred(STACK_T stack, NODE_T node);
NODE_T  Copy_to_scratch(NODE_T node);
void    Copy_node(NODE_T node1, NODE_T node2);
int     Allocate_lists(int max_depth, int max_children);
void    Free_lists(void);
int     Allocate_root(NODE_T* root_ptr);
void    Get_local_stack(STACK_T* local_stack_ptr);
int     Allocate_node(STACK_T stack, NODE_T* node_ptr);
void    Initialize(NODE_T node, STACK_T* stack_ptr);
void    Initialize_node(NODE_T node, NODE_T parent, 
            int sibling_rank, int seed);
int     Empty(STACK_T stack);
NODE_T  Pop(STACK_T stack);
void    Push(NODE_T node, STACK_T stack);
void    Print_node(NODE_T node, STACK_T stack, MPI_Comm comm);
void    Print_stack(char* title, STACK_T stack, MPI_Comm comm);
void    Print_stack_list(int count, STACK_T stack, MPI_Comm comm);

#endif
