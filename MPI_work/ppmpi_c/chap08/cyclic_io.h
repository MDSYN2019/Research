/* cyclic_io.h -- header file for cyclic array I/O functions
 *
 * See Chap 8, pp. 158 & ff in PPMPI
 */
#ifndef CYCLIC_IO_H
#define CYCLIC_IO_H
#include "mpi.h"
#include "cio.h"

#define ENTRIES_MAX 1024
#define LOCAL_ENTRIES_MAX 1024
#define PROC_MAX 64

typedef struct {
    MPI_Comm* comm;          /* Comm for collective ops */
#define Comm_ptr(array)       ((array)->comm)
#define Comm(array)           (*((array)->comm))

    int       p;             /* Size of array_comm      */
#define Comm_size(array)      ((array)->p)

    int       my_rank;       /* My rank in array_comm   */
#define Comm_rank(array)      ((array)->my_rank)

    int       global_order;  /* Global size of array    */
#define Order(array)          ((array)->global_order)

    int       padded_size;  /* Padded array size */
#define Padded_size(array)   ((array)->padded_size)

    float     entries[ENTRIES_MAX];
                            /* Elements of the array   */
#define Entries(array)        ((array)->entries)
#define Entry(array,i)        ((array)->entries[i])

    float     local_entries[LOCAL_ENTRIES_MAX]; 
                            /* Local elements of the array */
#define Local_entries(array)  ((array)->local_entries)
#define Local_entry(array,i)  ((array)->local_entries[i])

    int       local_size;    /* = n/p or n/p + 1        */
#define Local_size(array)        ((array)->local_size)

    int       stride;        /* Number of elements between two */
                             /* successive local entries */
#define Stride(array)         ((array)->stride)

    MPI_Datatype cyclic_mpi_t;
#define Type(array)           ((array)->cyclic_mpi_t)
} CYCLIC_ARRAY_STRUCT;

typedef CYCLIC_ARRAY_STRUCT* CYCLIC_ARRAY_T; 

void Initialize_params(
         MPI_Comm*      comm   /* in  */, 
         int            n      /* in  */, 
         CYCLIC_ARRAY_T  array  /* out */);

void Build_cyclic_type(
         MPI_Datatype* cyclic_mpi_t  /* out */,
         int           stride        /* in  */,
         int           array_size    /* in  */,
         int           p             /* in  */);
 
void Print_params(
         CYCLIC_ARRAY_T  array  /* in */);

void Print_entries(
         char*          title  /* in */,
         CYCLIC_ARRAY_T  array  /* in */); 

void Read_entries(
         char*          prompt  /* in */,
         CYCLIC_ARRAY_T  array   /* in */);
#endif
