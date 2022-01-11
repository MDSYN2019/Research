/* cyclic_io.c -- Functions for I/O of arrays using a cyclic distribution.
 *
 * See Chap 8, pp. 158 & ff in PPMPI
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "cio.h"
#include "cyclic_io.h"


/* 
 * Initialize all members except entries
 */
void Initialize_params(
         MPI_Comm*      comm   /* in  */,
         int            n      /* in  */,
         CYCLIC_ARRAY_T array  /* out */) {

    int   p;
    int   my_rank;
    int   q;
    int   quotient;
    int   remainder;

    Comm_ptr(array) = comm;

    MPI_Comm_size(*comm, &p);
    Comm_size(array) = p;

    MPI_Comm_rank(*comm, &my_rank);
    Comm_rank(array) = my_rank;

    Order(array) = n;

    quotient = n/p;
    remainder = n % p;

    if (remainder == 0)
        Padded_size(array) = n;
    else
        Padded_size(array) = p*(quotient+1);
    if (my_rank < remainder) {
        Local_size(array) = quotient+1;
    } else {
        Local_size(array) = quotient;
    }

    Stride(array) = p;

    Build_cyclic_type( &Type(array), Stride(array), 
        Padded_size(array), p);

}  /* Initialize_params */


/********************************************************/
void Build_cyclic_type(
         MPI_Datatype* cyclic_mpi_t  /* out */,
         int           stride        /* in  */,
         int           array_size    /* in  */,
         int           p             /* in  */) {

    MPI_Datatype  vector_mpi_t;
    int           blocksizes[2];
    MPI_Aint      displacements[2];
    MPI_Datatype  type_list[2];

    MPI_Type_vector(array_size/p, 1, stride, MPI_FLOAT,
        &vector_mpi_t);

    blocksizes[0] = blocksizes[1] = 1;
    displacements[0] = 0;
    displacements[1] = sizeof(float);
    type_list[0] = vector_mpi_t;
    type_list[1] = MPI_UB;

    MPI_Type_struct(2, blocksizes, displacements, type_list,
        cyclic_mpi_t);
    MPI_Type_commit(cyclic_mpi_t);
}  /* Build_cyclic_type */


/********************************************************/
void Print_params(
         CYCLIC_ARRAY_T array /* in */) {

    Cprintf(Comm(array),"p = ", "%d", Comm_size(array));
    Cprintf(Comm(array),"my_rank = ", "%d", Comm_rank(array));
    Cprintf(Comm(array),"order = ", "%d", Order(array));
    Cprintf(Comm(array),"padded size = ", "%d", 
        Padded_size(array));

    Cprintf(Comm(array),"my size","%d", Local_size(array)); 
    Cprintf(Comm(array),"stride","%d", Stride(array)); 

} /* Print_params */

/********************************************************/
/* Assumes that each process is using local_entries 
 *    member to store current contents of array.  If
 *    this is not the case, appropriate range of
 *    values from entries must be copied into 
 *    local_entries before call to MPI_Gather.
 */
void Print_entries(
         char*          title  /* in */,
         CYCLIC_ARRAY_T  array  /* in */) {

    int root;
    int q;
    int quotient;
    int remainder;
    int i, j, k;
    int send_size;

    Get_io_rank(Comm(array), &root);

    send_size = Padded_size(array)/Comm_size(array);
    MPI_Gather(Local_entries(array), send_size, MPI_FLOAT,
        Entries(array), 1, Type(array), root,
        Comm(array));

    if (Comm_rank(array) == root) {
        printf("%s\n", title);
        printf("    Processes\n");
        for (q = 0; q < Comm_size(array); q++)
            printf("%4d    ",q);
        printf("\n");
        for (q = 0; q < Comm_size(array); q++)
            printf("--------");
        printf("\n");

        quotient = Order(array)/Comm_size(array);
        remainder = Order(array) % Comm_size(array);
        k = 0;
        for (i = 0; i < quotient; i++) {
            for (j = 0; j < Comm_size(array); j++) {
                printf("%7.3f ",  Entry(array,k));
                k++;
            }
            printf("\n");
            fflush(stdout);
        }
        for (j = 0; j < remainder; j++) {
            printf("%7.3f ",  Entry(array,k));
            k++;
        }
        printf("\n");
        fflush(stdout);
    }
} /* Print_entries */


/********************************************************/
/* Reads values into local_entries member on each process.
 *     If values should go into entries member, it
 *     is necessary to add a loop to copy the values.
 */
void Read_entries(
         char*          prompt  /* in */,
         CYCLIC_ARRAY_T  array   /* in */) {

    int root;
    int i;
    int c;
    int recv_size;

    Get_io_rank(Comm(array), &root);

    if (Comm_rank(array) == root) {
        printf("%s\n",prompt);
        for (i = 0; i < Order(array); i++) {
            scanf("%f", &Entry(array,i));
        }
        /* Skip to end of line */
        while ((c = getchar()) != '\n');
        /* Fill padding with 0's */
        for (i = Order(array); i < Padded_size(array); i++)
            Entry(array,i) = 0.0;

    }

    recv_size = Padded_size(array)/Comm_size(array);
    MPI_Scatter(Entries(array), 1, Type(array),
        Local_entries(array), recv_size, MPI_FLOAT,
        root, Comm(array));

} /* Read_entries */
