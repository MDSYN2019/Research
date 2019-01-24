/* parallel_bitonic.c -- parallel bitonic sort of randomly generated list
 *     of integers
 *
 * Input:
 *     n: the global length of the list -- must be a power of 2.
 *
 * Output:
 *     The sorted list.
 *
 * Notes:
 *     1.  Assumes the number of processes p = 2^d and p divides n.
 *     2.  The lists are statically allocated -- size specified in MAX.
 *     3.  Keys are in the range 0 -- KEY_MAX-1.
 *     4.  Implementation can be made much more efficient by using
 *         pointers and avoiding re-copying lists in merges.
 *
 * See Chap 14, pp. 320 & ff. in PPMPI.
 */

#include <stdio.h>
#include <stdlib.h>
    /* Get rand and qsort */
#include "mpi.h"
#include "cio.h"

#define MAX 16384

#define LOW 0
#define HIGH 1

typedef int KEY_T;
#define KEY_MAX 32768
#define key_mpi_t MPI_INT

KEY_T temp_list[MAX]; /* buffer for keys received */
                      /* in Merge_split           */
KEY_T scratch_list[MAX]; /* temporary storage for */
                      /* merges                   */

void Generate_local_list(int list_size, KEY_T local_list[]);
void Print_list(char* title, int list_size, KEY_T local_list[],
         MPI_Comm io_comm);
void Local_sort(int list_size, KEY_T local_keys[]);
int Key_compare(const KEY_T* p, const KEY_T* q);
int log_base2(int x);
void Par_bitonic_sort_incr(int list_size, KEY_T local_list[],
        int proc_set_size, MPI_Comm comm);
void Par_bitonic_sort_decr(int list_size, KEY_T local_list[],
        int proc_set_size, MPI_Comm comm);
void Merge_split(int list_size, KEY_T local_list[],
        int which_keys, int partner, MPI_Comm  comm);
void Merge_list_low(int list_size, KEY_T list1[], KEY_T list2[]);
void Merge_list_high(int list_size, KEY_T list1[], KEY_T list2[]);


/********************************************************************/
main(int argc, char* argv[]) {
    int       list_size;         /* Local list size  */
    int       n;                 /* Global list size */
    KEY_T     local_list[MAX];
    int       proc_set_size;
    int       my_rank;
    int       p;
    unsigned  and_bit;
    MPI_Comm  io_comm;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);

    Cscanf(io_comm,"Enter the global list size","%d",&n);
    list_size = n/p;

    Generate_local_list(list_size, local_list);
/*
    Print_list("Before local sort", list_size, local_list, io_comm);
 */       
    Local_sort(list_size, local_list);
/*
    Print_list("After local sort", list_size, local_list, io_comm);
 */
 
    /* and_bit is a bitmask that, when "anded" with  */
    /* my_rank, tells us whether we're working on an */
    /* increasing or decreasing list                 */
    for (proc_set_size = 2, and_bit = 2;
         proc_set_size <= p;
         proc_set_size = proc_set_size*2, 
         and_bit = and_bit << 1)
        if ((my_rank & and_bit) == 0)
            Par_bitonic_sort_incr(list_size, 
                      local_list, proc_set_size, MPI_COMM_WORLD);
        else
            Par_bitonic_sort_decr(list_size, 
                      local_list, proc_set_size, MPI_COMM_WORLD);

    Print_list("After sort", list_size, local_list, io_comm);

    MPI_Finalize();
}  /* main */


/*********************************************************************/
void Generate_local_list(
         int    list_size     /* in  */, 
         KEY_T  local_list[]  /* out */) {
    int i;
    int my_rank;

    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    srand(my_rank);

    for (i = 0; i < list_size; i++)
        local_list[i] = rand() % KEY_MAX; 
}  /* Generate_local_list */ 


/*********************************************************************/
void Print_list(
         char*     title        /* in */,
         int       list_size    /* in */, 
         KEY_T     local_list[] /* in */,
         MPI_Comm  io_comm      /* in */) {
    int         i, q;
    int         p;
    int         my_rank;
    int         root;
    MPI_Status  status;

    MPI_Comm_size(io_comm, &p);
    MPI_Comm_rank(io_comm, &my_rank);
    Get_io_rank(io_comm, &root);

    if (my_rank == root) {
        printf("%s\n", title);
        for (q = 0; q < root; q++) {
            MPI_Recv(temp_list, list_size, key_mpi_t,
                q, 0, io_comm, &status);
            printf("Process %d > ",q);
            for (i = 0; i < list_size; i++)
                printf("%d ", temp_list[i]);
            printf("\n");
        }
        printf("Process %d > ", root);
        for (i = 0; i < list_size; i++)
            printf("%d ", local_list[i]);
        printf("\n");
        for (q = root+1; q < p; q++) {
            MPI_Recv(temp_list, list_size, key_mpi_t,
                q, 0, io_comm, &status);
            printf("Process %d > ",q);
            for (i = 0; i < list_size; i++)
                printf("%d ", temp_list[i]);
            printf("\n");
        }
    } else {
        MPI_Send(local_list, list_size, key_mpi_t,
            root, 0, io_comm);
    }
             
}  /* Print_list */


/*********************************************************************/
void Local_sort(
         int    list_size     /* in     */, 
         KEY_T  local_keys[]  /* in/out */) {
    
    qsort(local_keys, list_size, sizeof(KEY_T), 
        (int(*)(const void*, const void*))(Key_compare));
} /* Local_sort */


/*********************************************************************/
int Key_compare(const KEY_T* p, const KEY_T* q) {

    if (*p < *q)
        return -1;
    else if (*p == *q)
        return 0;
    else /* *p > *q */
        return 1;

}  /* Key_compare */


/********************************************************************/
int log_base2(int x) {
    int count = 0;

    while (x > 1) {
        x = x/2;
        count++;
    }

    return count;

}  /* log_base2 */


/********************************************************************/
void Par_bitonic_sort_incr(
        int       list_size      /* in     */, 
        KEY_T*    local_list     /* in/out */, 
        int       proc_set_size  /* in     */,
        MPI_Comm  comm           /* in     */ ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    int       my_rank;

    MPI_Comm_rank(comm, &my_rank);

    proc_set_dim = log_base2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);
    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank < partner)
            Merge_split(list_size, local_list, LOW,
                partner, comm);
        else
            Merge_split(list_size, local_list, HIGH,
                partner, comm);
        eor_bit = eor_bit >> 1;
    }
}  /* Par_bitonic_sort_incr */


/********************************************************************/
void Par_bitonic_sort_decr(
        int       list_size      /* in     */, 
        KEY_T*    local_list     /* in/out */, 
        int       proc_set_size  /* in     */,
        MPI_Comm  comm           /* in     */ ) {

    unsigned  eor_bit;
    int       proc_set_dim;
    int       stage;
    int       partner;
    int       my_rank;

    MPI_Comm_rank(comm, &my_rank);

    proc_set_dim = log_base2(proc_set_size);
    eor_bit = 1 << (proc_set_dim - 1);
    for (stage = 0; stage < proc_set_dim; stage++) {
        partner = my_rank ^ eor_bit;
        if (my_rank > partner)
            Merge_split(list_size, local_list, LOW,
                partner, comm);
        else
            Merge_split(list_size, local_list, HIGH,
                partner, comm);
        eor_bit = eor_bit >> 1;
    }

} /* Par_bitonic_sort_decr */


/********************************************************************/
void Merge_split(
        int       list_size     /* in     */, 
        KEY_T     local_list[]  /* in/out */, 
        int       which_keys    /* in     */,
        int       partner       /* in     */,
        MPI_Comm  comm          /* in     */ ) {

    MPI_Status status;
    
    /* key_mpi_t is an MPI (derived) type */
    MPI_Sendrecv(local_list, list_size, key_mpi_t, 
                 partner, 0, temp_list, list_size, 
                 key_mpi_t, partner, 0, comm, &status);
    if (which_keys == HIGH)
        Merge_list_high(list_size, local_list, 
                        temp_list);
    else
        Merge_list_low(list_size, local_list, 
                        temp_list);
} /* Merge_split */


/********************************************************************/
/* Merges the contents of the two lists. */
/* Returns the smaller keys in list1     */
void Merge_list_low(
        int    list_size  /* in     */,
        KEY_T  list1[]    /* in/out */,
        KEY_T  list2[]    /* in     */) {
    int  i;
    int  index1 = 0;
    int  index2 = 0;

    for (i = 0; i < list_size; i++)
        if (list1[index1] <= list2[index2]) {
            scratch_list[i] = list1[index1];
            index1++;
        } else {
            scratch_list[i] = list2[index2];
            index2++;
        }
    for (i = 0; i < list_size; i++)
        list1[i] = scratch_list[i];
        
}  /* Merge_list_low */


/********************************************************************/
/* Returns the larger keys in list 1.    */
void Merge_list_high(
        int    list_size  /* in     */,
        KEY_T  list1[]    /* in/out */,
        KEY_T  list2[]    /* in     */) {
    int  i;
    int  index1 = list_size - 1;
    int  index2 = list_size - 1;

    for (i = list_size - 1; i >= 0; i--)
        if (list1[index1] >= list2[index2]) {
            scratch_list[i] = list1[index1];
            index1--;
        } else {
            scratch_list[i] = list2[index2];
            index2--;
        }
    for (i = 0; i < list_size; i++)
        list1[i] = scratch_list[i];

}  /* Merge_list _high */

