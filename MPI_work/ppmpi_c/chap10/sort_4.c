/* sort_4.c -- level 4 version of sort program
 * 2.  Add fake Get_list_size 
 * 2.  Return scalars in Allocate_list
 * 2.  Add fake Get_local_keys
 * 3.  Definition of Redistribute keys.
 * 3.  Finish Allocate_list.
 * 3.  Add Insert.
 * 3.  Add Local_sort.
 * 3.  Add Print_list.
 * 4.  Add Find_alltoall_send_params.
 * 4.  Add Find_cutoff
 * 4.  Add Find_recv_displacements
 * 4.  Allow input for list size in Get_list_size
 *
 * Input: 
 *     list_size: global size of list to be sorted.
 *
 * Output: contents of list before and after sorting.
 *
 * See Chap 10, pp. 226 & ff, esp. pp. 236 & ff., in PPMPI.
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "cio.h"
#include "sort_4.h"

int       p;
int       my_rank;
MPI_Comm  io_comm;

/*********************************************************************/
main(int argc, char* argv[]) {
    LOCAL_LIST_T  local_keys;
    int           list_size;
    int           error;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm);
    Cache_io_rank(MPI_COMM_WORLD, io_comm);

    list_size = Get_list_size();

    /* Return negative if Allocate failed */
    error = Allocate_list(list_size, &local_keys);
    if (error < 0) {
        printf("Process %d > Can't allocate list\n", my_rank);
        printf("Process %d > Quitting\n", my_rank);
        MPI_Abort(MPI_COMM_WORLD, -1);
    }

    Get_local_keys(&local_keys);
    Print_list(io_comm, &local_keys);

    Redistribute_keys(&local_keys);

    Local_sort(&local_keys);
    Print_list(io_comm, &local_keys);

    MPI_Finalize();
} /* main */


/*********************************************************************/
int Get_list_size(void) {
    int size;

    Cscanf(io_comm,"How big is the list?","%d",&size);
    return size;
} /* Get_list_size */


/*********************************************************************/
/* Return value negative indicates failure */
int Allocate_list(
        int           list_size  /* in  */, 
        LOCAL_LIST_T* local_keys /* out */) {

    List_allocated_size(local_keys) = list_size/p;
    List_size(local_keys) = list_size/p;
    List(local_keys) = (KEY_T*) 
        malloc(List_allocated_size(local_keys)*sizeof(KEY_T));
    if (List(local_keys) == (KEY_T*) NULL)
        return -1;
    else
        return 0;
} /* Allocate_list */


/*********************************************************************/
void Get_local_keys(LOCAL_LIST_T* local_keys) {
    int i;

    /* Seed the generator */
    srand(my_rank);

    for (i = 0; i < List_size(local_keys);  i++)
        Insert_key(rand() % KEY_MOD, i, local_keys);
} /* Get_local_keys */


/*********************************************************************/
void Insert_key(KEY_T key, int i, 
    LOCAL_LIST_T* local_keys) {
    List_key(local_keys, i) = key;
} /* Insert_key */


/*********************************************************************/
void Redistribute_keys(
         LOCAL_LIST_T* local_keys  /* in/out */) {
    int new_list_size, i, error = 0;
    int* send_counts;
    int* send_displacements;
    int* recv_counts;
    int* recv_displacements;
    KEY_T* new_keys;
    
    /* Allocate space for the counts and displacements */
    send_counts = (int*) malloc(p*sizeof(int));
    send_displacements = (int*) malloc(p*sizeof(int));
    recv_counts = (int*) malloc(p*sizeof(int));
    recv_displacements = (int*) malloc(p*sizeof(int));

    Local_sort(local_keys);
    Find_alltoall_send_params(local_keys, 
        send_counts, send_displacements);

    /* Distribute the counts */
    MPI_Alltoall(send_counts, 1, MPI_INT, recv_counts,
        1, MPI_INT, MPI_COMM_WORLD);

    /* Allocate space for new list */
    new_list_size = recv_counts[0];
    for (i = 1; i < p; i++)
        new_list_size += recv_counts[i];
    new_keys = (KEY_T*) 
        malloc(new_list_size*sizeof(KEY_T));

    Find_recv_displacements(recv_counts, recv_displacements);

    /* Exchange the keys */
    MPI_Alltoallv(List(local_keys), send_counts, 
        send_displacements, key_mpi_t, new_keys, 
        recv_counts, recv_displacements, key_mpi_t,
        MPI_COMM_WORLD);

    /* Replace old list with new list */
    List_free(local_keys);
    List_allocated_size(local_keys) = new_list_size;
    List_size(local_keys) = new_list_size;
    List(local_keys) = new_keys;

    /* Free temporary storage */
    free(send_counts);
    free(send_displacements);
    free(recv_counts);
    free(recv_displacements);

} /* Redistribute_keys */


/******************************************************************/
void Find_alltoall_send_params(
         LOCAL_LIST_T* local_keys         /* in  */,
         int*          send_counts        /* out */, 
         int*          send_displacements /* out */) {
    KEY_T cutoff;
    int i, j;

    /* Take care of process 0 */
    j = 0;
    send_displacements[0] = 0;
    send_counts[0] = 0;
    cutoff = Find_cutoff(0);
    /* Key_compare > 0 if cutoff > key */
    while ((j < List_size(local_keys)) &&
           (Key_compare(&cutoff,&List_key(local_keys,j)) 
                > 0)) {
        send_counts[0]++;
        j++;
    }

    /* Now deal with the remaining processes */
    for (i = 1; i < p; i++) {
        send_displacements[i] = 
            send_displacements[i-1] + send_counts[i-1];
        send_counts[i] = 0;
        cutoff = Find_cutoff(i);
        /* Key_compare > 0 if cutoff > key */
        while ((j < List_size(local_keys)) &&
               (Key_compare(&cutoff,&List_key(local_keys,j)) 
                   > 0)) {
            send_counts[i]++;
            j++;
        }
    }
} /* Find_alltoall_send_params */


/*********************************************************************/
int Find_cutoff(int i) {
    return (i+1)*(KEY_MAX + 1)/p;
}  /* Find_cutoff */


/*********************************************************************/
void Find_recv_displacements(int recv_counts[],
    int recv_displacements[]){
    int i;

    recv_displacements[0] = 0;
    for (i = 1; i < p; i++)
        recv_displacements[i] =
            recv_displacements[i-1]+recv_counts[i-1];
} /* Find_recv_displacements */


/*********************************************************************/
void Local_sort(LOCAL_LIST_T* local_keys) {
    
    qsort(List(local_keys), List_size(local_keys), sizeof(KEY_T), 
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


/*********************************************************************/
void Print_list(MPI_Comm io_comm, LOCAL_LIST_T* local_keys) {
    char list_string[LIST_BUF_SIZE];
    char key_string[MAX_KEY_STRING];
    int  i;

    list_string[0] = '\0';
    for (i = 0; i < List_size(local_keys); i++) {
        sprintf(key_string,"%d ", List_key(local_keys,i));
        strcat(list_string, key_string);
    }
    Cprintf(io_comm,"Contents of the list", "%s", list_string);
} /* Print_list */
