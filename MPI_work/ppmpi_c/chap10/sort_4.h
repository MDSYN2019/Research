/* sort_4.h -- header file for sort_4.c
 * 2. New definition of LOCAL_LIST_T and member access macros
 * 2. Add prototype for Insert_key 
 * 3. Complete definition of LOCAL_LIST_T
 * 3. Add new list macros and def of key_mpi_t
 * 3. Add prototype for Key_compare
 * 3. Add macro for LIST_BUF_SIZE and MAX_KEY_STRING
 * 4. Add prototype for Find_cutoff
 */
#ifndef SORT_H
#define SORT_H

#define KEY_MIN 0
#define KEY_MAX 32767
#define KEY_MOD 32768
#define LIST_BUF_SIZE 128
#define MAX_KEY_STRING 10

typedef int KEY_T;
typedef struct {
    int allocated_size;
    int local_list_size;
    KEY_T* keys;
} LOCAL_LIST_T;

/* Assume list is a pointer to a struct of type
 * LOCAL_LIST_T */
#define List_size(list) ((list)->local_list_size)
#define List_allocated_size(list) ((list)->allocated_size)
#define List(list) ((list)->keys)
#define List_free(list) {free List(list);}
#define List_key(list,i) (*((list)->keys + i))

#define key_mpi_t MPI_INT

int Get_list_size(void);
int Allocate_list(int list_size,
    LOCAL_LIST_T* local_keys);
void Get_local_keys(LOCAL_LIST_T* local_keys);
void Insert_key(KEY_T key, int i,  
    LOCAL_LIST_T* local_keys);
void Redistribute_keys(LOCAL_LIST_T* local_keys);
void Find_alltoall_send_params(LOCAL_LIST_T* local_keys, 
        int send_counts[], int send_displacements[]);
int Find_cutoff(int i);
void Find_recv_displacements(int recv_counts[],
         int recv_displacements[]);
void Local_sort(LOCAL_LIST_T* local_keys);
int Key_compare(const KEY_T* p, const KEY_T* q);
void Print_list(MPI_Comm io_comm, LOCAL_LIST_T* local_keys);
#endif
