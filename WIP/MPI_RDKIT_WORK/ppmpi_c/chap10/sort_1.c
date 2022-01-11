/* sort_1.c -- level 1 version of sort program
 *
 * Input: none
 * Output:  messages indicating flow of control through program
 *
 * See Chap 10, pp. 226 & ff in PPMPI.
 */
#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "cio.h"
#include "sort_1.h"

int       p;
int       my_rank;
MPI_Comm  io_comm;

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

    Get_local_keys(&local_keys);
    Print_list(&local_keys);
    Redistribute_keys(&local_keys);
    Local_sort(&local_keys);
    Print_list(&local_keys);

    MPI_Finalize();
} /* main */

int Get_list_size() {
    Cprintf(io_comm,"","%s", "In Get_list_size");
    return 0;
} /* Get_list_size */

/* Return value negative indicates failure */
int Allocate_list(int list_size, 
    LOCAL_LIST_T* local_keys) {

    Cprintf(io_comm,"","%s", "In Allocate_key_list");
    return 0;
} /* Allocate_list */

void Get_local_keys(LOCAL_LIST_T* local_keys) {
    Cprintf(io_comm,"","%s", "In Get_local_keys");
} /* Get_local_keys */

void Redistribute_keys(LOCAL_LIST_T* local_keys) {
    Cprintf(io_comm,"","%s", "In Redistribute_keys");
} /* Redistribute_keys */

void Local_sort(LOCAL_LIST_T* local_keys) {
    Cprintf(io_comm,"","%s","In Local_sort");
} /* Local_sort */

void Print_list(LOCAL_LIST_T* local_keys) {
    Cprintf(io_comm,"","%s","In Print_list");
} /* Print_list */
