/* stats.h -- definitions and declarations for stats.c
 */

#ifndef STATS_H
#define STATS_H
#include "mpi.h"

#define MAX_TESTS 100

#define MEMBERS 9
#define INT_MEMBERS 6
#define DOUBLE_MEMBERS 3

typedef struct {
    int       nodes_expanded;  /* Nodes popped from stack */
    int       requests_sent;
    int       rejects_sent;
    int       work_sent;
    int       rejects_recd;
    int       work_recd;
    double    par_dfs_time;
    double    svc_req_time;
    double    work_rem_time;
} STATS_T;

static double  start_time;
extern STATS_T stats;
extern double  overhead_time;

#define Start_time(time) {start_time = MPI_Wtime();}
#define Finish_time(time) {double finish, diff; \
                           finish = MPI_Wtime(); \
                           diff = finish - start_time - overhead_time; \
                           stats.time = stats.time + diff;}
#define Incr_stat(member) {(stats.member)++;}

void Get_overhead(void);
void Build_stats_mpi_t(void);
void Set_stats_to_zero(STATS_T* stats);
void Update_totals(STATS_T* totals, STATS_T* new);
void Print_title(void);
void Print_ind_stats(int rank, STATS_T* stats);
void Print_stats(MPI_Comm io_comm);
#endif
