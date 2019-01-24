/* stats.c -- Functions for use in gathering performance information on
 *     parallel tree search.  
 *
 * Note:  All times are local!
 *
 * See Chap 14, pp. 328 & ff, in PPMPI for a discussion of parallel tree
 *     search.
 */

#include "stats.h"
#include "mpi.h"
#include <stdio.h>
#include "cio.h"

STATS_T      stats = {0, 0, 0, 0, 0, 0.0, 0.0, 0.0};
double       overhead_time;
MPI_Datatype stats_mpi_t;

/********************************************************************/
void Get_overhead() {
    int i;

    stats.par_dfs_time = 0.0;
    for (i = 0; i < MAX_TESTS; i++) {
        Start_time(par_dfs_time);
        Finish_time(par_dfs_time);
    }
    overhead_time = stats.par_dfs_time/MAX_TESTS;
    stats.par_dfs_time = 0.0;
}  /* Get_overhead */


/********************************************************************/
void Build_stats_mpi_t(void) {
         int           i;
         int           block_lengths[MEMBERS];
         MPI_Datatype  types[MEMBERS];
         MPI_Aint      displacements[MEMBERS];
         MPI_Aint      start;
         MPI_Aint      address;

    for (i = 0; i < MEMBERS; i++)
        block_lengths[i] = 1;

    for (i = 0; i < INT_MEMBERS; i++) {
	types[i] = MPI_INT;
    }
    for (i = INT_MEMBERS; i < MEMBERS; i++) {
        types[i] = MPI_DOUBLE;
    }

    MPI_Address(&(stats.nodes_expanded), &start);
    displacements[0] = (MPI_Aint) 0;
    MPI_Address(&(stats.requests_sent), &address);
    displacements[1] = address - start;
    MPI_Address(&(stats.rejects_sent), &address);
    displacements[2] = address - start;
    MPI_Address(&(stats.work_sent), &address);
    displacements[3] = address - start;
    MPI_Address(&(stats.rejects_recd), &address);
    displacements[4] = address - start;
    MPI_Address(&(stats.work_recd), &address);
    displacements[5] = address - start;
    MPI_Address(&(stats.par_dfs_time), &address);
    displacements[6] = address - start;
    MPI_Address(&(stats.svc_req_time), &address);
    displacements[7] = address - start;
    MPI_Address(&(stats.work_rem_time), &address);
    displacements[8] = address - start;

    MPI_Type_struct(MEMBERS, block_lengths, displacements, types,
        &stats_mpi_t);
    MPI_Type_commit(&stats_mpi_t);
}  /* Build_stats_mpi_t */


/********************************************************************/
void Set_stats_to_zero(
         STATS_T*  stats  /* out */) {
    stats->nodes_expanded = 0;
    stats->requests_sent = 0;
    stats->rejects_sent = 0;
    stats->work_sent = 0;
    stats->rejects_recd = 0;
    stats->work_recd = 0;
    stats->par_dfs_time = 0.0;
    stats->svc_req_time = 0.0;
    stats->work_rem_time = 0.0;
}  /* Set_stats_to_zero */


/********************************************************************/
void Update_totals(
         STATS_T*  totals  /* in/out */,
         STATS_T*  new     /* in     */) {
    totals->nodes_expanded += new->nodes_expanded ;
    totals->requests_sent += new->requests_sent ;
    totals->rejects_sent += new->rejects_sent ;
    totals->work_sent += new->work_sent ;
    totals->rejects_recd += new->rejects_recd ;
    totals->work_recd += new->work_recd ;
    if (totals->par_dfs_time < new->par_dfs_time)
        totals->par_dfs_time = new->par_dfs_time;
    if (totals->svc_req_time < new->svc_req_time)
        totals->svc_req_time = new->svc_req_time;
    if (totals->work_rem_time < new->work_rem_time)
        totals->work_rem_time = new->work_rem_time;
}  /* Update_totals */


/********************************************************************/
void Print_title(void) {
printf("                Performance Statistics\n");
printf("     (Totals are sums for counts and maxima for times)\n");
printf("      Nodes  Reqs  Rejs  Work  Rejs  Work    DFS     Svc     Wk_rm\n");
printf("Proc   Exp   Sent  Sent  Sent  Recd  Recd    Time    Time     Time\n");
printf("----  -----  ----  ----  ----  ----  ----    ----    ----    -----\n");
}  /* Print_title */


/********************************************************************/
void Print_ind_stats(
         int       rank   /* in */,
         STATS_T*  stats  /* in */) {

    if (rank < 0)
        printf("Totals ");
    else
        printf(" %2d    ", rank);

    printf("%3d    %2d    %2d    %2d    %2d    %2d   ",
        stats->nodes_expanded, stats->requests_sent,
        stats->rejects_sent, stats->work_sent,
        stats->rejects_recd, stats->work_recd);
    printf("%7.2e %7.2e %7.2e\n",
        stats->par_dfs_time, stats->svc_req_time,
        stats->work_rem_time);
}  /* Print_ind_stats */


/********************************************************************/
void Print_stats(
         MPI_Comm  io_comm  /* in */) {
     STATS_T     recd_stats;
     STATS_T     totals;
     int         p;
     int         io_rank;
     int         my_rank;
     int         q;
     MPI_Status  status;

     MPI_Comm_size(io_comm, &p);
     MPI_Comm_rank(io_comm, &my_rank);
     Get_io_rank(io_comm, &io_rank);

     Build_stats_mpi_t();
     Set_stats_to_zero(&totals);

     if (my_rank == io_rank) {
        printf("\n");
        Print_title();
        for (q = 0; q < io_rank; q++) {
            MPI_Recv(&recd_stats, 1, stats_mpi_t, q, 0, io_comm, 
                &status);
            Print_ind_stats(q, &recd_stats);
            Update_totals(&totals, &recd_stats);
        } /* for */
        Print_ind_stats(io_rank, &stats);
        Update_totals(&totals, &stats);
        for (q = io_rank+1; q < p; q++) {
            MPI_Recv(&recd_stats, 1, stats_mpi_t, q, 0, io_comm, 
                &status);
            Print_ind_stats(q, &recd_stats);
            Update_totals(&totals, &recd_stats);
        } /* for */
        Print_ind_stats(-1, &totals);
     } else {
         MPI_Send(&stats, 1, stats_mpi_t, 0, 0, io_comm);
     }
}  /* Print_stats */
