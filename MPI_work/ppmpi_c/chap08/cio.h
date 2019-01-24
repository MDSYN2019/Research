/* cio.h -- header file for use with cio.c -- basic collective I/O functions
 *
 * See Chap 8, pp. 142 & ff in PPMPI
 */
#ifndef IO_H
#define IO_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include "mpi.h"
#include "vsscanf.h"

/* All process ranks < HUGE */
#define HUGE 32768
#define NO_IO_ATTR -1

extern int IO_KEY;

int Cache_io_rank(
        MPI_Comm   orig_comm   /* in     */,
        MPI_Comm   io_comm     /* in/out */);

int Copy_attr(
        MPI_Comm   comm1       /* in     */,
        MPI_Comm   comm2       /* in/out */,
        int        KEY         /* in     */);

void Get_corresp_rank(
        MPI_Comm   comm1       /* in  */,
        int        rank1       /* in  */,
        MPI_Comm   comm2       /* in  */,
        int*       rank2_ptr   /* out */);

int Get_io_rank(
        MPI_Comm io_comm      /* in  */,
        int*     io_rank_ptr  /* out */);

int Cscanf(
        MPI_Comm  io_comm  /* in  */,
        char*     prompt   /* in  */,
        char*     format   /* in  */,
                  ...      /* out */);

int Cprintf(
        MPI_Comm  io_comm  /* in */,
        char*     title    /* in */,
        char*     format   /* in */,
                  ...      /* in */);

int Cerror_test(
        MPI_Comm  io_comm       /* in */,
        char*     routine_name  /* in */,
        int       error         /* in */);
#endif
/* End of io.h */
