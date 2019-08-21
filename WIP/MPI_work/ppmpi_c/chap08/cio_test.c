/* cio_test.c -- program for testing the functions in cio.c
 *
 * Compile with Makefile.cio
 *
 * Input: An int, a float, and a string, numerous times
 *
 * Output: Information on the communicators and input.
 *
 * See Chap 8, pp. 142 & ff in PPMPI
 */
#include <stdio.h>
#include "mpi.h"
#include "cio.h"

void Print_attr(
         char* source_comm_name, 
         char* new_comm_name,
         MPI_Comm new_comm);

main(int argc, char* argv[]) {
    MPI_Comm io_comm_world;
    MPI_Comm even_comm;
    MPI_Comm odd_comm;
    MPI_Comm duped_comm;
    int split_key;
    int p;
    int my_rank;
    int ival;
    float fval;
    char sval[100];
    int ret_val;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    printf("Process %d > IO_KEY = %d\n", my_rank, IO_KEY);
    MPI_Comm_dup(MPI_COMM_WORLD, &io_comm_world);
    Cache_io_rank(MPI_COMM_WORLD, io_comm_world);

    printf("Process %d > IO_KEY = %d\n", my_rank, IO_KEY);
    Print_attr("MPI_COMM_WORLD", "io_comm_world", io_comm_world);

    split_key = my_rank % 2;
    if (split_key == 0) {
        MPI_Comm_split(io_comm_world, split_key, my_rank, &even_comm);
        Cache_io_rank(io_comm_world, even_comm);
        Print_attr("io_comm_world", "even_comm", even_comm);
        Cache_io_rank(MPI_COMM_WORLD, even_comm);
        Print_attr("MPI_COMM_WORLD", "even_comm", even_comm);
        MPI_Comm_dup(even_comm, &duped_comm);
        Print_attr("even_comm", "duped_comm", duped_comm);
    } else {
        MPI_Comm_split(io_comm_world, split_key, my_rank, &odd_comm);
        Cache_io_rank(io_comm_world, odd_comm);
        Print_attr("io_comm_world", "odd_comm", odd_comm);
        Cache_io_rank(MPI_COMM_WORLD, odd_comm);
        Print_attr("MPI_COMM_WORLD", "odd_comm", odd_comm);
        MPI_Comm_dup(odd_comm, &duped_comm);
        Print_attr("odd_comm", "duped_comm", duped_comm);
    }
    fflush(stdout);
    
    ret_val = Cscanf(MPI_COMM_WORLD, "Enter an int, a float, and a string",
           "%d %f %s", &ival, &fval, &sval);
    ret_val = Cprintf(MPI_COMM_WORLD, "MPI_COMM_WORLD read:","%d %f %s", 
           ival, fval, sval);

    ret_val = Cscanf(io_comm_world, "Enter an int, a float, and a string",
           "%d %f %s", &ival, &fval, &sval);
    ret_val = Cprintf(io_comm_world, "io_comm_world read:","%d %f %s", 
           ival, fval, sval);
    if (split_key == 0) {
        ret_val = Cscanf(even_comm, "Enter an int, a float, and a string",
               "%d %f %s", &ival, &fval, &sval);
        ret_val = Cprintf(even_comm, "even_comm read:","%d %f %s", 
               ival, fval, sval);
        ret_val = Cscanf(duped_comm, "Enter an int, a float, and a string",
           "%d %f %s", &ival, &fval, &sval);
        ret_val = Cprintf(duped_comm, "duped_comm read:","%d %f %s", 
           ival, fval, sval);
    } else {
        ret_val = Cscanf(odd_comm, "Enter an int, a float, and a string",
           "%d %f %s", &ival, &fval, &sval);
        ival = 2*ival; fval = 2*fval;
        ret_val = Cprintf(odd_comm, "odd_comm read:","%d %f %s", 
           ival, fval, sval);
        ret_val = Cscanf(duped_comm, "Enter an int, a float, and a string",
           "%d %f %s", &ival, &fval, &sval);
        ival = 2*ival; fval = 2*fval;
        ret_val = Cprintf(duped_comm, "duped_comm read:","%d %f %s", 
           ival, fval, sval);
    }
    ret_val = Cerror_test(MPI_COMM_WORLD,"main MPI_COMM_WORLD",0);
    ret_val = Cerror_test(io_comm_world,"main io_comm_world",0);
    if (split_key == 0) {
        ret_val = Cerror_test(even_comm,"main even_comm",0);
        ret_val = Cerror_test(duped_comm,"main duped_comm",0);
    } else {
        ret_val = Cerror_test(odd_comm,"main odd_comm",0);
        ret_val = Cerror_test(duped_comm,"main duped_comm",0);
    }
    ret_val = Cerror_test(io_comm_world,"main io_comm_world", -1);

    MPI_Finalize();
} /* main */

void Print_attr(
         char*     source_comm_name,
         char*     new_comm_name,
         MPI_Comm  new_comm) {

    int*  io_rank_ptr;
    int   flag;
    int   my_rank;
  
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

    MPI_Attr_get(new_comm, IO_KEY, &io_rank_ptr, &flag);
    if (flag == 0) {
        printf("Process %d > No attribute associated to IO_KEY in %s from %s\n",
                my_rank, new_comm_name, source_comm_name);
    } else {
        printf("Process %d > io_rank = %d in %s (received from %s)\n", 
                my_rank, *io_rank_ptr, new_comm_name, source_comm_name);
    }
} /* Print_attr */
