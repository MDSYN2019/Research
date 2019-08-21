/* cio.c -- Collective functions for basic I/O operations
 *
 * See Chap 8, pp. 142 & ff in PPMPI
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "mpi.h"
#include "cio.h"

/* BUFSIZ is defined in stdio.h */
char         io_buf[BUFSIZ];  

/* Key identifying IO_Attribute */
int          IO_KEY = MPI_KEYVAL_INVALID;  

/* unused */
void*        extra_arg;

static int*  error_buf;
static int   error_bufsiz = 0;

/********************************************************/
/* Attempt to identify a process in io_comm that can be
 *     used for I/O.
 *
 * First see whether program defined io rank has been
 *     cached with either communicator.  If this fails
 *     try MPI defined io rank.
 *
 * Return values:
 *     1.  0: rank of I/O process cached with io_comm.
 *     2.  NO_IO_ATTR: couldn't find processor that could
 *         carry out I/O.  MPI_PROC_NULL cached with
 *         io_comm.
 *
 * Notes:
 *     1.  This is a collective operation, since function
 *         Copy_attr may use collective comm.
 *     2.  Only possible values cached are a valid process
 *         rank in comm2 or MPI_PROC_NULL.  (MPI_ANY_SOURCE
 *         won't be cached.)
 */
int Cache_io_rank(
        MPI_Comm   orig_comm       /* in     */,
        MPI_Comm   io_comm         /* in/out */) {

    int retval;   /* 0 or NO_IO_ATTR */

#ifdef DEBUG
    int my_rank;
    MPI_Comm_rank(orig_comm, &my_rank);

    printf("Process %d > In Cache_io_rank\n", my_rank);
#endif

    /* Check whether IO_KEY is defined.  If not, define */
    if (IO_KEY == MPI_KEYVAL_INVALID) {
        MPI_Keyval_create(MPI_DUP_FN, 
            MPI_NULL_DELETE_FN, &IO_KEY, extra_arg);
    } else if ((retval = Copy_attr(io_comm, io_comm, 
                IO_KEY)) != NO_IO_ATTR) {
        /* Value cached */
        return retval;
    } else if ((retval = Copy_attr(orig_comm, io_comm, 
                IO_KEY)) != NO_IO_ATTR) {
        /* Value cached */
        return retval;
    } 

#ifdef DEBUG
     printf("Process %d > In Cache_io_rank, give up on IO_KEY\n", 
             my_rank);
#endif
    /* Now see if we can find a value cached for MPI_IO */
    if ((retval = Copy_attr(orig_comm, io_comm, 
                MPI_IO)) != NO_IO_ATTR) {
        /* Value cached */
#ifdef DEBUG
       printf("Process %d > Copied attribute from orig to io\n",
               my_rank);
#endif
        return retval;
    } else if ((retval = Copy_attr(io_comm, io_comm, 
                MPI_IO)) != NO_IO_ATTR) {
        /* Value cached */
#ifdef DEBUG
        printf("Process %d > Copied attribute from io to io\n",
                my_rank);
#endif
        return retval;
    } 

#ifdef DEBUG
    printf("Process %d > In Cache_io_rank, returning at end\n", 
            my_rank);
#endif
    /* Couldn't find process that could carry out I/O */
    /* Copy_attr has cached MPI_PROC_NULL             */
    return NO_IO_ATTR;

}  /* Cache_io_rank */



/********************************************************/
/* Get attribute value associated with attribute key KEY
 *     in comm1, and cache with comm2 IO_KEY
 *
 * KEY can be either IO_KEY or MPI_IO.
 *
 * Return values:
 *     1.  0:  valid attribute successfully cached.
 *     2.  NO_IO_ATTR:  Couldn't find process that could
 *         carry out I/O.  MPI_PROC_NULL is cached with
 *         comm2.
 */ 
int Copy_attr(
        MPI_Comm   comm1   /* in     */,
        MPI_Comm   comm2   /* in/out */, 
        int        KEY     /* in     */) {

    int   io_rank;
    int   temp_rank;
    int*  io_rank_ptr;
    int   equal_comm;
    int   flag;

#ifdef DEBUG
    int my_rank;
    MPI_Comm_rank(comm1, &my_rank);
    printf("Process %d > In Copy_attr \n",
        my_rank);
#endif

    MPI_Attr_get(comm1, KEY, &io_rank_ptr, &flag);

#ifdef DEBUG
    if (flag == 0) {
        printf("Process %d > Attr_get failed\n", my_rank);
    } else if (*io_rank_ptr == MPI_ANY_SOURCE) {
        printf("Process %d > attr = MPI_ANY_SOURCE\n",
                my_rank);
        printf("Process %d > MPI_ANY_SOURCE = %d\n",
                my_rank, MPI_ANY_SOURCE);
    } else {
        printf("Process %d > attr = %d\n",
                my_rank, *io_rank_ptr);
    }
#endif

    if (flag == 0) {
        /* Attribute not cached with comm1 */
        io_rank_ptr = (int*) malloc(sizeof(int));
        *io_rank_ptr = MPI_PROC_NULL;
        MPI_Attr_put(comm2, IO_KEY, io_rank_ptr);
        return NO_IO_ATTR;
    } else if (*io_rank_ptr == MPI_PROC_NULL) {
        MPI_Attr_put(comm2, IO_KEY, io_rank_ptr);
        return NO_IO_ATTR;
    } else if (*io_rank_ptr == MPI_ANY_SOURCE) {
        /* Any process can carry out I/O.  Use */
        /* process 0                           */
#ifdef DEBUG
        printf("Process %d > Returning from Copy on MPI_ANY_SOURCE\n",
                my_rank);
#endif
        io_rank_ptr = (int*) malloc(sizeof(int));
        *io_rank_ptr = 0;
        MPI_Attr_put(comm2, IO_KEY, io_rank_ptr);
#ifdef DEBUG
        {
        int* temp_ptr;
        MPI_Attr_get(comm2, IO_KEY, &temp_ptr, &flag);
        if (flag == 0)
            printf("Process %d > In Copy, no value cached!\n", my_rank);
        else
            printf("Process %d > In Copy, cached io_rank = %d\n",
                    my_rank, *temp_ptr);
        }
#endif
        return 0;
    }

    /* Value in *io_rank_ptr is a valid process  */
    /* rank in comm1.  Action depends on whether */
    /* comm1 == comm2.                           */
    MPI_Comm_compare(comm1, comm2, &equal_comm);

    if (equal_comm == MPI_IDENT) {
        /* comm1 == comm2.  Valid value already */
        /* cached.  Do nothing.                 */
        return 0;
    } else {
        /* Check whether rank returned is valid */
        /* process rank in comm2                */
        Get_corresp_rank(comm1, *io_rank_ptr,
            comm2, &temp_rank);

        /* Different ranks may have been returned */
        /* on different processes.  Get min       */
        if (temp_rank == MPI_UNDEFINED)
            temp_rank = HUGE;
        MPI_Allreduce(&temp_rank, &io_rank, 1, MPI_INT,
            MPI_MIN, comm2);

        io_rank_ptr = (int*) malloc(sizeof(int));
        if (io_rank < HUGE) {
            *io_rank_ptr = io_rank;
            MPI_Attr_put(comm2, IO_KEY, io_rank_ptr);
            return 0;
        } else {
            /* No process got a valid rank in comm2 */
            /* from Get_corresp_rank                */
            *io_rank_ptr = MPI_PROC_NULL;
            MPI_Attr_put(comm2, IO_KEY, io_rank_ptr);
            return NO_IO_ATTR;
        }
    }
}  /* Copy_attr  */


/********************************************************/
/* Determines whether the process with rank rank1 in 
 *     comm1 is a valid rank in comm2.
 * If it is, it returns the rank in *rank2_ptr.  If it 
 *     isn't it returns MPI_UNDEFINED.
 *
 */
void Get_corresp_rank(
        MPI_Comm   comm1       /* in  */,
        int        rank1       /* in  */,
        MPI_Comm   comm2       /* in  */,
        int*       rank2_ptr   /* out */) {

    MPI_Group  group1;
    MPI_Group  group2;

    MPI_Comm_group(comm1, &group1);
    MPI_Comm_group(comm2, &group2);

    MPI_Group_translate_ranks(group1, 1, 
            &rank1, group2, rank2_ptr);

}  /* Get_corresp_rank */


/********************************************************/
/* Check whether IO_KEY is valid.  If it is, attempt to
 *     access it.  If it isn't attempt to define it from
 *     MPI_COMM_WORLD.
 *
 * Return values:
 *     1.  0: Valid rank returned.
 *     2.  NO_IO_ATTR:  Unable to find rank.
 */
int Get_io_rank(
        MPI_Comm io_comm      /* in  */,
        int*     io_rank_ptr  /* out */) {

    int*  temp_ptr;
    int   flag;

    if (IO_KEY == MPI_KEYVAL_INVALID) {
        MPI_Keyval_create(MPI_DUP_FN, 
            MPI_NULL_DELETE_FN, &IO_KEY, extra_arg);
    } else {
        MPI_Attr_get(io_comm, IO_KEY, &temp_ptr, &flag);
        if ((flag != 0) && (*temp_ptr != MPI_PROC_NULL)) {
            *io_rank_ptr = *temp_ptr;
            return 0;
        }
    }

    if (Copy_attr(MPI_COMM_WORLD, io_comm, MPI_IO) 
             == NO_IO_ATTR) {
        return NO_IO_ATTR;
    } else {
        MPI_Attr_get(io_comm, IO_KEY, &temp_ptr, &flag);
        *io_rank_ptr = *temp_ptr;
        return 0;
    }

}  /* Get_io_rank */


/********************************************************/
/* Prompt for input, read one line and broadcast.  
 *
 * Return values:
 *     1. 0:  input read
 *     2. NO_IO_ATTR:  no rank cached with IO_KEY.
 *
 * Notes:
 *     1. Prompt is significant only on IO_process 
 */
int Cscanf(
        MPI_Comm  io_comm  /* in  */,
        char*     prompt   /* in  */, 
        char*     format   /* in  */,
                  ...      /* out */) {

    va_list  args;
    int      my_io_rank;
    int      root;                            
/*
    int      count = 0;
    char     c;
    char*    ptr;
*/
    if (Get_io_rank(io_comm, &root) == NO_IO_ATTR)
        return NO_IO_ATTR;
    MPI_Comm_rank(io_comm, &my_io_rank);

    /* Read in data on root */
    if (my_io_rank == root) {
        printf("%s\n", prompt);
        gets(io_buf);
    }
/*
    ptr = io_buf;
    if (my_io_rank == root) {
        printf("%s\n", prompt);
        while ( (count < BUFSIZ - 1) &&
                ((c = getchar()) != EOF) &&
                (c != '\n') ) {
            *ptr = c;
            count = count + 1;
            ptr = ptr + 1;
        }
        *ptr = '\0';
    }  */ /* my_io_rank == root */

    /* Broadcast the input data */
    MPI_Bcast(io_buf, BUFSIZ, MPI_CHAR, root, io_comm);
            
    /* Copy the input data into the parameters */
    va_start(args, format);
    vsscanf(io_buf, format, args);
    va_end(args);
    
    return 0;
} /* Cscanf */


/********************************************************/
/* Prints data from all processes.  Format of data must
 *     be the same on each process.
 *
 * Return values:
 *     1.  0:  data printed
 *     2.  NO_IO_ATTR:  no rank cached with IO_KEY
 *
 * Notes:
 *     1.  Title is significant only on root.
 */
int Cprintf(
        MPI_Comm  io_comm  /* in */,
        char*     title    /* in */,
        char*     format   /* in */, 
                  ...      /* in */) {

    int         q;                            
    int         my_io_rank;
    int         io_p;
    int         root;
    MPI_Status  status;
    va_list     args;

    if (Get_io_rank(io_comm, &root) == NO_IO_ATTR)
        return NO_IO_ATTR;
    MPI_Comm_rank(io_comm, &my_io_rank);
    MPI_Comm_size(io_comm, &io_p);
    
    /* Send output data to io_process */
    if (my_io_rank != root) {
        /* Copy the output data into io_buf */
        va_start(args, format);
        vsprintf(io_buf, format, args);
        va_end(args);

        MPI_Send(io_buf, strlen(io_buf) + 1, MPI_CHAR, 
            root, 0, io_comm);
    } else { /* my_io_rank == root */
        printf("%s\n",title);
        fflush(stdout);
        for (q = 0; q < root; q++) {
            MPI_Recv(io_buf, BUFSIZ, MPI_CHAR, q,
                0, io_comm, &status);
            printf("Process %d > %s\n",q, io_buf);
            fflush(stdout);
        }

        /* Copy the output data into io_buf */
        va_start(args, format);
        vsprintf(io_buf, format, args);
        va_end(args);
        printf("Process %d > %s\n",root, io_buf);
        fflush(stdout);

        for (q = root+1; q < io_p; q++) {
            MPI_Recv(io_buf, BUFSIZ, MPI_CHAR, q,
                0, io_comm, &status);
            printf("Process %d > %s\n",q, io_buf);
            fflush(stdout);
        }
        printf("\n");
        fflush(stdout);
    }

    return 0;
} /* Cprintf */


/********************************************************/
/* Gathers error codes from all processes to all processes.
 *     If any error is negative, all processes abort.   
 *
 * Return values:
 *     1. 0:  no error detected.
 *     2. NO_IO_ATTR:  No valid rank cached with IO_KEY.
 *
 * Notes:
 *     1. "routine_name" only has significance on io_process.
 */ 
int Cerror_test(
        MPI_Comm  io_comm       /* in */,
        char*     routine_name  /* in */, 
        int       error         /* in */) {

    int q;
    int io_p;
    int error_count = 0;
    int io_process;
    int my_io_rank;

    if (Get_io_rank(io_comm, &io_process) == NO_IO_ATTR)
        return NO_IO_ATTR;
    MPI_Comm_size(io_comm, &io_p);
    MPI_Comm_rank(io_comm, &my_io_rank);

    if (error_bufsiz == 0) {
        error_buf = (int*) malloc(io_p*sizeof(int));
        error_bufsiz = io_p;
    } else if (error_bufsiz < io_p) {
        realloc(error_buf, io_p);
        error_bufsiz = io_p;
    }

    MPI_Allgather(&error, 1, MPI_INT, error_buf, 1, 
        MPI_INT, io_comm);
    for (q = 0; q < io_p; q++) {
        if (error_buf[q] < 0) {
            error_count++;
            if (my_io_rank == io_process) {
                fprintf(stderr,"Error in %s on process %d\n",
                    routine_name, q);
                fflush(stderr);
            }
        }
    }
    if (error_count > 0)
        MPI_Abort(MPI_COMM_WORLD, -1);

    return 0;
} /* Cerror_test */
