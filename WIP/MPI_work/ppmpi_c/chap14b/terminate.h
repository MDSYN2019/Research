/* terminate.h
 *
 * Definitions and declarations for distributed termination detection
 */
#ifndef TERMINATE_H
#define TERMINATE_H

#include "mpi.h"

#define RETURN_ENERGY_TAG 1000
#define HALF_ENERGY_TAG   2000
#define COMPLETE_TAG      3000

/* Number of primes stored in prime_list */
/*     Currently number of primes < 100  */
#define MAX_PRIMES 25

typedef struct {
    int   num_divisors;
#define Num_divisors(divs) ((divs).num_divisors)
    int   divisor_list[MAX_PRIMES];
#define Divisor_list(divs) ((divs).divisor_list)
#define Divisor(divs, i)   (*(((divs).divisor_list)+i))
} DIVISOR_T;

typedef struct{
    long numerator;
#define Numerator(y) ((y)->numerator)
    long denominator;
#define Denominator(y) ((y)->denominator)
} RATIONAL_T;

extern MPI_Datatype rational_mpi_t;

void  Print_divisors(void);
void  Setup_term_detect(void);
void  Build_rational_mpi_t(void);
void  Find_divisors(int x);
void  Reduce(RATIONAL_T* y);
void  Add(RATIONAL_T* x, RATIONAL_T* y);
void  Divide_by_2(RATIONAL_T* x);
int   Equal(RATIONAL_T* x, RATIONAL_T* y);

void  Return_energy(MPI_Comm comm);
void  Send_half_energy(int destination, MPI_Comm comm);
void  Recv_half_energy(int source, MPI_Comm comm);
void  Receive_returned_energy(MPI_Comm comm);
int   Search_complete(MPI_Comm comm);

#endif
