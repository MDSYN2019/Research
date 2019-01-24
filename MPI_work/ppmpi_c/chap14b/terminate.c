/* terminate.c -- functions for determining whether every process has
 *     exhausted its stack.
 *
 * Here's the basic idea.  At the start of the program, process 0 has
 * one unit of some indestructible quantity -- I called it energy.
 * When the data is distributed, the energy is divided into p
 * equal parts and also distributed among the processes.  So each
 * process has 1/p of the energy.  (I don't actually send the energy
 * out initially, since every process knows it will get 1/p units of
 * energy, each can just set its energy to 1/p).  Whenever
 * a process exhausts its local stack, it sends whatever energy
 * it has to process 0.  Whenever a process satisfies a request
 * for work, it splits its energy in two equal pieces, keeping half
 * for itself and sending half to the process receiving the work.
 * Each process keeps track of its available energy in "my_energy".
 * Process 0, also keeps track of returned energy in "returned_energy."
 * Since energy is never destroyed, the sum of all energy on all the
 * processes will always be 1.  When process 0 finds that
 * "returned_energy" is 1, no process (including itself) will have
 * any work left, and it can broadcast a termination message to all
 * processes.  The only catch here is that floating point
 * arithmetic isn't exact.  So there are functions for
 * rational arithmetic:  a fraction a/b is represented as a pair of
 * longs (a,b).  Dividing by 2 changes (a,b) to (a,2b).  Adding
 * (a,b) + (c,d) = (ad + bc,bd).  This can still cause problems,
 * but there are no security checks.  Future work, no doubt.
 *
 * See Chap 14, pp. 334 & ff, in PPMPI.
 */
#include "mpi.h"
#include "node_stack.h"
#include <stdio.h>
#include <stdlib.h>
#include "terminate.h"

MPI_Datatype rational_mpi_t;
extern int p;
extern int my_rank;

static DIVISOR_T divs;
static int prime_list[MAX_PRIMES] = 
    {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97};

static RATIONAL_T ONE = {1,1};
static RATIONAL_T my_energy; 
static RATIONAL_T returned_energy = {0,1};  /* significant only on proc 0 */



/*********************************************************************/
void Print_divisors(void) {
    int i;
    printf("Process %d > Divisors of %d are ", my_rank, 2*p);
    for (i = 0; i < Num_divisors(divs); i++)
        printf("%d ", Divisor(divs,i));
    printf("\n");
}  /* Print_divisors */


/*********************************************************************/
void Setup_term_detect(void) {

    Find_divisors(2*p);  /* Use 2*p to force inclusion of 2 in list */
                         /*     of divisors                         */
    my_energy.numerator = 1;
    my_energy.denominator = p;

    Build_rational_mpi_t();


}  /* Setup_term_detect */
        

/*********************************************************************/
void Build_rational_mpi_t(void) {
    int           count = 2;
    int           block_lengths[2];
    MPI_Datatype  types[2];
    MPI_Aint      displacements[2];
    MPI_Aint      start;
    MPI_Aint      address;

    block_lengths[0] = block_lengths[1] = 1;
    types[0] = types[1] = MPI_LONG;

    MPI_Address(&my_energy, &start);
    MPI_Address(&(my_energy.numerator), &address);
    displacements[0] = address - start;
    MPI_Address(&(my_energy.denominator), &address);
    displacements[1] = address - start;

    MPI_Type_struct(count, block_lengths, displacements, types,
        &rational_mpi_t);
    MPI_Type_commit(&rational_mpi_t);
}  /* Build_rational_mpi_t */


/*********************************************************************/
void Find_divisors(int x) {
    int quotient = x;
    int remainder;
    int prime_indx = 0;
    int divisor_indx = 0;
    int prime;

    while(quotient != 1) {
        prime = prime_list[prime_indx];
        remainder = quotient % prime;
        if (remainder == 0) {
            Divisor(divs,divisor_indx) = prime;
            divisor_indx++;
            do {
                quotient = quotient/prime;
                remainder = quotient % prime;
            } while (remainder == 0);
        }
        prime_indx++;
        if ((quotient != 1) && (prime_indx >= MAX_PRIMES)) {
            Print_divisors();
            fprintf(stderr, "x = %d has too many distinct divisors\n", x);
            fprintf(stderr, "Increase the size of prime_list\n");
            exit(-1);
        }
    }

    Num_divisors(divs) = divisor_indx;
        
} /* Find_divisors */


/*********************************************************************/
void Reduce(RATIONAL_T* y) {
    int   i;
    int   divisor;

    for (i = 0; i < Num_divisors(divs); i++) {
        divisor = Divisor(divs,i);
        while (((Numerator(y) % divisor) == 0) && 
               (Denominator(y) % divisor) == 0) {
            Numerator(y) = Numerator(y)/divisor;
            Denominator(y) = Denominator(y)/divisor;
        }           
    }
}  /* Reduce */


/*********************************************************************/
/* Add x = x + y */
void Add(RATIONAL_T* x, RATIONAL_T* y) {

    Numerator(x) = Denominator(y)*Numerator(x) + Denominator(x)*Numerator(y);
    Denominator(x) = Denominator(x)*Denominator(y);
    Reduce(x);
}  /* Add */


/*********************************************************************/
void Divide_by_2(
         RATIONAL_T*  x  /* in/out */) {

    if ((Numerator(x) % 2) == 0)
        Numerator(x) = Numerator(x)/2;
    else
        Denominator(x) = 2*Denominator(x);
}  /* Divide_by_2 */


/*********************************************************************/
int Equal(RATIONAL_T* x, RATIONAL_T* y) {

    if ( Numerator(x)*Denominator(y) == 
             Numerator(y)*Denominator(x) )
        return TRUE;
    else
        return FALSE;
}  /* Equal */


/*********************************************************************/
void Send_half_energy(
         int       destination  /* in */,
         MPI_Comm  comm         /* in */) {

    Divide_by_2(&my_energy);

    MPI_Send(&my_energy, 1, rational_mpi_t, destination,
        HALF_ENERGY_TAG, comm);
}  /* Send_half_energy */


/*********************************************************************/
void Recv_half_energy(
         int       source  /* in  */,
         MPI_Comm  comm    /* out */) {
     MPI_Status  status;
     
     MPI_Recv(&my_energy, 1, rational_mpi_t, source,
         HALF_ENERGY_TAG, comm, &status);
}  /* Recv_half_energy */


/*********************************************************************/
void Return_energy(MPI_Comm comm) {

    MPI_Send(&my_energy, 1, rational_mpi_t, 0, RETURN_ENERGY_TAG, comm);
    
}  /* Return_energy */


/*********************************************************************/
void Receive_returned_energy(
         MPI_Comm  comm  /* in */) {
    int         done = FALSE;
    int         source;
    int         energy_returned;
    RATIONAL_T  energy;
    MPI_Status  status;

    while (!done) {
        MPI_Iprobe(MPI_ANY_SOURCE, RETURN_ENERGY_TAG, comm,
            &energy_returned, &status);
        if (energy_returned) {
            source = status.MPI_SOURCE;
            MPI_Recv(&energy, 1, rational_mpi_t, source, 
                RETURN_ENERGY_TAG, comm, &status);
            Add(&returned_energy, &energy);
        } else {
            done = TRUE;
        }
    }  /* while */
}  /* Receive_returned_energy */


/*********************************************************************/
int Search_complete(
        MPI_Comm  comm  /* in */) {
    int         dest;
    int         x = 0;
    int         completed;
    MPI_Status  status;
    
    if (my_rank == 0) {
        Receive_returned_energy(comm);
        if (Equal(&returned_energy, &ONE)) {
            for (dest = 1; dest < p; dest++) {
                MPI_Send(&x, 1, MPI_INT, dest, COMPLETE_TAG, comm);
            }
            return TRUE;
        } else {
            return FALSE;
        }
    } else { /* my_rank != 0 */
        MPI_Iprobe(0, COMPLETE_TAG, comm, &completed, &status);
        if (completed) {
            MPI_Recv( &x, 1, MPI_INT, 0, COMPLETE_TAG, comm, &status);
            return TRUE;
        } else {
            return FALSE;
        }
    }

}  /* Search_complete */
