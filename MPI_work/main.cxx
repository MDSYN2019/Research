#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cassert>
#include <mpi.h>

// MPI headers

#include "MPI_IO.hpp"
#include "MPI_broadcast.hpp"
#include "MPI_functions.hpp"

std::vector<int> bb;

// Creates an array of random numbers. Each number has a value from 0 - 1
float *create_rand_nums(int num_elements) {
  float *rand_nums = (float *)malloc(sizeof(float) * num_elements);
  assert(rand_nums != NULL);
  int i;
  for (i = 0; i < num_elements; i++) {
    rand_nums[i] = (rand() / (float)RAND_MAX);
  }
  return rand_nums;
}


void Read_vector(double* local_a, int local_n, int n, std::string vec_name, int my_rank, MPI_Comm comm) {
  double* a = NULL;
  int i;
  if (my_rank == 0) {
    // TODO
  }

}

void Build_derived_type(
         float*         a_ptr           /* in   */,
         float*         b_ptr           /* in   */,
         int*           n_ptr           /* in   */,
         MPI_Datatype*  mesg_mpi_t_ptr  /* out  */);



void Build_derived_type(
         float*         a_ptr           /* in   */,
         float*         b_ptr           /* in   */,
         int*           n_ptr           /* in   */,
         MPI_Datatype*  mesg_mpi_t_ptr  /* out  */) {
                        /* pointer to new MPI type */

    /* The number of elements in each "block" of the   */
    /*     new type.  For us, 1 each.                  */
    int block_lengths[3];      

    /* Displacement of each element from start of new  */
    /*     type.  The "d_i's."                         */   
    /* MPI_Aint ("address int") is an MPI defined C    */
    /*     type.  Usually an int.                      */
    MPI_Aint displacements[3];

    /* MPI types of the elements.  The "t_i's."        */
    MPI_Datatype typelist[3];  
                              
    /* Use for calculating displacements               */
    MPI_Aint start_address; 
    MPI_Aint address;

    block_lengths[0] = block_lengths[1] 
                     = block_lengths[2] = 1;

    /* Build a derived datatype consisting of  */
    /* two floats and an int                   */
    typelist[0] = MPI_FLOAT;
    typelist[1] = MPI_FLOAT;
    typelist[2] = MPI_INT;

    /* First element, a, is at displacement 0      */
    displacements[0] = 0;

    /* Calculate other displacements relative to a */
    MPI_Address(a_ptr, &start_address);

    /* Find address of b and displacement from a   */
    MPI_Address(b_ptr, &address);
    displacements[1] = address - start_address;

    /* Find address of n and displacement from a   */
    MPI_Address(n_ptr, &address);
    displacements[2] = address - start_address;

    /* Build the derived datatype */
    MPI_Type_struct(3, block_lengths, displacements, 
        typelist, mesg_mpi_t_ptr);

    /* Commit it -- tell system we'll be using it for */
    /* communication.                                 */
    MPI_Type_commit(mesg_mpi_t_ptr);
} /* Build_derived_type */

void Get_data4(
         float*  a_ptr    /* out */, 
         float*  b_ptr    /* out */, 
         int*    n_ptr    /* out */,
         int     my_rank  /* in  */) {

    char  buffer[100];  /* Store data in buffer        */
    int   position;     /* Keep track of where data is */    
                        /*     in the buffer           */

    if (my_rank == 0){
        printf("Enter a, b, and n\n");
        scanf("%f %f %d", a_ptr, b_ptr, n_ptr);

        /* Now pack the data into buffer.  Position = 0 */
        /* says start at beginning of buffer.           */
        position = 0;  

        /* Position is in/out */
        MPI_Pack(a_ptr, 1, MPI_FLOAT, buffer, 100,
            &position, MPI_COMM_WORLD);
        /* Position has been incremented: it now refer- */
        /* ences the first free location in buffer.     */

        MPI_Pack(b_ptr, 1, MPI_FLOAT, buffer, 100,
            &position, MPI_COMM_WORLD);
        /* Position has been incremented again. */

        MPI_Pack(n_ptr, 1, MPI_INT, buffer, 100,
            &position, MPI_COMM_WORLD);
        /* Position has been incremented again. */

        /* Now broadcast contents of buffer */
        MPI_Bcast(buffer, 100, MPI_PACKED, 0,
            MPI_COMM_WORLD);
    } else {
        MPI_Bcast(buffer, 100, MPI_PACKED, 0,
            MPI_COMM_WORLD);

        /* Now unpack the contents of buffer */
        position = 0;
        MPI_Unpack(buffer, 100, &position, a_ptr, 1,
            MPI_FLOAT, MPI_COMM_WORLD);
	
        /* Once again position has been incremented: */
        /* it now references the beginning of b.     */

        MPI_Unpack(buffer, 100, &position, b_ptr, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(buffer, 100, &position, n_ptr, 1, MPI_INT, MPI_COMM_WORLD);
    }
} /* Get_data4 */



int main(int argc, char** argv) {
  
  MPI_BC aa;
  /*
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  
  if (rank == 0) {
    int value = 17;
    int result = MPI_Send(&value, 1, MPI_INT, 1, 0, MPI_COMM_WORLD);

    if (result == MPI_SUCCESS)
      std::cout << "Rank 0 OK!" << std::endl;
  } else if (rank == 1) {
    int value;
    int result = MPI_Recv(&value, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    if (result == MPI_SUCCESS && value == 17)
      std::cout << "Rank 1 OK!" << std::endl;
  }
  
  MPI_Finalize();
  */
  int         my_rank;   /* My process rank           */
  int         p;         /* The number of processes   */
  float       a;         /* Left endpoint             */
  float       b;         /* Right endpoint            */
  int         n;         /* Number of trapezoids      */
  float       h;         /* Trapezoid base length     */
  float       local_a;   /* Left endpoint my process  */
  float       local_b;   /* Right endpoint my process */
  int         local_n;   /* Number of trapezoids for  */
                           /* my calculation            */
  float       integral;  /* Integral over my interval */
  float       total;     /* Total integral            */
  int         source;    /* Process sending integral  */
  int         dest = 0;  /* All messages go to 0      */
  int         tag = 0;
  MPI_Status  status;
  
  void Get_data3(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank);
  float Trap(float local_a, float local_b, int local_n, float h);    /* Calculate local integral  */

    MPI_Init(&argc, &argv);    /* Let the system do what it needs to start up MPI */
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);    /* Get my process rank */
    MPI_Comm_size(MPI_COMM_WORLD, &p);     /* Find out how many processes are being used */
    Get_data3(&a, &b, &n, my_rank);

    h = (b-a)/n;    /* h is the same for all processes */
    local_n = n/p;  /* So is the number of trapezoids */

    /* Length of each process' interval of 
     * integration = local_n*h.  So my interval
     * starts at: */
    local_a = a + my_rank*local_n*h;
    local_b = local_a + local_n*h;
    integral = Trap(local_a, local_b, local_n, h);

    /* Add up the integrals calculated by each process */
    MPI_Reduce(&integral, &total, 1, MPI_FLOAT,
        MPI_SUM, 0, MPI_COMM_WORLD);

    /* Print the result */
    if (my_rank == 0) {
        printf("With n = %d trapezoids, our estimate\n", 
            n);
        printf("of the integral from %f to %f = %f\n", 
            a, b, total); 
    }

    /* Shut down MPI */
    MPI_Finalize();


  return 0;
}
