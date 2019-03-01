#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <array>
#include <cassert>
#include "mpi.h"

// gsl libraries

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include "MPI_IO.hpp"

template <int m1, int l1, int t1, int m2, nt l2, int t2>
Physical<m1+m2,l1+l2,t1+t2> operator*(Physical<m1, l1, t1> lhs, Physical<m2, l2, t2> rhs) {
  return Physical<m1+m2, l1+l2, t1+t2>::unit*lhs.value() * rhs.value();
}

template <class T> class placeHolder {
public:
  placeHolder() { create(); }
  explicit placeHolder(std::size_t n, const T& val = T()) {create(n, val);}  
  // iterators
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef T value_type;
  typedef std::ptrdiff_t difference_type;
  typedef T& reference;
  typedef const T& const_reference;
};

typedef struct {
  double a;
  double b;
  int c;
} myStruct;

float A[10][10];


/*

Other derived types

If the data to be transmitted consists of a subset of the entries in an array, 
we shouldn't need to provide such detailed information ince all the information
have the same basic type

If we wish to send the third column of A, this won't work, since A[0][2] ...
aren't stored in contiguous memory locations. 

However, we can use MPI_type_vector to create a derived datatype, since the 
displacement of succesive elements of the derived types is constant
- A[1][2] is displaced 10 floats beyond A[0][2], A[2][2] 

It is fairly expensive to build a derived datatype. So applications that
make use of derived datatypes typically use the types many times.


*/


class mpiMatrixMultiplication {

  /*

    A communicator is a collection of processes that can send 
    messages to each other. A topology is a structure imposed 
    on the processes in a communicator that allow the processes
    in different ways 


    c_{ij} = a_{i,0}b_{0,j} + a_{i,1}b_{1,j} + a_{i,n-1}b_{n-1,j}
   */
  
public:
  
  
private:

};

MPI_Datatype newMPIDT2;  

MPI_Type_vector(10, 1, 10, MPI_FLOAT, &newMPIDT2);
MPI_Type_commit(&newMPIDT2);

/* 
   newMPIDT2 can be used to send any column of A. If we want to send the jth column of A,
   j = 0, 1, 2, 3, ..9, we simply call the communcation routine with the first argument 
   &(A[0][j]). 

   Also note that in fact column_mpi_t can be used to send any column of any 10 x 10 
   matrix of floats, since the stride and element type will be the same. 

  => Type Matching 

  At this point it is natural to ask, what are the rules for matching MPI datatypes
  For example, suppose a program contains the following code --> (2)

  
*/

// (2)

if (my_rank == 0) {
  MPI_Send(message, send_count, send_mpi_t, 1, 0, MPI_COMM_WORLD);
 } else if (my_rank == 1) {
  MPI_Recv(message, recv_count, recv_mpi_t, 0, 0, MPI_COMM_WORLD, &status);
 }
   



// --

if (my_rank == 0) {
  MPI_Send(&A[0][2], 1, column_mpi_t, 1, 0, MPI_COMM_WORLD);
 } else {
  MPI_Recv(&(A[0][2]), 1, column_mpi_t, 0, 0, MPI_COMM_WORLD, &status);
 }

// TODO
MPI_Group group_world;
MPI_Group first_row_group;
MPI_Comm first_row_comm;

/*

It is fairly expensive to build a derived datatype. So applications that make use of derived 
datatypes typically use the types many times

 */

void MPIGetData(float*  a_ptr, float*  b_ptr, int* n_ptr , int my_rank) {
    char  buffer[100];  /* Store data in buffer        */
    int   position;     /* Keep track of where data is */                            /*     in the buffer           */
    if (my_rank == 0){
      std::cout << "Enter a, b, and n\n";
        scanf("%f %f %d", a_ptr, b_ptr, n_ptr);
        /* Now pack the data into buffer.  Position = 0 */
        /* says start at beginning of buffer.           */
        position = 0;  
        /* Position is in/out */
        MPI_Pack(a_ptr, 1, MPI_FLOAT, buffer, 100, &position, MPI_COMM_WORLD);
        /* Position has been incremented: it now refer- */
        /* ences the first free location in buffer.     */
        MPI_Pack(b_ptr, 1, MPI_FLOAT, buffer, 100, &position, MPI_COMM_WORLD);
        /* Position has been incremented again. */
        MPI_Pack(n_ptr, 1, MPI_INT, buffer, 100, &position, MPI_COMM_WORLD);
        /* Position has been incremented again. */
        /* Now broadcast contents of buffer */
        MPI_Bcast(buffer, 100, MPI_PACKED, 0, MPI_COMM_WORLD);
    } else {
      MPI_Bcast(buffer, 100, MPI_PACKED, 0, MPI_COMM_WORLD);
        /* Now unpack the contents of buffer */
        position = 0;
        MPI_Unpack(buffer, 100, &position, a_ptr, 1, MPI_FLOAT, MPI_COMM_WORLD);
        /* Once again position has been incremented: */
        /* it now references the beginning of b.     */
        MPI_Unpack(buffer, 100, &position, b_ptr, 1, MPI_FLOAT, MPI_COMM_WORLD);
        MPI_Unpack(buffer, 100, &position, n_ptr, 1, MPI_INT, MPI_COMM_WORLD);
    }
}


void Build_mpi_type (double* a_p,  double* b_p, int* n_p, myStruct* stct, MPI_Datatype* input_mpi_t_p) {
  int array_of_blocklengths[3] = {1,1,1};
  MPI_Datatype array_of_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint array_of_displacements[3] = {0};
  MPI_Aint a_addr, b_addr, n_addr;

  stct->.a = *a_p;
  stct->.b = *b_p;
  stct->c = *n_p;


  //(*stct).b = *b_p;
  //(*stct).c = *n_p;
  
  MPI_Get_address(&stct->a, &a_addr);
  MPI_Get_address(&stct->b, &b_addr);
  MPI_Get_address(&stct->c, &n_addr);
  
  array_of_displacements[1] = b_addr - a_addr;
  array_of_displacements[2] = n_addr - a_addr;
  MPI_Type_create_struct(3, array_of_blocklengths, array_of_displacements, array_of_types, input_mpi_t_p);
  MPI_Type_commit(input_mpi_t_p);
}

void Get_Input(int my_rank, int comm_sz,  double* a_p,  double* b_p, int* n_p, myStruct* myStruct) {
  MPI_Datatype input_mpi_t;  
  Build_mpi_type(a_p, b_p, n_p, myStruct, &input_mpi_t);
  if (my_rank == 0) {
    MPI_Bcast(myStruct, 1, input_mpi_t, 0, MPI_COMM_WORLD);
  }
  MPI_Type_free(&input_mpi_t);
}

std::vector<int> bb;
myStruct exampleDatatype;
double a, b;
int c;

int main(int argc, char** argv) {
  int my_rank, comm_sz;
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  std::cout << "Enter a, b and n \n";
  scanf("%lf %lf %d", &a, &b, &c);
  Get_Input(0, comm_sz, &a, &b, &c, &exampleDatatype);

  for (int q = 1; q < comm_sz; q++) {
    std::cout << q << " " << exampleDatatype.a << " " << exampleDatatype.b << " " << exampleDatatype.c << std::endl;
  }

  MPI_Finalize();


  int idfdf;
   double xi, yi, x[10], y[10];

  std::cout << "#m=0,S=17\n";

  for (int idfdf = 0; idfdf < 10; idfdf++)
    {
      x[idfdf] = idfdf + 0.5 * sin (idfdf);
      y[idfdf] = idfdf + cos (idfdf * idfdf);
      std::cout << "%g %g\n" <<  x[idfdf] << y[idfdf];
    }

  std::cout << "#m=1,S=0\n";

  
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *spline = gsl_spline_alloc (gsl_interp_cspline, 10);
  gsl_spline_init (spline, x, y, 10);

  for (xi = x[0]; xi < x[9]; xi += 0.01) {
    yi = gsl_spline_eval (spline, xi, acc);
    std::cout << "%g %g\n" <<  xi <<  yi;
  }

  gsl_spline_free (spline);
  gsl_interp_accel_free (acc);
  
  return 0;
}
