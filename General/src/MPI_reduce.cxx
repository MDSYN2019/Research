#include <iostream>
#include <sstream>
#include <cstdio>
#include <cassert>
#include <vector>
#include <mutex>
#include <string>
#include <mpi.h>

class MPI_sorting_methods {
public:
  MPI_sorting_methods();
  ~MPI_sorting_methods();
  void Bubble_sort(int a[], int n) {
    for (list_length = n; list_length >= 2; list_length--) {
      for (i = 0; i < list_length - 1; i++) {
	if (a[i] > a[i+1]) {
	  temp = a[i];
	  a[i] = a[i+1];
	  a[i+1] = temp;
	}
      }
    }
  }
private:
  int list_length, i, temp;
};

void Get_input(int my_rank, int comm_sz, double * a_p, double * b_p, int * n_p ) {
  int dest;
  if (my_rank == 0) {
    std::cout << "Enter a, b and n " << std::endl;
    // scanf
    for (dest = 1; dest < comm_sz; dest++) {
      MPI_Send(a_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD); // send all values from not 0 to 0 
      MPI_Send(b_p, 1, MPI_DOUBLE, dest, 0, MPI_COMM_WORLD); // ditto 
      MPI_Send(n_p, 1, MPI_INT, dest, 0, MPI_COMM_WORLD);  // ditto 
    }
  }
  else {
    MPI_Recv(a_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // if process is 0, receive 
    MPI_Recv(b_p, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // ditto
    MPI_Recv(n_p, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // ditto again
  }
}

void Vector_sum(double x[], double y[], double z[], int n) {
  int i;
  for (i = 0; i < n; i++) {
    z[i] = x[i] + y[i];
  }
}

void Parallel_vector_sum(double local_x[], double local_y[], double local_z[], int local_n) {
  int local_i;
  for (local_i = 0; local_i < local_n; local_i++) {
    local_z[local_i] = local_x[local_i] + local_y[local_i]; 
  }
}

void Read_vector(double local_a[], int local_n, int n, char vec_name[], int my_rank, MPI_Comm comm) {
  double * a = NULL;
  int i;
  if (my_rank == 0) {
    a = new double [n];
    std::cout << "Enter the vector %s \n" << vec_name;
    std::stringstream ss;
    
    for (std::size_t i = 0; i < n; i++) {
      std::scanf("%lf ", &a[i]);
      MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0, comm);
      delete[] a;
    }
  }
  
  else {
    MPI_Scatter(a, local_n, MPI_DOUBLE, local_a, local_n, MPI_DOUBLE, 0, comm);
  }
}

void Print_vector (double local_b[], int local_n, int n, char title[], int my_rank*, MPI_Comm comm) {
  double b* = NULL;
  int i;
  if (my_rank == 0) {
    b = malloc (n * sizeof(double));
    MPI_Gather(local_b, local_n, MPI_DOUBLE, b, local_n, MPI_DOUBLE, 0, comm);    
  }

}



/*
  Now suppose we want to test our vector addition function. It would be convenient to be able to read the dimension
of the vectors and then read in the vectors x and y.
 */

double NewtonInterpolant(double x, int npts, double *xpts, double * newton_coeff) {
  double sum = 0.0, xval;

  for (int i = 0; i < npts; i++) {
    xval = 1.0; 
  }
}

int main(int argc, char ** argv) {
  int my_rank, comm_sz, n = 1024, local_n;
  double a = 0.0, b = 3.0; //  the limits of the function

  // Local variables are whose contents are significant only on the process that's using them.
  
  double h, local_a, local_b;
  double local_int, total_int;
  int source;
  
  MPI_Status status; // variable to contain status information
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz); // processes in the communicator 
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // process ranks 


  h = (b - a) / n;
  local_n = n / comm_sz; // dividing the number of trapezoid for each node
  
  local_a = a + my_rank*local_n*h;
  local_b = local_a + local_n*h;
  local_int = Trap(local_a, local_b, local_n, h);

  /*
  if (my_rank != 0) {
    MPI_Send(&local_int, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD); // send to rank 0 the node with the rank 0 to the MPI_COMM_WORLD 'pool'
  } else {

    total_int = local_int;

    for (source = 1; source < comm_sz; source++) {
      // Make sure to add up what the total_int at node 0 must recieve. 
      MPI_Recv(&local_int, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE); // receive nodes to local_int, of source MPI_COMM_WORLD with the tag 0  -> At this point, the total_int has received the values of each of the nodes local_int 
      total_int += local_int;
    }
  }
  */
  
  MPI_Reduce(&local_int, &total_int, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  if (my_rank == 0) {
    std::cout << "With n = %d trapezoids, our estimate " << n << std::endl;
    std::cout << "of the integral " << a << " " <<  b << " " << total_int << std::endl;
  }
  
  MPI_Finalize();

  return 0;  
}
