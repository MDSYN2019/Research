/*
  One obvious problem with our program is its lack of generality. The function, 
  f(x), and the input data, a, b, n, are hardwired. So if we want to change any 
  of these, we must edit and recompule the program. 
 */

#ifndef __MPI_IO__
#define __MPI_IO__

#include <iostream>
#include <iomanip>
#include <vector>
#include "mpi.h"

class MPI_input {
public:
  MPI_input();
  MPI_input(float*, float*, int*, int, int);
  virtual ~MPI_input();
  void MPI_start();
  void Get_data();
private:
  int* n_ptr;
  int p;
  int source = 0;
  int dest;
  int tag;

  float* a_ptr;
  float* b_ptr;
  int my_rank;
  
  MPI_Status status;
};

#endif 
