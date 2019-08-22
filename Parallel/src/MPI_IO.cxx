#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>

#include "mpi.h"

// cppunit tests

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

// Class template for this ..class. Redundant statement.

#include "MPI_IO.hpp"


// Non-object methods

/*
  An alterntive appraoch to grouping data is provided by the MPI functions 
  MPI_pack and MPI_Unpack. MPI_Pack allows one to explicitiy store noncontiguous
  (translation - non-array continuous data) in contiguous memory locations,
  and MPI_Unpack can be used to copy data from a contiguous buffer from 
  a contiguous buffer into noncontiguous memory locations
 */


MPIInput::MPI_input() {
} // constructor 

MPIInput::MPI_input(int mr, int pe) {
  my_rank = mr;
  p = pe;
}

void MPIInput::MPIStart() { 
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);  
}


void MPIInput::I_send() {
  int power_2_stage;

  // 2^stage = 1 << stage
  power_2_stage = 1 << stage;

  if (my_rank < power_2_stage) {
    *dest_ptr = my_rank + power_2_stage;

    if (*dest_ptr >= p) return 0;
    else return 1;  
  }  
}


void MPIInput::bubbleSort(int a[], int n) {
   // --  bubble sort variables --
  int temp;
  for (int list_length = n; list_length >= 2; list_length--) {
    for (int i = 0; i < list_length - 1 ; i++) {
      if (a[i] > a[i+1]) {
	temp = a[i];
	a[i] = a[i+1];
	a[i+1] = temp;
      }
    }
  }
}


void MPIInput::getData() {

  if (my_rank == 0) {
    std::cout << "Enter a, b and n \n";
    scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
    for (int dest = 1; dest < p; dest++) {
      tag = 0;
      MPI_Send(a_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 1;
      MPI_Send(b_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 2;
      MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }
  } else {
    tag = 0;
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 1;
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 2;
    MPI_Recv(n_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
  }
}

//! Class destructor for MPI
/*! 
  Destructor 
*/
MPIInput::~MPI_input() {
  MPI_Finalize();
} 
