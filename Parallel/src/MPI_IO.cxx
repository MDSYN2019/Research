// Last updated: 26/08/2019

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <iterator>
#include <numeric>
#include <cassert>


#include "mpi.h"

// cppunit tests

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

// Class template for this ..class. Redundant statement.

#include "MPI_IO.hpp"

// Boost libraries to check MPI datatype

#include <boost/mpi/datatype.hpp>

// Boost libraries to read input


// Non-object methods

/*
  An alterntive appraoch to grouping data is provided by the MPI functions 
  MPI_pack and MPI_Unpack. MPI_Pack allows one to explicitiy store noncontiguous
  (translation - non-array continuous data) in contiguous memory locations,
  and MPI_Unpack can be used to copy data from a contiguous buffer from 
  a contiguous buffer into noncontiguous memory locations
 */



//! MPIInput class - constructor
/*! Sets up the processor ranks and size for use in other functions 

  At the moment, we do not have a default constructor 
*/


MPIInput::MPIInput(int mr, int pe) {
  my_rank = mr;
  p = pe;

  // Error checking - need to make sure there are valid in
  if (mr == NULL || me == NULL) {
    // Halt the program 
  }

  // Initiate MPI setup, including the number of processes and the ranks of processes.
  // Generally speaking, 0 is the master process
  
 
  MPI_Init(NULL, NULL); // 

  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // Get the process rank
  MPI_Comm_size(MPI_COMM_WORLD, &p); // allocate the number of processes being used
  
}


//! MPIInput class - MPIStart method
/*! Sets up the processor ranks and size for use later - could really be incoporated into the consturcotr */

void MPIInput::MPIStart() { 
}

void MPIInput::bubbleSort(int a[], int n) {
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


void MPIInput::getDataPack(float* a_ptr, float* b_ptr, int* b_ptr) {

  int position; /*!< Placeholder */
  char buffer[100]; /*!< Keeping track of where the data is */


  if (my_rank == 0) {
    std::cout << "Enter a, b, and n \n";
    scanf("%f %f %d", a_ptr, b_ptr, n_ptr); /*< Read data */
  
    //! MPI_pack 
    /*!
      MPI_pack allows one to explicitly store uncontiguous data in 
      contiguous memory locations, and MPI_unpack can be used to copy
      data from a contiguous buffer into noncontiguous memory locations.      
     */
   
    position = 0; /*!< Now pack the data into buffer. Positon = 0 says start at beginning of buffer */

    MPI_Pack(a_ptr, 1, MPI_FLOAT, buffer, 100, &position, MPI_COMM_WORLD);
    MPI_Pack(b_ptr, 1, MPI_FLOAT, buffer, 100, &position, MPI_COMM_WORLD);
    MPI_Pack(n_ptr, 1, MPI_INT, buffer, 100, &position, MPI_COMM_WORLD);

    MPI_Bcast(buffer, 100, MPI_PACKED, 0, MPI_COMM_WORLD); /*<! Now broadcast contents of the buffers */
  } else {

    
    //! MPI_unpack 
    /*!
      MPI_unpack can be used to copy data from a contiguous buffer into noncontiguous memory locations.
      i.e. it is just the complete opposite of pack
     */

    // Now unpack the contents of the buffer
    MPI_Bcast(buffer, 100, MPI_PACKED, 0, MPI_COMM_WORLD);
    MPI_Unpack(buffer, 100, &position, a_ptr, 1, MP_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, 100, &position, b_ptr, 1, MPI_FLOAT, MPI_COMM_WORLD);
    MPI_Unpack(buffer, 100, &position, n_ptr, 1, MPI_INT, MPI_COMM_WORLD);
  }
  
}



/*
Functino getData

Reads in the user input a, b and n 

Input parameters:
-----------------

1. int my_rank: rank of the current processor 

2. int p: number of processors

Output parameters:
------------------

1. float* a_ptr: pointer to left endpoint a 

2. float* b_ptr: pointer to right endpoint b 

3. int* n_ptr: pointer to the number  

IO on parallel systems
----------------------

Last updated: 26/08/2019

Converting the example in the book to a vector creation getdata and to send it to the nodes

*/

void MPIInput::getData(int* start_inp, int* end_inp, int* n_ptr) { 

  if (my_rank == 0) {
    
    std::cout << "Enter a, b and n \n";  
    std::cout << "Please ensure that a is smaller than b";
    scanf("%lf %lf %d", start, end, n_ptr); // I should really change this into sstream for strict c++, but I dont think that is necessary
    assert(end > start); // Program will stop if this is not accurate
    
    /* Ensure that a C++ style vector can be made */
    std::vector<uint32_t> unintVec;
    std::vector<int> intVec;
    
    for (int i = start; i <= end; i++) {
      intVec.push_back(i);
    }
    int vecSize = intVec.size();

    if (my_rank == 0) {
    for (int dest = 1; dest < p; dest++) { // Looping over the number of processes
      tag = 0;
      MPI_Send(a_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 1;
      MPI_Send(b_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
      tag = 2;
      MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
      // Send vectors
      tag = 3;
      MPI_Send(&intVec[0], vecSize, MPI_INT, dest, tag, MPI_COMM_WORLD); // This should point to the address of the first 
                                                                               // element of the vector 
    }
  } else {
    tag = 0; // The purpose of the tag
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 1;
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
    tag = 2;
    MPI_Recv(n_ptr, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
    tag = 3;
    MPI_Recv(&intVec[0], vecSize, MPI_INT, source, tag, MPI_COMM_WORLD, &status);
    }
}

//! Class destructor for MPI
/*! 
  
Destructor 
----------
  MPI_Finalize() - What does it do?
*/

  MPIInput::~MPI_input() {
  MPI_Finalize();
} 
