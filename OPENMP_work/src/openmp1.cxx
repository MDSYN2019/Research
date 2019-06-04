/*

|---------------------------------|
|OpenMP - "Open multi processing" |
|---------------------------------|

OpenMP (Open multi-processing), on the other hand, sometimes allows the programmer 
to simply state that a block fo code should be exectuted in parallel, and the precise 
determination of the takes and which thread should execute them is left to the compiler and 
the run-time system.

OpenMP is a way to program shared memory devices. This means that the parellism occurs where every parallel thread has access 
to all your data 

Example:

You can think of it as: parallelism can appen during the execution of a specific for loop by splittig up the loop
among the different threads.

OpenMP provides hat's known as a "directives-based" shared memory API. In C and C++, this means that there are speciial preprocessor 
instructions known as pragmas. Pragmas are typcailyl added to a system to allow behaviours that aren't part of the basic C specification

*/

#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <string>
#include <memory>

/* Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP is defined. */

#ifndef _OPENMP_
#include <omp.h>

#endif

// Main header to include 
#include "openmp1.hpp"

/* CPPunit tests */

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

/*

For every class, we need to define these properly 

T::T() - one or more constructors, perhaps with arguments

T::~T() - the destructor 

T::T(const T&) - the copy constructor

T::operator=() - the assignment operator

*/

int RandomNumber() {
  return (std::rand() % 100);
}

template <class T> void Vec<T>::create() {
  data = avail = limit = 0;
}

template <class T> void Vec<T>::create(size_type n, const T& val) {
  data = alloc.allocate(n);
  limit = avail  = data + n;
  std::uninitialized_fill(data, limit, val);
}

template <class T>
void Vec<T>::create(const_iterator i, const_iterator j) {
  data = alloc.allocate(j - i);
  limit = avail = std::uninitialized_copy(i,j,data); // allocate additional space at the end of the vector
}

template <class T> void Vec<T>::uncreate() {

  if (data) { // If data exists
    iterator it = avail;
    while (it != data)
      alloc.destroy(--it);
    // Return all the space that was deallocated
    alloc.deallocate(data, limit-data); 
  }
  
  data = limit = avail = 0;
}

template <class T> void Vec<T>::grow() {
  // when growing, allocate twice as much space as currently in use
  size_type new_size = max(2 * (limit - data), ptrdiff_t(1));
  // allocate new space and copy existing elements to the new space
  iterator new_data = alloc.allocate(new_size);
  //iterator new_avail
}
/* 
OMP::OMP(int N) {
  thread_count = N;
}

OMP& OMP::OMP::operator(const OMP& ref) {
  //
}

OMP::~OMP() {
}

void OMP::add(int a) {
  val += a;
}

int OMP::Linear_search(int key, int* A, int n) {
  int i;
  // thread count is global
# pragma omp parallel for num_threads(thread_count)
  for (int i = 0; i < n; i++) {
    if (A[i] == key) {
      return i;
    }
    else {
      return -1;
    }
  }
}

void OMP::pi() {

  double factor = 1.0;
  double sum = 0.0;
  
# pragma omp parallel for num_threads(thread_count)	\
  reduction(+:sum)
  for (k = 0; k < n; k++) {
    sum += factor / (2*k + 1);
    factor = -factor;
  }
}

void OMP::Compute_trapezium() {

  double h, x, my_result;
  double local_a, local_b;
  int i, local_n;

  my_rank = omp_get_thread_num();
  thread_count = omp_get_num_threads();

  h = (b - a) / n; // each slice of the trapezium
  local_a = a + my_rank * local_n * h;
  local_b = local_a + local_n * h;

  my_result = ((float)local_a + (float)(local_b)) / 2.0;
  for (int i = 1; i <= local_n - 1; i++) {
    x = local_a + i * h;
    my_result += (float)x;
  }
  my_result = my_result * h;
  // TODO
  // #pragma omp critical
  //
  if (n & thread_count != 0) {
    std::cerr << error_type << " n must be evenly divisible by thread_count " << std::end;
    exit(0); // Exit the program 
  }
}

void OMP::addup() {
  global_result = 0.0;
#pragma omp parallel num_threads(thread_count)
  reduction(+:val) // In OpenMP it may be possible to spcift that th result of a reduction is a reduction variable.
                   // To do this, a reduction clause can be added to a parlllel directive

    this->add(3); // Use the function that has already been allocated onto this class 
#pragma omp critical
  global_result = val;
}
*/

