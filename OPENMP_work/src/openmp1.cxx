/*
OpenMP - Open multi processing

Pthreads requires the programmer to explicitly epecifiy the behaviour of each thread. 

OpenMP (Open multi-processing), on the other hand, sometimes allows the programmer 
to simply state that a block fo code should be exectuted in parallel, and the precise 
determination of the takss and which thread should execute them is left to the compiler and 
the run-time system

OpenMP is a way to program shared memory devices. This means that the parellism occurs where every parallel thread has access 
to all your data 

Example:
You can think of it as: parallelism can appen during the execution of a specific for loop by splittig up the loop
among the different threads.

OpenMP provides hat's known as a "directives-based" shared memory API. In C and C++, this means that there are speciial preprocessor 
instructions known as pragmas. Pragmas are typcailyl added to a system to allow behaviours that aren't part of the basic C specification

*/

#include <iostream>
#include <cstdlib>
#include <cstdio>
//#include <Eigen/Dense>

/* Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP is defined. */

#ifndef _OPENMP
#include <omp.h>
#endif

#include "openmp1.h"

/* cppunit tests */                                                                                                                                      
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

/*
  Data dependencies

  If a for loop fails to satisfy one of the rules outlined in the 
  preceding section, the compiler will simply reject it.

  For example, suppose we try to compile a program with the following 
  linear search function

  OpenMP will only parallelize for loops that are in the canonical form

  Canonical form: 
  
  -  The variable index must have an integer or pointer type (cannot be a float!!)  
  
  -  The expressions start, end, and incr must have a compatible type. For example, 
     if index is a pointer, then incr must have an integer type.

  -  The expressions start, end and incr must not change during the execution for the loop
  
  -  During execution of the loop, the variable index can only be modified by the 
     'increment expression' in the for statement.

 */

// Estimating pi
// serial pi estimator 


/*
double factor = 1.0;
double sum = 0.0;
for (int k = 0; k < n; k++) {
  sum += factor / (2*k + 1);
  factor -= factor;
 }
pi_approx = 4.0 * sum;


// parallel OpenMP pi estimator
double factor = 1.0;
double sum = 0.0;

# pragma omp parallel for  num_threads(thread_count)	\
  reduction(+:sum)
for (int k = 0; k < n; k++) {
  sum += factor(2*k + 1);
  factor -= factor;
 }

int Linear_search(int key, int A[], int n) {
  int i;
  // thread count is global
# pragma omp parallel for num_threads(thread_count)
  for (int i = 0; i < n; i++) {
    if (A[i] == key) {
      return i;
    }
    return -1;
    
  }
}

template <class T> class Vec {
public:
  // Interface
  // Two constructors defined here -
  // We need to provide typedefs for the const and nonconst iterator types.
  
  Vec() {create();} // Default constructor 
  explicit Vec(std::size_t n, const T& = T()) { create(n, val);} // When we say that a constructor is explicit, we're saying that the compiler
                                                                 // will use specifically invokes the constructor

  typedef T* iterator;
  typedef const T* const_iterator
  typedef size_t size_type;
  
private:
  // Implemntation
  T* data;
  T* limit;
};

*/

OMP::OMP(int N) {
  thread_count = N;
}
OMP::~OMP() {
}
void OMP::add(int a) {
  val += a;
}
void OMP::addup() {
  global_result = 0.0;
#pragma omp parallel num_threads(thread_count)
  this->add(3);
#pragma omp critical
  global_result = val;
}

