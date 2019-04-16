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

// Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP is defined.

#ifndef _OPENMP
#include <omp.h>
#endif

#include "openmp1.h"

// cppunit tests                                                                                                                                      
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

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


OMP::OMP() {
# pragma omp parallel num_threads(thread_count) \ // OpenMP directive - the program should start some threads. Each thread that's forked
  reduction(+: global_result) // Code specified that global result is a reduction variable
  global_result += Local_trap(double a, double b, int n);
}
  
  OMP::~OMP() {
  }
  
  void OMP::OMP_reduce() {
    h = (b-a)/n;
    approx = ((float)a + (float)b)/2.0;
# pragma omp parallel for num_threads(thread_count) \
  reduction (+: approx)
    for (int i = 1; i <= n -1; i++) {
      approx += (float)(a+ i*h);
    }
    pprox = h * approx;
    
  }
