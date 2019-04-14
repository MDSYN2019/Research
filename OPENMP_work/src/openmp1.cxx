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

void Trap(double a, double b, int n, double * global_result_p) {
  double h, x, my_result;
  double local_a, local_b;
  int i, local_n;

  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();

  h = (b-a)/n;
  local_n = n/thread_count;
  local_a = a + my_rank * local_n * h;
  local_b = loycal_a + local_n * h;
  my_result = ((float)local_a + (float)local_b) / 2.0;
  for (int i = 1; i <= local_n - 1; i++) {
    x = local_a + i * h;
    my_result += float(x);
  }

  my_result = my_result * h;
# pragma omp critical
  *global_result_p += my_result;
}

/*
  
Trapezium rule 

1. Weidentified two types of tasks:
a. Computation of the areas of individual trapezouids and 
b. Adding the areas of trapezoids

2. There is no communication among the tasks in the first collection, but each task in the first 
   collection communicated with the task in task 1 b

OpenMP consists of a library of functions and macros, so we usually need to iunclude a header 
file with prototypes and macro definitions. 


*/

void Hello(void);
int main(int argc, char *argv[]) {
  int thread_count = strtol(argv[1], NULL, 10);
  
# pragma omp parallel num_threads(thread_count) // OpenMP directive - the program should start some threads. Each thread that's forked
                                                // should execute the Hello function, and when the threads return from the call to Hello,
                                            // hey should be terminate
  
  // The first part is the parallel directive, and as you might have guessed it   

  /* Recollect that thread is short for thread of executon. The name is meant to suggest
     a sequence of statements executed by a program. 
     
     Threads are typically started or forked by a process, and they share most of the resources of the process
     that starts them - for example, access to stdin and stdout - but each thread has its own stack and program 
     counter.

     When a thread completes execution it joins the process that started it.

     It should be noted tha there may be system0defined limitations on the number of 
     threads that a program can start. The OpenMP standard deosnt guarantee that this will
     actually start thread_count threads. 
     
   */
  
  Hello();
 return 0;
}


OMP::OMP() {
# pragma omp parallel num_threads(thread_count) // OpenMP directive - the program should start some threads. Each thread that's forked                                                // should execute the Hello function, and when the threads return from the call to Hello,
}
  
  OMP::~OMP() {
  }

