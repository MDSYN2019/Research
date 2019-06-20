#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <string>
#include <memory>

#ifndef _OPENMP_
#include <omp.h>
#endif

// Main header to include 

#include "openmp1.hpp"

// QAT headers

//#include "Argument.h"

/* CPPunit tests */

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>



//int RandomNumber() {
//  return (std::rand() % 100);
//}

// --------

// Is in the hpp file

//template <class T>
//Vec<T>::Vec(size_type n, const T& t = T() ) {
// create(n,t);
//}

/*
template <class T>
OMP<T>::OMP(int N) {
  thread_count = N;
  my_rank = omp_get_thread_num();
  n = omp_get_num_threads();
} // Constructor 

template <class T>
OMP<T>::OMP(const OMP& ref) {  
} // Reference allocator 

template <class T>
OMP<T>::~OMP() {
} // Destructor 

template <class T>
void OMP<T>::add(T a) {
  val += a;
}

template <class T>
void OMP<T>::addup() {
  global_result = 0.0;
#pragma omp parallel num_threads(thread_count) reduction(+:val) // In OpenMP it may be possible to spcift that th result of a reduction is a reduction variable.
  this->add(3); // Use the function that has already been allocated onto this class 
#pragma omp critical
  global_result = val;
}

template <class T>
int OMP<T>::Linear_search(int key, int* A, int n) {
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

template <class T>
void OMP<T>::pi() {
  double factor = 1.0;
  double sum = 0.0;
  // Computing Pi with openMP
  
# pragma omp parallel for num_threads(thread_count) reduction(+:sum) private(factor)
  for (int k = 0; k < n; k++) {
    if (k % 2 == 0) {

      factor = 1.0;

    } else {

      factor = -1.0;
    }
    
    sum += factor / (2*k + 1);
  }
}

*/
