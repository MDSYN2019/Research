#ifndef _OPENMP1_
#define _OPENMP1_

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <memory>
//#include <Eigen/Dense>
#include <omp.h>

/*
  
As we noted earlier, we'll usually specifiy the number of threads on the command line, so we'll modify our parllalel directive ws with the 
num_threads clause. A clause in OpenMP is just some text that modifies a directive. The num_thereads clause 
can be aded to a parallel directive. 

It allows the programmer to speucfy the number of threads that should be execute the following block;

# prgama omp parlallel nun)threads(thread)count)

/*
What actually happens when the program gets to the parallel directive? Prior to the parllel directive, the program is using a 
single thread, the process started when the program started execution.
*/

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

template <class T> class Vec {
 public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;

  Vec() { create(); }
  explicit Vec(size_type n, const T& t = T() ) {
    create(n,t);
  }

  Vec(const Vec& v) {create(v.begin(), v.end()); }
  Vec& operator=(const Vec&);
  ~Vec() {uncreate();}

  const T& operator[] (size_type i) const {return data[i];}

  void push_back(const T& t) {
    if (avail == limit) {
      grow();
      unchecked_append(t);
    }
  }
  
  size_type size() const {return avail - data;}
  iterator begin() {return data;}
  const_iterator begin() const {return data;}
  iterator end() {return avail;}
  const_iterator end() const {return avail;}
  // Testing positions
  
 private:
  iterator data;
  iterator avail;
  iterator limit;
  // facilities for memory allocation
  std::allocator<T> alloc;
  // allocate and initlize the underlying array
  void create();
  void create(size_type, const T&);
  void create(const_iterator, const_iterator);
  // destroy the elements in the array and free the memory
  void uncreate();
  // support functions for push_back
  void grow();
  void unchecked_append(const T&);
};

class OMP {
public:
  // Constructors and destructors 
  OMP(int);
  //  OMP(const OMP& OMPCopy); // Copy constructor 
  // OMP& operator=(const OMP& ref); // self-assignment operator
  ~OMP(); // Destructor  
  void add(int);
  void addup();
  void pi();

  int Linear_search(int, int*, int n);
 private:
  int val;
  int thread_count;
  int global_result;
  int my_rank; // get current rank
  int n;
};

#endif
  
