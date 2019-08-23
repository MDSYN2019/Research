
//! Openmp explanation 

/*!

As we noted earlier, we'll usually specify the number of threads on the command line, so we'll modify parallel directive was with the 
num_threads clause. A clause in OpenMP is just some text that modifies a directive. The num_threaads clause, it allows a programmer 
to specify the number of threads that should be executed in the following block:

# pragma omp parallel num_thrads (thread_count)

What actually happens when the program gets to the parallel directive? Prior to the parallel directive, 
the program is using a single thrad, the process started when tbe program started execution.

*/


#ifndef GUARD_VEC_H
#define GUARD_VEC_H

#include <iostream>	 /*!< std::cout, std::endl */
#include <cstddef> 	 /*!< std::size_t */
#include <algorithm>     /*!< std::max */
#include <memory>	 /*!< std::allocator, std::uninitialized_fill, std::uninitialized_copy */
#include <cstdlib>
#include <cstdio>
#include <memory>
#include <type_traits>

//#include <Eigen/Dense>
#include <omp.h>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>



template <class T> class OMP {
public:
  OMP(int);   // Constructors and destructors 
  explicit OMP(const OMP& OMPCopy); // Copy constructor 
  OMP& operator=(const OMP& ref); // self-assignment operator  
  ~OMP(); // Destructor  
  void add(T);
  void addup();
  void pi();
 private:
  int val;
  int thread_count;
  int global_result;
  int my_rank; // get current rank
  int n;
  std::allocator<T> alloc;
};



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


/*

For every class, we need to define these properly 

T::T() - one or more constructors, perhaps with arguments

T::~T() - the destructor 

T::T(const T&) - the copy constructor

T::operator=() - the assignment operator

*/


// Inherits from the testcase 
template <class T> class Vec : public CppUnit::TestCase {
 public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef T& reference;
  typedef const T& const_reference;

  Vec() {create();}
  explicit Vec(size_type n, const T& t = T()) {
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
  void runTest();

 private:
  iterator data;
  iterator avail;
  iterator limit;
  // facilities for memory allocation
  std::allocator<T> alloc;
  // allocate and initlize the underlying array
  void create(); // checked 
  void create(size_type, const T&); // checked 
  void create(const_iterator, const_iterator); // checked 
  // destroy the elements in the array and free the memory
  void uncreate(); // checked 
  // support functions for push_back
  void grow();
  void unchecked_append(const T&);
};

template <class T>
Vec<T>& Vec<T>::operator=(const Vec& rhs) {
  if (&rhs != this) {
    uncreate();
  }
  create(rhs.begin(), rhs.end());
  return *this;
}

template <class T> void Vec<T>::create() {
  data = avail = limit = 0;
}

// Is in the hpp file
template <class T> void Vec<T>::create(size_type n, const T& val) {
  data = alloc.allocate(n);
  limit = avail  = data + n;
  std::uninitialized_fill(data, limit, val);
}

// Is in the hpp file 
template <class T>
void Vec<T>::create(const_iterator i, const_iterator j) {
  data = alloc.allocate(j - i);
  limit = avail = std::uninitialized_copy(i,j,data); // allocate additional space at the end of the vector
}

//---------

// Is in the hpp file 
template <class T>
void Vec<T>::uncreate() {
  if (data) { // If data exists
    iterator it = avail;
    while (it != data){
      alloc.destroy(--it);
    }
    // Return all the space that was deallocated
    alloc.deallocate(data, limit-data); 
  }
  
  data = limit = avail = 0;
}

template <class T> void Vec<T>::unchecked_append(const T& val)
{
    alloc.construct(avail++, val);
}


template <class T> void Vec<T>::grow()
{
    // when growing, allocate twice as much as space as currently in use
    size_type new_size = std::max(2*(limit-data), ptrdiff_t(1));
    // allocate new space and copy existing elements to the new space
    iterator new_data = alloc.allocate(new_size);
    iterator new_avail = std::uninitialized_copy(data, avail, new_data);
    // return the old space
    uncreate();
    // reset pointers to point to the newly allocated space
    data = new_data;
    avail = new_avail;
    limit = data + new_size;
}

template <class T> void Vec<T>::runTest() {
  int A = 1;
  CPPUNIT_ASSERT(A == 1);  
}

/*

  Using inheritance and dynamic binding

  The ability to build our own datatypes is fundamental to the idea of object orientated programming (OOP)

  Now we want to look at the other key compononet of OOP - which is inheritance and dynamic binding

 */


#endif /* GUARD_VEC_H */
