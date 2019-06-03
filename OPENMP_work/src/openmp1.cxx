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
#include "openmp1.h"
#include "openmp2.h"

/* CPPunit tests */

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



/*

For every class, we need to define these properly 

T::T() - one or more constructors, perhaps with arguments

T::~T() - the destructor 

T::T(const T&) - the copy constructor

T::operator=() - the assignment operator

*/

// Estimating pi
// serial pi estimator 

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

template <class ForwardIterator, class Generator>
void generate(ForwardIterator first, ForwardIterator last, Generator gen) {
  while (first != last) {
    *first = gen();
    ++first;
  }

}

int RandomNumber() {
  return (std::rand() % 100);
}

template <class T> class Vec {
public:  
  typedef T* iterator; // Defines a type parameter T pointer which acts a a inner iterator in the vector function 
  typedef const T* const_iterator; // Ditto as above, but a constant iterator
  typedef size_t size_type;
  typedef T value_type;
  typedef std::ptrdiff_t difference_type; // What does ptrdiff do?
  typedef T& reference;
  typedef const T& const_reference;
  
  // We form the name of an overloaded operator by appending the operator to the word operator
  Vec() {create();} // Need to replace create here 

  Vec(const Vec& v) {std::generate (v.begin(), v.end(), RandomNumber);} // copy constructor - we allocate new space for the copy vector

  Vec& operator=(const Vec& ref) {  // Defines what it means to assign one value of the class type to another.
    if (&ref != this) {
      uncreate();
      create(rhs.begin(), rhs.end());
    }
  }
  
  explicit Vec(size_type n, const T& val = T()) { create(n, val);} // We are saying to the compiler that we will use the constuctor
                                                                   // only in contexts in which the user expresly invokes the constructor 
  // Prevents weird default constructor allocation 
  // new operations: size and index
  size_type size() const {return avail - data; } // const allocated onto the return value to ensure the return data is not mutable
  

  
  T& operator[](size_type i) {return data[i]; }
  const T& operator[] (size_type i) const {return data[i]; }

  // new functions to return iterators
  iterator begin() { return data;}
  const_iterator begin() const {return data;} // return a constant type
  
  iterator end() {return avail;}
  const_iterator end() const {return avail;}

  void push_back(const T& val) {
    if (avail == limit) // get space if needed 
      grow();
    unchecked_append(val); // append the new element
  }
 
private:
  // Using the iterator type created above
  iterator data;
  iterator limit;
  iterator limit;

  // Facilities for memory allocation
  allocator<T> alloc; // object to handle memory allocation

  void create();
  void create(size_type, const T&);
  void create(const_iterator, const_iterator);

  // destroy the elements in the array and free the memory
  void uncreate();

  // support functions for push_back
  void grow();
  void unchecked_append(const T&);
};

template <class T> void Vec<T>::create() {
  data = avail = limit = 0;
}

template <class T> void Vec<T>::create(size_type n, const T& val) {
  data = alloc.allocate(n);
  limit = avail = data + n;
  uninitialized_fill(data, limit, val);
  
}

template <class T>
void Vec<T>::create(const_iterator i, const_iterator j) {
  data = alloc.allocate(j - 1);
  limit = avail = uninitialized_copy(i, j, data);
}



template <class T>
Vec<T>& Vec<T>::operator=(const Vec& rhs) {
  // check for self-assignment
}

class ComplexNumberTest : public CppUnit::TestCase { 
public: 
  ComplexNumberTest( std::string name ) : CppUnit::TestCase( name ) {} 
  void runTest() {
    CPPUNIT_ASSERT( Complex (10, 1) == Complex (10, 1) );
    CPPUNIT_ASSERT( !(Complex (1, 1) == Complex (2, 2)) );
  }
};


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

struct Student_info {
  std::string name;
  double midterm, final;
  std::vector<double> homework;
};


