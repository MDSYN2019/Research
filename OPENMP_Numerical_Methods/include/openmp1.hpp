
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


#ifndef GUARD_VEC_H
#define GUARD_VEC_H

#include <iostream>	 // std::cout, std::endl
#include <cstddef> 	 // std::size_t
#include <algorithm>     // std::max
#include <memory>	 // std::allocator, std::uninitialized_fill, std::uninitialized_copy
#include <cstdlib>
#include <cstdio>
#include <memory>
//#include <Eigen/Dense>
#include <omp.h>


#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


template <class T> class Vec {
public:
    // member types
    typedef T* iterator;
    typedef const T* const_iterator;
    typedef std::size_t size_type;
    typedef T value_type;

    // constructors
    Vec() {
        std::cout << "calling default constructor" << std::endl;
        create();
    }

    explicit Vec(size_type n, const T& t = T()) {
        std::cout << "calling the explicit constructor" << std::endl;
        create(n, t);
    }

    // copy constructor,
    Vec(const Vec& v) {
        std::cout << "calling copy constructor" << std::endl;
        create(v.begin(), v.end());
    }

    //assignment operator
    Vec& operator=(const Vec&);

    // destructor
    ~Vec() {
        std::cout << "calling destructor" << std::endl;
        uncreate();
    }

    // indexing operator
    const T& operator[](size_type i) const {
      std::cout << "calling operation[]" << std::endl;
        return data[i];
    }

    // push_back function
    void push_back(const T& t){
        if(avail == limit)
            grow();
        unchecked_append(t);
    }

    // size function
    size_type size() const { return avail - data; }

    // begin(), end() function
    iterator begin() { return data; }
    const_iterator begin() const { return data; }
    iterator end() { return avail; }
    const_iterator end() const { return avail; }

private:
    iterator data;  // first element in the Vec
    iterator avail; // (one past) the last element in the Vec
    iterator limit; // (one past) the allocated memory

    // facilities for memory allocation
    std::allocator<T> alloc; // object to handle memory allocation

    // allocate and initialize the underlying array
    void create();
    void create(size_type, const T&);
    void create(const_iterator, const_iterator);

    // destroy the elements in the array and free the memory
    void uncreate();

    // support functions for push_back
    void grow();
    void unchecked_append(const T&);
};

// initialize data members to nullptr
template <class T> void Vec<T>::create()
{
    data = avail = limit = nullptr;
}

// create and initialize data members with a size and a value
template <class T> void Vec<T>::create(size_type n, const T& val)
{
    data = alloc.allocate(n);
    limit = avail = data + n;
    std::uninitialized_fill(data, limit, val);
}

// create and initialize data members by copying values from an input sequence
template <class T> void Vec<T>::create(const_iterator i, const_iterator j)
{
    data = alloc.allocate(j - i);
    limit = avail = std::uninitialized_copy(i, j, data);
}

// destruct the class object using destroy and deallocate functions
template <class T> void Vec<T>::uncreate()
{
    if(!data){
        // destroy the elements in reverse order
        iterator it = avail;
        while(it != data)
            alloc.destroy(--it);
        alloc.deallocate(data, limit - data);
    }
    // reset pointers to indicate that Vec is empty again
    data = limit = avail = 0;
}

// assign values from right-hand operand to the left-hand operand
template <class T>
Vec<T>& Vec<T>::operator= (const Vec& rhs)
{
    std::cout << "calling operator= function" << std::endl;

    // check for self-assignment
    if(&rhs != this)
    {
        // free the array in the left-hand side
        uncreate();

        // copy elements from the right-hand to the left-hand side
        create(rhs.begin(), rhs.end());
    }
    return *this;
}

// reallocate storage to hold more elements
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

// add new element at the end of the vector
template <class T> void Vec<T>::unchecked_append(const T& val)
{
    alloc.construct(avail++, val);
}


template <class T> class OMP {
public:
  // Constructors and destructors 
  OMP(int);
  explicit OMP(const OMP& OMPCopy); // Copy constructor 
  OMP& operator=(const OMP& ref); // self-assignment operator
  ~OMP(); // Destructor  
  void add(T);
  void addup();
  void pi();
  int Linear_search(int, int*, int n);
 private:
  int val;
  int thread_count;
  int global_result;
  int my_rank; // get current rank
  int n;
  std::allocator<T> alloc;
};




#endif /* GUARD_VEC_H */

/*
template <class T> class Vec {
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
*/
  
