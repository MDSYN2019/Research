#ifndef __EIGEN__LIB__
#define __EIGEN__LIB__

#include <cstddef>
#include <iostream>
#include <complex>
#include <string>
#include <sstream>
#include <vector>


namespace Skiena {
  template <class T> class linked_list {
  public:
    list *search_list(list *l, item_type x);
    linked_list(); // Constructor
    ~linked_list(); // Destructor 
  private:
    double item;
    typedef struct list {
      T val;
      struct list *next;
    };
  };
  
  linked_list::linked_list() {
  }
  
  linked_list::*search_list(list *l, item_type x) {
    if (l == NULL) {
      return (NULL);
    }
    if (l->item == x) {
      return (l);
    } else {
      return (search_list(l->next, x));
    }
  }
    
    
}

class Schrodinger {
public:
  // Constructors
  Schrodinger();
  Schrodinger(int* , std::string);
  void inputY();
private:
  // These two values should be the same in this case
  int nrows; // number of rows in the matrix
  int ncolumns; // number of columns in the matrix 
  double v = 0.5;
  double k = 0.2;
  std::string usage;
  Eigen::VectorXcd Y(nrows);
  
};

Schrodinger::Schrodinger(int* argc, std::string* argv) {
  usage = std::string("Usage: ") + argv[0] 
    std::cout << "" << std::endl; 
}

Schrodinger::inputY() {
  if (nrows.size() == 0 || ncolumns.size() ==0):
    throw std::runtime_error("Need to input the right number of columns/rows for the matrix");
}


class complex {
private:
  double real, img;
public:
  complex(double real = 0, double image = 0);
  complex operator+(const complex&) const; // make a overloaded function with the + operation which takes in the address of another complex number to return a complex number again.
};

// define contructor
complex::complex (double r, double i) {
  real = r;
  imag = i;
}

complex complex::operator+ (const complex& c) const {
  complex result;
  result.real = (this->real + c.real);
  result.imag = (this->imag + c.imag);
  return result; // This should be a const, as we defined the overloaded operator to return a const 
}

template <class T> class Vec {
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef std::ptrdiff_T difference_type;
  typedef T& reference;
  typedef const T& const_reference;
  /*
    What about self-assignment? It is possible that a user might wind up assigning an object to itself. As we shall see, it iss crucial 
    that assignment operators deal correctly with self-assignment
    
   */
  template <class T>
  Vec<T>& Vec<T>::operator=(const Vec& rhs) {
   // check for self-assignemnt
    if (&rhs != this) {
      // free the array in the left hand side
      uncreate();
      create (rhs.begin(), rhs.end());
    }
    return *this;
  }
    
    
  &Vec operator=(const Vec&);
  
  Vec() {create();}
  // we'll assume that one of our utiity functions will handle the allocation and copy so that the copy constucorr can forward its wworkd to taht function
  Vec(const Vec& v){(create(v.begin(), v.end();} // copy constructor
		     
		     explicit Vec(std::size_t n, const T& val = T()) {create(n,val);} // this explicit constructor takes a size_type and a value - this will allocate neough mempry of type T of number n, and initialize it with the values val 

  size_type size() const {return limit - data;}
  T& operator[] (size_type i) { return data[i]}
  const T& operator[](size_type i) const {return data[i];}

  iterator begin() {return data;}
  const_iterator begin() const {return data;}

  iterator_end() {return limit;}
  const_iterator end() {return limit;}
  
private:
  iterator data; // first the first element of the data
  iterator limit;
};

#endif
