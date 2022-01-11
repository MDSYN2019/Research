#ifndef __STR__
#define __STR__


/*
Making class objects act like values 
------------------------------------


When we copy an object of such a type, the original and copy 
have the same value but are otherwise independent.

The class author can control conversions and related operations on class objects, thereby 
providing classes whose objects behave even more similairly to objects of built-in types.

The standard library string class is a good example of such a type 
because of its rich set of operators for automatic conversions


We will focus on writing a simplified version of sting, str, where 
we will focus on the operators and conversions that let us write expressions 
involving strings


We do not need to worry much about the implementation details of our Str class,
because we did most of the work when we implemented to Vec class.


 */


#include <vector>
#include <algorithm>

// Our Vec class

template <class T> class Vec {
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  
private:
  
}

// Our class delegates the work of managing its data to the Vec class
// that we wrote. That class is almost good enough to support our Str,
// it lacks only a clear function.


class Str {
  typedef Vec<char>::size_type size_type;

  // default constructor; create an empty str
  Str(size_type n, char c): data(n, c) {}

  // Create a Str from a null-terminated array of chars

  Str(const char* cp) {
    std::copy(cp, cp + std::strlen(cp), std::back_inserter(data));
  }


  // Create a Str from the range denoted by iterators b and e
  template <class In> Str(In b, In e){
    std::copy(b, e, std::back_inserter(data));
  }
};


#endif

