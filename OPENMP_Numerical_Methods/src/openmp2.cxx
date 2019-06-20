#ifndef __NEW__
#define __NEW__ 


#include <exception>
#include <iostream>
#include <deque>
#include <thread>
#include <vector>
#include <cstring>
#include "openmp1.hpp"

/*
template <typename T, typename Container = std::deque<T> >
class stack {
public:
private:
};


template <typename T, typename Container = std::vector<T> >
class AA {
public:
private:
};

*/

class Str {
public:
  typedef Vec<char>::size_type size_type;
  Str() {}

  // Create a Str containing n copies of c
  Str(size_type n, char c): data(n, c) {}

  // Create a Str from a null-terminated array of char
  Str(const char* cp) {
    std::copy(cp, cp + std::strlen(cp), std::back_inserter(data));
  }
  template <class In> Str(In b, In e) {
    std::copy(b, e, std::back_inserter(data));
  }
private:
  Vec<char> data;
};


class Complex { 
  friend bool operator==(const Complex& a, const Complex& b);
  double real, imaginary;
  
public:
  Complex( double r, double i = 0 ) 
    : real(r) , imaginary(i) 
  {}
};

bool operator==( const Complex &a, const Complex &b )
{ 
  return a.real == b.real  &&  a.imaginary == b.imaginary; 
}


#endif 
