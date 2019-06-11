#ifndef __NEW__
#define __NEW__ 

#include <exception>
#include <iostream>
#include <deque>
#include <thread>

//template <typename T, typename Container = std::deque<T> >
//class stack {
//public:
//private:
//};


class Complex { 
  friend bool operator ==(const Complex& a, const Complex& b);
  double real, imaginary;
public:
  Complex( double r, double i = 0 ) 
    : real(r) , imaginary(i) 
  {}
};

bool operator ==( const Complex &a, const Complex &b )
{ 
  return a.real == b.real  &&  a.imaginary == b.imaginary; 
}


#endif 
