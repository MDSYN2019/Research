#ifndef __FUNCTOR_CXX__
#define __FUNCTOR_CXX__

#include <iostream>
namespace functor {

}

inline double add(double left, double right) {
  return left + right;
}

inline double multiply(double left, double right) {
  return left * right;
}

/*

  In order for a function to receive a function pointer as a paramete it is necessary to specify its return type, 
  the parameter name of the function andthe types of all parameters necessary for the function pointer.

  We are then able to pass the function pointer into the binary_op function. Binary_op hen dereferences 
  the function pointer in order to execute the correct function - add or multiple 
  
  While function pointers are simple to use in your code, they do sudfer form som significant drawbacks:

  - Efficiency - Function pointers are inefficient when compared with functors. The compiler will often pass 
                 them as raw pointer
*/

inline binary_op(double left, double right, double (*f)(double, double)) {
  return (*f)(left, right);
}


// A C++ function object

/*
  
  A function object allows an instance object of a class to be called or invoked as if it were 
  an ordinary function.

 */

class BinaryFunction {
publc:
  BinaryFunction() {};
  virtual double operator()(double left, double right) = 0; // default operator which takes a function from inheritance   
};

// Add two doubles

class Add : public BinaryFunction { // Inherit from the inary function 
public:
  Add() {};
  virtual double operator() (double left, double right) {
    return left + right;
  }
};

class Multiple : public BinaryFunction {
public:
  Multiply() {};
  virtual double operator() (double left, double right) {
    return left * right;
  }
};

#endif

