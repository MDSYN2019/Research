#ifndef __NEW__
#define __NEW__ 

#include <exception>
#include <iostream>

int add(int i, int j) {
    return i + j;
}

/*
PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring
    m.def("add", &add, "A function which adds two numbers");
}
*/

/*
class Core {
public:
  Core() : midterm
private:
};
*/

class Complex { 
  friend bool operator ==(const Complex& a, const Complex& b);
  double real, imaginary;
public:
  Complex( double r, double i = 0 ) 
    : real(r)
        , imaginary(i) 
  {
  }
};

bool operator ==( const Complex &a, const Complex &b )
{ 
  return a.real == b.real  &&  a.imaginary == b.imaginary; 
}


#endif 
