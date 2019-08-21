#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <string>
#include <complex>

/* Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP is defined. */

#ifndef _OPENMP
#include <omp.h>
#endif

#include "openmp1.hpp"

/* cppunit tests */                                                                                                                                      
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

/*

A numeric progression is a seuqence of numbers, where the value of each number depends on one or more of the 
previous value.
 
*/


/*
Our class delegates the work of managed ing ts data to the vec class that we write in chapter 11. The class is almost as good enough to 


 */


class Str {
public:
  typedef Vec<char>::size_type size_type;
  // default constructor, create an empty str
  Str () {}
  // create a Str containing n copies of c
  Str (size_type n, char c): data(n, c) {};

  // create a Str from a null-terminated array of char
  Str(const char* cp) {
    std::copy(cp, cp + std::strlen(cp), std::back_inserter(data));
  }
  // create a Str from the range denoted by iterators b and e

  template <class In> Str(In b, In e) {
    std::copy(b, e, std::back_inserter(data));
  }
  
};


typedef std::complex<double> Complex;

class Progression {
public:
  Progression (long f = 0) : first(f), cur(f) {} // allocate values of first and cur as f
  virtual ~Progression() {};
  void printProgression(int n);
protected:
  virtual long firstValue();
  virtual long nextValue();
protected:
  long first;
  long cur;
};

class ArithProgression : public Progression {
public:
  ArithProgression(long i = 1);
protected:
  virtual long nextValue();
protected:
  long inc;
};

template <class T>
class TemplateUnderTest {
  T *t_;
public:
  TemplateUnderTest(T *t) : t_(t) {} // allocate to t_ the value of t 
  void SomeMethod() {
    t->DoSomething();
    t->DoSomeOtherThing();
  }
};

