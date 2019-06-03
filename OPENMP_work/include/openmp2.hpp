#include <iostream>
#include <cstdlib>
#include <cstdio>

/* Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP is defined. */

#ifndef _OPENMP
#include <omp.h>
#endif

#include "openmp1.h"

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

