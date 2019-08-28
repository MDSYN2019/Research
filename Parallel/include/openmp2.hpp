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

/*!
----------------
Custom Str class
----------------

A numeric progression is a seuqence of numbers, where the value of each number depends on one or more of the 
previous value.


Objects of built-in types generally behave like values: Whenever we copy an object of 
such a type, the original and copy have the same value but are otherwise indepedent. 

For most of the built-in types, the language also defines a rich set of operators and provides 
automatic conversions between logically similar types. For example, if we add 
an int and a double, the compiler automatically converts the int into a double

When we define our own classes, we control the extent to which the resulting objects behave like values.
By defining copying and assigning appropriately, the class author an arrange for objects of that class 
to act like values - that is, the class author can arrange for each object to have state that is independent
of any other object. 

Our Vec and Student_info classes are examples of types that act like values

We shall see that the class author an also control conversions and related operations on class objects, thereby 
providing classes whose objects behave even more similarly to objects of built-in types. 

*/

//! A constructor with a character pointer input
/*! 
  Defining a Str class that lets us create objects that behave approproximately as we would like.
 */

class Str {
public:
  typedef Vec<char>::size_type size_type;
  Str () {}   /*!< default constructor, create an empty str */  
  Str (size_type n, char c): data(n, c) {};  /*!<  create a Str containing n copies of c */


  //! Last two custom constructors
  /*! 

    The last two consturctors are similar to each oter. Their constructor initializers 
    are empty, which means that data is implicitly initializes as an empty Vec. Each 
    constructor asks copy to append the supplied characters  to the initally 
    empty data.


   */
  Str(const char* cp) {  
    std::copy(cp, cp + std::strlen(cp), std::back_inserter(data); // Copy input string into the vector
  } 
  //! create a Str from the range denoted by iterators b and e
  /*!
    The most interesting constructor is the ifnal one, which takes two iterators 
    and creates a new Str that contains a copy of the characters in the given 
    sequence. 

    Like the 


   */    
  template <class In> Str(In b, In e) {
    std::copy(b, e, std::back_inserter(data));
  }
private:
  Vec<char> data;
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

