#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <string>
#include <complex>

/* Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP is defined. */

#ifndef _OPENMP2_
#include <omp.h>
#endif

#include "openmp1.hpp"

/* cppunit tests */                                                                                                                                      
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

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

    What is interesting about this constructor is that it is itself a template 
    function. Because it is a template, it effectively defines a family of constrctors
    that can be instantiated for different types of iterators.

    This constructor could be used to create a Str from an array of characters.

  */    
      template <class In> Str(In b, In e) {
      std::copy(b, e, std::back_inserter(data));
    }

    //! Str operations

    /*! 
      If we think about the kind of core we've used that used strings, we can see that 
      we used several operators:

      cin >> s // use the input operator to read a string 
      cout << s  // use the output operator to write a string  
      s[i] // use the index operator to access a character 
      s1 + s2 // use the addition oeprator to concatenate two strings

      All these are binary operators, so that if we define as funcitons, 
      each function will have two parameters, one of which may be implicit
      if the function is a member 
      
    */

    //! 

    /*!  
      The index operators just forward their work to the corresponding 
      Vec operations. It is worth noting that, as we did for class Vec,
      we define two version of the index operator.
    */
    
    char& operator[] (const size_type i) {return data[i];}
    const char& operator[] (const size_type i) const {return data[i];}
    std::istream& operator>> (std::istream&, Str&);    
  private:
    Vec<char> data;
  };
}

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

// TODO 
template <class T>
class TemplateUnderTest {
  T *t_; // Pointer of type T - template class
public:
  TemplateUnderTest(T *t) : t_(t) {} // allocate to t_ the value of t  - I think this was called a value constructor 
  void SomeMethod() {
    t->DoSomething();
    t->DoSomeOtherThing();
  }
};

#endif 
