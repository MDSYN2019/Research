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

  Str(const char* cp) {  
    std::copy(cp, cp + std::strlen(cp), std::back_inserter(data); // Copy input string into the vector
  }
      
      template <class In> Str(In b, In e) {
      std::copy(b, e, std::back_inserter(data));
    }

    /*

      There are two versions of the index operator. One that can operate 
      on const objects and another which cannot.

      By returning a reference to the character, the nonconst version
      gives access to the character.

      The const version returns a reference to a const char, 
      thereby preventing the user from writing the underlying 
      character 


     */

    char& operator[] (const size_type i) {return data[i];} // The index operator -> Take a index i and return .. 
    const char& operator[] (const size_type i) const {return data[i];}

    // Input/Output Operators

    /*
      What about the other functions? 
      -------------------------------

      The most interesting problem in defining these functions 
      is deciding whether these operations should be members 
      of the Str class. 
      
    */
    
    std::istream& operator>>(std::istream&, Str&);    
    std::ostream& operator<<(std::ostream&, const Str&); 

    
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
