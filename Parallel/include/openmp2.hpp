#ifndef __openmp2__
#define __openmp2__

#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <sstream>
#include <string>
#include <complex>

/* Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP2_ is defined. */

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


/* Tpyedefs */

typedef std::complex<double> Complex;

//! A constructor with a character pointer input

/*! 
  Defining a Str class that lets us create objects that behave approproximately as we would like.
*/

class Str {
public:
  typedef Vec<char>::size_type size_type;
  Str() {}   /*!< default constructor, create an empty str */  
  Str(size_type n, char c): data(n, c) {};  /*!<  create a Str containing n copies of c */

  Str(const char* cp) {  
    std::copy(cp, cp + std::strlen(cp), std::back_inserter(data)); // Copy input string into the vector
  }
      
  template <class In> Str(In b, In e) {
    std::copy(b, e, std::back_inserter(data));
  }

  // Returns a Str
  Str& operator+=(const Str& s) {
    std::copy(s.data.begin(), s.data.end(), std::back_inserter(data));
    return *this; // Return the references Str - dereferenced 
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
    -------------------------------
    What about the other functions? 
    -------------------------------
    The most interesting problem in defining these functions 
    is deciding whether these operations should be members 
    of the Str class. 
    
  */
  friend std::istream& operator>> (std::istream&, Str&);
  friend std::ostream& operator<< (std::ostream&, const Str&);
  size_type size() const {return data.size(); } 

private:
  Vec<char> data;
};

/*
  Other binary operators
  ---------------------

  What remains of our work on the Str class is to implement the + operator.

  Before we can do so, we must make several decisions: Should the operator 
  be a member? 
  
  For now, lets make some initial guesses about the answers. First, we know that 
  we want to be able to concatenate values that are of type Str.

  
 */


std::istream& operator>>(std::istream& is, Str& s) {
      
  /*
    
    The input operator isn't much harder to write than the output operator. It needs to 
    read and remember characters from the input stream. Each time we call the input 
    operator, it should read and discard any leading whitespace, and then read and rememeber 
    characters until it hits whitespace or end-of-life 
    
    
    Each time we call the input operator, it should read 
    and discard any leading whitespace, and then and remember 
    characters until it hits whitespace or the end of file 
    
  */
  
  //s.data.clear();  TODO
  char c;
  
  // isspace - check whether c is a whitespace character 
  while (is.get(c) && isspace(c)) {
    // Nothing to do, except testing the condition
  }
  // If there is still something to read, do so until next whitespace character 
  if (is) {
    do s.data.push_back(c);
    while (is.get(c) && !isspace(c));
    // if we read whitespace, then put it back on the stream
    if (is) {
      is.unget(); 
    }
  }
  return is;
}

std::ostream& operator<<(std::ostream& os, const Str& s) {
  /*
    output operator - iterate through the Str, writing a 
    single character at a time
  */
  for (Str::size_type i = 0; i != s.size(); ++i) {
    os << s[i];
  }
  return os;
}


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

#endif 
