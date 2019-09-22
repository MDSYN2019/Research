#ifndef __NEW__
#define __NEW__ 


#include <exception>
#include <iostream>
#include <deque>
#include <thread>
#include <vector>
#include <sstream>
#include <cstring>
#include "openmp1.hpp"

//! Making class objects act like values

/*!

Objects of built-in types generally behave like values: whatever we copy 
an object of such a type, the original and the copy have the same value 
but are otherwise independent. Subsequent changes to one object do not affect the other. We can create 
objects of these types.

--------------------------
Class A -> copy -> Class B
--------------------------

For most of the built-in types, the language also defines a rich set of operators
and provides automatic conversions between logically similar types. For example, 
if we add an int and a double, the compiler automatically converts the int 
into a double.

Whene we define our own classes, we control the extent to which the resulting objects 
behave like values.

 */


//! A Str class - a lightweight version of the string container
/*!
  Our class delegates the work of managing its data to the Vec class 
  that we wrote. 

  The Str class has four constructors, each of which arranges to create data 
  as an appropriately initiated Vec object.

  The default constructor for Str implicitly invokes the Vec default 
  constructor to create an empty Str. Note that because our class 
  has other contructors, we must explicitly define the default constructor,
  even though it does exactly what the synthesized defauly constructor 
  would have done.

 */

Str::ostream& operator<<(std::ostream& os, const Str& s) {
  /*
    output operator - iterate through the Str, writing a 
    single character at a time
  */
  for (Str::size_type i = 0; i != s.size(); ++i) {
    os << s[i];
  }
  return os;
}


/*
  Friends 
  -------
 */

Str::istream& operator>>(std::istream& is, const Str& s) {

  /*

    The input operator isn't much harder to write than the output operator. It needs to 
    read and remember characters from the input stream. Each time we call the input 
    operator, it should read and discard any leading whitespace, and then read and rememeber 
    characters until it hits whitespace or end-of-life 


    Each time we call the input operator, it should read 
    and discard any leading whitespace, and then and remember 
    characters until it hits whitespace or the end of file 

  */

  s.data.clear();
  
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

#endif 
