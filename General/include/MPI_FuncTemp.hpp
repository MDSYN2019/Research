#ifndef __FUNCTEMP__
#define __FUNCTEMP__

/*
Date: 16/08/2019

When the program executer he function call insruction the CPU 
stp

inline functions - one of the most important feature of C++
-----------------------------------------------------------


 */

// maximum of two int values
inline int const& max(int const& a, int const& b) { // inline reduces computation overhead for simple functions
  return a<b?b:a;
}

// maximum of two value of any type
template <typename T>
inline T const& max (T const& a, T const& b) {
  return a<b?b:a;
}

// maximum of three values of any type
template <typename T>
inline T const& max (T const& a, T const& b, T const& c) {
  return max(max(a,b), c);
  
}

#endif 
