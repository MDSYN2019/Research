#ifndef __temp__
#define __temp__



/*

Function templates provide a functional behaviour that can be called for different types.

In other words, a function template represents a family of functions

 */

template <typename T>
inline T const& max (T const& a, T const& b) {
  return a < b?b:a;
}

#endif
