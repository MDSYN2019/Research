#pragma once 


// maximum of two int values

inline int const& max(int const& a, int const& b) {
  return a<b?b:a;
}

// Overloading the function

// maximum of two values of any type
template <typename T>
inline T const& max(T const& a, T const& b) {
  // if a < b then use b else use a
  return a < b? b:a;
  
}


// maximum of three values of any type
template <typename T>
inline T const& max(T const& a, T const& b, T const& c) {
  return max(max(a,b), c);
}
