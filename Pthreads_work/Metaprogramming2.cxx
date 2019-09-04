

/*

Metaprogramming

https://www.youtube.com/watch?v=lrziylOWBT4

- If we have a input of a string to a certain function, we could 
  check if that string has certain properties at compile time.

- 

- Motivation and uses for TMP

- Type traits and type checking

- Member detection

- Traits classes and algorithm selection

- If time permits: compile-time computation


Motivation: 

Algorithm Selection:

- Some STL algorithms have a more eficient or practical implementation 
for some restricted subset of types

- std::copy for pointers to POD elements

- std::distance for random access iterators

<type_traits>

=> A large collection of metafunctions that check properties of types
   
   - Some of them can't be implemented without compiler support
   - Some of them can be easily implemented if necessary (or for practise)



*/


#include <iostream>
#include <iterator>
#include <vector>

//---------------------------------------
// Reminder: Class Template Specialization
//---------------------------------------


template <typename T> struct less {
  bool operator()(T a, T b) const {return a < b; }
};

template <typename T> struct less {
  bool operator()(T a, T b) const {return a < b; }
};


int main(void) {

  std::vector<int> v{3, 1, 4};
  std::cout << std::distance(v.begin(), v.end()) << "\n";
  
  
  // std::distance

  

  return 0;
}
