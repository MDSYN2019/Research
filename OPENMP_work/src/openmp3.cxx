
#include <iostream>
#include <iomanip>


class Str {
public:
  typedef Vec<char>::size_type size_type;
  // default constructor, create an empty Str
  Str() {}

  // Create a Str containing n copies of c
  Str(size_type n, char c) : data(n, c) { }

  
 

};
