
#include <iostream>
#include <string>
#include "max.hpp"
#include "classtemp.hpp"

int main(void) {
  int i = 42;
  std::cout << ::max(7,i) << std::endl;

  double f1 = 3.4;
  double f2 = -6.7;

  ::max(7, 42, 68);
  ::max(7.0, 42.0);
  

  return 0;
}
