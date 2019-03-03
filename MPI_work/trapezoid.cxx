/*  
  Implementation of the trapezium rule  - A form of numerical integration 
*/

#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <sstream>
#include <cassert>
#include <cmath>
#include <array>
#include <cassert>
#include "mpi.h"

#include "trapezoid.hpp"

Trap::Trap() {}
Trap::~Trap() {}

void Trap::read() {
  std::cout << "Please input Values:" << std::endl;
  std::cin >> a >> b >> n;
  std::cout << a << b << n << std::endl;
}

void Trap::computeTrapezium() {
  h = (b - a)/n;
  integral = ((float)(a) + (float)(b))/2.0;
  x = a;
  for (unsigned int i = 1; i <= n -1; i++) {
    x = x + h;
    integral = integral + (float)(x);
  }
  
  integral = integral * h;
}


