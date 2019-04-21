#include <iostream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>

#include "openmp1.h"
#include "syn_dbg.hpp"

void test_debug() {
  debug("A message");
  // passing in arguments like printf
  debug("I am %d years old ", 37);
}

OMP AA(5);

int main () {
  AA.addup();
  test_debug();
  return 0;
  
}
