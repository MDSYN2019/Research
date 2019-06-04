#include <iostream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>

#include "openmp1.hpp"
#include "syn_dbg.hpp"

void test_debug() {
  debug("A message");
  // passing in arguments like printf
  debug("I am %d years old ", 37);
}

//OMP AA(5);

int main () {
  int x = 5;

# pragma omp parallel num_threads(thread_count) private(x)
  {
    int my_rank = omp_get_thread_num();
    std::cout << "Thread %d ";
    x = 2 * my_rank + 2;
  }
  // AA.addup();
  // test_debug();
  return 0;
}
