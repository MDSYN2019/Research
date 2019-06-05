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
  for (int phase = 0; phase < n; phase++) {
    if (phase % 2 == 0) {
      for (int i = 1; i < n; i += 2) {
	if (a[i-1] > a[i]) {
	  Swap(&a[i-1], &a[i]);
	}
      }
    }
  }

  for (int phase = 0; phase < n; phase++) {
    if (phase % 2 == 0) {
# pragma omp parallel for num_threads(thread_count) default(none) shared(a, n) private(i, tmp)
      for (int i = 1; i < n; i += 2) {
	if (a[i-1] > a[i]) {
	  tmp = a[i-1];
	  a[i-1] = a[i];
	  a[i] = tmp;
	}
      }
    } else {
# pragma omp parallel for num_threads(thread_count) default(none) shared(a, n) private(i, tmp)

      for (int i = 1; i < n-1; i += 2) {

	if (a[i] ) {

	}

      }
    }
    
  }

  // AA.addup();
  // test_debug();
  return 0;
}
