#include <iostream>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>

#include "openmp1.h"

OMP AA(5);

int main () {
  AA.addup();
  return 0;
  
}
