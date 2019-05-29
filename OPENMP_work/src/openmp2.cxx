#ifndef __openmp2__
#include "openmp2.hpp"

void Progression::printProgression(int n) {
  std::cout << firstValue();
  for (int i = 2; i <= n; i++) {
    std::cout << ' ' << nextValue();
  }
  std::cout << std::endl;
}

long Progression::firstValue() {
  cur = first;
  return cur;
}
// # pragma omp parallel num_threads(thread_count)

long Progression::nextValue() {
  return ++cur;
}


ArithProgression::ArithProgression(long i) : Progression(), inc(i) { }
// Constructor 
long ArithProgression::nextValue() {
  cur += inc;
  return cur;
}



#endif
