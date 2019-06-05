#include <iostream>
#include <omp.h>

double output = 0.0;
double finaloutput = 0.0;

void Add(double* output) {
  for (int i = 0; i <= 10; i++) {
    *output += i;
  }
}

int main() {

  int nProcessors = omp_get_max_threads();
  std::cout << nProcessors << std::endl;  
  omp_set_num_threads(nProcessors); // Set the number of threads

  //std::cout << "There are " <<  omp_get_num_threads() << " threads" << std::endl;
  //int my_rank = omp_get_thread_num();
  //int thread_count = omp_get_num_threads();
  //std::cout << thread_count << std::endl;
  
# pragma omp parallel num_threads(nProcessors) 
  Add(&output);

# pragma omp critical
  finaloutput += output;
  std::cout << finaloutput << std::endl;
  
  return 0;
  
}
