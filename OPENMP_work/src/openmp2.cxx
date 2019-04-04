#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <omp.h>

/*

  When the block of code is completed - in our example, when the threads return from the call to Hell, there's an implicit barrier. 
  
  This means that a thread that has completed the block of code will wait for other threads in the team to complete the block.  

  ------------------
  Scope of Variables 
  ------------------


  In serial programming, the scope of a variable consists of those parts of a program in which the 
  variable can be used. For example, a variable declared at the beginning of a C function has a 
  "function-wide" sceop, that is, it can be only be accessed in the body of the function. 

  On the other hand, a variable declared at the beginning of a .c file but outside any function has a "file-wide" scope, tat is, any function in the file 




 */


void Trap(double a, double b, int n, double* global_result_p);

int main(int argc, char* argv[]) {
  double global_result = 0.0;
  double a, b;
  int n;
  int thread_count;
  Trap(a, b, n, &global_result);
  
  
  return 0;
}

void Trap(double a, double b, int n, double* global_result_p) {
  double h, x, my_result;
  double local_a, local_b;
  int i, local_n;

  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();

  h = (b -a)/n;
  local_n = n / thread_count; 
  local_a = a + my_rank * local_n * h;
  local_b = local_a + local_n * h;

  my_result = ((float)local_a + (float)local_b)/2.0;

  for (unsighed int i = 1; i <= local_n - 1; i++) {
    x = local_a + i * h;
    my_result += (float)x;
  }

  my_result = my_result * h;
  # pramga omp critical
  *global_result_p += my_result;
}
