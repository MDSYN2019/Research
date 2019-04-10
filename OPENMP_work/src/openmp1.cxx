#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <omp.h>

/*
  OpenMP provides what's known as a "directives-based" sahred memory API. In C and C++, this means that ther are special preprocessor instrutions
  known as pragmas. Pragmas are typically added to a system to allow behaviours that aren't part of the basic C specification. 

  Compilers that dont support the pragmas are free to ignore them. This allows a program that uses the pragmas to run on platforms that
  don't support them.

 */


/*
  
Trapezium rule 

1. Weidentified two types of tasks:
a. Computation of the areas of individual trapezouids and 
b. Adding the areas of trapezoids

2. There is no communication among the tasks in the first collection, but each task in the first 
   collection communicated with the task in task 1 b




 */


/*
OpenMP consists of a library of functions and macros, so we usually need to iunclude a header 
file with prototypes and macro definitions. 
*/

void Hello(void);

int main(int argc, char *argv[]) {
  int thread_count = strtol(argv[1], NULL, 10);
  
# pragma omp parallel num_threads(thread_count) // OpenMP directive - the program should start some threads. Each thread that's forked
                                                // should execute the Hello function, and when the threads return from the call to Hello,
                                            // hey should be terminate
  
  // The first part is the parallel directive, and as you might have guessed it   

  /* Recollect that thread is short for thread of executon. The name is meant to suggest
     a sequence of statements executed by a program. 
     
     Threads are typically started or forked by a process, and they share most of the resources of the process
     that starts them - for example, access to stdin and stdout - but each thread has its own stack and program 
     counter.

     When a thread completes execution it joins the process that started it.

     It should be noted tha there may be system0defined limitations on the number of 
     threads that a program can start. The OpenMP standard deosnt guarantee that this will
     actually start thread_count threads. 
     
   */
  
  Hello();
 return 0;
}


void Hello(void) {
  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();
  std::cout << "Hello from thread " << my_rank  << " of " <<  thread_count << std::endl;
}

void Trap(double a, double b, int n, double * global_result_p) {

  double h, x, my_result;
  double local_a, local_b;
  int i, local_n;

  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();
}

class OMP {
public:
  OMP();
  ~OMP();
  
  int input;
  
private:
  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();
  
};
  



