#ifndef _OPENMP1
#define _OPENMP1

#include <iostream>
#include <cstdlib>
#include <cstdio>
//#include <Eigen/Dense>
#include <omp.h>

/*
  
As we noted earlier, we'll usually specifiy the number of threads on the command line, so we'll modify our parllalel directive ws with the 
num_threads clause. A clause in OpenMP is just some text that modifies a directive. The num_thereads clause 
can be aded to a parallel directive. 

It allows the programmer to speucfy the number of threads that should be execute the following block;

# prgama omp parlallel nun)threads(thread)count)

/*
What actually happens when the program gets to the parallel directive? Prior to the parllel directive, the program is using a 
single thread, the process started when the program started execution.
*/


class OMP {
public:
  OMP(int);
  ~OMP();  
  OMP& operator=(const OMP& ref); // self-assignment operator
  void addup();
  void add(int);
  int Linear_search(int, int*, int n);
  void Compute_trapezium(); 
  void pi();
 private:
  int val;
  int my_rank; // 
  int thread_count;
  int global_result;
  int my_rank; // get current rank
};

#endif
  
