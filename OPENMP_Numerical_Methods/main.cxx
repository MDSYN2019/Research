#include <iostream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <omp.h>


// Eigen Library

#include "Eigen/Dense"

using namespace Eigen;

//#include <eigen3/Eigen/Core>

#include "openmp1.hpp"
#include "syn_dbg.hpp"

// QAT headers

#include "Variable.h"
#include "Argument.h"

void test_debug() {
  debug("A message");
  // passing in arguments like printf
  debug("I am %d years old ", 37);
}

double f (int i) {
  int j, start = i * (i + 1) / 2, finish = start +  i;
  double return_val = 0.0;

  for (j = start; j <= finish; j++) {
    return_val += sin(j);
  }
  return return_val;
}

//OMP AA(5);
using namespace Genfun;

int main () {
  Genfun::Variable X;
  std::cout << X(3.14) << std::endl;
  
  Matrix2d a;
  a << 1, 2,
       3, 4;
  MatrixXd b(2,2);
  b << 2, 3,
       1, 4;
  std::cout << "a + b =\n" << a + b << std::endl;

  // Computing the inverse and the determinant

  Matrix3f A;

  A << 1, 2, 1,
    2, 1, 0,
    -1, 1, 2;
  
  std::cout << "Here is the matrix A:\n" << A << std::endl;
  std::cout << "The determinant of A is " << A.determinant() << std::endl;
  std::cout << "The inverse of A is:\n" << A.inverse() << std::endl;
 
  /*
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
	if (a[i] > a[i+1]) {
	  tmp = a[i+1];
	  a[i+1] = a[i];
	  a[i] = tmp;
	}
      }
    } 
  }
  */
  // AA.addup();
  // test_debug();
  return 0;
}
