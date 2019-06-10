#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <string>
#include <memory>
#include <vector>

#ifndef _OPENMP_
#include <omp.h>

#endif

// Main header to include 

#include "openmp1.hpp"
#include "openmp_LA.hpp"

// QAT headers

#include "Argument.h"

/* CPPunit tests */

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

// Eigen Library

#include "Eigen/Dense"
#include "Eigen/LU"

typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Matrix<double, 2, 2> Matrix2x2;


SYN_Mat<T>::SYN_Mat() {} // default constructor
SYN_Mat<T>::~SYN_Mat() {}

SYN_Mat<T>::thomas_algorithm(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<double>&) {
  c_star[0] = c[0] / b[0];
  d_star[0] = d[0] / b[0];

  for (int i = 1; i < N; i++) {
    double m = 1.0/ (b[i] - a[i] * c_star[i-1]);
    c_star[i] = c[i] * m;
    d_star[i] = (d[i] - a[i] * d_star[i-1]) * m;
  }

  // This is the reverse sweep, used to update he solution vector f 
  for (int i = N-1; i--> 0;) {
    f[i] = d_star[i] - c_star[i] * d[i+1];
  }
}



