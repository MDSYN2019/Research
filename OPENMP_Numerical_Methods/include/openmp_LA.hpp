#ifndef __MATRIX_CPP__
#define __MATRIX_CPP__


#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include <string>
#include <memory>
#include <vector>
#include <cmath>

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

template <typename T> class SYN_Mat {
public:
  SYN_Mat(int);
  SYN_Mat(const SYN_Mat<T>& alloc);
  // Operator overloading
  SYN_Mat<T>& operator=(const SYN_Mat<T>& alloc);

  virtual ~SYN_Mat(); // virtual destructor
  
  // Matrix mathematical operations

  SYN_Mat<T>& operator=(SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator+(SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator-(SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator-=(SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator*(SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator*=(SYN_Mat<T>& rhs);
  SYN_Mat<T> transpose();

  // Matrix/scalar operations

  SYN_Mat<T> operator+(const T& rhs);
  SYN_Mat<T> operator-(const T& rhs);
  SYN_Mat<T> operator*(const T& rhs);
  SYN_Mat<T> operator/(const T& rhs);

  void thomas_algorithm(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<double>&);
  
private:
  int N;
  std::vector<double> c_star(N, 0.0);
  std::vector<double> d_star(N, 0.0);
  unsigned rows;
  unsigned cols;
};
#endif
