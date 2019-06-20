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


#include "openmp_LA.hpp"

// 1

/*

Our task with the source file is to implement all of the methods outlines in the header file. In particular 
we need to implement methods for the following:

 */
template <typename T>
SYN_Mat<T>::SYN_Mat(unsigned _rows, unsigned _cols, const T& _initial) {
  mat.resize(_rows);

  for (unsigned i = 0; i < mat.size() ; i++) {
    mat[i].resize(_cols, _initial);
  }

  rows = _rows;
  cols = _cols;
}

template <typename T>
SYN_Mat<T>::SYN_Mat(const SYN_Mat<T>& alloc) {
  mat = alloc.mat;
  rows = alloc.get_rows();
  cols = alloc.get_cols();
}

// There is no dynamic memory allocation, we don't need to do anything. We can let the compuler handle the destruction of the individual type
// members
template <typename T>
SYN_Mat<T>::~SYN_Mat() {}

// Defining operators when interacting on SYN_Mat classes of another instantiation

template <typename T>
SYN_Mat<T>& SYN_Mat<T>::operator=(const SYN_Mat<T>& rhs) {

  if (&rhs == this) {
    return *this;
  }
  
  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols(); 
  mat.resize(new_rows);

  for (unsigned i = 0; i < mat.size(); i++) {
    mat[i].resize(new_cols);
  }

  for (unsigned i = 0; i < new_rows; i++) {
    for (unsigned j = 0; j < new_cols; j++) {
      mat[i][j] = rhs(i,j);
    }
  }

  rows = new_rows;
  cols = new_cols;

  return *this;
}


/*
Mathematical Operations Implementation

The next part of the implmentation concerns the methods overloading the binary operators that allow 
matrix algebra such as addition, subtraction and militplication. There are two types of operators to be overleaded here.

*/

template <typename T>
SYN_Mat<T> SYN_Mat<T>::operator+(const SYN_Mat<T>& rhs) {
  SYN_Mat result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[i][j] + rhs(i,j);
    }
  }

  return result;
}

template <typename T>
SYN_Mat<T>& SYN_Mat<T>::operator+=(const SYN_Mat<T>& rhs) {

  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      this->mat[i][j] += rhs(i,j);
    }
  }
  return *this;
}

template <typename T>
SYN_Mat<T>& SYN_Mat<T>::operator-(const SYN_Mat<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  SYN_Matresult(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[i][j] - rhs(i,j);
    }
  }
  return result; 
}

template <typename T>
SYN_Mat<T>& SYN_Mat<T>::operator-=(const SYN_Mat<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      this->mat[i][j] -= rhs(i,j);
    }
  }
  return *this;
}

template <typename T>
SYN_Mat<T>& SYN_Mat<T>::operator*(const SYN_Mat<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  SYN_Mat result(rows, cols, 0.0);
  
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
       for (unsigned k = 0; k < cols; k++) {
	 result(i,j) += this->mat[i][k] * rhs(k,j);
       }
    }
  }
  return result;
}


template <typename T>
SYN_Mat<T>& SYN_Mat<T>::operator*=(const SYN_Mat<T>& rhs) {
  SYN_Mat result = (*this) * this;
  (*this) = result;
  return *this;
}


// calculate a transpose of this matrix
template <typename T>
SYN_Mat<T> SYN_Mat<T>::transpose() {
  SYN_Mat result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[j][i];
    }
  }
  return result;
}


