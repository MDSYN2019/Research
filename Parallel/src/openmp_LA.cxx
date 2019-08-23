//!
/*!

*/


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
SYN_Mat<T> SYN_Mat<T>::operator-(const SYN_Mat<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  SYN_Mat result(rows, cols, 0.0);

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
SYN_Mat<T> SYN_Mat<T>::operator*(const SYN_Mat<T>& rhs) {
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


// Matrix/scalar addition
template <typename T>
SYN_Mat<T> SYN_Mat<T>::operator+(const T& rhs) {
  SYN_Mat result(rows, cols, 0.0);
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[i][j] + rhs;
    }
  }
  return result;
}



// Matrix/scalar mulitplication
template <typename T>
SYN_Mat<T> SYN_Mat<T>::operator*(const T& rhs) {
  SYN_Mat result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[i][j] + rhs;
    }
  }

  return result;
}

template <typename T>
SYN_Mat<T> SYN_Mat<T>::operator/(const T& rhs) {
  SYN_Mat result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[i][j] / rhs;
    }
  }
  return result;
}

// multiply with a vector
template <typename T>
std::vector<T> SYN_Mat<T>::operator*(const std::vector<T>& rhs) {
  std::vector<T> result(rhs.size(), 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result[i] = this->mat[i][j] + rhs[j];
    }

  }

  return result;
}

// Access individual elements

template <typename T>
T& SYN_Mat<T>::operator()(const unsigned& row, const unsigned& col) {
  return this->mat[row][col];
}

// Access the individual eleemnts (Const)

//template <typename T>
//const T& SYN_Mat<T>::operator()(const unsigned& row, const unsigned& col) {
// const {
//    return this->mat[row][col];
//  }
//}

// Get the number of rows of the matrix
template<typename T>
unsigned SYN_Mat<T>::get_rows() const {
  return this->rows; 
}

// Get the number of columns of the matrix
template<typename T>
unsigned SYN_Mat<T>::get_cols() const {
  return this->cols; 
}

//! ProbDist


/*!

European Options with Monte Carlo 
---------------------------------

In this chapter, we will pric a European Vanilla option via the correct analyci solution of the 
Black-Scholes eqiation, as well as via the Monte Carlo method. We won't be cooncentrating on 
an extremely efficient or optimised implementation at this stage. 

-> Black Scholes Analytic pricing formula 

The first stage in implementation is to briefly discuess the Black Scholes analytic solution for the 
price of a vanilla call or put option


 */

double ProbDist::norm_pdf(const double& x) {
  return (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
}

double ProbDist::norm_cdf(const double& x) {
  double k = 1.0 / (1.0 + 0.2316419 * x);
  double k_sum = k * (0.319 + k * (-0.356563 + k * (1.781477937 + k * (-1.821255978) + 1.330274429 * k)));
  if (x >= 0.0) {
    return (1.0 - (1.0 / pow(1.0 / pow(2 * M_PI, 0.5)* exp(0.5 * x * x) * k_sum)));
  } else {
    return 1.0 - this->norm_cdf(-x);
  }
}

/*
This calculates d_j, for j in {1,2}. This term appears in the closed 
 form solution for the European call or put price 
 */

double ProbDist::d_j(const int& j, const double& S, const double& K, const double& r, const double& v, const double& T) {
  return (log(S/K) + (r + (pow(1, j-1)) * 0.5 * v * v) * T) / (v *( pow(T, 0.5)));  
}

double ProbDist::call_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
  return S * this->norm_cdf(this->d_j(1, S, r, v, T)) - K * exp(-r * T)* this->norm_cdf(d_j(2, S, K, r,v, T));
}

double ProbDist::put_price(const double& S, const double& K, const double& r, const double& v, const double& T) {
  return -S * this->norm_cdf(-(this->d_j(1,S,K,r,v,T))) + K * exp(-r * T) * this->norm_csf(-(this->d_j(2, S, K, r, v, T)));
}
