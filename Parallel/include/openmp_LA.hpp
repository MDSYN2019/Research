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

// Typedefs

typedef Eigen::Matrix<double, 4, 4> Matrix4x4;
typedef Eigen::Matrix<double, 3, 3> Matrix3x3;
typedef Eigen::Matrix<double, 2, 2> Matrix2x2;

template <typename T> class SYN_Mat : public CppUnit::TestCase {
public:
  SYN_Mat(unsigned _rows, unsigned _cols, const T& _initial); // 1 
  SYN_Mat(const SYN_Mat<T>& alloc);
  // Operator overloading
  SYN_Mat<T>& operator=(const SYN_Mat<T>& alloc);
  virtual ~SYN_Mat(); // virtual destructor
  
  // Matrix mathematical operations

  SYN_Mat<T> operator+(const SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator+=(const SYN_Mat<T>& rhs);
 
  SYN_Mat<T> operator-(const SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator-=(const SYN_Mat<T>& rhs);

  SYN_Mat<T> operator*(const SYN_Mat<T>& rhs);
  SYN_Mat<T>& operator*=(const SYN_Mat<T>& rhs);

  SYN_Mat<T> transpose();

  // Matrix/scalar operations
  SYN_Mat<T> operator+(const T& rhs);
  SYN_Mat<T> operator-(const T& rhs);
  SYN_Mat<T> operator*(const T& rhs);
  SYN_Mat<T> operator/(const T& rhs);


  // Matrix Vector operations
  std::vector<T> operator*(const std::vector<T>& rhs);
  std::vector<T> diag_vec();
  
  T& operator()(const unsigned& row, const unsigned& col);
  const T& operator()(const unsigned& row, const unsigned& col) const;

  // Access the row and column sizes

  unsigned get_rows() const;
  unsigned get_cols() const;
  
  
  /*
    Testing functions
   */

  void test1();
  void test2();
  
private:
  Matrix4x4 l_4, u_4, p_4; // 4 by 4 
  Matrix3x3 l_3, u_3, p_3; // 3 by 3 
  Matrix2x2 l_2, u_2, p_2; // 2 by 2   
  int N;
  //  std::vector<double> c_star(N, 0.0);
  // std::vector<double> d_star(N, 0.0);
  std::vector<std::vector<T> > mat;
  unsigned rows;
  unsigned cols;

  // double S; // Option price
  // double K; // Strike price 
  //double r; // Risk-free price 
  // double v; // Volatility of the underlying (20 %) 
  //  double T; // One year until expiry
};

/*
class ProbDist {
public:
  // Constructor/Destructor section
  ProbDist();
  ProbDist(const ProbDist<T>& alloc);
  // Operator overloading
  ProbDist& operator=(const ProbDist& alloc);
  virtual ~ProbDist(); // virtual destructor
  
  // Functions

  double norm_pdf(const double& x);
  double norm_cdf(const double& x);
  double d_j(const int&, const double&, const double&, const double&, const double&, const double&);
  double call_price(const double&, const double&, const double&, const double&, const double&);
  double put_price(const double&, const double&, const double&, const double&, const double&); 
};
*/

/*
void thomas_algorithm(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<double>&);
void cholesky_decomposition();
*/

// --------------------------------------//


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

template <typename T>
const T& SYN_Mat<T>::operator()(const unsigned& row, const unsigned& col)  const {
    return this->mat[row][col];
}


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



#endif


