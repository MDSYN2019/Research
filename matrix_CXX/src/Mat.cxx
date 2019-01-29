#ifndef __SIMPLE_MATRIX_CPP
#define __SIMPLE_MATRIX_CPP

#include <vector>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "Mat.hpp"

// Default constructor

template <typename Type>
SimpleMatrix<Type>::SimpleMatrix() {
  // No need for implementation, as the vector "mat"
  // will create the necessary storage  
}
// Constructor with row/col specification and default values
template <typename Type>
SimpleMatrix<Type>::SimpleMatrix(const int& rows, const int& cols, const Type& val) {
  for (int i = 0; i < rows; i++) {
    std::vector<Type> col_vec (cols, val);
    mat.push_back(col_vec);
  }
}
// Copy contructor
template <typename Type>
SimpleMatrix<Type>::SimpleMatrix(const SimpleMatrix<Type>::operator=(const SimpleMatrix<Type>& _rhs)) {
  if (this == &_rhs) return *this; // If we point the = to itsself, return the pointer to this 
  return *this; // else return the pointer to this function
}
// Destructor
template <typename Type>
SimpleMatrix<Type>::~SimpleMatrix() {
  // No need for implementation, as there is no
  // manual dynamic memory allocation
}

// Matrix access method, via copying
template <typename Type>
SimpleMatrix<Type> SimpleMatrix<Type>::get_mat() const {
  return mat;
}
// Matrix access method, via rows and columns index
template <typename Type>
Type& SimpleMatrix<Type>::value(const int& row, const int& col) {
  return mat[row][col];
}
// QSMatrix

/*
  Our task with the source file is to implement all of the methods outlines in the header file.
  In particular we need to implement methods for the following:

  - Constructors, destructor and assignment operator
  - Matrix mathematial methods : Additon, subtraction, multiplication and the transpos 
  - Matrix/ 
 */

// Parameter construction
template<typename T>
QSMatrix<T>::QSMatrix(unsigned _rows, unsigned _cols, const T& _inital) {
  mat.resize(_rows);
  for (unsigned i = 0; i < mat.size(); i++) {
    mat[i].resize(_cols, _initial);
  }
  rows = _rows;
  cols = _cols;
}

// Copy constructor
template<typename T>
QSMatrix<T>::QSMatrix(const QSMatrix<T>& rhs) {
  mat = rhs.mat;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// (Virtual) destructor
template<typename T>
QSMatrix<T>::~QSMatrix() {}

// Assignment operator
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator=(const QSMatrix<T>& rhs) {
  if (&rhs == this)
    return *this;
}

template<typename T>
QSMatrix<T>& QSMatrix<T>::operator+(const QSMatrix<T>& rhs) {
  QSMatrix result(rows, cols, 0.0);
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; i < cols; i++) {
      result(i,j) = this->mat[i][j] + rhs(i,j);
    }
  }
  return result;
}

// Cumulative addition of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator+=(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      this->mat[i][j] += rhs(i,j);
    }
  }
  return *this;
}

// Left multiplication of this matrix and another

template <typename T>
QSMatrix<T> QSMatrix<T>::operator*(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  QSMatrix result(rows, cols, 0.0);

  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; i < cols; j++) {
       for (unsigned k = 0; k < cols; k++) {
	 result(i,j) += this->mat[i][k] * rhs(i,j);
       }
    }
  }
  
  return result;
}

// Cumulative left multiplication of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator*=(const QSMatrix<T>& rhs) {
  QSMatrix result = (*this) * rhs;
  (*this) = result;
  return *this;
}

// Matrix/scalar addition 
template <typename T>
QSMatrix<T> QSMatrix<T>::operator+(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[i][j] + rhs;
    }
  }
  return result;
}

// Matrix/scalar subtraction
template <typename T>
QSMatrix<T> QSMatrix<T>::operator-(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0 j < cols; j++) {
      result(i,j) = this->mat[i][j] - rhs;
    }
  }
  
  return result;
}

// Matrix/scalar multiplication
template <typename T>
QSMatrix<T> QSMatrix<T>::operator*(const T& rhs) {
  QSMatrix result(rows, cols, 0.0);
  for (unsigned i = 0; i < rows; i++) {
    for (unsigned j = 0; j < cols; j++) {
      result(i,j) = this->mat[i][j] * rhs;
    }
  }
  return result;
}


#endif
