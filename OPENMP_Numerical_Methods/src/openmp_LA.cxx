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

//template constructor 

ProbDist::norm_pdf() {

}
