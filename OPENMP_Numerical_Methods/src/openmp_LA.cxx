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


SYN_Mat<T>::SYN_Mat() {} // default constructor
SYN_Mat<T>::~SYN_Mat() {}

SYN_Mat<T>::thomas_algorithm(const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, const std::vector<double>&, std::vector<double>&) {}



