#include <iostream>
#include <iomanip> // has the setprecision
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <omp.h>

// Eigen Library

#include "Eigen/Dense"

// cppunit libraries

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

//#include <eigen3/Eigen/Core>

#include "openmp1.hpp"
#include "openmp_LA.hpp"
#include "openmp_dynamicbindingandinheritance.hpp"


// Non-standard namespaces

using namespace Eigen;


//#include "syn_dbg.hpp"


// QAT headers

//#include "Variable.h"
//#include "Argument.h"


double f (int i) {
  int j, start = i * (i + 1) / 2, finish = start +  i;
  double return_val = 0.0;
  for (j = start; j <= finish; j++) {
    return_val += sin(j);
  }
  return return_val;
}

/*
void A() {

  
  Matrix2f A;
  A << 1, 2, 2, 3;

   SelfAdjointEigenSolver<Matrix2f> eigensolver(A);
   if (eigensolver.info() != Success) abort();

   std::cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
   std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
        << "corresponding to these eigenvalues:\n"
	     << eigensolver.eigenvectors() << std::endl;

   SYN_Mat<double> mat1(10, 10, 10);
   SYN_Mat<double> mat2(10, 10, 10);
   SYN_Mat<double> mat3 = mat1 + mat2;

   for (int i = 0; i < mat3.get_rows(); i++) {
     for (int j = 0; j < mat3.get_rows(); j++) {
       std::cout << mat3(i,j) << ", "; 
     }
     std::cout << std::endl;
   } 
}

*/

//OMP AA(5);
//using namespace Genfun;

// Matrix2d a;
//  a << 1, 2,
//       3, 4;
//  MatrixXd b(2,2);
//  b << 2, 3,
//       1, 4;
//  std::cout << "a + b =\n" << a + b << std::endl;

  // Computing the inverse and the determinant

//  Matrix3f A;
//  A << 1, 2, 1,
//    2, 1, 0,
//    -1, 1, 2;
  
// std::cout << "Here is the matrix A:\n" << A << std::endl;
//std::cout << "The determinant of A is " << A.determinant() << std::endl;
// std::cout << "The inverse of A is:\n" << A.inverse() << std::endl;
//MatrixXf p(3,3);
//p << 12, -51, 4,
//     6 , 167 ,-68,
//     -4, 24 , -41;

// A local class - one where there is a class inside the function


int main (void) {
 
  std::vector<Core*> students;
  
  Core* record;
  char ch;
  std::string::size_type maxlen = 0;
  
  //std::string::size_type maxlen = 0;

  // Read and store the data
  // U = undergrad
  // G = Grad

  
  while (std::cin >> ch) {

    if (ch == 'U') {
      record = new Core(); // allocate a Core object  - pointer to a class on the heap 
    } else {
      record = new Grad(); // allocat a Grad object - pointer to a class on the heap 
    }
    record->read(std::cin);          // `virtual' call
    maxlen = std::max(maxlen, record->name().size());// dereference
    students.push_back(record);
  }

  std::sort(students.begin(), students.end(), compare_Core_ptrs);

  // Write the names and grades

  for (std::vector<Core*>::size_type i = 0; i != students.size(); ++i) {
    // students[i] is a pointer that we dereference to call the functions
    std::cout << students[i]->name();
    try {
      double final_grade = students[i]->grade();

      std::streamsize prec = std::cout.precision();
      std::cout << std::setprecision(3) << final_grade << std::setprecision(prec) << std::endl;
	
    } catch (std::domain_error e) {
      std::cout << e.what() << std::endl;
    }
    delete students[i];
  }
  
  return 0;
}
