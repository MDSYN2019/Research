#include <iostream>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <omp.h>

// Eigen Library

#include "Eigen/Dense"
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

//#include <eigen3/Eigen/Core>

#include "openmp1.hpp"
#include "openmp_LA.hpp"
#include "openmp_dynamicbindingandinheritance.hpp"

//#include "syn_dbg.hpp"


// QAT headers

//#include "Variable.h"
//#include "Argument.h"

using namespace Eigen;

double f (int i) {
  int j, start = i * (i + 1) / 2, finish = start +  i;
  double return_val = 0.0;
  for (j = start; j <= finish; j++) {
    return_val += sin(j);
  }
  return return_val;
}


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


void local_function1() {
  class Test {
    static void method() {
      std::cout << "Local class method() " << std::endl;      
    }
  };
}

int main (void) {

  //std::vector<Core> students;
  // Core records;
  //std::string::size_type maxlen = 0;

  // read and store the data

  //  while (records.read(std::cin)) {
  // maxlen = std::max(maxlen, records.name().size()); // get the maxiumum size of the students name
  //  students.push_back(records); // push back the record
  //}

  // Not sure if I need the std::sort
  
  // std::sort(student.begin(), students.end(), compare);
  
  // Write the names and grades

  // for (std::vector<Core>::size_type i = 0; i != students.size(); ++i) {
  //  std::cout << students[i].name() << std::string(maxlen + 1 - students[i].name.size(), ' ');
  //  }

  /*
  try {
   double final_grade = students[i].grade(); // Core::grade
   std::streamsize prec = std::cout.precision();
   std::cout << std::setprecision(3) << final_grade << std::setprecision(prec) << std::endl;
  } catch (std::domain_error e) {
    std::cout << e.what() << std::endl;
  } // This block of code will need to be tested
  */
  
  // Checking the utility of the custom built Vec

  Vec<int> EE;
  EE.push_back(1);

  std::cout << EE[0] << std::endl;

  EE.runTest();
  
  //EE.push_back(3);
  //std::cout << EE[0] << std::endl;
  
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
    
   return 0;
}

/// TODO

   //   Eigen::HouseHolderQR<MatrixXf> qr(p);
   //Eigen::MatrixXf q = qr.householderQ();
   //std::cout << "Q Matrix :\n" << q << std::endl << std::endl ;
  
   int x = 5;

   //# pragma omp parallel num_threads(thread_count) private(x)
     //{
   // int my_rank = omp_get_thread_num();
   // std::cout << "Thread %d ";
   // x = 2 * my_rank + 2;
   // }
  
 
  /*
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
