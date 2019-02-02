//#include "Mat.h"
#include <iostream>
#include <string>
#include <Eigen/Dense>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

class Complex {
  friend bool operator==(const Complex& a, const Complex& b);
  double real, imaginary;
public:
  Complex(double r, double i = 0) : real(r), imaginary(i) {
  }
};

bool operator==(const Complex &a, const Complex &b) {
  return true;
}

class ComplexNumberTest : public CppUnit::TestCase { 
public: 
  ComplexNumberTest( std::string name ) : CppUnit::TestCase( name ) {}
  
  void runTest() {
    CPPUNIT_ASSERT( Complex (10, 1) == Complex (10, 1) );
    CPPUNIT_ASSERT( !(Complex (1, 1) == Complex (2, 2)) );
  }
};

ComplexNumberTest A("df");

int main(int argc, char **argv) {

  A.runTest();
  
  Eigen::MatrixXd p(2,2);

  p(0,0) = 3;
  p(1,0) = 3;
  p(0,1) = 3;
  p(1,1) = 3;
  
  std::cout << p.sum() << std::endl;
  /*
  QSMatrix<double> mat1(10,10,1.0);
  QSMatrix<double> mat2(10,10,2.0);
  QSMatrix<double> mat3 = mat1 + mat2;

  for (int i = 0; i < mat3.get_rows(); i++) {
    for (int j = 0; k < mat3.get_cols(); j++) {
      std::cout << mat3(i,j) << ",";    
    }
    std::cout << std::endl;
  }


  Eigen::Matrix3d p;
  */
  
  return 0;
  
}
