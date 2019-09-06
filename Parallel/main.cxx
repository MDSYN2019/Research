//! Generic C++ Headers

#include <Eigen/Core>
#include <iostream>
#include <vector>
#include <complex>
#include <limits>

//! MPI Headers

#include "MPI_broadcast.hpp"
#include "openmp_LA.hpp"
#include "MPI_IO.hpp"

//! Accelerated C++ Headers

#include "Core.h"

//! Statistics Headers

#include "statistics.h"

//! QAT Headers

#include "QatGenericFunctions/Variable.h"
#include "QatGenericFunctions/Sin.h"
#include "QatGenericFunctions/NormalDistribution.h"
#include "QatGenericFunctions/Parameter.h"


// ------------------------
// Example - root finding
// -----------------------

double newtonRaphson(double x, Genfun::GENFUNCTION P) {
  double x1 = x;
  while (1) {
    double deltaX = -P(x) / P.prime() (x);
    x += deltaX;
    if (float(x1) == float(x)) break;
    x1 = x;
  }
  
  return x;
}

/*
double f(const std::vector<double> &a, double x) {
  unsigned int N = (a.size() - 1) / 2;
  // allocate and zero the arrays
  std::vector<double> e(2 * N +1, 0.0);
  std::vector<double> q(2*N, 0.0);
  std::vector<double> eold(2*N+1, 0.0);
  std::vector<double> quold(2*N+1, 0.0);
  std::vector<double> cf(2*N, 0.0);

  for (int j = 0; j <= 2*N-1; j++) {
    qold[j] = a[j+1]/a[j];
    cf[0] = a[0];
    cf[1] = -qold[0];
  }

  // Evaluate the coefficients  
}
*/

int main(void) {

  Genfun::Variable X;
  std::cout << X(3.14) << std::endl;
  Genfun::Sin sin;
  
  Genfun::GENFUNCTION f = 1 + 2 * X + X * X * X;
  Genfun::GENFUNCTION FF = (X-1) * (X-2) * (X-3) * (X - M_PI) * (X-4);
 
  std::cout << 1 << " " << f(1) <<  " " << f.prime()(1) << std::endl;
  std::cout << 1 << " " << FF(1) <<  " " << FF.prime()(1) << std::endl;


  // Root Finding

  //  const Genfun::AbsFunction * f = &FF;

  //  for (int i = 0; i < 5; i++) {
  // double x = newtonRaphson(-1.0, *f);
  //  Genfun::GENFUNCTION F1 = (*f) / (X - x);
    //  if (f != &FF) delete f;
    //f = F1.clone();
  //}
  
  // MPIInput MPIobj(3,3);
  // MPIobj.getData(&A, &B, &C);
  // Call MPI

  int A, B, C;

  
  // Furst we create the parameter list
  double S = 100.0;
  double K = 100.0;
  double r = 0.05;
  double v = 0.2;
  double T = 1.0;

  ProbDist probdist;
  double call, put;

  call = probdist.call_price(S, K, r, v, T);
  put = probdist.put_price(S, K, r, v, T);

  std::cout << "call_price: " << call << std::endl; 
  std::cout << "Put_price: " << put << std::endl; 

  std::vector<Core> students;
  Core record;
  std::string::size_type maxlen = 0;

  // read and store the data

  while (record.read(std::cin)) {
    maxlen = std::max(maxlen, record.name().size());
    students.push_back(record);
  }

  std::sort(students.begin(), students.end(), compare);

  // write the names and grades

  

  /*
  
    Exceptions provide a way to react to exception circumstances, like runtime erros, in programs
    by transferring control to special functions called handlers
    
    To catch exceptions, a portion of code is placed under exception inspection. This is done 
    by enclosing that portion of code in a try-block. When an exceptional circumtance arises within
    that block, an exception is thrown that transfers the control to the exception handler. 
    
    An exception is thrown by using the throw keyword from inside the try block. Exception 
      handlers are declared with the keyword catch, which must be placed immediately after the try block
      
  */
  
  
  for (std::vector<Core>::size_type i = 0; i != students.size(); ++i) {
    std::cout << students[i].name() << std::string(maxlen + 1 - students[i].name().size(),  ' ');

    try {
      double final_grade = students[i].grade(); // Core::grade
      std::streamsize prec = std::cout.precision();
      std::cout << std::setprecision(3) << final_grade << std::setprecision(prec) << std::endl;
    } catch (std::domain_error e) { // catches the domain error 
      std::cout << e.what() << std::endl; 
    }
  }

  //  Genfun::Variable X;

  /*
    Managing memory (almost) manually
    ---------------------------------

    When we built our student_info handle class, we combined 
    two separable abstractions. Not only was that class an interface
    to the operations on student records.
    
    What we'd like is to be able to define a class that is similar to 
    student_info, but that is strictly an interface class. Such interface
    classes are common in C++, especially when they interface to an 
    interface hierachy. 


    
    
   */

  
  return 0;
}
