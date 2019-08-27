#include <Eigen/Core>
#include <iostream>
#include <vector>
#include "statistics.h"


//! MPI headers
#include "MPI_broadcast.hpp"
#include "openmp_LA.hpp"
#include "MPI_IO.hpp"

// Accelerated C++
#include "Core.h"

int main(void) {


  
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
  
  return 0;
}
