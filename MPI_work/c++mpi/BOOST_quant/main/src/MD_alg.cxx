/*
Detailed instructions written here.
*/

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <string>

// Using Boost libraries

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/tuple/tuple.hpp>

// Custom headers

#include "MD_object.hpp"

  
MCalg::MCalg() {}

MCalg::MCalg(int a) {
  num = a; // Number of hits for Pi 
}

double MCalg::randist() {

  std::uniform_real_distribution<double> unif(-1,1);
  std::default_random_engine re;

  // boost::random::mt19937 randgen;
  // boost::random::uniform_int_distribution<> dist(-1, 1);

  return unif(re);
}

void MCalg::direct_pi() {  

  for (int i = 0; i < num; i++) {
    
    x = MCalg::randist(); 
    y = MCalg::randist();

    std::cout << x << std::endl;
    std::cout << y << std::endl;
    
    if ( ((x*x) + (y*y)) < 1 ) {
      N_hits = N_hits + 1;
      
    }
  }
}

void MCalg::print_pi() {

  std::cout << N_hits << std::endl;
}

void MCalg::print_xy() {
  std::cout << x << std::endl;
  std::cout << y << std::endl;
}

MCalg MCalg::operator+(MCalg aso) {  

  MCalg brandNew;
  
  brandNew.num = this->num + aso.num;

  return(brandNew); // returns the brandnew sally object 
}
