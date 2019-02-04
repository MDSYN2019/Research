 /*
 *
 *                This source code is part of
 *                    ******************
 *                    ***   SYNMD  ***
 *                    ******************
 *                 molecular modeling library
 *
 * Copyright (c) 2018, Sang Young Noh
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of Artistic License:
 *
*/

/*

-> Parameter input with completeness and consistency checks

-> Runtime array allocation, with array sizes determined by the actual system size

-> Initialization of variables 

-> The main loop which cycles through the force comptuations and trajectory integration, and performs data collection at specified intervals

*/

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cmath>

#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <array>

// Using Boost libraries

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/tuple/tuple.hpp>

// Custom headers

#include "MD_object.hpp"
#include "MD_libraries_structures.hpp"

MD_structures::MD_structures() {} // We only have one constructor so we don't need to specifically define anything
MD_structures::MD_structures(int N) {
  nMols = N;
}

// We only have one constructor so we don't need to specifically define anything
void MD_structures::allocMat() {}
void MD_structures::printArray() {}
void MD_structures::LeapfrogStep(int part) {
  if (part == 1) { 
    for (int n = 0; n < nMols; n++) {
      // Add vectors here
      positionArray.at(n).rv = (0.5 * deltaT * positionArray.at(n).ra);
      positionArray.at(n).r = (0.5 * deltaT * positionArray.at(n).rv);
    }
  }
  else {
    for (int n = 0; n < nMols; n++) {
      positionArray.at(n).rv = (0.5 * deltaT * positionArray.at(n).ra); 
    }
  }
} // Iterator

void MD_structures::computeForce() {
  //  std::vector<Mol> positionArray;
  for (int i = 0; i <= nMols; i++) {    
    mol.r(0) = 0.0;     // Mol mol;
    mol.r(1) = 0.0;     // Mol mol;
    positionArray.push_back(mol); // push back inital type of data, then we initialize it 
  }
  
  for (int j1 = 0; j1 < nMols - 1; j1++) {
    for (int j2 = j1+1; j2 < nMols; j2++) {
      dr(0) = positionArray.at(j1).r(0) - positionArray.at(j2).r(0); 
      dr(1) = positionArray.at(j1).r(1) - positionArray.at(j2).r(1);  
      // x-coordinates block - minimum image convention   
      if (dr(0) >= 0.5 * region(0)) {
	dr(0) -= region(0);
      } else if (dr(0) >= 0.5 * region(0)) {
	dr(0) += region(0);
      }
      // y-coordinates block - minimum image convention 
      if (dr(1) >= 0.5 * region(1)) {
	dr(1) -= region(1);
      } else if (dr(0) >= 0.5 * region(1)) {
	dr(1) += region(1);
      }
      rr = (dr(0) * dr(0)) + dr(1) * dr(1);
      if (rr < rrCut) {	
	rri = 1.0/rr; // 
	rri3 = rri * rri * rri; // 
	fcVal = 48. * rri3 * (rri3 - 0.5) * rri; //
	positionArray.at(j1).ra(0) += fcVal * dr(0);
	positionArray.at(j1).ra(1) += fcVal * dr(1);
	positionArray.at(j2).ra(0) += fcVal * dr(0);
	positionArray.at(j2).ra(1) += fcVal * dr(1);
	uSum += 4.0 * rri3 * (rri3 - 1.0) + 1.0;
      }      
    }
  }
}

void MD_structures::initCoords() {

  gap(0) = region(0) / initUcell(0);
  gap(1) = region(1) / initUcell(1);
  
  n = 0;

  for (int ny = 0; ny <= initUcell(1); ny++) {
    for (int nx = 0; ny <= initUcell(0); nx++) {
      c(0) = nx + 0.5;
      c(1) = ny + 0.5;
      c = c * gap;
      c = c + -0.5 * region;          
    }   
  }
}

void MD_structures::initVels() {

  int n;
  // Initialized the VSum as 0 
  vSum(0) = 0.0;
  vSum(1) = 0.0;

  for (int ny = 0; ny <= initUcell(1); ny++) {
    for (int nx = 0; ny <= initUcell(0); nx++) {
      mol.r(0);
      mol.r(1);
      
    }
  }
}

void MD_structures::SetParams()
{
  // Page 46 
  // TODO

  rCut = pow(2., 1./6);
  region(0) = s1 * vector(0);
  region(0) = s2 * vector(0);
  
  nMol = placeholder;
  valMag = placeholder;
}

/*
// Calculating Pi through monte-carlo algorithm

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

    if (((x*x) + (y*y)) < 1 ) {
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


*/
