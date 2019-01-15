#pragma once
#include <vector>
#include <map>
#include <random>


#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

/*
We use the __mt19937 with the default seed as a source of 
randomness. The numbers produced will be the same 
every time the program is run. 

One common method to change this is to seed the current time 
 */

using namespace std;

class MCalg {
public:
  MCalg();
  MCalg(int);
  //~MCalg();
  
  MCalg operator+(MCalg); 
  double randist(); // distribution function 
  void direct_pi(); // approximtion of pi
  void print_pi();
  void print_xy();
private:
  int num;
  int N_hits = 0;
  double x;
  double y;
};
