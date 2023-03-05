#ifndef __analysis__
#define __analysis__

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// boost libraries
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string.hpp>

class analysis {
public:
  analysis(std::string); // constructor 
  void read_file();
  void generate_parameter();
  void store_scaling_parameters();
  void test_parameters();
  ~analysis(); // destructor 
private:
  // private variables - mostly vectors 
  std::string filename; 
  std::vector<std::vector<std::string> > mass_store;
  std::vector<std::vector<std::string> > pair_coeff_store;
  std::vector<std::vector<std::string> > bond_coeff_store;
  std::vector<std::vector<std::string> > angle_coeff_store;
  std::vector<std::string> placeholder;
  std::vector<std::vector<std::string> > placeholder_II;
  // iterators
  std::vector<std::string>::iterator i;
  std::vector<std::string>::iterator stringiter;
  std::vector<std::vector<std::string>>::iterator row;
};

#endif
