#include "analysis.h"
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/split.hpp>
#include <iostream>
#include <typeinfo>
#include <vector>

using boost::is_any_of;

std::vector<std::vector<std::string> > recorded_file;

analysis::analysis(std::string input_file) {
  /*
    constructor 
  */
  filename = input_file;
}

analysis::~analysis() {
  std::cout << "class instance destroyed" << std::endl;
}

void analysis::read_file() {
  /*
    read the file of input_file 
  */
  std::ifstream infile(filename);
  for (std::string line; getline(infile, line);) {
    // ignore all lines that start with hash 

    if (line[0] == '#') {
      continue;
    }    

    else {
      boost::algorithm::split(placeholder, line, is_any_of("\t "), boost::token_compress_on);

      // store mass values 

      if (placeholder[0] ==  "mass") {
	mass_store.push_back(placeholder);
      }

      // store pair coefficient values 

      else if (placeholder[0] == "pair_coeff") {
	pair_coeff_store.push_back(placeholder);	
      }

      // store bond coefficient values 

      else if (placeholder[0] == "bond_coeff") {
	bond_coeff_store.push_back(placeholder);
      }
      
      // store angle coefficient values 

      else if (placeholder[0] == "angle_coeff") {
	angle_coeff_store.push_back(placeholder);
      }   
      placeholder.clear();
    }
  }
}
    
void analysis::generate_parameter() {
  /*
    produce a working parameter file that 
    has been modified with the correct parameters 

    the ones we have the keep concious about at the moment are:

    Mimic - Real 
    ------------
    13-7 
    12-6 
    9-3 
    10-4
    TODO - need to change the conditional input so that the parameters can easily be modified 
  */  

  for (row = pair_coeff_store.begin(); row != pair_coeff_store.end(); row++) {   
    if (row->at(1) == "7" && row->at(2) == "7") { 
      std::cout << " " << row->at(0) << " " << row->at(1) << " " << row->at(2)
                << " " << row->at(3) << " " << row->at(4) << " " << row->at(5) <<  " "  << row->at(6) << std::endl;
    }
  }
}

