/*
---------------------------------------------------------------------------------
| Jarzynski-Equality Algorithm based on multiple implementations/corrections    |
|                                                                               | 
| VERSION: 0.0.2                                                                |  
---------------------------------------------------------------------------------

This work directly implements work that was published in the following publications:
 
1. "Nonequilibrium equality for free energy differences", Phys Rev Lett, 78, 2690 

2. "Free energy reconstruction from nonequilibrium single-molecule pulling experiments", PNAS, 3658 - 3661, 98, 7, 2001

3. "Free energy calculation from steered molecular dynamics simulations using Jarzynski's equality", Journal of Chemical Physics, 119, 6, 2003
 
4. "Bias and error in estimates of equilibrium free-energy differences from nonequilibrium measurements", PNAS, 12564 - 12569, 100, 22, 2003 


The idea of the steered molecular dynamics / Jarzynski equality is to sample the work values and to assess the free energy change. As work is path-dependent and not a state function, whilst free energy is a equilibrium state function, the main novelty of the equality is that an equilirium property can be extracted from a non-equilbrium experiment.

Here, the experiment was carried out using LAMMPS (https://lammps.sandia.gov/). The simulation of a toluene going through a lipid bilayer was carried out at various velocities, ones that are far out of the `equilibrium' regime (i.e. very fast pulling simulations), and those closer to the `equilbrium' regime (slower simulations).  Each output file prints out the index, reaction coordinate, force and work values respectively.  

Here, there are two main free energy interpreters - the raw Jarzynski equality as seen from reference 1, and the taylor series improved evaluator seen in reference 3. 
 
This code has the serial code for estimating the Jarzynski estimate. I have used primarily C++ stl-based libraries, with some Boost code as well.

Hence, this code will require the boost libraries and in the future will be using the GSL libraries for managing the statistical aspects of the algorithm.

Abbreviation: JE - Jarzynski Equalily
*/


#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <cassert>
#include <numeric>
#include <functional>
#include <exception>

// STL libraries

#include <string>
#include <vector>
#include <list>
#include <array>
#include <map>
#include <algorithm>
#include <stdexcept>
#include <experimental/filesystem>

#include "JE_compute.hpp"

int main(int argc, char *argv[]) {
  std::string filenameString = argv[1]; // Pull file names
  std::string pullDirectory = argv[2]; // Directory for where the pull files are in the computer
  int numberOfFiles = atoi(argv[3]); // Number of pull filenameString
  std::cout << filenameString << " " << pullDirectory << " " << numberOfFiles << " " << argc << " " << std::endl;
  
  try {
    if (argc <= 0 || argc > 4) {
      throw std::invalid_argument("received the wrong amount of input parameters");
    }
  } catch(const std::invalid_argument& e) {
    std::cout << "Usage: filename directory numberoffiles" << std::endl;
    std::cout << "Aborting Program, due to error:" << e.what() << " "  << std::endl;
    exit(1);
  }

  
  JarzynskiFreeEnergy FE; // Class instance  

  int index = 0;
  std::string filename;
  std::string path;
  std::vector<std::string> traj;

  FE.resetIndex();  

  for (int i = 0; i <= numberOfFiles; i++) {     
    index = i + 1; // Filename index
    auto b = std::to_string(index);
    filename = filenameString + "." + b; // Concatenate string  
    traj.push_back(filename); // store the filesname    

  }
  
  // Loop through all the filesnames and work through them 

  for (std::vector<std::string>::const_iterator i = traj.begin(); i != traj.end(); ++i) {

    path = pullDirectory + "/" + std::string(*i);
    std::cout << "Reading pull files from directory:" << path << std::endl;
    FE.read(path);

  }   

  FE.vecProcess(); // Compute the JE interpreters

  return 0;
}
