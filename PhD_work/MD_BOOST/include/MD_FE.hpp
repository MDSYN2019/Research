/*
-------------------------------------------------
Jarzynski-Equality Algorithm based on the  Hummer-Szabo implementation/correction
VERSION: 1.0.1
-------------------------------------------------
This work directly implements work that was published in PNAS March 27, 2001. 98 (7) 3658-3661:

"Free energy reconstruction from nonequilibrium single-molecule pulling e
xperiments"

LINK: 'http://www.pnas.org/content/98/7/3658' 

Hummer and Szabo implemented 

Here, I've included the raw Jarzynski equality interpretations along with the HS implementation. 

This program uses the SMD module in LAMMPS - the input file format needs to follow:
- - - - - -

*/


#ifndef __MD_FE__
#define __MD_FE__

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <cmath>
#include <cassert>

// STL libraries

#include <string>
#include <vector>
#include <list>
#include <array>
#include <map>
#include <algorithm>

// Boost libraries

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

typedef boost::tuple<int, double, double, double> tuple;
typedef std::vector< boost::tuple<int, double, double, double> > tuple_list;

class Jarzynski_Free_energy {
public:
  Jarzynski_Free_energy(); // Default Constructor
  ~Jarzynski_Free_energy(); // Destructor   
  void readInput();
  void setVar(); // Set the spring constant and others
  void calculate();
  void setIndex();
  void vecProcess();

  // Free energy calculation from work distribution
  
  double JE_raw(double);
  double JE_taylor(double,int);
  double JE_alternative(double);

  double read(std::string input);
private:
  /* Private variables */
  // integers 

  int numbins = 0; // The number of bins we want to divide the trajectory with  
  int index = 0;
  int nlines = 0;    
  int numberOfTrajectories = 0;  
  int numBins = 0;
  // doubles
  
  double boltzmann = 0.0019872041; // units are kcal mol^-1 //
  double spring_const = 0.0; // what units - Kcal/mole-Angstrom 

  double f; // force
  double w; // work
  double z; // z coordinate 
  double U; // second column
  int number; // number 
 
  // For the purpose of IO of files
  
  std::string filename;
  std::string line;
 
  // ifstream
  std::ifstream myfile;

  // STL containers
  std::vector<double> bin; // lower interval of a bin - distances
  std::vector<double> timeBin; 
  std::vector<double> HSwork;
  std::vector<double> normalwork;
  std::vector<double> work1;
  std::vector<double> pos1;
  std::vector<double> G0;
  
  std::vector<int> num;
  std::vector<int> cz; // z coordinates
  std::vector<double> Un;
  std::vector<double> force;
  std::vector<double> work;
  std::vector<std::string> vecOfFiles;

  tuple_list coordinateBin;

  // STL iterators
  
  typedef std::vector<double>::iterator diter;
  // A const iterator is an iterator that points to const value (like a const T* pointer). Dereferencing it returns
  // it returns a reference to a constant value (const T&) and prevents the modification of the referenced value.    
  typedef std::vector<std::string>::const_iterator striter;
  
  // Boost containers
  
  tuple E_tuple;  

};



#endif
