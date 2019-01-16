/*
-------------------------------------------------
Jarzynski-Equality Algorithm based on the  Hummer-Szabo implementation/correction

VERSION: 1.0.1
-------------------------------------------------
This work directly implements work that was published in PNAS March 27, 2001. 98 (7) 3658-3661:

"Free energy reconstruction from nonequilibrium single-molecule pulling experiments"

LINK: 'http://www.pnas.org/content/98/7/3658' 

Hummer and Szabo implemented 

Here, I've included the raw Jarzynski equality interpretations along with the HS implementation. 

- - - - - -
*/

#include <iostream>
#include <fstream>
#include <sstream>
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

#include "MD_FE.hpp"

/* Including the GSL libraries */

//#include <gsl/gsl_math.h>
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_histogram2d.h>
//#include <gsl/gsl_statistics_int.h> 

typedef boost::tuple<int, double, double, double> tuple;
    
// Constructors and Destructors
Jarzynski_Free_energy::Jarzynski_Free_energy() {} // Default Constructor
// Jarzynski_Free_energy(int , int, int); // Input Constructor 
Jarzynski_Free_energy:: ~Jarzynski_Free_energy() {} // Destructor 

void Jarzynski_Free_energy::readInput() {}
void Jarzynski_Free_energy::vecProcess() {  

  int index;

  // Iterators to work with
  typedef std::vector<int>::iterator intIter;
  typedef std::vector<double>::iterator dubIter;

  std::cout << "Which Free energy calculator would you like to use?" << std::endl;
  std::cout << "Insert Option:" << std::endl;
  // getline(std::cin, index);
 
  for (intIter i = num.begin(); i != num.end(); ++i) {
    std::cout << *i << std::endl;
  }
  
  double max_z = *max_element(cz.begin(), cz.end()); // Define minimum z coordinate
  double min_z = *min_element(cz.begin(), cz.end()); // Define maximum z coordinate     
  for (int i = min_z; i != max_z; i++) {

    tuple temp{i, i - 0.5, i + 0.5, 0.0};

    coordinateBin.push_back(temp);  

  }
  
  for (tuple_list::const_iterator i = coordinateBin.begin(); i != coordinateBin.end(); ++i) {
    std::cout << "Bin Center: " << i->get<0>() << " " << i->get<1>()  << " " << i->get<2>() << " " << i->get<3>() << std::endl;
  }
  //  switch(index) {
  // case 1:
  //}
  // Store the value into bins  
}

void Jarzynski_Free_energy::setVar() {}  // Set the spring constant and others
void Jarzynski_Free_energy::calculate() {}
void Jarzynski_Free_energy::setIndex(){
    
  bin.clear(); // lower interval of a bin - distances
  timeBin.clear(); 
  HSwork.clear();
  normalwork.clear();
  work1.clear();
  pos1.clear();
  G0.clear();
  
  num.clear();
  cz.clear(); // z coordinates
  Un.clear();
  force.clear();
  work.clear();
}

double Jarzynski_Free_energy::JE_raw(double W) {
  /*
   Raw Jarzynski equality iterpretation
  */
  double E;
  E = exp(-beta*W);
  return E;
}

double Jarzynski_Free_energy::JE_taylor(double W, int N) {
  /*
    Second order taylor expansion interpretation
  */
  double E;
  E = -beta*(W/N) + ((beta*beta)/2)*((W*W)/N - (W/N)*(W/N));
  return E;
}

double Jarzynski_Free_energy::read(std::string input) {
  myfile.open(input, std::ifstream::in);
  //  assert(myfile); // Ensure that the file is open
  // If there is a error with this input file, then print out an error
  std::cout << myfile.is_open() << std::endl;
  if (myfile.is_open() == 0) {    
    std::cout << "Cannot open files" << std::endl;     
      exit(1);
  }
  
  else {	
    std::cout << "Opening file:" << " " << input << std::endl;	
  }
  
  //while (getline(myfile, line)) {
  // /// myfile >> number >> z >> U >> f >> w;
  //  nlines++;
  // std::cout << line << " " << std::endl;
  // }
  while (myfile >> number >> z >> U >> f >> w) {
    std::cout << " " <<  number << " "  << " " <<  z << " " << " " <<  U << " " <<  f << " " << w << std::endl; 
    num.push_back(number); // index
    cz.push_back(z); // z coordinates
    Un.push_back(U); // .. 
    force.push_back(f); // .. 
    work.push_back(w); // ..
    nlines++;
  }

  myfile.close();
  std::cout << "Number of lines:" << " " << nlines << " " << std::endl;
  return nlines;
}

int main(int argc, char *argv[]) {  
  std::string path;
  Jarzynski_Free_energy FE; // Class instance
  
  //FE.vecProcess();
  int nlines = 0;
  int numberOfTrajectories = 40;
  int index = 0;
  
  typedef std::vector<std::string>::const_iterator striter; 
  std::string filename = "pull";
  std::string f;
  std::string bb;
  std::vector<std::string> traj;
  std::string line;

  for (int i = 0; i < numberOfTrajectories; i++) {     
    index = i + 1;
    auto b = std::to_string(index);
    f = filename + "." + b; 
    traj.push_back(f);
  }
  
  for (std::vector<std::string>::const_iterator i = traj.begin(); i != traj.end(); ++i) {
    path = std::string("/home/noh/Desktop/Program/GIT/MD_BOOST/pullfiles/pull") +  "/" + std::string(*i);
    FE.read(path);
  }
  
  return 0;
}

    /*
      template<class InputIt, class T>
      InputIt find(InputIt first, InputIt last, const T& value)
      {
      for (; first != last; ++first) {
      if (*first == value) {
	return first;
      }
    }
    return last;
  }

    */


