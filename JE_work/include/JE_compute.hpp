#ifndef __JE__
#define __JE__

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

/* STL iterators -  typedefs for iterators, and typedef for vectors and tuples */

typedef boost::tuple<int, double, double, double> tuple; 
typedef std::vector<boost::tuple<int, double, double, double> > tupleList; 
typedef std::vector<int>::iterator intIter; // Integer iterator 
typedef std::vector<double>::iterator doubleIter; // Double iterator
typedef std::vector<std::string>::const_iterator stringIter; // String iterator

class JarzynskiFreeEnergy {  
public:
  JarzynskiFreeEnergy(); // Default Constructor  
  JarzynskiFreeEnergy(int, double); // Default Constructor with parameters  
  ~JarzynskiFreeEnergy(); // Destructor 

  void vecProcess();
  void resetIndex();
  void read(std::string input);  
  // Free energy calculation functions from work distribution
  
  double JERaw(std::vector<double> *JEVector); // Raw JE interpreter
  double JETaylor(std::vector<double> *JEVector); // Taylor Series (Second term) JE interpreter  
  double JEprocessVector(int, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector);
  double alpha(double, double, double);  

private:  

  /*

    Scientific constants used in this work 
  
  */

  double BOLTZMANN = 0.0019872041; // units for the boltzmann constant are in kcal mol^-1 //
  double Temperature; // temperature
  int numberOfPullFiles;
  
  int index = 0;
  int nLines = 0;    
  int numberOfTrajectories = 0;  
  double force; // force
  double work; // work
  double z; // z coordinate 
  double bilayerCOM; // Bilayer COM z coordinate
  int number; // number
  
  // For the purpose of IO of files

  std::string filename;
  std::string line;
  std::ifstream myfile;  // ifstream

  // STL containers
  std::vector<int> lineNumberVector; // index vector
  std::vector<double> coordinateZVector; // z coordinates vector
  std::vector<double> bilayerCOMVector; // bilayer COM vector
  std::vector<double> forceVector; // force vector
  std::vector<double> workVector; // work vector   
  std::vector<double> JERawVector; // Storing the work for the raw JE interpreter 
  std::vector<double> JETaylorVector; // Storing the work ffor the taylor series JE interpreter
  tupleList JERawCoordinateBin; // vector for storing raw JE tuples 
  tupleList JETaylorCoordinateBin; // vecto for storing taylor series JE tuples
};


#endif 
