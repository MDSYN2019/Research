
  /*!
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

    
  This code has the serial code for estimating the Jarzynski estimate. I have used primarily C++ stl-based libraries, with some Boost code as well. Hence, this code will require the boost libraries and in the future will be using the GSL libraries for managing the statistical aspects of the algorithm.
    
  Abbreviation: JE - Jarzynski Equality

  */
#ifndef __JE__
#define __JE__

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

//!
/*! STL iterators -  typedefs for iterators, and typedef for vectors and tuples */  
typedef boost::tuple<int, double, double, double> tuple; /*!< Boost-constructed tuple to store the index, center point in reaction coordinate, lower bound, upper bound and the free energy value */ 
typedef std::vector<boost::tuple<int, double, double, double> > tupleList; /*!< Vector of the boost::tuple */  
typedef std::vector<int>::iterator intIter; /*!< Integer iterator */ 
typedef std::vector<double>::iterator doubleIter; /*!< Double iterator */ 
typedef std::vector<std::string>::const_iterator stringIter; /*!< String iterator */

//! 
/*! The main class */  

class JarzynskiFreeEnergy {  
public:
  JarzynskiFreeEnergy(); /*!< Default Constructor */ 
  JarzynskiFreeEnergy(int, double); /*!< Default Constructor with parameters */ 
  ~JarzynskiFreeEnergy(); /*!< Destructor */ 

  void vecProcess();
  void resetIndex();
  void read(std::string input);  

  /*! Free energy calculation functions from work distribution */
  double JERaw(std::vector<double> *JEVector); /*< Raw JE interpreter */
  double JETaylor(std::vector<double> *JEVector); /*< Taylor Series (Second term) JE interpreter */  
  double JEprocessVector(int, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector); /*< functor which can take a JE algorithm then processing it accordingly */
  double alpha(double, double, double); /*< */   
  friend void MPI_setup();
  friend void MPI_vec_send();
  friend void MPI_parameter_send();
private:
  //! /
  /*! Scientific constants used in this work */
  const double BOLTZMANN = 0.0019872041; /*< units for the boltzmann constant are in kcal mol^-1 */
  double Temperature; /*< Temperature */
  int numberOfPullFiles; /*< The number of work files to compute the free energy values with */
  
  int index = 0;
  int nLines = 0;    
  int numberOfTrajectories = 0;  
  double force; /*< Force */
  double work; /*< Work */
  double z; /*< z coordinate */ 
  double bilayerCOM; /*< Bilayer COM z coordinate */
  int number; /*< Number */
  //! 
  /*! For the purpose of IO of files */
  std::string filename;
  std::string line;
  std::ifstream myfile;  
  //!
  /*! STL containers  */
  std::vector<int> lineNumberVector; /*< Index vector */
  std::vector<double> coordinateZVector; /*< z coordinates vector */
  std::vector<double> bilayerCOMVector; /*< Bilayer COM vector */
  std::vector<double> forceVector; /*< Force vector */
  std::vector<double> workVector; /*< Work vector */   
  std::vector<double> JERawVector; /*< Storing the work for the raw JE interpreter */ 
  std::vector<double> JETaylorVector; /*< Storing the work ffor the taylor series JE interpreter */
  tupleList JERawCoordinateBin; /*< Vector for storing raw JE tuples */
  tupleList JETaylorCoordinateBin; /*< Vector for storing taylor series JE tuples */
};

#endif 
