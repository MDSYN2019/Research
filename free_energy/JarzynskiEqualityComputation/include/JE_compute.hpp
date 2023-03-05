  /*!
  ---------------------------------------------------------------------------------
  | Jarzynski-Equality Algorithm based on multiple implementations/corrections    |
  |                                                                               | 
  | VERSION: 0.0.3                                                                |  
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
#include "mpi.h"

// cppunit tests                                                                                                                                      
//#include <cppunit/extensions/TestFactoryRegistry.h>
//#include <cppunit/CompilerOutputter.h>
//#include <cppunit/ui/text/TestRunner.h>
//#include <cppunit/TestFixture.h>
//#include <cppunit/extensions/HelperMacros.h>

//!
/*! STL iterators -  typedefs for iterators, and typedef for vectors and tuples */  
typedef boost::tuple<int, double> indexDoubletuple; // .. 
typedef boost::tuple<int, double, double, double, double> tuple; /*!< Boost-constructed tuple to store the index, center point in reaction coordinate, lower bound, upper bound and the free energy value */ 
typedef std::vector<boost::tuple<int, double, double, double, double> > tupleList; /*!< Vector of the boost::tuple */  
typedef std::vector<int>::iterator intIter; /*!< Integer iterator */ 
typedef std::vector<double>::iterator doubleIter; /*!< Double iterator */ 
typedef std::vector<std::string>::const_iterator stringIter; /*!< String iterator */

// structs cannot be externed 

typedef struct {
  double BM;
  double T; 
} parameterData;

typedef struct {
double val;
double err;
} FE;

//! 

/*! The main class */  

class JarzynskiFreeEnergy {  
public:
  JarzynskiFreeEnergy(); /*!< Default Constructor */ 
  JarzynskiFreeEnergy(int, double); /*!< Default Constructor with parameters */ 
  JarzynskiFreeEnergy(const JarzynskiFreeEnergy& J); /*< Copy constructor */
  ~JarzynskiFreeEnergy(); /*!< Destructor */ 
  void vecProcess();
  void resetIndex();
  void read(std::string input);  
  /*! Free energy calculation functions from work distribution */
  double JERaw(std::vector<double> *JEVector); /*< Raw JE interpreter */
  double JETaylor(std::vector<double> *JEVector); /*< Taylor Series (Second term) JE interpreter */  
  double JEalpha(std::vector<double> *JEVector); /*< Taylor Series (Second term) JE interpreter */  
  double JERawErr(std::vector<double> *JEVector); /*< Raw JE interpreter */
  double JETaylorErr(std::vector<double> *JEVector); /*< Taylor Series (Second term) JE interpreter */  
  double JEalphaErr(std::vector<double> *JEVector); /*< Taylor Series (Second term) JE interpreter */  

  double JEprocessVector(int, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector); /*< functor which can take a JE algorithm then processing it accordingly */
  double alpha(double, double, double); /*< */   

  // friend functions to take care of the MPI implementation
  // friend class MPI_setup (const JarzynskiFreeEnergy&);

  //friend class JEunitTest;
private:
  /*! Scientific constants used in this work */
  double BOLTZMANN = 0.0019872041; /*< units for the boltzmann constant are in kcal mol^-1 */
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
  std::vector<double> JEAlphaVector; /*< Storing the work ffor the taylor series JE interpreter */

  tupleList JERawCoordinateBin; /*< Vector for storing raw JE tuples */
  tupleList JETaylorCoordinateBin; /*< Vector for storing taylor series JE tuples */
  tupleList JEAlphaCoordinateBin; /*< Vector for storing taylor series JE tuples */
};

/*
class JarzynskiFreeEnergyTest : public CppUnit:TestCase, public JarzynskiFreeEnergy  { // Multiple inheritance - might be an issue 
public:
  // Might need some thinking regarding the constructors 
  void vectorTest();
  // Do I need any private variables for this class?
};
*/


// Reminder

/*
    Public mode: If we derive a sub class from a public base class. Then the public member of the base class will become public in the derived class and protected members of the base class will become protected in derived class.

    Protected mode: If we derive a sub class from a Protected base class. Then both public member and protected members of the base class will become protected in derived class.

    Private mode: If we derive a sub class from a Private base class. Then both public member and protected members of the base class will become Private in derived class.

// In short, usually, using public should be sufficient! 

 */

class MPI_setup : public JarzynskiFreeEnergy { // Make sure we inherit from the JarzynskiFreeEnergy 
public:
  MPI_setup();
  //  MPI_setup(int*, int*);
  MPI_setup(const MPI_setup& alloc);
  MPI_setup& operator=(const MPI_setup& alloc);
  virtual ~MPI_setup();
  
  void MPI_vec_send();
  void MPI_parameter_struct_constructor(MPI_Datatype*);
  void MPI_data_send(JarzynskiFreeEnergy*);
  void MPI_parameter_broadcast();  // Might not need to broadcasr 
  void MPI_divide_vector(std::vector<double>*);
  friend class JarzynskiFreeEnergy; // MPI cl ass inherits from the Jarzynski equality class 
private:
  int my_rank, comm_sz; // MPI address and total p size 
  parameterData parameters;
  std::vector<int> lineNumberVectorSplit; /*< Index vector */
  std::vector<double> coordinateZVectorSplit; /*< z coordinates vector */
  std::vector<double> bilayerCOMVectorSplit; /*< Bilayer COM vector */
  std::vector<double> forceVectorSplit; /*< Force vector */
  std::vector<double> workVectorSplit; /*< Work vector */   
  std::vector<double> JERawVectorSplit; /*< Storing the work for the raw JE interpreter */ 
  std::vector<double> JETaylorVectorSplit; /*< Storing the work ffor the taylor series JE interpreter */
  double low, high; // we get the highest and lowest values and divide
  MPI_Datatype VectorMPI, VectorMPI2;
  JarzynskiFreeEnergy sample; /*!< Instance of the class to get the function method */  

};

#endif 
