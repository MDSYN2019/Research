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

/*! STL libraries */

#include <string>
#include <vector>
#include <list>
#include <array>
#include <map>
#include <algorithm>
#include <experimental/filesystem>

/*! Boost libraries */

#include "boost/tuple/tuple.hpp"
#include "boost/tuple/tuple_io.hpp"

/*! GSL libraries */

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_statistics_int.h> 

/*! Import from header */

#include "JE_compute.hpp"

/*! MPI */

#include "mpi.h"

static inline double computeSquare (double x) { return x*x;} // function for squaring the elements in a vector

JarzynskiFreeEnergy::JarzynskiFreeEnergy() {} 
JarzynskiFreeEnergy::~JarzynskiFreeEnergy() {}

double JarzynskiFreeEnergy::JEprocessVector(int position, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector) { 
  /*!< Compute free energies with algorithm plugged in from JarzynskiFreeEnergy::*f, with the work values from std::vector<double> *VectorInput  */ 
  int work_index = 0;
  doubleIter diterator; /*!< Integer iterator */ 
  JarzynskiFreeEnergy sample; /*!< Instance of the class to get the function method */  
  JEVector->clear(); /*!< Empty vector before doing mean calculations - as we are working with a pointer to a vector, and not creating a copy each time, this is required */
  for (diterator = coordinateZVector.begin(); diterator <= coordinateZVector.end(); ++diterator, ++work_index) { 
    if (*diterator > position - 0.5 && *diterator < position + 0.5) { /*!< If the work values are within an angstrom range, add to vector */  
      JEVector->push_back((workVector[work_index])); /*!< store work values */ 
    }
  }
  return (sample.*f)(JEVector); /*!< Assess the free energy with the function as pointed to by sample.f* */ 
}

void JarzynskiFreeEnergy::vecProcess() {
  /*! This function uses the boost::tuple and create a vector of tuples; it stores the free energy values with the first index storing the coordinate value where we base the search. The second and third elements store the minimum and maximum range from which the work distribution was analyzed. Finally, the last element stores the free energy values. 
*/
  double maxZ = *max_element(coordinateZVector.begin(), coordinateZVector.end()); //!< Define minimum z coordinate                                
  double minZ = *min_element(coordinateZVector.begin(), coordinateZVector.end()); //!< Define maximum z coordinate          
  //! We want to accumulate the values via the bins: */                                                                                                     
  for (int index = 0; index < maxZ; ++index) {

    tuple JERawVal{index, index - 0.5, index + 0.5, JEprocessVector(index, &JarzynskiFreeEnergy::JERaw, &JERawVector)}; //!< Make bins for storing the work values between i - 0.5 and i + 0.5 - i.e. a 1 angstrom interval, using the raw JE interpreter  
    tuple JETaylorVal{index, index - 0.5, index + 0.5, JEprocessVector(index, &JarzynskiFreeEnergy::JETaylor, &JETaylorVector)}; //!< Make bins for storing the work values between i - 0.5 and i + 0.5 - i.e. a 1 angstrom interval, using the taylor series JE interpreter 
    
    JERawCoordinateBin.push_back(JERawVal); //!< Store free energy values from JERaw algorithm  
    JETaylorCoordinateBin.push_back(JETaylorVal); //!< Push back values in each 
  }

  std::cout << " RAW JE interpreter results" << std::endl;
  std::cout << " Columns: Coordinate, Coordinate range min, Coordinate range max, Free Energy (Kcal mol)" << std::endl;

  for (tupleList::const_iterator index = JERawCoordinateBin.begin(); index != JERawCoordinateBin.end(); ++index) {
    std::cout << "Bin Center: " << index->get<0>() << " " << index->get<1>()  << " " << index->get<2>() << " " << std::fixed << std::setprecision(5) << index->get<3>() << std::endl; //! Print out to 5 decimal places 
  }

  std::cout << " Taylor Series JE interpreter results" << std::endl;
  std::cout << " Columns: Coordinate, Coordinate range min, Coordinate range max, Free Energy (Kcal mol)" << std::endl;

  for (tupleList::const_iterator index = JETaylorCoordinateBin.begin(); index != JETaylorCoordinateBin.end(); ++index) {
    std::cout << "Bin Center: " << index->get<0>() << " " << index->get<1>()  << " " << index->get<2>() << " " << std::fixed << std::setprecision(5) << index->get<3>() << std::endl; //! Print out to 5 decimal places 
  }
}

void JarzynskiFreeEnergy::resetIndex(){
  /** Clear all bins before reading in anything */
  lineNumberVector.clear(); /**< Clear file line index */
  coordinateZVector.clear(); /**< Clear z coordinates */
  bilayerCOMVector.clear();  /**< Clear COM of the bilayer */
  forceVector.clear(); /**< Clear Force values */
  workVector.clear(); /**< Work values */
}

double JarzynskiFreeEnergy::JERaw(std::vector<double> *JEVector) {
  /** The Raw Jarzynski Equality computer */
  std::vector<double> RawVector; /**< Vector to copy the value into, as to not change the values of the elements inside the vector pointer */
  double G; /*!< Free energy (Gibbs) */
  double Beta = 1 / (-BOLTZMANN * 303); /**< Boltzmann >Factor, at 303K */
  doubleIter workIterator; /**< iterator for vector */
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVector.push_back(exp(*workIterator * Beta));
  }
  
  G = log(std::accumulate(RawVector.begin(), RawVector.end(), 0.0)) * Beta; /*!< compute the raw JE */
  return G; 
}

double JarzynskiFreeEnergy::JETaylor(std::vector<double> *JEVector) {
  //! The Taylor Series Jarzynski Equality computer 
  std::vector<double> RawVector; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 
  std::vector<double> squaredRawVector; /*!< Vector to store the squared work values. */  
  double G; /*!< Free Energy */
  double Beta = 1 / (-BOLTZMANN * 303);  /*!< Boltzmann Factor */
  doubleIter workIterator; /*!< double iterator */

  /*! Store the raw work values from the JEVector */
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVector.push_back(*workIterator);
  }  
  /*! Store the squared work values */

  squaredRawVector.resize(RawVector.size());
  std::transform(RawVector.begin(), RawVector.end(), squaredRawVector.begin(), computeSquare);
  double AvWork = std::accumulate(RawVector.begin(), RawVector.end(), 0.0)/ RawVector.size(); /*!< Average work */   
  double squaredAvWork = std::accumulate(squaredRawVector.begin(), squaredRawVector.end(), 0.0)/ squaredRawVector.size(); /*!< Average squared work */
  G = AvWork - (RawVector.size() / RawVector.size() -1) * (Beta / 2) * (squaredAvWork - AvWork* AvWork); /*!< Taylor series interpreter */
  return G;
}

void JarzynskiFreeEnergy::read(std::string input) {
  /*! IO - Input the column data values as computed from LAMMPS */
  std::ifstream myFile;
  myFile.open(input, std::ifstream::in);
  std::cout << myFile.is_open() << std::endl;
  if (myFile.is_open() == 0) {    
    std::cout << "Cannot open files" << std::endl;     
    exit(1);
  } else {	
    std::cout << "Successfully read. Opening file:" << " " << input << std::endl;	
    std::cout << std::endl;
  }
  
  while (myFile >> number >> z >> bilayerCOM >> force >> work) {
    lineNumberVector.push_back(number); /*!< Line index */
    coordinateZVector.push_back(z); /*!< z coordinates of the Nanoparticle/molecule */
    bilayerCOMVector.push_back(bilayerCOM); /*!< z coordinate of the bilayer COM */ 
    forceVector.push_back(force); /*!< Force values */
    workVector.push_back(work); /*!< Work values */
    nLines++; /*!< Add index for next line */ 
  }

  myFile.close();
  std::cout << "The number of lines in file:" << " " << nLines << " " << std::endl;
  std::cout << std::endl;
}

// MPI class 
MPI_setup::MPI_setup() { // Default constructor for MPI
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
}

void MPI_setup::MPI_parameter_struct_constructor(MPI_Datatype* input_mpi_t_p) {
  parameterData parameters;
  parameters.BM = 0.0019872041; /*< units for the boltzmann constant are in kcal mol^-1 */; 
  parameters.T = 303;

  // Define parameters for storing the variables 

  int array_of_blocklengths[2] = {1,1};
  MPI_Datatype array_of_types[2] = {MPI_DOUBLE, MPI_DOUBLE};
  MPI_Aint array_of_displacements[2] = {0};
  MPI_Aint BM_addr, T_addr;	     

  MPI_Get_address(&parameters.BM, &BM_addr);
  MPI_Get_address(&parameters.T, &T_addr);

  array_of_displacements[1] = T_addr - BM_addr;
  MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements, array_of_types, input_mpi_t_p);
  MPI_Type_commit(input_mpi_t_p);
}

/*
void MPI_setup::MPI_parameter_broadcast() {
}
*/
void MPI_setup::MPI_data_send(JarzynskiFreeEnergy* serialClass) {
  // double max_z = *max_element(serialClass->coordinateZVector.begin(), serialClass->coordinateZVector.end()); //!< Define minimum z coordinate     
  //  double min_z = *min_element(serialClass->coordinateZVector.begin(), serialClass->coordinateZVector.end()); //!< Define maximum z coordinate
  // Make a copy of the coordinate and workvector vector, so that we do not touch the original vector
  workVectorSplit.assign(serialClass->workVector.begin(), serialClass->workVector.end());
  coordinateZVectorSplit.assign(serialClass->coordinateZVector.begin(), serialClass->coordinateZVector.end());  
  //  int subvec_size = workVecctorSplit.size() / p;
  // std::cout << "For this MPI code, we are splitting the vector into " << subvec_size << "chunks"  << std::endl;
  if (my_rank == 0) {
    // This is a bit tricky - I will need to make sure that the vector is divided equally
  } 
}
