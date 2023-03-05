/*!
---------------------------------------------------------------------------------
| Jarzynski-Equality Algorithm based on multiple implementations/corrections    |
|                                                                               | 
| VERSION: 0.0.3                                                                |  
---------------------------------------------------------------------------------

Initialzation of the template defined in the header 

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

/*! cppunit tests */    

//#include <cppunit/extensions/TestFactoryRegistry.h>
//#include <cppunit/CompilerOutputter.h>
//#include <cppunit/ui/text/TestRunner.h>
//#include <cppunit/TestFixture.h>
//#include <cppunit/extensions/HelperMacros.h>

// Miscellaneous functions

static inline double computeSquare (double x) {return x*x;} // function for squaring the elements in a vector

void duplicate_remove(std::vector<double> &v) {
  auto end = v.end();
  for (auto it = v.begin(); it != end; ++it) {
    end = std::remove(it + 1, end, *it);
  }
  
  v.erase(end, v.end());
}

// ---
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

    tuple JERawVal{index, index - 0.5, index + 0.5, JEprocessVector(index, &JarzynskiFreeEnergy::JERaw, &JERawVector), JEprocessVector(index, &JarzynskiFreeEnergy::JERawErr, &JERawVector)}; //!< Make bins for storing the work values between i - 0.5 and i + 0.5 - i.e. a 1 angstrom interval, using the raw JE interpreter  
    tuple JETaylorVal{index, index - 0.5, index + 0.5, JEprocessVector(index, &JarzynskiFreeEnergy::JETaylor, &JETaylorVector), JEprocessVector(index, &JarzynskiFreeEnergy::JETaylorErr, &JETaylorVector)}; //!< Make bins for storing the work values between i - 0.5 and i + 0.5 - i.e. a 1 angstrom interval, using the taylor series JE interpreter 
    tuple JEAlphaVal{index, index - 0.5, index + 0.5, JEprocessVector(index, &JarzynskiFreeEnergy::JEalpha, &JEAlphaVector), JEprocessVector(index, &JarzynskiFreeEnergy::JEalphaErr, &JEAlphaVector)}; //!< Make bins for storing the work values between i - 0.5 and i + 0.5 - i.e. a 1 angstrom interval, using the taylor series JE interpreter 

    JERawCoordinateBin.push_back(JERawVal); //!< Store free energy values from JERaw algorithm  
    JETaylorCoordinateBin.push_back(JETaylorVal); //!< Push back values in each 
    JEAlphaCoordinateBin.push_back(JEAlphaVal); //!< Push back values in each 
  }

  std::cout << " RAW JE interpreter results" << std::endl;
  std::cout << " Columns: Coordinate, Coordinate range min, Coordinate range max, Free Energy (Kcal mol)" << std::endl;

  for (tupleList::const_iterator index = JERawCoordinateBin.begin(); index != JERawCoordinateBin.end(); ++index) {

    std::cout << index->get<0>() << " " << index->get<1>()  << " " << index->get<2>() << " " << std::fixed << std::setprecision(5) << index->get<3>() <<  " " << index->get<4>() << std::endl; //! Print out to 5 decimal places 

  }

  std::cout << " Taylor Series JE interpreter results" << std::endl;
  std::cout << " Columns: Coordinate, Coordinate range min, Coordinate range max, Free Energy (Kcal mol)" << std::endl;
  for (tupleList::const_iterator index = JETaylorCoordinateBin.begin(); index != JETaylorCoordinateBin.end(); ++index) {

    std::cout  << index->get<0>() << " " << index->get<1>()  << " " << index->get<2>() << " " << std::fixed << std::setprecision(5) << index->get<3>() << " " << index->get<4>() << std::endl; //! Print out to 5 decimal places 

  }

  std::cout << " Alpha Series JE interpreter results" << std::endl;
  std::cout << " Columns: Coordinate, Coordinate range min, Coordinate range max, Free Energy (Kcal mol)" << std::endl;
  for (tupleList::const_iterator index = JEAlphaCoordinateBin.begin(); index != JEAlphaCoordinateBin.end(); ++index) {
    std::cout  << index->get<0>() << " " << index->get<1>()  << " " << index->get<2>() << " " << std::fixed << std::setprecision(5) << index->get<3>() << " " <<  index->get<4>() << std::endl; //! Print out to 5 decimal places 
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
  std::vector<double> VarVector; /**< Vector to copy the value into, as to not change the values of the elements inside the vector pointer */

  FE result;
  double G; /*!< Free energy (Gibbs) */
  double Beta = 1 / (-BOLTZMANN * 303); /**< Boltzmann >Factor, at 303K */
  doubleIter workIterator; /**< iterator for vector */
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVector.push_back(exp(*workIterator * Beta));
  }

  for (int i = 0; i <= RawVector.size(); i++) {
    VarVector.push_back((RawVector[i]/RawVector.size()) * ( 1 / Beta));
  }

  double Varsum = std::accumulate(VarVector.begin(), VarVector.end(), 0.0);
  double Varmean = Varsum / VarVector.size();
  double Varsq_sum = std::inner_product(VarVector.begin(), VarVector.end(), VarVector.begin(), 0.0);
  double stdev = std::sqrt(Varsq_sum / VarVector.size() - Varmean * Varmean);
 
  G = log(std::accumulate(RawVector.begin(), RawVector.end(), 0.0) / RawVector.size()) * ( 1 / Beta); /*!< compute the raw JE */

  result.val = G;
  result.err = stdev;

  return G; 
}


double JarzynskiFreeEnergy::JERawErr(std::vector<double> *JEVector) {
  /** The Raw Jarzynski Equality computer */
  std::vector<double> RawVector; /**< Vector to copy the value into, as to not change the values of the elements inside the vector pointer */
  std::vector<double> VarVector; /**< Vector to copy the value into, as to not change the values of the elements inside the vector pointer */

  FE result;
  double G; /*!< Free energy (Gibbs) */
  double Beta = 1 / (-BOLTZMANN * 303); /**< Boltzmann >Factor, at 303K */
  doubleIter workIterator; /**< iterator for vector */
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVector.push_back(exp(*workIterator * Beta));
  }

  for (int i = 0; i <= RawVector.size(); i++) {
    VarVector.push_back((RawVector[i]/RawVector.size()) * ( 1 / Beta));
  }

  double Varsum = std::accumulate(VarVector.begin(), VarVector.end(), 0.0);
  double Varmean = Varsum / VarVector.size();
  double Varsq_sum = std::inner_product(VarVector.begin(), VarVector.end(), VarVector.begin(), 0.0);
  double stdev = std::sqrt(Varsq_sum / VarVector.size() - Varmean * Varmean);
 
  G = log(std::accumulate(RawVector.begin(), RawVector.end(), 0.0) / RawVector.size()) * ( 1 / Beta); /*!< compute the raw JE */

  return stdev; 
}


double JarzynskiFreeEnergy::JETaylor(std::vector<double> *JEVector) {
  //! The Taylor Series Jarzynski Equality computer 
  std::vector<double> RawVector; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 
  std::vector<double> squaredRawVector; /*!< Vector to store the squared work values. */  
  std::vector<double> VarVector; /**< Vector to copy the value into, as to not change the values of the elements inside the vector pointer */
  FE result;
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
  double squaredAvWork = std::accumulate(squaredRawVector.begin(), squaredRawVector.end(), 0.0)/ squaredRawVector.size(); 

  for (int i = 0; i <= RawVector.size(); i++) {
    VarVector.push_back((RawVector[i]));
  }

  double Varsum = std::accumulate(VarVector.begin(), VarVector.end(), 0.0);
  double Varmean = Varsum / VarVector.size();
  double Varsq_sum = std::inner_product(VarVector.begin(), VarVector.end(), VarVector.begin(), 0.0);
  double stdev = std::sqrt(Varsq_sum / VarVector.size() - Varmean * Varmean);
 
  G = AvWork - ((RawVector.size() / RawVector.size() - 1) * (Beta / 2) * ((squaredAvWork - (AvWork* AvWork)))); /*!< Taylor series interpreter */

  result.val = G;
  result.err = stdev;
  
  return G;
}


double JarzynskiFreeEnergy::JETaylorErr(std::vector<double> *JEVector) {
  //! The Taylor Series Jarzynski Equality computer 
  std::vector<double> RawVector; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 
  std::vector<double> squaredRawVector; /*!< Vector to store the squared work values. */  
  std::vector<double> VarVector; /**< Vector to copy the value into, as to not change the values of the elements inside the vector pointer */
  FE result;
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
  double squaredAvWork = std::accumulate(squaredRawVector.begin(), squaredRawVector.end(), 0.0)/ squaredRawVector.size(); 

  for (int i = 0; i <= RawVector.size(); i++) {
    VarVector.push_back((RawVector[i]));
  }

  double Varsum = std::accumulate(VarVector.begin(), VarVector.end(), 0.0);
  double Varmean = Varsum / VarVector.size();
  double Varsq_sum = std::inner_product(VarVector.begin(), VarVector.end(), VarVector.begin(), 0.0);
  double stdev = std::sqrt(Varsq_sum / VarVector.size() - Varmean * Varmean);
 
  G = AvWork - ((RawVector.size() / RawVector.size() - 1) * (Beta / 2) * ((squaredAvWork - (AvWork* AvWork)))); /*!< Taylor series interpreter */

  result.val = G;
  result.err = stdev;
  
  return stdev;
}


double JarzynskiFreeEnergy::JEalpha(std::vector<double> *JEVector) {

  double G; /*!< Free Energy */
  double Beta = 1 / (-BOLTZMANN * 303);  /*!< Boltzmann Factor */
  double Wdiss;
  double alpha;
  double B;
  
  std::vector<double> RawVectorJE; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 
  std::vector<double> RawVector; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 
  std::vector<double> VarVector; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 

  FE result;
  doubleIter workIterator; /*!< double iterator */

    /*! Store the raw work values from the JEVector */
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVectorJE.push_back(exp(*workIterator * Beta));
    RawVector.push_back(*workIterator);
  }
  
  double sum = std::accumulate(RawVector.begin(), RawVector.end(), 0.0);
  double mean = sum / RawVector.size();
  double sq_sum = std::inner_product(RawVector.begin(), RawVector.end(), RawVector.begin(), 0.0);

  Wdiss = 0.5 * Beta * std::sqrt(sq_sum / RawVector.size() - mean * mean);
  alpha = (log(15.0 * Beta * Wdiss)/log(15.0 * exp(2.0 * Beta * Wdiss) - 1.0));
  B = Wdiss/(pow(10.0 , alpha));

  for (int i = 0; i <= RawVector.size(); i++) {
    VarVector.push_back((RawVectorJE[i]) * (1/Beta) - B);
  }

  double Varsum = std::accumulate(VarVector.begin(), VarVector.end(), 0.0);
  double Varmean = Varsum / VarVector.size();
  double Varsq_sum = std::inner_product(VarVector.begin(), VarVector.end(), VarVector.begin(), 0.0);
  double stdev = std::sqrt(Varsq_sum / VarVector.size() - Varmean * Varmean);
  
  G = log(std::accumulate(RawVectorJE.begin(), RawVectorJE.end(), 0.0) / RawVectorJE.size()) * ( 1 / Beta) - B;
  result.val = G;
  result.err = stdev;
  return G;
}


double JarzynskiFreeEnergy::JEalphaErr(std::vector<double> *JEVector) {
  double G; /*!< Free Energy */
  double Beta = 1 / (-BOLTZMANN * 303);  /*!< Boltzmann Factor */
  double Wdiss;
  double alpha;
  double B;
  
  std::vector<double> RawVectorJE; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 
  std::vector<double> RawVector; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 
  std::vector<double> VarVector; /*!< Vector to copy the work value into, as to not change the values of the elements inside the vector pointer */ 

  FE result;
  doubleIter workIterator; /*!< double iterator */

    /*! Store the raw work values from the JEVector */
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVectorJE.push_back(exp(*workIterator * Beta));
    RawVector.push_back(*workIterator);
  }
  
  double sum = std::accumulate(RawVector.begin(), RawVector.end(), 0.0);
  double mean = sum / RawVector.size();
  double sq_sum = std::inner_product(RawVector.begin(), RawVector.end(), RawVector.begin(), 0.0);

  Wdiss = 0.5 * Beta * std::sqrt(sq_sum / RawVector.size() - mean * mean);
  alpha = (log(15.0 * Beta * Wdiss)/log(15.0 * exp(2.0 * Beta * Wdiss) - 1.0));
  B = Wdiss/(pow(10.0 , alpha));

  for (int i = 0; i <= RawVector.size(); i++) {
    VarVector.push_back((RawVectorJE[i]) * (1/Beta) - B);
  }

  double Varsum = std::accumulate(VarVector.begin(), VarVector.end(), 0.0);
  double Varmean = Varsum / VarVector.size();
  double Varsq_sum = std::inner_product(VarVector.begin(), VarVector.end(), VarVector.begin(), 0.0);
  double stdev = std::sqrt(Varsq_sum / VarVector.size() - Varmean * Varmean);
  
  G = log(std::accumulate(RawVectorJE.begin(), RawVectorJE.end(), 0.0) / RawVectorJE.size()) * ( 1 / Beta) - B;

  result.val = G;
  result.err = stdev;

  return stdev;
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

/*!<  MPI class - Sending the datatypes and vectors  */

//MPI_setup::MPI_setup() { // Default constructor for MPI
//  MPI_Init(NULL, NULL);
//  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank); // Allocate current rank to my_rank
//  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);  // allocate total number of nodes to p 
//
//  // Set up basic scientific constants
//  parameters.BM = 0.0019872041; // units for the boltzmann constant are in kcal mol^-1 ; 
//  parameters.T = 303;
//
//}
//
//
//void MPI_setup::MPI_parameter_struct_constructor(MPI_Datatype* input_mpi_t_p) {
//  parameters.BM = 0.0019872041; // units for the boltzmann constant are in kcal mol^-1 ; 
//  parameters.T = 303;
//  // Define parameters for storing the variables 
//
//  int array_of_blocklengths[2] = {1,1};
//  MPI_Datatype array_of_types[2] = {MPI_DOUBLE, MPI_DOUBLE};
//  MPI_Aint array_of_displacements[2] = {0};
//  MPI_Aint BM_addr, T_addr;	     
//
//  MPI_Get_address(&parameters.BM, &BM_addr);
//  MPI_Get_address(&parameters.T, &T_addr);
//
//  array_of_displacements[1] = T_addr - BM_addr;
//  MPI_Type_create_struct(2, array_of_blocklengths, array_of_displacements, array_of_types, input_mpi_t_p);
//  MPI_Type_commit(input_mpi_t_p);
//}
//
//void MPI_setup::MPI_data_bcast(JarzynskiFreeEnergy* serialClass) {
//  workVectorSplit.assign(serialClass->workVector.begin(), serialClass->workVector.end()); // Assign all workvector values read from the serial class to here 
//  coordinateZVectorSplit.assign(serialClass->coordinateZVector.begin(), serialClass->coordinateZVector.end()); // Assign all coordinate ZVector to here   
//  double maxZ = *max_element(serialClass->coordinateZVector.begin(), serialClass->coordinateZVector.end()); //!< Define minimum z coordinate                                
//  double minZ = *min_element(serialClass->coordinateZVector.begin(), serialClass->coordinateZVector.end()); //!< Define maximum z coordinate          
//  MPI_Bcast(&workVectorSplit[0], workVectorSplit.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//  MPI_Bcast(&coordinateZVectorSplit[0], coordinateZVectorSplit.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
//}
//
//
//// Using a function ppinter
//
//void MPI_setup::MPI_divide_vector(int position, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector) {
//
//  int work_index = 0;  
//  JarzynskiFreeEnergy sample; /*!< Instance of the class to get the function method */  
//  doubleIter diterator; /*!< Integer iterator */ 
//  std::vector<double> coordinateVector;
//  std::vector<double> positionVector;
//
//  // Largly copied from the serial code - May need to implement this through inheritance 
//  
//  // 1. Divide the workvector into managable chunks and then redsitribute each work bin into multiple nodes
//  
//  for (diterator = coordinateZVector.begin(); diterator <= coordinateZVector.end(); ++diterator, ++work_index) { 
//    if (*diterator > position - 0.5 && *diterator < position + 0.5) { /*!< If the work values are within an angstrom range, add to vector */  
//      JEVector->push_back((workVector[work_index])); /*!< store work values */ 
//      // Not pointer 
//      coordinateVector.push_back(*diterator);
//      positionVector.push_back(position);
//    }
//  }
//
//  // Remove duplicate values in positionVector
//
//  duplicate_remove(positionVector);
//  
//  // Get rough estimate of the 
//
//  std::cout << "Each node will receive a MPI_Type vector of element size: " << JEVector << std::sendl; 
//  double* vecPointer = JEVector.data(); // Make it slighly easier to allocate data onto the MPI_dervied datatype later on   
//  // Builidng a derived type, based on the first vector size  
//  MPI_Type_vector(JEVector.size(), 1, JEVector.size(), VectorMPI, &VectorMPI2); // As it is ex
//  MPI_Type_commit(&VectorMPI2);
//
//  std::cout << "The total number of nodes allocated for this work is .. " << std::endl;
//    
//  // TODO - need to make sure that 
//
//  if (my_rank == 0) {
//    MPI_send(&JEVector[0], 1, VectorMPI2, 1, 0, MPI_COMM_WORLD);  
//  } else {
//
//    for (int source = 1; source < comm_sz; source++) {
//      MPI_recv(&JEVector[0], 1, VectorMPI2, 0, 0, MPI_COMM_WORLD, &status); // If rank is not 0, get divide the values evenly
//    }
//  }
//
//											   
//  
//
//}
