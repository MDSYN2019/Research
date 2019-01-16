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


The idea of the steered molecular dynamics / Jarzynski equality is to sample the work values and to assess the free energy change. As work is path-dependent 
and not a state function, whilst free energy is a equilibrium state function, the main novelty of the equality is that an equilirium property can be extracted
from a non-equilbrium experiment.

Here, the experiment was carried out using LAMMPS (https://lammps.sandia.gov/). The simulation of a toluene going through a lipid bilayer was carried out at various velocities, 
ones that are far out of the `equilibrium' regime (i.e. very fast pulling simulations), and those closer to the `equilbrium' regime (slower simulations).  Each output file
prints out the index, reaction coordinate, force and work values respectively.  

Here, there are two main free energy interpreters - the raw Jarzynski equality as seen from reference 1, and the taylor series improved evaluator seen in reference 3. 
 
This code has the serial code for estimating the Jarzynski estimate. I have used primarily C++ stl-based libraries, with some Boost code as well.

Hence, this code will require the boost libraries and in the future will be using the GSL libraries for managing the statistical aspects of the algorithm.

Abbreviation: JE - Jarzynski Equalily

*/

/*

---------------------------------  
|           Headers             |
--------------------------------- 
 
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

// STL libraries

#include <string>
#include <vector>
#include <list>
#include <array>
#include <map>
#include <algorithm>
#include <experimental/filesystem>

/* Boost libraries */

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>

/* GSL libraries */

#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram2d.h>
#include <gsl/gsl_statistics_int.h> 

/* MPI libraries */

#include <mpi.h>

/* STL iterators -  typedefs for iterators, and typedef for vectors and tuples */

typedef boost::tuple<int, double, double, double> tuple; 
typedef std::vector<boost::tuple<int, double, double, double> > tupleList; 
typedef std::vector<int>::iterator intIter; // Integer iterator 
typedef std::vector<double>::iterator doubleIter; // Double iterator
typedef std::vector<std::string>::const_iterator stringIter; // String iterator

static inline double computeSquare (double x) { return x*x; } // function for squaring the elements in a vector

void get_data(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank, int p) {
  //TODO
}

void Check_for_error(
      int       local_ok   /* in */,
      char      fname[]    /* in */,
      char      message[]  /* in */,
      MPI_Comm  comm       /* in */) {
   int ok;

   MPI_Allreduce(&local_ok, &ok, 1, MPI_INT, MPI_MIN, comm);
   if (ok == 0) {
      int my_rank;
      MPI_Comm_rank(comm, &my_rank);
      if (my_rank == 0) {
         fprintf(stderr, "Proc %d > In %s, %s\n", my_rank, fname,
               message);
         fflush(stderr);
      }
      MPI_Finalize();
      exit(-1);
   }
}  

class MPIJarzynskiFreeEnergy {
public:
  MPIJarzynskiFreeEnergy(); // default constructor 
  MPIJarzynskiFreeEnergy(int, int); // default constructor 
  ~MPIJarzynskiFreeEnergy(); // Destructor 
  double MPIJEprocessVector(int, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector); // Acquire work
  void MPIJE_Bcast(int, int, std::vector<double>); // TODO 

private:
  int comm_sz;
  int my_rank; // rank of process
  int p; // number of processes
  int source;  // rank of sender
  int dest;
  int tag = 0;
  std::string S;
  MPI_Status status;
  tupleList JERawCoordinateBin; // vector for storing raw JE tuples 
};

MPI_JarzynskiFreeEnergy::MPI_JarzynskiFreeEnergy() {
  // default constructor
  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
}
void MPI_JarzynskiFreeEnergy::MPIJE_Bcast () {
}
double MPI_JarzynskiFreeEnergy::MPI_JEprocessVector(int position, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector) { 

  int work_index = 0;
  doubleIter diterator; // Integer iterator
  JarzynskiFreeEnergy sample; // Instance of the class to get the function method
  
  JEVector->clear(); // Empty vector before doing mean calculations - as we are working with a pointer to a vector, and not creating a copy each time, this is required

  for (diterator = coordinateZVector.begin(); diterator <= coordinateZVector.end(); ++diterator, ++work_index) { // loop over the work values, count up the work_index 
    if (*diterator > position - 0.5 && *diterator < position + 0.5) { // If the work values are within an angstrom range .. 
      JEVector->push_back((workVector[work_index])); // store work values
    }
  }  
  return JEVector; // Assess the free energy with the function as pointed to by sample.f* 
}


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
 
  // ifstream

  std::ifstream myfile;

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

JarzynskiFreeEnergy::JarzynskiFreeEnergy() {
} // Default Constructor

JarzynskiFreeEnergy::~JarzynskiFreeEnergy() {
} // Default Destructor

double JarzynskiFreeEnergy::JEprocessVector(int position, double (JarzynskiFreeEnergy::*f) (std::vector<double> *VectorInput), std::vector<double> *JEVector) { 

  int work_index = 0;
  doubleIter diterator; // Integer iterator
  JarzynskiFreeEnergy sample; // Instance of the class to get the function method
  
  JEVector->clear(); // Empty vector before doing mean calculations - as we are working with a pointer to a vector, and not creating a copy each time, this is required

  for (diterator = coordinateZVector.begin(); diterator <= coordinateZVector.end(); ++diterator, ++work_index) { // loop over the work values, count up the work_index 
    if (*diterator > position - 0.5 && *diterator < position + 0.5) { // If the work values are within an angstrom range .. 
      JEVector->push_back((workVector[work_index])); // store work values
    }
  }  
  return (sample.*f)(JEVector); // Assess the free energy with the function as pointed to by sample.f* 
}

void JarzynskiFreeEnergy::vecProcess() {

  /*
  
    This function uses the boost::tuple and create a vector of tuples; it stores the free energy values with the first index storing the coordinate value where 
    we base the search. The second and third elements store the minimum and maximum range from which the work distribution was analyzed. Finally, the last element 
    stores the free energy values.
   
  */
    
  double max_z = *max_element(coordinateZVector.begin(), coordinateZVector.end()); // Define minimum z coordinate                                    
  double min_z = *min_element(coordinateZVector.begin(), coordinateZVector.end()); // Define maximum z coordinate                                                                    

  /* We want to accumulate the values via the bins: */                                                                                                     
  for (int i = 0; i < max_z; ++i) {
    /* Create a new tuple for each bin. */
    tuple JERawVal{i, i - 0.5, i + 0.5, JEprocessVector(i, &JarzynskiFreeEnergy::JERaw, &JERawVector)}; // Make bins for storing the work values between i - 0.5 and i + 0.5 - i.e. a 1 angstrom interval, using the raw JE interpreter
    tuple JETaylorVal{i, i - 0.5, i + 0.5, JEprocessVector(i, &JarzynskiFreeEnergy::JETaylor, &JETaylorVector)}; // Make bins for storing the work values between i - 0.5 and i + 0.5 - i.e. a 1 angstrom interval, using the taylor series JE interpreter

    JERawCoordinateBin.push_back(JERawVal); // Store 
    JETaylorCoordinateBin.push_back(JETaylorVal); // Push back values in each
  }

  std::cout << " RAW JE interpreter results" << std::endl;
  std::cout << " Columns: Coordinate, Coordinate range min, Coordinate range max, Free Energy (Kcal mol)" << std::endl;
  
  for (tupleList::const_iterator i = JERawCoordinateBin.begin(); i != JERawCoordinateBin.end(); ++i) {
    std::cout << "Bin Center: " << i->get<0>() << " " << i->get<1>()  << " " << i->get<2>() << " " << std::fixed << std::setprecision(5) << i->get<3>() << std::endl; // Print out to 5 decimal places 
  }

  std::cout << " Taylor Series JE interpreter results" << std::endl;
  std::cout << " Columns: Coordinate, Coordinate range min, Coordinate range max, Free Energy (Kcal mol)" << std::endl;
  
  for (tupleList::const_iterator i = JETaylorCoordinateBin.begin(); i != JETaylorCoordinateBin.end(); ++i) {
    std::cout << "Bin Center: " << i->get<0>() << " " << i->get<1>()  << " " << i->get<2>() << " " << std::fixed << std::setprecision(5) << i->get<3>() << std::endl; // Print out to 5 decimal places
  }
}

void JarzynskiFreeEnergy::resetIndex(){

  /* Clear all bins before reading in anything */
  
  lineNumberVector.clear(); // file line index
  coordinateZVector.clear(); // z coordinates
  bilayerCOMVector.clear();  // COM of the bilayer
  forceVector.clear(); // force values
  workVector.clear(); // work values
}

double JarzynskiFreeEnergy::JERaw(std::vector<double> *JEVector) {
  
  std::vector<double> RawVector; // Vector to copy the value into, as to not change the values of the elements inside the vector pointer 
  double G; // Free energy (Gibbs)

  double Beta = 1.0 / (BOLTZMANN * 303); // Boltzmann Factor
  doubleIter workIterator; // iterator 

  // Storing the exp(-Beta * work) in the new vector
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVector.push_back(exp( (*workIterator) * -Beta ));
  }

  G = log (std::accumulate(RawVector.begin(), RawVector.end(), 0.0) / RawVector.size() ) / -Beta; // compute the raw JE 
  return G; 
}

double JarzynskiFreeEnergy::JETaylor(std::vector<double> *JEVector) {

  std::vector<double> RawVector; // Vector to copy the work value into, as to not change the values of the elements inside the vector pointer 
  std::vector<double> squaredRawVector; // Vector to store the squared work values.  
  double G;

  double Beta = 1 / (-BOLTZMANN * 303);  // Boltzmann Factor
  doubleIter workIterator; // double iterator

  // Store the raw work values from the JEVector  
  for (workIterator = JEVector->begin(); workIterator != JEVector->end(); ++workIterator) {
    RawVector.push_back(*workIterator);
  }
  
  // Store the squared work values
  squaredRawVector.resize(RawVector.size());
  std::transform(RawVector.begin(), RawVector.end(), squaredRawVector.begin(), computeSquare);
  
  double AvWork = std::accumulate(RawVector.begin(), RawVector.end(), 0.0)/ RawVector.size(); // Average work   
  double squaredAvWork = std::accumulate(squaredRawVector.begin(), squaredRawVector.end(), 0.0)/ squaredRawVector.size(); // Average squared work

  G = AvWork - (RawVector.size() / RawVector.size() -1) * (Beta / 2) * (squaredAvWork - AvWork* AvWork); // taylor series interpreter
  return G;
}

void JarzynskiFreeEnergy::read(std::string input) {

  /*
    IO - Input the column data values as computed from LAMMPS 
  */
  
  myfile.open(input, std::ifstream::in);
  std::cout << myfile.is_open() << std::endl;
  
  if (myfile.is_open() == 0) {    
    std::cout << "Cannot open files" << std::endl;     
    exit(1);
  }
  else {	
    std::cout << "Successfully read. Opening file:" << " " << input << std::endl;	
  }
    while (myfile >> number >> z >> bilayerCOM >> force >> work) {

    lineNumberVector.push_back(number); // Line index
    coordinateZVector.push_back(z); // z coordinates of the nanoparticle
    bilayerCOMVector.push_back(bilayerCOM); // z coordinate of the bilayer COM
    forceVector.push_back(force); // Force values
    workVector.push_back(work); // Work values
    nLines++; // Add index for next line 

    }
    myfile.close();
    std::cout << "The number of lines in file:" << " " << nLines << " " << std::endl;
    //  return nlines;
}

int main(int argc, char *argv[]) {  
  
  std::string filenameString = argv[1]; // Pull file names
  std::string pullDirectory = argv[2]; // Directory for where the pull files are in the computer
  int numberOfFiles = atoi(argv[3]); // Number of pull files
  
  JarzynskiFreeEnergy FE; // Class instance  
  int index = 0;
  
  std::string filename;
  std::string path;
  std::vector<std::string> traj;

  FE.resetIndex();  
  for (int i = 0; i < numberOfFiles; i++) {     
    index = i + 1; // Filename index
    auto b = std::to_string(index);
    filename = filenameString + "." + b; // Concatenate string  
    traj.push_back(filename); // store the filesname    
  }

  // Loop through all the filesnames and work through them 
  for (std::vector<std::string>::const_iterator i = traj.begin(); i != traj.end(); ++i) {
    path = pullDirectory + "/" + std::string(*i);
    std::cout << "Reading pull files from directory: " << path << std::endl;
    FE.read(path);
  } 
  
  FE.vecProcess(); // Compute the JE interpreters
  return 0;
}

