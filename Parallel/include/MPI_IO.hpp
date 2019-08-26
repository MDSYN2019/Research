/*


 */

#ifndef __MPI_IO__
#define __MPI_IO__

#include <iostream>
#include <iomanip>
#include <vector>
#include <mpi.h>

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>



class MPIInput : public CppUnit::TestCase { // Inherit from Cppunittest

public:
  //  MPIInput(); /*!< Default constructor - at the moment disabled*/  
  MPIInput(int, int); /*!< int int constructor */
  void getData();
  void bubbleSort();
  void oddEvenSort();
  virtual ~MPIInput(); 
  
private:
  int source = 0; /*!< process sending integral */ 
  int dest = 0; /*!< All messages go to 0 */
  int tag = 0; /*!< The tag int is used for marking MPI messages */
  int my_rank;
  int p;
  
  int* start;
  int* end;
  int* n_ptr;

  MPI_Status status;

};

#endif 
