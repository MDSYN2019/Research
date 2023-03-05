/*!
---------------------------------------------------------------------------------
| Crooks Fluctuation Algorithm based on multiple implementations/corrections    |
|                                                                               | 
| VERSION: 0.0.1                                                                |  
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

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


class CrooksFluc:
  
