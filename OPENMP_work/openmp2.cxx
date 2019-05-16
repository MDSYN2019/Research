#include <iostream>
#include <cstdlib>
#include <cstdio>
//#include <Eigen/Dense>

/* Instead of just calling the OpenMP functions, e can first check whetehr _OPENMP is defined. */

#ifndef _OPENMP
#include <omp.h>
#endif

#include "openmp1.h"

/* cppunit tests */                                                                                                                                      
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


/*

A numeric progression is a seuqence of numbers, where the value of each number depends on one or more of the 
previous value.
 
*/



class Progression {

};
