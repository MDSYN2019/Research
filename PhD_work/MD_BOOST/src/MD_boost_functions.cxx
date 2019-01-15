// Only including the boost libraries
// GSL libraries
#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

// Standard C++ libraries
#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <istream>
#include "MD_boostlib.hpp"

class boost_example {  

public:
  void testMacro();
};

/*

One of the nice features of C++ is that you can give special meanings to operators, when they are used with user-defined classes. This is called operator overloading

 */


// Implementation of the Gauss-Jordan elimination
void gaussJ(MatDoub_IO &a, MatDoub_IO)
