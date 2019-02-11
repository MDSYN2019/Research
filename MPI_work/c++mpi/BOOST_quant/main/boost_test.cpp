/** 
 * 
 *
 *
 **/

#include <stdio.h>

// Standard C++ libraries - including STL libraries

#include <map>
#include <string>
#include <algorithm>
#include <iostream>
#include <vector>
#include <istream>

// Boost libraries

#include <boost/current_function.hpp>
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/detail/lightweight_test.hpp>

// GSL libraries

#include <gsl/gsl_math.h>
#include <gsl/gsl_sf_bessel.h>

using std::map;  
using std::cin;  
using std::endl; 
using std::cout;
using std::string;
using std::vector;

class boost_example {   
public:
  void testMacro(); 
}; // Closes boost_example

void boost_example::testMacro() {

  vector<double> myVec(10);
  BOOST_FOREACH(double& x, myVec) x = 10.0;
  BOOST_FOREACH(double& x, myVec)  std::cout << x << std::endl;  

};



int main (void) {

  
  double x = 5.0;
  double y = gsl_sf_bessel_J0(x);
  cout << "J0(%g) = %.18e" <<  x <<  y;
  
  // Counting words, using associative containers   
  string s;
  map<string, int> counters; // store each word and an associated counter
  while (cin >> s) {
    ++counters[s];
  }
  for (map<string, int>::const_iterator it = counters.begin(); it != counters.end(); ++it) {
    //cout << it->first << "\t" << it->second << endl;
  }

  // using boost_foreach for looping over the map container 
  boost_example a;
  a.testMacro();
  
  // Find all the lines that refer to each word in the input
  
  return 0;
}
