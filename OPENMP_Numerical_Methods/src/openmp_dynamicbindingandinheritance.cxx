
#include <iostream>
#include <iomanip>
#include <istream>
#include <vector>
#include <map>
#include <fstream>


#include "openmp_dynamicbindingandinheritance.hpp"


/*

This code contains two more important differences. First note that 
the call to ::grade putting :: in front of a name insists on using 
a version of that name that is not a member of anything.

 */


/*

To complete our classes, we need to implement four constructors. Default
constructor and the constructor that takes an istream, once for each class.

 */


bool compare(const Core& c1, const Core& c2) {
  return c1.name() < c2.name();
}

/*

Our compare function does the right thing when we call it 
with a Grad object as an argument because the name function 
is shared by both Grad and Cor objects

As a first cut at solving this problem, we'd write 
a function that is similar to compare

*/

bool compare_grades(const Core& c1, const Core& c2) {
  return c1.grade() < c2.grade();
}

std::string Core::name() const {
  return n;
}

//double Core::grade() const {
//  return ::grade(midterm, final, homework);
//}

std::istream& Core::read_common(std::istream& in) {
  // read and store the student's name and exam grades
  in >> n >> midterm >> final;
  return in;
}

std::istream& Core::read(std::istream& in) {
  read_common(in);
  read_hw(in, homework);
  return in;
}

// inherited read method 
std::istream& Grad::read(std::istream& in) {
  /*
    Note that in the definition of Grad::read, we can refer to 
    elements from the base class without any special notation, 
    because these elements are also members of Grad. 
  */

  read_common(in);
  in >> thesis;
  read_hw(in, homework);
  return in; 
}


