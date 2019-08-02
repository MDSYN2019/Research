
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

bool compare(const Core& c1, const Core& c2) {
  return c1.name() < c2.name();
}

/*

To complete our classes, we need to implement four constructors. Default
constructor and the constructor that takes an istream, once for each class.

 */

bool compare_grades(const Core& c1, const Core& c2) {
  return c1.grade() < c2.grade();
}



// Filling in the details for the Core method

Core::Core() {
  midterm(0), final(0);
} // default constructor

Core::Core(std::istream& is) {
  read(is);
}

std::string Core::name() const {
  return n;
}

double Core::grade() const {
  return midterm; // 
}

std::istream& Core::read_common(std::istream& in) {
  // read and store the student's name and exam grades
  in >> n >> midterm >> final;
  return in;
}

virtual std::istream& Core::read(std::istream& in) {
  read_common(in);
  read_hw(in, homework);
  return in;
}


// Grad method

double Grad::grade() const {
  return std::min(Core::grade(), thesis);
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


