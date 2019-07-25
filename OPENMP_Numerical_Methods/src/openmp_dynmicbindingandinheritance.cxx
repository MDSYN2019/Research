
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include "openmp_dynamicbindingandinheritance.hpp"


/*

This code contains two more important differences. First note that 
the call to ::grade putting :: in front of a name insists on using 
a version of that name that is not a member of anything.


 */


std::string Core::name() const {
  return n;
}

double Core::grade() const {
  return ::grade(midterm, final, homework);
}

istream& Core::read_common(istream& in) {
  // read and store the student's name and exam grades
  in >> n >> midterm >> final;
  return in;
}

istream& Core::read(istream& in) {
  read_common(in);
  read_hw(in, homework);
  return in;
}

istream& Grad::read(in) {
  read_common(in);
  in >> thesis;
  read_hw(in, homework);
  return in; 
}

// TODO




// --
istream& Student_info::read(istream& in)  {

  in >> name >> midterm >> final;
  read_hw(in, homework);
  return in;
  
}

double Student_info::grade() const 
