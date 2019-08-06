#ifndef _CORE_H_
#define _CORE_H_

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

// Include for the virtual function, which needs to be described for the non-inherited class explicitly

#include "grade.hpp"
#include "median.hpp"

std::istream& read_hw(std::istream& in, std::vector<double>& hw);

class Core {
  friend class Student_info;
public:
  // Constructors/destructors
  Core(): midterm(0), final(0) { }
  Core(std::istream& is) { read(is); }
  virtual ~Core() {}
  
  std::string name() const; // implemented in openmp_....cxx
  virtual std::istream& read(std::istream& in);
  // Virtual because we need this to be dynamically bound because of the inherited method having an identical name 
  virtual double grade() const {return ::grade(midterm, final, homework);}  // By using a virtual implementation, the function will now determine which function to run (the original or inherited version) binpsecting each object

protected: // protection label allows inherited objects to use the variables/functions
  std::istream& read_common(std::istream&);
  double midterm, final;
  std::vector<double> homework; 
  std::string n;
  // TODO
  //  virtual Core* clone() const {return new Core(*this);}
private:
  //std::filebuf fb;
};

class Grad: public Core { // inherit from core
public:
  Grad(): thesis(0) { }
  Grad(std::istream& is) { read(is); }
  double grade() const { return std::min(Core::grade(), thesis); }
private:
  double thesis;  
};


/*
  TODO
 */

class Student_info {
  // constructor and copy control
public:
  Student_info(): cp(0) {}
  Student_info(std::istream& is) : cp(0) {read(is);}
  Student_info(const Student_info&);

  
private:
  Core* cp;
};
#endif
