#ifndef _CORE_H
#define _CORE_H
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <istream>


std::istream& read_hw(std::istream& in, std::vector<double>& hw) {
  if (in) {
    hw.clear();
    double x;
    while (in >> x) {
      hw.push_back(x);
      // Clear the stream so that input will work for the next student
      in.clear();
    } 
  }  
  return in;
}

class Core {
public:
  Core() {};
  // ~Core();
  Core(std::istream&);
  // virtual destructor
  //virtual ~Core();
  std::string name() const; // implemented in openmp_....cxx
  std::istream& read(std::istream&); // Virtual because we need this to be dynamically bound because of the inherited method having an identical name 
  double grade() const; // By using a virtual implementation, the function will now determine which function to run (the original or inherited version) binpsecting each object
  
protected: // protection label allows inherited objects to use the variables/functions
  std::istream& read_common(std::istream&);
  double midterm, final;
  std::vector<double> homework; 
private:
  std::string n;
  std::filebuf fb;
};

class Grad: public Core { // inherit from core
public:
  Grad();
  Grad(std::istream&);
  double grade() const;
  std::istream& read(std::istream&);  
private:
  double thesis;  
};


/*
Our users ahve to remember to allocate space for the records as they read them.


 */

#endif
