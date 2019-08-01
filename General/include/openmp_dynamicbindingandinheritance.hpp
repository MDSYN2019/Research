#ifndef _CORE_H
#define _CORE_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <istream>

/*

Inheritence and constructors

Before we write the constructors for Core and Grad, we need to understand how the implementation
creates object of a derived type. As with any class type, the implementation 
begins by allocating space for the object. The fact that the object is of a derived type adds 
an extra step to the construction process in order to construct the base-class 
part of object. 

Derived objects are constructed by:


-> Allocating space for the entire object (base-class members as well as as derived members)

-> Calling the base-class constructor to initialize the base-class parts of the object

-> Iniitlaizing the members of the derived class as directed by constructor initializer

-> Executing the body of the derived class constructor 


 */


/*
Dynamic binding

The run-time selection of the virtual function to exceut is relveant only when the function is called throguh a reference
or pointer

 */


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
  Core():  midterm(0), final(0) {};
  Core(std::istream&);
  std::string name() const; // implemented in openmp_....cxx
  virtual std::istream& read(std::istream&); // Virtual because we need this to be dynamically bound because of the inherited method having an identical name 
  virtual double grade() const; // By using a virtual implementation, the function will now determine which function to run (the original or inherited version) binpsecting each object
  
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


#endif
