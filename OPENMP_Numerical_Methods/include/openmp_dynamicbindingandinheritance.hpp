#ifndef GUARD_VEC_H
#define GUARD_VEC_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <istream>

/*
struct Student_info {
  std::string name;
  double midterm, final;
  std::vector<double> homework;  

  // member functions
  std::istream& read(std::istream&); // Added
  double grade() const; // added  
};
*/

//bool compare (const Student_info&, const Student_info&);

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
  std::istream& read(std::istream&);
  double grade() const;
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

bool compare (const Core& c1, const Core& c2) {
  return c1.name() < c2.name();  
}

/*

Operations


To complete our classes, we need to implement four 
constructors: the defaukt constructor and the constructor that 
takes an istream, once for each class.

We must also implment six operations: the name and read_common 
operations in class Core, and the read and grae functions for both classes. 

Before writing our code, we need to think about how student records
will be structured. As befroe, we'll want to accomodate a variable 
number of homework assignments, so those grades must come at the end of each record.

 */

#endif
