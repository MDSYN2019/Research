
#include <iostream>
#include <iomanip>
#include <istream>
#include <vector>
#include <map>
#include <fstream>


#include "grade.hpp"
#include "median.hpp"
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

Core::Core(std::istream& is) {
  read(is);
}

std::string Core::name() const {
  return n;
}

double Core::grade() const {
  return ::grade(midterm, final, homework); // 
}

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


/*
  A simple handle class
  --------------------
  Although the approach that we have just seen is straightfoerwad, it has problems: 
  the program has sacquired a lot of extra complexity related to managing pointers

  Our users have to remember to allocate space for the records as they read them,
  and to remember to free that space when they no longer need the data

  The code is constantly dereferencing the pointers to get at the underlying objects
  

  -- 

  What we'd like to do is to find a way to preserve the good properties of our 
  simpler programs 

  Our code became cluttered when we realized that we needed to be able to deal 
  with objects whose type we could not know until run time. We knew that 
  each object would be either a Core, or something derived from Core. 

  Our solution used pointers, because we could allocate a pointer to Core and then
  make that pointer point to either a Core or a Grad object.

  The trouble with our solution is that it imposed error-prone bookkeeping 
  on our users. We can't eliminate that bookkeeping, but we can hide it 
  from our users by writing a new class that will encapsulate the pointer to Core 
  
  Key point : Encapsulate the pointer to Core 

*/

/*
class Student_info {
  // constructor and copy control
  Student_info(): cp(0) {}
  Student_info(std::istream is) : cp (0) { read(is)}
  Student_info(const Student_info&);
  Student_info& operator=(const Student_info&);
  ~Student_info() {delete cp;}
  // operations
  std::istream read&(std::istream&);

  std::string name() const {
    if (cp) return cp->name;
    else throw std::runtime_error("uninitialilized student");
  }

  double grade() const {
    if (cp) return cp->grade();
    else throw std::runtime_error("unintialized student");
  }

  static bool compare (const Student_info& s1, const Student_info& s2) {
    return s1.name() < s2.name();
  }

private:
  Core* cp;
};
*/
