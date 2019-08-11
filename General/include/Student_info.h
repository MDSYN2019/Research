#ifndef GUARD_Student_info_h
#define GUARD_Student_info_h

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
#include "openmp_dynamicbindingandinheritance.hpp"


/*
  Our code became cluttered when we realized that we neeeded to be able to deal with 
  objects whose type we could not know until run time.

  We knew that each object would be either a Core, or something derived from core. 
  Our solution used pointers, bcause we could allocate a pointer to Core and then make that 
  pointer point to either a Core or a Grad object. The troiiuble with our soliton is that is 
  imposed err-rprobe bookkeeping on our users. We can't eliminiate that bookkeeping, but we 
  can hide it from our users writing a new class
  
 */

class Student_info {
 public:
 Student_info(): cp(0) {}
 Student_info(std::istream& is) : cp(0) {read(is);}
  Student_info(const Student_info&);
  Student_info& operator=(const Student_info&);
  ~Student_info() { delete cp;}
  // operations

  std::istream& read(std::istream&);

  std::string name() const {
    if (cp) {
      return cp->name();
    }
    else {
      throw std::runtime_error("blah Blah blah");
    }		       
  }

  double grade() const {
    if (cp) { return cp->grade();
    } else {
      throw std::runtime_error("uninitialized student");
    }
  }
  
  static bool compare (const Student_info& s1, const Student_info& s2) {
    return s1.name() < s2.name();
  }
  
 private:
  Core* cp;
};

/*
  The idea her eis that a Student info can represent either a Core 
or a grad.


Each Student_info object will hold a pointer, called cp, that points
to an object that has eiehter type Core or a type derived from Core.

In the read function, we'll allocate the object to which cp points

Therefore, both constructors initi\lize cp to 0, indiciating that the 
student_info is as yet unbund. In the constructor that takes an istream,
we call the Student_info::read function. That function will allocate 
a new object a new object of the appropriate type, and will give that object the value that read from the indicated istream





 */


#endif

