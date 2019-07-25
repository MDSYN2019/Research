#ifndef GUARD_VEC_H
#define GUARD_VEC_H

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

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

class Student_info {
public:
  std::string name() const {
    return n;
      }

  bool valid() const {return !homework.empty();} // return true when empty, false if not empty
  std::istream& read(std::istream&);
  double grade() const;
  
private:
  std::string n;
  double midterm, final;
  std::vector<double> homework;
};

bool compare (const Student_info&, const Student_info&);


class Core {
public:
  Core();
  Core(std::istream&);
  std::string name() const; // implemented in openmp_....cxx
  std::istream& read(std::istream&);
  double grade() const;
protected:
  std::istream& read_common(std::istream&);
  double midterm, final;
  std::vector<double> homework;
private:
  std::string n;
};


class Grad: public Core {
public:
  Grad();
  Grad(std::istream&);
  double grade() const;
  std::istream& read(std::istream&);
private:
  double thesis;
};


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
