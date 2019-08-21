#include <iostream>
#include <stdexcept>

#include "openmp_dynamicbindingandinheritance.hpp"
#include "Student_info.h"



/*

The read function has three responsibilities: It must free the object, if any, to which 
this handle is bound. It must then decide what kind of obejct we are about to read, and 
it must allocate the right kind of obecjt, which it can initialize from the stream it was 
given

 */

std::istream& Student_info::read(std::istream& is) {

  delete cp;

  char ch;
  is >> ch;

  if (ch == 'U') cp = new Core(is);
  else if (ch == 'G') cp = new Grad(is);
  else if (ch == 'P') cp = new PassFail(is);
  else if (ch == 'A') cp = new Audit(is);
  else throw std::runtime_error("read invalid student type");

  return is;
}

Student_info::Student_info(const Student_info& s) : cp(0) {
  if (s.cp) {
    cp = s.cp->clone();    
  }
}





Student_info& Student_info::operator=(const Student_info& s) {
  if (&s != this) {
    delete cp;
    if (s.cp) cp = s.cp->clone();
    else cp = 0;
  }

  return *this;
}

