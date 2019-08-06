#include <stdexcept>
#include <vector>

#include "grade.hpp"
#include "median.hpp"

// compute a student's overall grade from midterm and final exam grades and homework grade

double grade(double midterm, double final, double homework) {
  return 0.2 * midterm + 0.4 * final + 0.4 * homework;
}


/*
  Compute a student's overall grade from midterm and final exam grades and vector of 
  homework grades. 

  This function does not copy irts aeguemnt because median does so for us
*/

double grade(double midterm, double final, const std::vector<double>& hw) {
  if (hw.size() == 0) throw std::domain_error("Student has done no homrwork");
  return grade(midterm, final, median(hw));
}
