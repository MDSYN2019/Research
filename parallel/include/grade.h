#pragma once

#include <vector>

// Overloaded functions - multiple interpretations of the same functions


double grade(double, double, double);
double grade(double, double, const std::vector<double>&);
