#ifndef _UNITS_
#define _UNITS_

// This file contains conversion factors from SI units to other common systems


namespace units {
  // Length

  static const double meter = 1.0;
  static const double cm = 0.01;
  static const double mm = 0.001;
  static const double micron = 1.0E-6;
  static const double nm = 1.0E-9;
  static const double feet = 0.3048;
  static const km  = 1000.0;
  static const mile = 1609.34;

};


/*
Parallelizing the trapezoidal rule

1. Partition the problem solution into tasks

2. Identify the communication channels between the tasks

3. Aggregate the tasks into composite tasks

4. Map the composite tasks to cores


In the partitioning phase, we usually try to identify as many tasks as possible. For the trapezoidal rule, 
we might identify two types of tasks - one type is finding the area of a single trapezoid, and the other 
is computing the sum of these areas.

1
 */
