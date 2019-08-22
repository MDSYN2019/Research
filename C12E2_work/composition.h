//! Composition study for the C12E2 budding paper
/*! 
Placeholder
 */

#ifndef __composition__
#define __composition__

#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <utility>
#include <cerrno>
#include <cstdlib> 
#include <cmath>
#include <string>

int numberofSS = 100; /*The number of screenshots in the dump file*/
const int numberOfPolymers = 1000; // The number of polymers of each type - C12E2 or mimic 
const int numberofatoms = 71313; // Total number of beads in the simulation
const int indexCG = 7;
const int boxdim = 3;

class compute {
 public:
  compute();
  compute(std::string);
  void storeFile();
  void printVectorElements();
  void AllocationVec();
 private:
  // Vectors to store trajectory values 
  std::vector<double> xco; 
  std::vector<double> yco; 
  std::vector<double> zco; 
  std::vector<int> a;
  std::vector<int> b;
  std::vector< std::vector<double>> xcoTotal; 
  std::vector< std::vector<double>> ycoTotal; 
  std::vector< std::vector<double>> zcoTotal; 
  std::vector< std::vector<int>> aTotal;
  std::vector< std::vector<int>> bTotal;
  //double xco[numberofatoms],yco[numberofatoms],zco[numberofatoms]; 
  /* ditto for a and b, which represent index and atomtype respectively */
  //int a[numberofatoms],b[numberofatoms];
  // double xco[numberofatoms], yco[numberofatoms], zco[numberofatoms]; 
  FILE *ipf; /* input file */  
  // COM vectors - C12E2
  std::vector<double> headGroupC12_E2xCOM; // COM C12E2 X
  std::vector<double> headGroupC12_E2yCOM; // COM C12E2 Y
  std::vector<double> headGroupC12_E2zCOM; // COM C12E2 Z
  // COM vectors - C12E2-M
  std::vector<double> headGroupMIMICxCOM; // COM MIMIC X
  std::vector<double> headGroupMIMICyCOM; // COM MIMIC Y
  std::vector<double> headGroupMIMICzCOM; // COM MIMIC Z
  // -----------------------------------------------------------------------//
  // Defining normals (A7 -- A6_1 -- A6_2 -- A3_1 -- A3_2 -- A3_3 -- A4     //
  // -----------------------------------------------------------------------//
  int A7, A6_1, A6_2, A3_1, A3_2, A3_3, A4; // First Batch
  int A7_2, A6_1_2, A6_2_2, A3_1_2, A3_2_2, A3_3_2, A4_2; // Second Batch
  int A7_3, A6_1_3, A6_2_3, A3_1_3, A3_2_3, A3_3_3, A4_3; // Third Batch 
  int A7_4, A6_1_4, A6_2_4, A3_1_4, A3_2_4, A3_3_4, A4_4; // Fourth Batch

  // ----------------------------------------------------------------------//
  // Defining mimics (A13 -- A12_1 -- A12_2 -- A9_1 -- A9_2 -- A9_3 -- A10 //
  // ----------------------------------------------------------------------//  

  int A13, A12_1, A12_2, A9_1, A9_2, A9_3, A10; // First Batch
  int A13_2, A12_1_2, A12_2_2, A9_1_2, A9_2_2, A9_3_2, A10_2; // Second Batch
  int A13_3, A12_1_3, A12_2_3, A9_1_3, A9_2_3, A9_3_3, A10_3; // Third Batch
  int A13_4, A12_1_4, A12_2_4, A9_1_4, A9_2_4, A9_3_4, A10_4; // Fourth Batch
  int MimicCounter = 0;    
  int PolymerCounter = 0;    
  int atomtype; 
  int index, l, n;
  int nlines = numberofatoms + 9;      
  double tophead = 0;
  double downhead = 0;
  double mimictophead = 0;
  double mimicdownhead = 0;  
  double x,y,z; /*coordinates for the atoms in the box*/
  double box1;
  double box2; 
  double boxlength[boxdim];
  char line[100];  
  struct C12E2_skeleton C12E2_struct[numberOfPolymers];                                                                                              
  struct C12E2_skeleton C12E2M_struct[numberOfPolymers];    
};

#endif
