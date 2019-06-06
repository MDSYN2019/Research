#ifndef __MD__libraries__structures__
#define __MD__libraries__structures__

/*
*/

#include <Eigen/Dense>
#include <cmath>
#include <vector>

// Using the Eigen library
using Eigen::MatrixXd;

typedef struct {
  MatrixXd r = MatrixXd(1,1);
  MatrixXd rv = MatrixXd(1,1);
  MatrixXd ra = MatrixXd(1,1);
} Mol;

typedef struct {
  int val, sum, sum2;
} Prop;

class MD_structures {
private:  
  MatrixXd dr = MatrixXd(1,1); // You must initialize the size of the container before 
  MatrixXd vector = MatrixXd(1,1); // Ditto ;   
  MatrixXd VecI = MatrixXd(1,1); // Ditto ;   
  MatrixXd initUcell = MatrixXd(1,1);
  MatrixXd region, vSum = MatrixXd(1,1);
  MatrixXd c, gap = MatrixXd(1,1);
  int j1, j2, n, nx, ny;
  int nMols;
  int deltaT;
  int rCut, temperature, timeNow, uSum, velMag;
  int moreCycles, nMol, stepAvg, stepCount, stepEquil, stepLimit;
  double fcVal, rr, rrCut, rri, rri3;
  Mol mol;
  Prop kinEnergy, pressure, totEnergy;

  // STL libraries  
  std::vector<Mol> positionArray;
public:
  MD_structures(); // constuctor
  MD_structures(int); // constuctor with int input
  ~MD_structures(); // destructor
  // PBC functions
  void Wrap(int, int);
  inline void allocMat();
  inline void printArray();
  void computeForce();
  void LeapfrogStep(int);
  void initCoords();
  void initVels();
  void SetParams();
  void EvalProps(); // TODO
  void setNum(int N);
};

// What are the global variabls used by this program?

Mol *mol;
MatrixXd region, vSum;
MatrixXd initUcell;
Prop kinEnergy, pressure, totEnergy;

#endif 
