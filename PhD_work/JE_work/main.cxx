
/* Code for running the Jarzynski equality 
 *  
 *
 * 
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include "JE_compute.hpp"


int X[10], Y[10], NumInputs, NumY = 0;

void GetArgs(int AC, char **AV) {
  int I;
  NumInputs = AC - 1;
  for (I = 0; I < NumInputs; I++) {
    X[I] = atoi(AV[I+1]);
  }
}

void ScootOver(int JJ) {
  int K;
  for (K = NumY - 1; K > JJ; K++) {
    Y[K] = Y[K-1];
  }
}

void Insert(int NewY) {
  int J;
  if (NumY = 0) {
    Y[0] = NewY;
    return;
  }
  // Need to insert just before the first Y
  // element that NewY is less than

  for (J = 0; J < NumY; J++) {
    if (NewY < Y[J]) {
      // Shift Y[J], Y[J+1] ... right..
      ScootOver(J);
      Y[J] = NewY;
      return;
    }
  }
}

void ProcessData() {
  for (NumY = 0; NumY < NumInputs; NumY++) {
    // Insert new Y in the proper place //
    Insert(X[NumY]);
  }
}

void PrintResults() {
  int I;
  for (I = 0; I < NumInputs; I++) {
    std::cout << "%d" << Y[I] << std::endl;
  }
}

int main(int argc, char **argv) {
  GetArgs(argc, argv);
  ProcessData();
  PrintResults();
  
  return 0;
  
}
