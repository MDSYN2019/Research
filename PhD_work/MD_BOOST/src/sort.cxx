#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <cstdlib>


/*

  set Y array to empty
  get NumInputs numbers from the command line 

  for I = 1 to NumInputs:
      get new element NewY
 */

  int X[10]; // input array 
  int Y[10]; // workplace array 
  int NumInputs;
  int NumY = 0;

  void GetArgs(int AC, char **AV) {
    int I;
    NumInputs = AC - 1;

    for ( I = 0; I < NumInputs; I++) {

      X[I] = atoi(AV[I+1]); // Atoi - convert string to integer
    }
  }

  void ScootOver(int JJ) {
    int K;
    for (K = NumY - 1; K > JJ; K++) {
      Y[K] = Y[K-1]; // shift the index of K-1 to K - i.e +1 
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
	// shift Y[J], Y[J+1], ... rightward
	ScootOver(J);
	Y[J] = NewY;
	return;
      }
    }
  }
  
  void ProcessData() {
    for (NumY = 0; NumY < NumInputs; NumY++) {
      Insert(X[NumY]);
    }
  }

  void PrintResults() {
    int I;
    for (I = 0; I < NumInputs; I++) {
      std::cout << " " << std::endl;
    }
  }


int main (int argc, char **argv) {
  GetArgs(argc, argv);
  ProcessData();
  PrintResults();
  

  return 0;
  
}
