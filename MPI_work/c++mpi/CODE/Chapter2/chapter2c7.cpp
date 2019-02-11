
#include <iostream.h>

float FloatMachineEps();
double DoubleMachineEps();

int main(int argc, char ** argv){
  float fep;
  double dep;
  
  fep = FloatMachineEps();
  dep = DoubleMachineEps();
  
  cout << "Machine epsilon for single precision floating point numbers is: " << fep << endl;;
  
  cout << "Machine epsilon for double precision floating point numbers is: " << dep << endl;
  
  
}

float FloatMachineEps(){
  float  fmachine_e, ftest;
  fmachine_e = 1.0;
  
  ftest = 1.0 + fmachine_e;
  while(1.0 != ftest){
    fmachine_e = fmachine_e/2.0;
    ftest = 1.0 + fmachine_e;
  }

  return fmachine_e;
}


double DoubleMachineEps(){
  double dmachine_e, dtest;
  dmachine_e = 1.0;
  
  dtest = 1.0 + dmachine_e;
  while(1.0 != dtest){
    dmachine_e = dmachine_e/2.0;
    dtest = 1.0 + dmachine_e;
  }

  return dmachine_e;
}
