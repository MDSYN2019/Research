

#include <iostream.h>
#include <iomanip.h>
#include "SCchapter4.h"

double f(double x);

int main(int  argc, char * argv[]){
  int i,m;
  const int levels = 10;
  double value1,value2,value3;
  double exact = .2; 
  double xleft =  0.0, xright = 1.0;


  cout << setprecision(6) << setiosflags(ios::scientific);
  cout << "Level\tMidpoint\tTrapezoidal\tSimpsons\n";
  for(i=0;i<levels;i++){
    m = (int) pow(2,i);
    value1 = MidpointRule(i, xleft, xright,f); 
    value2 = TrapezoidRule(i, xleft, xright,f); 
    value3 = SimpsonsRule(m, xleft, xright,f); 
    cout << i << "\t" << fabs(value1-exact) << "\t" <<
      fabs(value2-exact) << "\t" <<
      fabs(value3-exact) << "\t" << endl;
  }
}

double f(double x){
  double value = x*x*x*x;
  return value;
}

