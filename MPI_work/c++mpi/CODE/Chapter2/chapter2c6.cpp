
#include <iostream.h>

int main(int argc, char ** argv){
  int i,xn;
  int max_iterations = 1000;
  int initial_guess = 100;
 
  xn = initial_guess;
  
  for(i=0;i<max_iterations;i=i+1){
    
    cout << i << " " << xn << endl;
    
    if(xn == 1)
      break;
    
    if(xn%2==0)
      xn = xn/2;
    else
      xn = 3*xn+1;
  }
}



