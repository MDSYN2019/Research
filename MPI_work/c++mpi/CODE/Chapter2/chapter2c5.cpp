
#include <iostream.h>

int main(int argc, char ** argv){
  int i,xn;
  int initial_guess = 100;


  xn = initial_guess;
  i = 0;

  while(xn != 1){

    cout << i << " " << xn << endl;

    if(xn%2==0)
      xn = xn/2;
    else
      xn = 3*xn+1;

    i=i+1;
  }
}



