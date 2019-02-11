
#include <iostream.h>

float f(float x);

int main(int argc, char ** argv){
   float w;
   w = f(2);
   cout << "The value of f(2) is: " << w << endl;
}

float f(float x){
   float y;
   y = x*x*x - x*x + 2;
   return y;
}


