#include<iostream.h>


int main(int argc, char **argv){
  int sum;

  sum = 0;

  for(int i=1;i<=1000;i=i+1)
     sum = sum + i;

  cout << "The sum from 1 to 1000 is: " << sum << endl;	
}
