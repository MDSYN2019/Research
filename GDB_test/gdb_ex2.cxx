#include <iostream>

long factorial(int n);

int main() {

  int n(0);

  std::cin >> n;
  
  long val = factorial(n);

  std::cout << val;

  std::cin.get();

  
  return 0;
  
}

long factorial(int n) {

  long result(1);

  while (n--) {
    result *= n;
  }
  
  return result;
}
