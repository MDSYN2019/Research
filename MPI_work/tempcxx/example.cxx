#include <iostream>

template <typename T>
inline T max(T a, T b) {
    return a > b ? a : b;
}

int main() {

  std::cout << max(3,7) << std::endl;
  return 0;
}
