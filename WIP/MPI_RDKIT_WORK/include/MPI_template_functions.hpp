/*
  Last Updated: 16/11/2021
  ------------------------


*/

#include <vector>
#include <cassert>

template<typename T>
T max (T a, T b) {
  //
  return b < a ? a : b;
}

template<typename T>
class Stack {
private:
  std::vector<T> elems;
public:
  void push(T const& elem);
  void pop();
  T const& pop() const;
  bool empty() const {
    return elems.empty(); 
  }
};
  
