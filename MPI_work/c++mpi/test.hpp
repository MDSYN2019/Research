#ifndef __VECA__
#define __VECA__

#include <iostream>
#include <iomanip>
#include <string>

/*

*/

template <typename E>
class Stack {
  public:
  int size() const;
  bool empty() const;
  const E& top() const throw(StackEmpty);
  void push(const E& e);
  void pop() throw (StackEmpty);
};

class StackEmpty : public RuntimeException {
public:
  StackEmpty(const std::string& err) : RuntimeException(err){}
};

template <class T> class Vec {

public:
  Vec();
  Vec(T, T);
  void initiate(T, T);
  void out();
  ~Vec();

private:
  T data;
  T limit;
};

template <typename T>
Vec<T>::Vec() {}

template <typename T>
void Vec<T>::initiate(T a, T b) {
  data = a;
  limit = b;
}

template <typename T>
void Vec<T>::out() {
  std::cout << data << " "  << " " <<limit << std::endl;
}

template <typename T>
Vec<T>::~Vec() {}

#endif
