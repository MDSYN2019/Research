#pragma once

/*
Default template arguments

For class templates you can also define defualt values for template parameters. These 
values are called default template arguments.


  
 */


#include <vector>
#include <stdexcept>

template <typename T, typename CONT = std::vector<T> > // second input - default template arguments
class Stack {
private:
  CONT elems; // elements
public:
  void push(T const&); // push element 
  void pop(); // pop element 
  T top() const; // return top element 
  bool empty() const { // return whether the stack is empty 
    return elems.empty();
  }
};

template <typename T, typename CONT>
void Stack<T, CONT>::push(T const& elem) {
  elems.push_back(elem);
}
