#pragma once

/*
Class Templates

Similar to functions, calsses can also be parameterized with one or more types. Container classe,s which 
are used to manage elements of a certain type, are a typical example of this feature. By using class templates, you 
can implement such container classes hile the element typei sis sitl open. IN this chapter we use a stack 
as a na example of a class template 

 */

#include <vector>
#include <stdexcept>


/*
Just on the return const value :

A function becomes const when const keyword is sued i th a functions delcaration. THe idea of const funcitons is not to allow them
to modift the object on which the are called. IT is recommended pracice to make as amny functions const as possible so that accidental changes to obejcts are avoided


 */


template <typename T>
class Stack {
private:
  std::vector<T> elems;
public:
  virtual void push(T const&); // push element
  void pop(); // pop element 
  T top() const; // return top element 

  bool empty() const { // 

  }
  // I will TRY to add the virtual example as seen in accelerated c++ 
};
