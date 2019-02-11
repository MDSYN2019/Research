#include <iostream>
#include <vector>
#include <string>
#include <typeinfo>

#include "MD_friend_function.hpp"

// boost

using namespace boost; 

// Generic class 
template <class T>
Rectangle<T> duplicate (const Rectangle<T>& param)  {
  Rectangle<T> res;
  res.width = param.width*2; // accessing the private members 
  res.height = param.height*2; // ditto 
  return res; 
}


