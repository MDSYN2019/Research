#pragma once

#include <iostream>
#include <cassert>
#include <boost/tuple/tuple.hpp>
#include <string>

using namespace boost; 

template <class T>
class Rectangle {
protected:
  T width, height;
public:
  void set_values (T a, T b) {
    width = a;
    height = b;
  }

  Rectangle() {} // constructor 
  Rectangle(T x, T y) : width(x), height(y) {}; // automatucally assign x to width, and y to height
  
  T area() {return width * height;}

  template<class S>
  friend Rectangle<S> duplicate (const Rectangle<S>&); // defining friend class that has access to the 
};



/* Inheritance  */

template <class T>

class square : public Rectangle<T> {

public:
  T aarea() {    
    return this->width * this->height;
  }
};

enum class Suit {Diamonds, Hearts, Clubs, Spades};
void PlayCard (Suit suit) {
  switch(suit) {
  case Suit::Diamonds : std::cout << "Diamonds \n"; break;
  case Suit::Hearts : std::cout << "Hearts \n"; break;
  case Suit::Clubs : std::cout << "Clubs \n"; break;
  case Suit::Spades : std::cout << "Spades \n"; break;	  
  }
}

class FriendF {
public:  
  // constructor
  FriendF(){FriendVar = 0;}
  ~FriendF();
private:
  int FriendVar;
  friend void FFriend(FriendF &sfo); // stinkfirstobject = sfo 
}; 

void FFriend(FriendF &sfo) {
  sfo.FriendVar = 99;
  std::cout << sfo.FriendVar << std::endl;}

FriendF::~FriendF() {
};

