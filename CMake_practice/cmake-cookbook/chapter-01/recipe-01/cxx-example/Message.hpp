#pragma once // Non-standard but widely supported preprocessor directive designed to cause the current source file
             // to be included
#include <iosfwd>
#include <string>

class Message {
public:
  Message(const std::string &m) : message_(m){}
  friend std::ostream &operator<< (std::ostream &os, Message &obj) {
  }
private:
  std::string message_;
  std::ostream &printObject(std::ostream &os);
};
