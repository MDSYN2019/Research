#include <iostream>
#include <string>

template <class T> class Vec {
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef std::ptrdiff_t difference_type;
  typedef T& reference;
  typedef const T& const_reference;

  Vec() {create();}
  explicit Vec(size_type n, const T& val = T()) { create(n, val);}
private:
  iterator data;
  iterator limit;   
};

class Box {

public:
  explicit Box(int width, int length, int height) : m_width(width), m_length(length), m_height(height) {}
  
  int Volume() {
    return m_width * m_length * m_height;
  }
  void p() {
    std::cout << m_width;
  }

  void set_m_width(int new_width) {
    Box::m_width = new_width;
  }

  Box operator+() 
private:
  int m_width;
  int m_length;
  int m_height;
};

/*
class StorageBox : public Box {
public:
  StorageBox(int width, int length, int height, const std::string label&, 
  
};
*/

int main() {

  Box A(2,3,4);
  A.p();
  A.set_m_width(3);
  
  return 0;
  
}
