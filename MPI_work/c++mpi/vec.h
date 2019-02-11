#ifndef __VEC__
#define __VEC__

template <class T> class Vec {
 public:
  typedef T* iterator;
  typedef const T* const iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef std::ptrdiff_t difference_type;
  
};

#endif
