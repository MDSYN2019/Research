#ifndef __MD__LAMDA__
#define __MD__LAMDA__

#include <vector>
#include <iostream>
#include <string>

/*

The reduce functions takes two required arguments and one additional, optional argument.

1. The reduction argument used to reduce a pair of arguments to a single result - i.e. sum takes two arguments 
   and returns the sum of the two arguments

2. The vector of values to be redued

3. An (optional) initial value that is used as the first value for the reduction

 */


/*

The benefit of map/reduce over a separate map and then reduce, is that we don't need to create an intermediate vector of results 
from the map. Results from the map can be immediately fed into the reduce as and when they are avaliable.

 */


template<class MAPFUNC, class REDFUNC, class ..ARGS>
auto mapReduce(MAPFUNC mapfunc, REDFUNC redfunc, const std::vector<ARGS>& ..args) {

  typedef typename std::result_of<MAPFUNC(ARGS...)>::type RETURN_TYPE;

  int nvals = detail::get_min_container_size(args...);
  if (nvals == 0) {
    return RETURN_TYPE();
  }

  RETURN_TYPE result = mapfunc(args[0]...);

  if (nvals == 1) {
    return result;
  }

  for (size_t i = 1; i < nvals; ++i) {    
    result = redfunc( result, mapfunc(args[i]...));
  }
  return result;
}


template<class FUNC, class T> // Take a template that takes a function and the value type

T reduce(FUNC func, const std::vector<T> &values) { 

  if (values.empty())  {
    return T();
  }

  else {

    T result = values[0];

    for (size_t i = 1; i < values.size(); ++i) {
      result = func(result, values[i]);
    }
    
    return result;
  }
}

template<class FUNC, class T>
auto map(FUNC func, const std::vector<T> &arg1, const std::vector<T> &arg2) {

  int nvalues = std::min(arg1.size(), arg2.size());
  auto result  = std::vector<T>(nvalues);

  for (int i = 0; i < nvalues; ++i ) {
    result[i] = func(arg1[i], arg2[i]);
  }

  return result; 
}

template<class T>
void print_vector(const std::vector<T> &values) {
  std::cout << "[";
  for (const T &value : values) {
    std::cout << " " << value;
  }
  std::cout << " ]" << std::endl;
}

#endif
