#include <functional>
#include "lib_mpi.hpp"
#include <vector>
#include <boost/current_function.hpp>
#include <boost/foreach.hpp>
#include <boost/static_assert.hpp>
#include <boost/detail/lcast_precision.hpp>

using namespace part1;

template<class FUNC, class ARG1, class ARG2>
auto call_function(FUNC func, ARG1 arg1, ARG2 arg2) { // function, argument1, argument2

  auto result = func(arg1, arg2);
  return result;

}

int sum(int x, int y) {
  return x + y;
}

int main (int argc, char **argv) {

  int result = sum(3,7);
  std::cout << result << std::endl;
  return 0;
}
