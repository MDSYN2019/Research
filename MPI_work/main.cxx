#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <numeric>
#include "mpi.h"
#include "MPI_broadcast.hpp"

MPI_BC AA;

int main () {
  AA.add_vector();
  AA.broadcast_vector();
  return 0;
}
