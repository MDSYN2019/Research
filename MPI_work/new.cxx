#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include <array>
#include <cassert> 
#include "mpi.h"
#include "new.hpp"

template <class ItemType>
InitiateVectorMethod::InitiateVectorMethod() {
  MPI_Init(NULL, NULL); 
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);                               
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
};

template <class ItemType>
InitiateVectorMethod::setup() {

};

template <class ItemType>
InitiateVectorMethod::trait() {
  struct mpi_type_traits {

    typedef ItemType element_type;
    typedef ItemType* element_addr_type;

    static inline MPI_Datatype get_type(ItemType&& raw);
    static inline size_t get_size(ItemType& raw);
    static inline element_addr_type get_addr(ItemType& raw);

  }
};

