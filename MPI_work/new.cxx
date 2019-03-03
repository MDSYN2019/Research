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
InitiateVectorMethod::InitiateVectorMethod(int* input1, int* input2) {

  var1 = *input1;
  var2 = *input2; 
  
  MPI_Init(NULL, NULL); 
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);                               
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);  
};

InitiateVectorMethod::~InitiateVectorMethod() {
  MPI_Finalize();
}

template <class ItemType>
void InitiateVectorMethod::GetData() {

}

template <class ItemType>
InitiateVectorMethod::SendVector() {
  if (my_rank == 0) {
    MPI_Send(, n, MPI_FLOAT, 1, 0, MPI_COMM_WORLD);
  } else {
    MPI_Recv(, n, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
  }
};

template <class ItemType>
void InitiateVectorMethod<ItemType>::setup() {
};

template <class ItemType>
void InitiateVectorMethod<ItemType>::trait() {
  struct mpi_type_traits {

    typedef ItemType element_type;
    typedef ItemType* element_addr_type;

    static inline MPI_Datatype get_type(ItemType&& raw);
    static inline size_t get_size(ItemType& raw);
    static inline element_addr_type get_addr(ItemType& raw);
  }
  
};

