#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <algorithm>
#include <iterator>
#include <mpi.h>
#include "MPI_broadcast.hpp"

/*
template <class T> class Vec {
public:
  typedef T* iterator;
  typedef const T* const_iterator;
  typedef size_t size_type;
  typedef T value_type;
  typedef std::ptrdiff_T difference_type;
  typedef T& reference;
  typedef const T& const_reference;

  template <class T>
  Vec<T>& Vec<T>::operator=(const Vec& rhs) {
   // check for self-assignment
    if (&rhs != this) {
      // free the array in the left hand side
      uncreate();
      create (rhs.begin(), rhs.end());
    }
    return *this;
  }
  
  &Vec operator=(const Vec&);
  Vec() {create();}
  // we'll assume that one of our utiity functions will handle the allocation and copy so that the copy constucorr can forward its wworkd to taht function
  Vec(const Vec& v){(create(v.begin(), v.end()));} // copy constructor
  explicit Vec(std::size_t n, const T& val = T()) {create(n,val);} // this explicit constructor takes a size_type and a value - this will allocate neough mempry of type T of number n, and initialize it with the values val 
  size_type size() const {return limit - data;}
  
  T& operator[] (size_type i) { return data[i]}
  const T& operator[](size_type i) const {return data[i];}
  iterator begin() {return data;}
  const_iterator begin() const {return data;}
  iterator_end() {return limit;}
  const_iterator end() {return limit;}
private:
  iterator data; // first the first element of the data
  iterator limit;
};
*/

std::map<std::string, std::string> typeConvDict; // TODO

enum MPI_TYPE {MPI_CHAR, MPI_SHORT, MPI_INT, MPI_LONG, MPI_LONG_LONG, MPI_UNSIGNED_CHAR, MPI_UNSIGNED_SHORT, MPI_UNSIGNED, MPI_UNSIGNED_LONG, MPI_FLOAT, MPI_DOUBLE, MPI_LONG_DOUBLE, MPI_BYTE, MPI_PACKED};

void my_bcast(void* data, int count, MPI_Datatype datatype, int root,
              MPI_Comm communicator) {
  int world_rank;
  MPI_Comm_rank(communicator, &world_rank);
  int world_size;
  MPI_Comm_size(communicator, &world_size);

  if (world_rank == root) {
    // If we are the root process, send our data to everyone
    int i;
    for (i = 0; i < world_size; i++) {
      if (i != world_rank) {
        MPI_Send(data, count, datatype, i, 0, communicator);
      }
    }
  } else {
    // If we are a receiver process, receive the data from the root
    MPI_Recv(data, count, datatype, root, 0, communicator, MPI_STATUS_IGNORE);
  }
}

MPI_BC::MPI_BC() {
  MPI_Init(NULL, NULL);
  MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
} // constructor 

MPI_BC::~MPI_BC() {
  MPI_Finalize();
} // destructor 

void MPI_BC::InitializeVec(int lenOfVec) : vectorOfBlockLengths(lenofVec), MPItype(lenOfVec), MPIDatatype(lenOfVec), MPIdisplacements(lenOfVec) {}

void MPI_BC::parallelAllocateVec(double* aa, double* bb, int lenOfVec, std::vector<int>& vecpart, MPI_Datatype* input_mpi_t_p) {
  std::iota(MPItype.begin(), MPItype.end(), 1); // Vector allocation of types
  std::iota(MPIDatatype.begin(), MPIDatatype.end(), MPI_INT); // Vector allocation of MPI_INt
  std::iota(MPIdisplacements.begin(), MPIdisplacements.end(), sizeof(int));  // vector allocation of the size of the vector 
  
  MPI_Get_address(&vecpart[0], aint);
  MPI_Type_create_struct(lenOfVec, MPItype, MPIdisplacements, MPItype, input_mpi_t_p);
  MPI_type_commit(input_mpi_t_p);
}

void MPI_BC::buildMpiType(double* a_p, double* b_p, int* n_p, MPI_Datatype* input_mpi_t_p) {
  /*
    
    A derived datatype can bbe used to represent any collection 
    of data items by storing both the types of items and their 
    relative locations in memory. 

    If a function that sends data knows the types and the relative 
    locations in memory of a collection of data items, 
    
  */
  
  int array_of_blocklengths[3] = {1,1,1};

  MPI_Datatype array_of_types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_INT};
  MPI_Aint a_addr, b_addr, n_addr;
  MPI_Aint array_of_displacements[3] = {0};  

  MPI_Get_address(a_p, &a_addr);
  MPI_Get_address(b_p, &b_addr);
  MPI_Get_address(n_p, &n_addr);

  array_of_displacements[1] = b_addr - a_addr;
  array_of_displacements[2] = n_addr - a_addr;
  MPI_Type_create_struct();
  MPI_Type_commmit(input_mpi_t_p);
}
// Build MPI type

void MPI_BC::Get_input(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) { // input, input, input, output, output
  if (my_rank == 0) {
    std::cout << "Enter a, b and n \n";
    scanf("%lf %lf %d", a_p, b_p, n_p); 
  }
  MPI_Bcast(a_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(b_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(n_p, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  
  // Get Input
}

void MPI_BC::Get_input2(int my_rank, int comm_sz, double* a_p, double* b_p, int* n_p) { // input, input, input, output, output

  MPI_Datatype input_mpi_t;
  buildMpiType(a_p, b_p, n_p, &input_mpi_t); // TODO
  // Get Input
  if (my_rank == 0) {

  }
}

void MPI_BC::Send(float a, float b, int n, int dest) {
    MPI_Send(&a, 1, MPI_FLOAT, dest, 0, MPI_COMM_WORLD);
    MPI_Send(&b, 1, MPI_FLOAT, dest, 1, MPI_COMM_WORLD);
    MPI_Send(&n, 1, MPI_INT, dest, 2, MPI_COMM_WORLD);
} /* Send */

void MPI_BC::Receive(float* a_ptr, float* b_ptr, int* n_ptr, int     source) {
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, 0, MPI_COMM_WORLD, &status);
    MPI_Recv(b_ptr, 1, MPI_FLOAT, source, 1, MPI_COMM_WORLD, &status);
    MPI_Recv(n_ptr, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
} /* Receive */

void MPI_BC::GetData() {
  if (my_rank == 0) {
    std::cout << "Enter a, b, and n \n";
    scanf("%lf %lf %d", a_ptr, b_ptr, n_ptr);
  }
 
  
}
