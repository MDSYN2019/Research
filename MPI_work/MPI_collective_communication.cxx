#ifndef CC
#define CC

#include <iostream>
#include <string>

/*
  Parallel Trapezoidal Rule

  Input: None
  Output: Estimate of the integral from a to b of f(x) using the trapezoidal rule and n trapezoids

*/

extern int cat;

class Person {
public:
  Person();
  Person(const Person& anotherPerson); // copy constructor 
  Person& operator=(const Person& anotherPerson);
  virtual ~Person();
  std::string getName();

  std::string getName() const;
  std::string idNum() const;
private:
  std::string name;
  std::string idNum;
};

// inheriting from person

class Student : public Person {
  // All methods from Person, in addition to all the private variables in Person are now avaliable to this example as well
private:
  std::string major;
  int gradYear;
};
  
  
class T_rule {
public:

  float Trap(float, float, int, float);
  float f(float);

  void get_data(float*, float*, int*, int, int);
  void processTRMPI();
  void setupMPI();  

private:
  int my_rank;
  int p;
  float a = 0.0;
  float b = 1.0;
  int n = 1024;
  float h;
  float local_a;
  float local_b;
  int local_n;

  float integral;
  float total;
  int source;
  int dest = 0;
  int tag = 0;
 
  // -- //  
  MPI_Status status;
};

float T_rule::Trap(float local_a, float local_b, int local_n, float h) {
  float integral;
  float x;
  int i;  
  integral = (local_a + local_b) / 2.0;
  x = local_a;

  for (int i = 1; i <= local_n-1; i++) {

    x = x + h;

    integral = integral + f(x);
  }
}

void T_rule::processTRMPI() {

  h = (b - a)/n; // h is the same for all processes
  local_n = n/p; // So is the number of trapezoids
  
  // Length of each process's interval of integration =  n * h, so my intervaql starts at:
  
  local_a = a + my_rank * local_n*h;
  local_b = local_a + local_n * h;
  integral = Trap(local_a, local_b, local_n, h);

  if (my_rank == 0) {    
    total = integral;

    for (int source  = 1; source < p; source++) {
      
      MPI_Recv(&integral, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);
      total = total + integral; 
    }
  } else {
    MPI_Send(&integral, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);
  }
}

float T_rule::f(float x) { 
  float return_val;

  // Do something with x and return return_val

  return return_val;
}

/*
  One obvious problem with our program is its lack of generality. The function,
  f(x), and the input data, a, b and n, are hardwired.  
*/

void T_rule::get_data(float* a_ptr, float* b_ptr, int* n_ptr, int my_rank, int p) {

  int source = 0 ;
  int tag;

  if (my_rank == 0) {
    std::cout << "Enter a, b, and n " << std::endl;

    scanf("%f %f %d", a_ptr, b_ptr, n_ptr);

    for (dest = 1; dest < p; dest++) {

      tag = 0;

      MPI_Send(a_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);

      tag = 1;

      MPI_Send(b_ptr, 1, MPI_FLOAT, dest, tag, MPI_COMM_WORLD);

      tag = 2;

      MPI_Send(n_ptr, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
    }

  } else {
    tag = 0;
    MPI_Recv(a_ptr, 1, MPI_FLOAT, source, tag, MPI_COMM_WORLD, &status);    
  }
}

void T_rule::setupMPI() {
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
}

class Counter {
public:
  Counter();
  int getCount();
  void increaseBy(int x);  
private:
  int count;
};

int main() {
  return 0;
  
}

#endif
