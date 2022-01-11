
#include <iostream>
#include <thread>

/* CPPunit tests */

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


class Distance {
private:
  int feet;
  int inches;
public:
  
  // required constructors
  Distance() {}
  Distance (int f, int i) {
    feet = f;
    inches = i;
  }
  // overload function call

  // overload function call

  Distance operator() (int a, int b, int c) { // Operator is () - i.e. we call a class() function that does something 
    Distance D;                               // which in this case, is simply makeing another Distance objec 
    // just put a random calculation
    D.feet = a + c + 100;
    D.inches = b + c + 100;
    return D;
  }

  // Method to display distance

  void displayDistance() {
    std::cout << "F: " << feet << "I: "<< inches << std::endl;
  }
};

struct func {
  int& i;
  
  func(int& i_): i(i_){}

  void operator()() {

    for (unsigned j = 0; j < 1000; j++) {
      std::cout << j << std::endl;
    }
  }
};

void f() { // Function that runs func

  int some_local_state = 0;

  std::thread t(func(some_local_state)); // Start a thread running the func
  try {
    // Anything we do here is done in this current thread
    std::cout << "Doing" << std::endk;
  }
  

}
// RAII - Resource Acquisition Initialization idiom (RAII) 


class thread_guard {
  std::thread& t;
public:
  explicit thread_guard(std::thread& t_): t(t_) {}
  
  ~thread_guard() { // a class with a join in its destructor 
    if (t.joinable()) {
      t.join();
    }
    thread_guard(thread_guard const&)=delete; // 
    thread_guard& operator=(thread_guard const&)=delete;
    
  }
};


/*

Running Threads in the Background
--------------------------------

Detarched threads are often called daemon threads 

*/

void counting() {
  for (int i = 0; i < 1000; i++) {
    std::cout << i << std::endl; 
  }
}

int main(void) {
  int placeholder = 3;
  func A(placeholder);
  std::cout << A.i  <<  std::endl;
  //  A();

  std::thread t(counting);
  t.join();
  
  //  std::thread t(counting);
  //t.join();

}
