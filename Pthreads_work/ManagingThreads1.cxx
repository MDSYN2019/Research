
#include <iostream>
#include <thread>

/* CPPunit tests */

#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/CompilerOutputter.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>


/*
Example of a function call operator
 */

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

void hello() {

  std::cout << "Hello concurrent world" << std::endl;
  
}



struct func {
  
  int& i;

  func (int& i_): i(i_) {}

  void operator()() { // operator callable by brackets () that takes in no parameters

    for (unsigned j =  0; j < 10000; ++j) {
      std::cout << i << std::endl;     
    }    

  }
};

void oops() {
  int some_local_state = 0;
  std::thread my_thread(func(some_local_state));
}



void f() {
  int some_local_state = 0;

  std::thread t(func(some_local_state)); // open thread with some local state that does 
 
  try {
    do_something_in_current_thread();
  }

  catch (..)
    {
      t.join();
      throw; // If there is an error, the thread finishes
    }
  t.join(); // finish the thread after computatin has finished
  
}

/*

Running threads in the background
---------------------------------

Detached threads are often called demon threads after the unix concept of daemon 
process that runs in the background without any explitiy uder interface. Such threads 


*/


// using RAII to wait for a thread to complete

class thread_guard {
  std::thread& t;

public:
  explicit thread_guard(std::thread& t_): t(t_) {}
  ~thread_guard() {
    if (t.joinable()) {
      t.join();
    }
    thread_guard(thread_guard const&)=delete;
    thread_guard& operator=(thread_guard const&)=delete;
  }
};

int main(void) {

  oops();
  std::thread t(hello);
  t.join();

}
