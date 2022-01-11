#ifndef __T__
#define __T__

/*

Launching a thread
------------------

Threads are started by constructing a std::thread object that specifies the task to run on that 
thread. In the simplest case, that task is just a plain, ordinary void returning function that takes 
no parameters.

It doesnt matter what the thread is going to do, or where its launched from, but starting a thread
using the C++ thread library always boils down to constructing a std::thread object.

void do_some_work();



 */

struct func {
  int& i;
  func (int& i_):i(i_){}
  void operator()() { // operator callable by brackets () that takes in no parameters
    for (unsigned j =  0; j < 10000; ++j) {
      std::cout << j << std::endl;     
    }
    
  }
};


void oops() {
  int some_local_state = 0;
  std::thread my_thread(func(some_local_state));

}

/*
  Waiting for a thread to complete

  If you need for a thread to complete, this can be done by calling join()
  on the assoicated std::thread instance.

  Inserting a call to my_thread.join() before the closing brace of the function
  body would therefore be sufficient to ensure that the thread was finished before the function 
  was exited.

  

 */

#endif 
