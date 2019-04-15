This includes some notes I made in the process of learning openMP (open multiple processing)

Wh\t actually happens when the program gets to the parallel directive? Prior to the parallel directive, the program is using a single thread, the process
started when the program started execution.

When the program reaches the parallel directive, the original thread continrues executing and (thread_count - 1) additonal threads are started.

In opneMP parlaance, the collection of threads executing the parallel block.

1. The original thread is called the master

2. The additional threads are called the slaves

Each thread executes the block following the directive, so in our example, each thread calls the hello function.

When the block of code is completed - in our example, when the threads return
from the call to hello, theres an implciit barrier. THis means that a thread has completed the block of code will wait for all the other threads in the team to complete the block.
 In the hello example, a thread that has completed the call to hello will wait for all the other threads in the team to return. when all the thread have ceompleted the block ,the slave thread will terminate and the master thread will continue executing the code that follows the block.


Example: The trapwezoidal rukle.

The trapezoidal rule for extimating the area under the curve. Recall that if
y = f(x) is a reasonably nice funciton and a < b are real numbers, then we can
estimate the area between the graph of f(x).

// Input: a, b, and n
h = (b-a) / n
approx = (f(a) + f(b))/2.0

Steps

1. We identified two types of tasks:
a. Computation of the areas of individual trapezoids, and
b. Adding the areas of trapezoids

2. There is no communication among the tasks in the first collection,
   but each task in the first collection communicates task 1(b).

3. We assumed that there would be many more trapezoids than cores, so
   we aggregated tasks by assigning a contiguous block to each thread.
   Effectively, this paritioned the interval [a,b] into larger subintervals,
   and each thread simply applied the serial trapezoidal rule to its
   subinterval.


We aren't quite done, however, since we still need to add up the thread's results. An obvious solution is to use a shared variable for the sums of all the thread's results, and each thread can add its (private) result into the shared variable.

The actual sequence of events might well be different but unlessone thread
finishes the computation global_result += my_result before the other
starts, the result will be incorrect. Recall that this is an exmaple of a
race condition: Multiple threads are attempting to access a shared resource,
at least one of the accesses is an update, and the accesses can result in an error.

We therefore need some mechanism to make sure that once one thread has started execting global_result += my_reuslt, no other thread can staet eecuting this code until the first thread has finished. 


=> Scope of variables

In serial programming, the scope of a variable consists of those parts of a program in which the variable can be used. For example, a variable declared at the
beginning of a C function has a function-wide scop,e that is it can only be accessed in the body of a function.

In openMP, the scope of a variable refers to the set of threads that can access
the variable in a parallel block.

1. A variable that can be accessed by all the threads in the team has shared scope.

2. A variable that can only be accessed by a single thread has private scope.

In the "hello World" program, the variables used by each thread (my_rank and thread_count) were declared in the Hello function, which is called inside the
parallel block.
	 Consequently, the variables used by each thread are allocated
	 from the thread's (private) stack, and hence all of the variables
	 have private scope.

=> The Reduction Clause



