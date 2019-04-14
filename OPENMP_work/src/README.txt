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
   we aggregated tasks by assigning a contiguous block to each thread
   
