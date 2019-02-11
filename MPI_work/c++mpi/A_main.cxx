#include <iostream>
#include <iomanip>
#include "test.hpp"

/*

Stacks and Queues
-----------------

We use the term container to denote a data structure that permmits storage and retrieval of data
independent of content. 

Containers are distinguished by the particular retrieval order they support. In the two 
most important types of containers, this retrieval order depends on the insertion order. 

- Stacks - Support retrieval by last-in, first out order. Stacks are simple to implement 
and very efficint. For this reason, stacks are probably the right container to use when 
retrieval order doesn't matter at all.

- Queues - Support retrieval of first in, first out. This is surely the fairest way to control
           waiting times for services. You want the container holding jobs to be processed 
           in FIFO (first in first out) order to minimize the maximum time spent waiting. 
           Many computing 

Stacks and Queues can be effectively implemented using either arrays or linked lists. The 
key issue is whether an upper bound on the size of the container is known in advance, 
thus permitting the use of a statically-allocated array.

Dictionaries
------------

The dictionary data type permits access to data items by content. You stick an item 
into a dictionary so you can find it when you need it.

The primary operations of a dictionary support are:

- search - Given a search key k, return a pointer to the element in dictionary D whose key 
           value is k, if one exists
- 

 */


int main () {

  Vec<int> A;
  A.initiate(1,2);
  A.out();
  
  return 0;
  

}
