#include <iostream>

typedef struct list {
  double item;
  struct list *next;
} list;

list *search_list(list *l, double x) {
  if (l == NULL) return NULL;

  if (l->item == x) {

    return l; 

  } else {
    return (search_list(l->next, x));
  }
}


// Insertion into a linked list

void insertion_list(list **l, double x) {
  list *p;
  p = new list [sizeof(list)]; // allocate a new list of size list on the heap
  p->item = x;
  p->next = *l;
  *l = p;
}

// Deletion from the list

/*
Deletion from a linked list is somewhat more complicated. First, we must find a 
pointer to the predecessor of the item to be deleted. 
 */

//list *predecessor_list(list *l, item_type x) {
//
//}

int main() {

  list A;
  list B;
  A.next = &B;
  list *C;
  list **D;
 
  A.item = 2.0;
  B.item = 3.0;
  
  std::cout << A.item << " " << B.item << std::endl;
  C = search_list(&A, 3.0);
  std::cout << C->item << std::endl;

  insertion_list(D,3.0);
  return 0;
  
}
