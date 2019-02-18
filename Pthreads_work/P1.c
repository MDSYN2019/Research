#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

int thread_count;

void* Hello(void* rank); // Thread function


int main() {

  long thread;
  pthread_t* thread_handles;

  thread_count = strtol(argv[1], NULL, 10);
  
}
