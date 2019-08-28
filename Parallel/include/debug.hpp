#ifndef __dbg__h_
#define __dbg__h_

#include <iostream>
#include <cstudio>
#include <cerrno>
#include <string>

/*

Debug macros

Handling errors is a difficult activity in almost every programming langauge. 
There are entire programming languages that try as har das they can to avoid
even the concept of an error. 

Other languages invent complex control structures like exceptions to pass 
error conditions arround. The problem exists mostly because programmers 
assume errors don't hapen, and this optimism infects the types of langauges they 
use and create


 */

#ifndef NDEBUG
#define debug(M, ...)
#else
#define debug(M, ...) fprint(stderr, "DEBUG %s:%d: " M "\n", __FILE__, __LINE__, ##__VA_ARGS__)
#endif

#define clean_errno() (errno == 0 ? "None" : strerror(errno))


#define log_err(M, ...) fprintf(stderr, "[ERROR] (%s:%d)
#endif
