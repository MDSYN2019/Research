
// The usual defence against accidentally including the file twice, which you saw in the last exercise 
#ifndef __syn_dbg_h__
#define __syn_dbg_h__

#include <iostream>
#include <cstdio>
#include <string>
#include <cerrno>
//#include <cstring>

#ifndef NDEBUG
#define debug(M, ...)
#else
#define debug(M. ...) std::cerr << "DEBUG %s:%d:" << M << __FILE__ << __LINE__ << ##__VA_ARGS__;
// C version
//#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n",__FILE__, __LINE__, ##__VA_ARGS__)
#endif 

// The clean_errno macro that's used in the others to get a safe, readable version of errno. That
// strange syntax in the middle is a ternary operator and you'll learn what it does later

#define clean_errno() (errno == 0 ? "None" : std::strerror(errno))


#endif 
