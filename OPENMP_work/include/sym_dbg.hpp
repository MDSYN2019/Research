#ifndef __syn_dbg_h__
#define __syn_dbg_h__

#include <iostream>
#include <cstdio>
#include <string>
#include <cerrno>
#include <cstring>

#ifndef NDEBUG
#define debug(M, ...)
#else
#define debug(M. ...) std::cerr << "DEBUG %s:%d:" << M << __FILE__ << __LINE__ << ##__VA_ARGS__;
// C version
//#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n",__FILE__, __LINE__, ##__VA_ARGS__)
#endif 

#define clean_errno() (errno == 0 ? "None" : std::strerror(errno) << "\n";


#endif 
