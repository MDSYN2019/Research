
/*

Adopted from Zed Shaw's macros

Handling errors is a difficult actiivity in almost every programming langauge. There are entire programming languages 
that try as hard as they can to avoid even the concept of an error. Otjher languages invent complex control 
structures like exceptions to pass error conditions around.

C (C++) tackles the problem by returning error codes and setting a global errno value that you can check 

 */


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
/*
The alternative #define that translates any use of debug ("format", arg1, arg2) into fprintf 
call to stderr. Many C programmers on't know this, but you can create macros taht actually work like 
printf and take variable arguemnts. Some C compilers 

## __VA
 */

// C version
//#define debug(M, ...) fprintf(stderr, "DEBUG %s:%d: " M "\n",__FILE__, __LINE__, ##__VA_ARGS__)
#endif 

// The clean_errno macro that's used in the others to get a safe, readable version of errno. That
// strange syntax in the middle is a ternary operator and you'll learn what it does later

#define clean_errno() (errno == 0 ? "None" : std::strerror(errno))

#define log_err(M, ...) fprintf(stderr, "[WARN] (%s:%d : errorno: %s)" M "\n" __FILE__, __LINE__, clean_errno(), ##__VAR_ARGS__) 

#define log_warn(M, ...)

#endif 
