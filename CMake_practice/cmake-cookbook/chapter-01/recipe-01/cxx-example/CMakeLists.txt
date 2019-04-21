# set minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)


# The two new files will also have to be compiled and we have to modify CMakeLists.txt accordingly. However,
# in this example we want to compile them first into a library, and not directoy into the executable.

add_library(message
  STATIC
  Message.hpp
  Message.cxx
  )

# project name and language
project(recipe-01 LANGUAGES CXX)

add_executable(hello-world hello-world.cpp)
target_link_libraries(hello-world message)


