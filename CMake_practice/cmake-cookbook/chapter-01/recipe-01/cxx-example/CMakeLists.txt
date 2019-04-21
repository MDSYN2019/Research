# set minimum cmake version
cmake_minimum_required(VERSION 3.5 FATAL_ERROR)

# project name and language
project(recipe-01 LANGUAGES CXX)

add_library(message-objs
  OBJECT
  Message.hpp
  Message.cxx
  )

# This is only needed for older compilers
# but doesn't hurt either to have it
set_target_properties(message-objs
  PROPERTIES
  POSITION_INDEPENDENT_CODE 1
  )


add_library(message-shared
  SHARED
  $<TARGET_OBJECTS:message-objs>
  )


add_library(message-static
  STATIC
  $<TARGET_OBJECTS:message-objs>
  )

add_executable(hello-world hello-world.cpp)

target_link_libraries(hello-world message-static)


