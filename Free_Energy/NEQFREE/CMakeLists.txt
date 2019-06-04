# Version number 
cmake_minimum_required (VERSION 1.0)
project (MD_JE_interpretation)

include(CheckCXXCompilerFlag)

# Ensure C++11 is invoked
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

# option for using MPI functions

option(USE_MPI "Use MPI functions for the compilation of this code" ON)

# Add the conditionals for USE_MPI

if (USE_MPI)
  add_subdirectory(MPI)
  list(APPEND MPILIBS MPI)
  list(APPEND MPI_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/MPI")
endif(USE_MPI)

if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
elseif(COMPILER_SUPPORTS_CXX14)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

# Include directories for headers and CXX files

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/src")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories (${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

# GSL libraries linking

set(CMAKE_GSL_LINKER_FLAGS "${CMAKE_GSL_LINKER_FLAGS} -lgsl")
set(CMAKE_GSL_LINKER_FLAGS_CBLAS "${CMAKE_GSL_LINKER_FLAGS_CBLAS} -lgslcblas")
set(MATH "${MATH} -lm")

# Find package boost 
find_package(Boost) 
if(Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif()

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

# Configuring MPI libraries
find_package(MPI REQUIRED)
if (MPI_FOUND)
  include_directories(${MPI_INCLUDE_PATH})
endif()

# Find Eigen libraries

find_package (Eigen3 3.3 REQUIRED NO_MODULE)
if (Eigen3_FOUND) 
  include_directories ("${EIGEN3_INCLUDE_DIRS}")
# Specifically add the src library definition for linking
endif()

add_subdirectory (src)
set (EXTRA_LIBS ${EXTRA_LIBS} src)

include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/MPI")

add_executable (JE main.cxx)
target_link_libraries (JE  ${EXTRA_LIBS} ${BOOST_LIBRARIES} ${Boost_LIBRARIES} ${MATH} ${EIGEN3_INCLUDE_DIR} ${MPI_LIBRARIES} ${MPI_INCLUDES} Eigen3::Eigen)

# Compile flags with MPI
if(MPI_COMPILE_FLAGS)
  set_target_properties(JE PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(JE PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
