cmake_minimum_required (VERSION 1.0)
project (OMP_C++)
# Enable Fortran
enable_language(Fortran)

if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none")
    set(bounds "-fbounds-check")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-stand f08 -free -implicitnone")
    set(bounds "-check bounds")
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "PGI")
    set(dialect "-Mfreeform -Mdclchk -Mstandard -Mallocatable=03")
    set(bounds "-C")
endif()

# Ensure C++11 is invoked
include(CheckCXXCompilerFlag)

FIND_LIBRARY(CPPUNIT_LIBRARY NAMES cppunit)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${DevIL_SOURCE_DIR}/cmake/Modules")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

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
#include_directories (${MPI_CXX_INCLUDE_PATH})

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

#find_package(MPI REQUIRED)
#include_directories(${MPI_INCLUDE_PATH})
#find_package(GSL)
#include_directories(${GSL_INCLUDE_DIR})

find_package(Boost) 

if(Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif()

find_package(MPI REQUIRED)
add_subdirectory (src)

set (EXTRA_LIBS ${EXTRA_LIBS} src)
set(MATH "${MATH} -lm")
set(CPPUNIT "${CPPUNIT} -lcppunit")


# include_directories ("${EIGEN3_INCLUDE_DIR}")
# Specifically add the src library definition for linking

# GSL libraries linking
#set(CMAKE_GSL_LINKER_LIB "${CMAKE_GSL_LINKER_LIB} -L/home/noh/anaconda3/pkgs/gsl-2.1-2/lib/")
#set(CMAKE_GSL_LINKER_FLAGS "${CMAKE_GSL_LINKER_FLAGS} -lgsl")
#set(CMAKE_GSL_LINKER_FLAGS_CBLAS "${CMAKE_GSL_LINKER_FLAGS_CBLAS} -lgslcblas")
# add the executable

include_directories (${Boost_INCLUDE_DIRS})
#target_compile_options(MPI_EX PRIVATE ${MPI_CXX_COMPILE_FLAGS})
add_executable (OMP_EX ${UNIT_TEST} main.cxx)
target_include_directories(OMP_EX PRIVATE)

ADD_LIBRARY(gsllib 
    main.cxx
)

target_link_libraries(gsllib -L/home/noh/anaconda3/pkgs/gsl-2.1-2/lib/)
target_link_libraries(gsllib libgsl.a)
target_link_libraries(gsllib libgslcblas.a)
target_link_libraries(OMP_EX ${EXTRA_LIBS} ${Boost_LIBRARIES} ${MATH} ${CPPUNIT} ${MPI_CXX_LIBRARIES} ${CMAKE_GSL_LINKER_LIB}  gsllib)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(OMP_EX  OpenMP::OpenMP_CXX)
endif()


