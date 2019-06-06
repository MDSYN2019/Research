cmake_minimum_required (VERSION 1.0)
project (OMP_C++)

# Enable Fortran
#enable_language(Fortran)
# Enable eigen (C++)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

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
include_directories ("${EIGEN3_INCLUDE_DIR}") # Eigen libaries

# Boost directories 
include_directories (${Boost_INCLUDE_DIRS})

# QAT directory
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/QAT/include")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/QAT/src")

#add_subdirectory (QAT)
#set (QAT ${QAT} "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/QAT/include")
#add_subdirectory (QATsrc)
#set (QATsrc ${QATsrc} "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/QAT/src")

add_subdirectory (src)
set (EXTRA_LIBS ${EXTRA_LIBS} src)
set (MATH "${MATH} -lm")
set (CPPUNIT "${CPPUNIT} -lcppunit")

# Include GSL directories
include_directories (${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(Boost) 

if(Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif()

find_package(MPI REQUIRED)

set(QATCOMP AbsFunction.cpp )
add_executable (OMP_EX ${UNIT_TEST} ${QATCOMP} main.cxx)
target_include_directories(OMP_EX PRIVATE)

add_library(gsllib main.cxx)
target_link_libraries(gsllib -L/home/noh/anaconda3/pkgs/gsl-2.1-2/lib/)
target_link_libraries(gsllib libgsl.a)
target_link_libraries(gsllib libgslcblas.a)
target_link_libraries(OMP_EX ${EXTRA_LIBS} ${Boost_LIBRARIES} ${MATH} ${CPPUNIT} ${MPI_CXX_LIBRARIES} ${CMAKE_GSL_LINKER_LIB}  gsllib)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(OMP_EX  OpenMP::OpenMP_CXX)
endif()

