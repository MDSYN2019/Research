cmake_minimum_required (VERSION 1.1)
project (OMP_C++)
include(CheckCXXCompilerFlag)



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

if(WITH_QATZIP)
  find_package(qatzip REQUIRED)
  set(HAVE_QATZIP ${QATZIP_FOUND})
endif(WITH_QATZIP)

# TODO - how is this associated with virtual functions?
set (CMAKE_AUTOMOC ON)
set (CMAKE_INCLUDE_CURRENT_DIR ON)

# Ensure C++11 is invokedinclude(CheckCXXCompilerFlag)

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
include_directories (${MPI_CXX_INCLUDE_PATH})


find_package(Eigen3 REQUIRED)
include_directories ("${EIGEN3_INCLUDE_DIR}") # Eigen libaries

find_package(pybind11 REQUIRED)  # or add_subdirectory(pybind11)

# Boost directories 
include_directories (${Boost_INCLUDE_DIRS})
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})


# QAT directory
include_directories ("${PROJECT_SOURCE_DIR}/thirdparty/include")

if (USE_QAT)
  include_directories ("${PROJECT_CURRENT_SOURCE_DIR}/thirdparty/QATsrc")
  add_subdirectory(QATsrc)
  set(QAT_LIBS ${QAT_LIBS} QATsrc)
endif (USE_QAT)

option (USE_QAT "Use QAT implementation" OFF) 

#add_subdirectory (QAT)
#set (QAT ${QAT} "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/QAT/include")
#add_subdirectory (QATsrc)
#set (QATsrc ${QATsrc} "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/QAT/src")

add_subdirectory (src)
set(EXTRA_SRC ${EXTRA_SRC} src)
set(EXTRA_LIBS ${EXTRA_LIBS} include)

set(MATH ${MATH} "-lm")
set(CPPUNIT ${CPPUNIT} "-lcppunit")
set(LBOOST ${LBOOST} "-lboost_system")

set(IN ${IN} "-Iinclude")

# Include GSL directories
include_directories (${GSL_INCLUDE_DIRS} ${GSLCBLAS_INCLUDE_DIRS})

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(Boost REQUIRED filesystem program_options regex timer chrono) # Finding the chrono library .. Ugh, I do hate it

if(Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif()

find_package(MPI REQUIRED)

#set(QATHEAD ${QAT}/*.h)

## XXL

# print some info (this can be removed)
message(STATUS "STXXL_CXX_FLAGS: ${STXXL_CXX_FLAGS}")
message(STATUS "STXXL_INCLUDE_DIRS: ${STXXL_INCLUDE_DIRS}")
message(STATUS "STXXL_LIBRARIES: ${STXXL_LIBRARIES}")

# apply CXXFLAGS to our configuration
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${STXXL_CXX_FLAGS}")

# add STXXL include directory
include_directories(${STXXL_INCLUDE_DIRS})

#set (INC "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp")
#file(GLOB source_files ${INC}) # Linking src and include 

add_executable (OMP_EX ${QAT_LIBS} main.cxx)
add_library(gsllib main.cxx)
target_link_libraries(gsllib -L/home/noh/anaconda3/pkgs/gsl-2.1-2/lib/)
target_link_libraries(gsllib libgsl.a)
target_link_libraries(gsllib libgslcblas.a)
target_link_libraries(OMP_EX ${Boost_LIBRARIES} ${MATH} ${CPPUNIT} ${MPI_CXX_LIBRARIES} ${CMAKE_GSL_LINKER_LIB} ${QAT_LIBS} ${STXXL_LIBRARIES} ${LBOOST} pybind11::module gsllib)

find_package(OpenMP)

if(OpenMP_CXX_FOUND)
    target_link_libraries(OMP_EX  OpenMP::OpenMP_CXX)
endif()

