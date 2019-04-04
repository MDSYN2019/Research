cmake_minimum_required (VERSION 1.0)
project (MPI_C++)
#find_package( OpenCV REQUIRED COMPONENTS core)
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Ensure C++11 is invoked

include(CheckCXXCompilerFlag)

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

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})
#find_package(GSL)
#include_directories(${GSL_INCLUDE_DIR})

find_package(Boost) 
find_package(CPPUNIT REQUIRED)
include_directories(${CPPUNIT_INCLUDE_DIR})
if(CPPUNIT_FOUND)
	set(CPPUNIT_LIBRARIES ${CPPUNIT_LIBRARY} ${CMAKE_DL_LIBS})
	set(CPPUNIT_INCLUDE_DIRS "${CPPUNIT_INCLUDE_DIR}")
	mark_as_advanced(CPPUNIT_ROOT_DIR)
endif()

#Include the "unit-tests" directory
include_directories("unit-tests")

if(Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif()

find_package(MPI REQUIRED)
add_subdirectory (src)

set (EXTRA_LIBS ${EXTRA_LIBS} src)
set(MATH "${MATH} -lm")

# include_directories ("${EIGEN3_INCLUDE_DIR}")
# Specifically add the src library definition for linking

# GSL libraries linking
#set(CMAKE_GSL_LINKER_LIB "${CMAKE_GSL_LINKER_LIB} -L/home/noh/anaconda3/pkgs/gsl-2.1-2/lib/")
#set(CMAKE_GSL_LINKER_FLAGS "${CMAKE_GSL_LINKER_FLAGS} -lgsl")
#set(CMAKE_GSL_LINKER_FLAGS_CBLAS "${CMAKE_GSL_LINKER_FLAGS_CBLAS} -lgslcblas")
# add the executable

include_directories (${Boost_INCLUDE_DIRS})
#target_compile_options(MPI_EX PRIVATE ${MPI_CXX_COMPILE_FLAGS})
add_executable (MPI_EX main.cxx)
target_include_directories(MPI_EX PRIVATE ${MPI_CXX_INCLUDE_PATH})

ADD_LIBRARY(gsllib 
    main.cxx
)

target_link_libraries(gsllib -L/home/noh/anaconda3/pkgs/gsl-2.1-2/lib/)
target_link_libraries(gsllib libgsl.a)
target_link_libraries(gsllib libgslcblas.a)
target_link_libraries(MPI_EX ${EXTRA_LIBS} ${Boost_LIBRARIES} ${MATH} ${MPI_CXX_LIBRARIES} ${CMAKE_GSL_LINKER_LIB} gsllib)

if(MPI_COMPILE_FLAGS)
  set_target_properties(MPI_EX PROPERTIES
    COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
  set_target_properties(MPI_EX PROPERTIES
    LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
