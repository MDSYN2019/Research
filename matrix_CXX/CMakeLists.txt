cmake_minimum_required (VERSION 1.0)
project (quant)
#find_package( OpenCV REQUIRED COMPONENTS core)
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

# Ensure C++11 is invoked

include(CheckCXXCompilerFlag)

CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

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

set(Boost_USE_STATIC_LIBS ON) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost) 
if(Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS})
endif()


find_library( CPPUNIT_LIBRARY_DEBUG NAMES cppunit cppunit_dll cppunitd cppunitd_dll
            PATHS   ${FOO_PREFIX}/lib
                    /usr/lib
                    /usr/lib64
                    /usr/local/lib
                    /usr/local/lib64
            PATH_SUFFIXES debug )

find_library( CPPUNIT_LIBRARY_RELEASE NAMES cppunit cppunit_dll
            PATHS   ${FOO_PREFIX}/lib
                    /usr/lib
                    /usr/lib64
                    /usr/local/lib
                    /usr/local/lib64
		    PATH_SUFFIXES release )


if(CPPUNIT_LIBRARY_DEBUG AND NOT CPPUNIT_LIBRARY_RELEASE)
    set(CPPUNIT_LIBRARY_RELEASE ${CPPUNIT_LIBRARY_DEBUG})
endif(CPPUNIT_LIBRARY_DEBUG AND NOT CPPUNIT_LIBRARY_RELEASE)

set( CPPUNIT_LIBRARY debug     ${CPPUNIT_LIBRARY_DEBUG}
  optimized ${CPPUNIT_LIBRARY_RELEASE} )





include_directories ("${EIGEN3_INCLUDE_DIR}")
# Specifically add the src library definition for linking

add_subdirectory (src)
set (EXTRA_LIBS ${EXTRA_LIBS} src)

# add the executable
if(Boost_FOUND)
  include_directories (${Boost_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
  add_executable (quant main.cxx)
  target_link_libraries (quant ${EXTRA_LIBS} ${BOOST_LIBRARIES} ${Boost_LIBRARIES} ${MATH} ${CMAKE_GSL_LINKER_FLAGS} ${CMAKE_GSL_LINKER_FLAGS_CBLAS}  ${CPPUNIT_LIBRARY})
endif()

