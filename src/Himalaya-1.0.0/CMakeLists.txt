cmake_minimum_required(VERSION 3.5)
project(Himalaya CXX Fortran)

# Locate Eigen, unless the caller has already said where it is
# (-DEIGEN3_INCLUDE_DIR=...; SOFTSUSY's build passes the directory its own
# configure found, so that both use the same Eigen).  Not fatal if this
# fails: the include path may instead be supplied through CMAKE_CXX_FLAGS,
# and the compiler will complain if it is missing.
if(NOT EIGEN3_INCLUDE_DIR)
  find_package(Eigen3 QUIET)
endif()

# if needed, set your Fortran and C++ compiler
#set (CMAKE_Fortran_COMPILER gfortran)
#set (CMAKE_CXX_COMPILER clang++-3.8)

# set the build type to release. If one needs debugging symbols, compile
# with "CMAKE_BUILD_TYPE Debug"
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set version number
set(Himalaya_VERSION_MAJOR 1)
set(Himalaya_VERSION_MINOR 0)
set(Himalaya_VERSION_RELEASE 0)

# set the CXX flags
# Eigen >= 3.4.90/5.0 requires at least C++14; honour a -std= given by the
# caller, otherwise default to C++14 (the Himalaya sources are C++11-clean).
if(NOT CMAKE_CXX_FLAGS MATCHES "-std=")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
endif()
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# set paths
set(INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source/include)

#include c++ header files
include_directories(${INCLUDE_PATH})
include_directories(${INCLUDE_PATH}/hierarchies)
if(EIGEN3_INCLUDE_DIR)
  include_directories(${EIGEN3_INCLUDE_DIR})
elseif(TARGET Eigen3::Eigen)
  get_target_property(_eigen_inc Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
  include_directories(${_eigen_inc})
endif()

# setting the correct version to all headers
configure_file (
  "${INCLUDE_PATH}/version.hpp.in"
  "${INCLUDE_PATH}/version.hpp"
  )

# himalaya files
file(GLOB_RECURSE SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/*.cpp)
set(LIBSOURCES ${SOURCES})

# mh2l files
file(GLOB_RECURSE FSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/include/mh2l/*f)

# build the himalaya and DSZ library
file(GLOB_RECURSE HEADERS ${CMAKE_CURRENT_SOURCE_DIR}/source/*.hpp)
file(GLOB LIBSOURCES ${SOURCES} ${HEADERS})
list(REMOVE_ITEM LIBSOURCES ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)
add_library(DSZ STATIC ${FSOURCES})
add_library(Himalaya STATIC ${LIBSOURCES})

# The example program is not needed by SOFTSUSY; build it only on request
# (cmake -DHIMALAYA_BUILD_EXAMPLE=ON .)
option(HIMALAYA_BUILD_EXAMPLE "Build the Himalaya example program" OFF)
if(HIMALAYA_BUILD_EXAMPLE)
  # set the executable and specify the c++ version
  add_executable(example ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)
  target_compile_features(example PRIVATE cxx_range_for)

  # link mh2l and himalaya to the executable
  target_link_libraries(example Himalaya DSZ)
endif()

# install
install(TARGETS Himalaya ARCHIVE DESTINATION lib)
