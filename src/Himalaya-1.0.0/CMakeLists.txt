cmake_minimum_required(VERSION 3.0)
project(Himalaya CXX Fortran)

# require an Eigen installation
find_package(Eigen3 REQUIRED)

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
set(CMAKE_CXX_FLAGS "-fPIC -std=c++11 -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

# set paths
set(INCLUDE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/source/include)

#include c++ header files
include_directories(${INCLUDE_PATH})
include_directories(${INCLUDE_PATH}/hierarchies)
include_directories(${EIGEN3_INCLUDE_DIR})

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

# set the executable and specify the c++ version
add_executable(example ${CMAKE_CURRENT_SOURCE_DIR}/source/example.cpp)
target_compile_features(example PRIVATE cxx_range_for)

# link mh2l and himalaya to the executable
target_link_libraries(example Himalaya DSZ)

# install
install(TARGETS Himalaya ARCHIVE DESTINATION lib)
