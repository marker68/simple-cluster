# You should have cmake >= 2.8 installed on your system.
# For Ubuntu, just run 'sudo apt-get install cmake'
cmake_minimum_required(VERSION 2.8)

# Project's name
project(simple_cluster)
 
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
 
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
 
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
 
include_directories("${PROJECT_INCLUDE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")

# Check for OpenCV, just for comparing with our methods
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(OpenCV)

# This project uses Google Testing Framework and
# we added GoogleTest as a sub-project. The version
# of this distribution is 1.5.0.
add_subdirectory(test/googletest)
# For googletest
include_directories("${PROJECT_SOURCE_DIR}/test/googletest/include")

set(PROJECT_SRCS 
${PROJECT_SOURCE_DIR}/src/k-means.cpp 
${PROJECT_SOURCE_DIR}/src/utilities.cpp
)

# Some of the features need C++ 11
list(APPEND CMAKE_CXX_FLAGS "-Wall -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs -lm")

# Create a .so file
add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS})
add_dependencies(${PROJECT_NAME} gtest_main)

# Make some tests
set(CMAKE_CXX_FLAGS "-L${CMAKE_BINARY_DIR} ${CMAKE_CXX_FLAGS} ${CMAKE_BINARY_DIR}/libgtest_main.a -lgtest -pthread")
add_executable(test_kmeans ${PROJECT_SOURCE_DIR}/test/test_kmeans.cpp ${PROJECT_SRCS})
add_dependencies(test_kmeans gtest_main)