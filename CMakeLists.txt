# You should have cmake >= 2.8 installed on your system.
# For Ubuntu, just run 'sudo apt-get install cmake'
cmake_minimum_required(VERSION 2.8)

# Project's name
project (simple-k-means)
 
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
 
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})
 
set(PROJECT_LIB_DIR ${PROJECT_SOURCE_DIR}/lib)
 
include_directories("${PROJECT_LIB_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")


# This project uses Google Testing Framework and
# we added GoogleTest as a sub-project. The version
# of this distribution is 1.5.0.
add_subdirectory(test/googletest)

set(PROJECT_SRCS 
${PROJECT_SOURCE_DIR}/src/k-means.cpp
${PROJECT_SOURCE_DIR}/src/util.cpp
)
add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS})
