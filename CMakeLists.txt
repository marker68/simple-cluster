# SIMPLE CLUSTERS: A simple library for clustering works.
# Copyright (C) 2014 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# You should have cmake >= 2.8 installed on your system.
# For Ubuntu, just run 'sudo apt-get install cmake'
cmake_minimum_required(VERSION 2.8)

# Project's name
project(simplecluster)
 
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
${PROJECT_SOURCE_DIR}/src/kd-tree.cpp 
)

# Some of the features need C++ 11 and math library
find_package(OpenMP)
if (OPENMP_FOUND)
    list (APPEND CMAKE_CXX_FLAGS 
    "-Wall -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs ${OpenMP_CXX_FLAGS}")
else()
    list(APPEND CMAKE_CXX_FLAGS 
    "-Wall -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")
endif()
# Create a shared library file
add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS})
add_dependencies(${PROJECT_NAME} gtest_main)

# Make some tests
set(TEST_LIBS_FLAGS ${CMAKE_BINARY_DIR}/libgtest_main.a m gtest pthread)
add_executable(test_utilities 
    ${PROJECT_SOURCE_DIR}/test/test_utilities.cpp 
    ${PROJECT_SOURCE_DIR}/src/utilities.cpp)
target_link_libraries(test_utilities ${TEST_LIBS_FLAGS})
set_target_properties(test_utilities PROPERTIES APPEND_STRING 
    PROPERTY COMPILE_FLAGS " -L${CMAKE_BINARY_DIR} ${CMAKE_CXX_FLAGS}")
add_dependencies(test_utilities gtest_main)
add_executable(test_kdtree 
    ${PROJECT_SOURCE_DIR}/test/test_kdtree.cpp 
    ${PROJECT_SOURCE_DIR}/src/utilities.cpp 
    ${PROJECT_SOURCE_DIR}/src/kd-tree.cpp)
target_link_libraries(test_kdtree ${TEST_LIBS_FLAGS})
set_target_properties(test_kdtree PROPERTIES APPEND_STRING 
    PROPERTY COMPILE_FLAGS " -L${CMAKE_BINARY_DIR} ${CMAKE_CXX_FLAGS}")
add_dependencies(test_kdtree test_utilities gtest_main)

# Only build this example when found OpenCV
if(OpenCV_FOUND AND NOT OpenCV_VERSION VERSION_LESS 2.4.0)
    # OpenCV paths
    set(OPENCV_INCLUDES "")
    foreach(dir ${OpenCV_INCLUDE_DIRS})
        set(OPENCV_INCLUDES "-I${dir} ${OPENCV_INCLUDES}")
    endforeach(dir)
    add_executable(test_kmeans 
        ${PROJECT_SOURCE_DIR}/test/test_kmeans.cpp 
        ${PROJECT_SRCS} )
    target_link_libraries(test_kmeans ${TEST_LIBS_FLAGS} ${OpenCV_LIBRARIES})
    set_target_properties(test_kmeans PROPERTIES APPEND_STRING 
        PROPERTY COMPILE_FLAGS "${OPENCV_INCLUDES} -L${CMAKE_BINARY_DIR} ${CMAKE_CXX_FLAGS}")
    add_dependencies(test_kmeans test_kdtree gtest_main)
else(OpenCV_FOUND AND NOT OpenCV_VERSION VERSION_LESS 2.4.0)
    message(WARNING "You do not have OpenCV >= 2.4.0! Your OpenCV version is ${OpenCV_VERSION}")
endif(OpenCV_FOUND AND NOT OpenCV_VERSION VERSION_LESS 2.4.0)