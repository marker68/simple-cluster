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

# PROJECT'S NAME
project(simplecluster)

# CHECK FOR SUBMODULES AND LIBS
# Check for OpenCV, just for comparing with our methods
set (CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
find_package(OpenCV)
find_package(OpenMP REQUIRED)
# This project uses Google Testing Framework and
# we added GoogleTest as a sub-project. The version
# of this distribution is 1.7.0.
add_subdirectory(test/googletest)

# THE LOCATION OF OUTPUT BINARIES
if(MSVC)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin/Debug)
else()
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
endif()
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
 
# INCLUDE DIRECTORIES
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories("${PROJECT_INCLUDE_DIR}")
include_directories("${PROJECT_SOURCE_DIR}")
link_directories("${CMAKE_BINARY_DIR}")
# For googletest
include_directories("${PROJECT_SOURCE_DIR}/test/googletest/include")

# SOURCE FILES
set(PROJECT_SRCS 
${PROJECT_SOURCE_DIR}/src/k-means.cpp 
${PROJECT_SOURCE_DIR}/src/utilities.cpp
${PROJECT_SOURCE_DIR}/src/kd-tree.cpp 
)

# COMPLATION FLAGS
if(NOT MSVC)
	if (OPENMP_FOUND)
    	list (APPEND CMAKE_CXX_FLAGS 
    	"-Wall -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs ${OpenMP_CXX_FLAGS}")
	else()
    	list(APPEND CMAKE_CXX_FLAGS 
    	"-Wall -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")
	endif()
endif()

# OUTPUT LIBRARIES
add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS})
if(MSVC)
	set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/MTd ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(${PROJECT_NAME} gtest_main)

# OUTPUT EXECUTABLES
if(NOT MSVC)
set(TEST_LIBS_FLAGS ${CMAKE_BINARY_DIR}/libgtest_main.a m gtest pthread)
else()
set(TEST_LIBS_FLAGS gtest)
endif()
add_executable(test_utilities 
    ${PROJECT_SOURCE_DIR}/test/test_utilities.cpp 
    ${PROJECT_SOURCE_DIR}/src/utilities.cpp)
target_link_libraries(test_utilities debug ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_utilities PROPERTIES COMPILE_FLAGS "/MTd")
endif()
add_dependencies(test_utilities gtest_main)
add_executable(test_kdtree 
    ${PROJECT_SOURCE_DIR}/test/test_kdtree.cpp 
    ${PROJECT_SOURCE_DIR}/src/utilities.cpp 
    ${PROJECT_SOURCE_DIR}/src/kd-tree.cpp)
target_link_libraries(test_kdtree debug ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_kdtree PROPERTIES COMPILE_FLAGS "/MTd")
endif()
add_dependencies(test_kdtree test_utilities gtest_main)

# Only build this example when found OpenCV
if(OpenCV_FOUND AND NOT OpenCV_VERSION VERSION_LESS 2.4.0)
    # OpenCV paths
    set(OPENCV_INCLUDES "")
    foreach(dir ${OpenCV_INCLUDE_DIRS})
        include_directories(dir)
    endforeach()
    add_executable(test_kmeans 
        ${PROJECT_SOURCE_DIR}/test/test_kmeans.cpp 
        ${PROJECT_SRCS} )
    target_link_libraries(test_kmeans debug ${TEST_LIBS_FLAGS} ${OpenCV_LIBRARIES})
    if(MSVC)
        set_target_properties(test_kmeans PROPERTIES COMPILE_FLAGS "/MTd ${OpenMP_CXX_FLAGS}")
    endif()
    add_dependencies(test_kmeans test_kdtree gtest_main)
else()
    message(WARNING "You do not have OpenCV >= 2.4.0! Your OpenCV version is ${OpenCV_VERSION}")
endif()
