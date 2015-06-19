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
project(simplecluster CXX)

# THE LOCATION OF OUTPUT BINARIES
set(CMAKE_LIBRARY_DIR ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_LIBRARY_DIR})
 
# INCLUDE DIRECTORIES
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR})
# For googletest
include_directories(${CMAKE_LIBRARY_DIR}/googletest/include)

# CHECK FOR SUBMODULES AND LIBS
# Check for OpenCV, just for comparing with our methods
set (CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
#find_package(OpenCV)
find_package(OpenMP)

# This project uses Google Testing Framework and
# we added GoogleTest as a sub-project. The version
# of this distribution is 1.7.0.
add_subdirectory(lib/googletest)

# SOURCE FILES
file(GLOB PROJECT_SRCS
    ${PROJECT_SOURCE_DIR}/src/*.cpp
)

# COMPLATION FLAGS
if(NOT MSVC)
	if (OPENMP_FOUND)
    	list (APPEND CMAKE_CXX_FLAGS 
    	"-fopenmp-simd -Wopenmp-simd -O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs ${OpenMP_CXX_FLAGS}")
	else()
    	list(APPEND CMAKE_CXX_FLAGS 
    	"-O3 -std=c++0x ${CMAKE_CXX_FLAGS} -g -ftest-coverage -fprofile-arcs")
	endif()
endif()

if(NOT MSVC)
    set(TEST_LIBS_FLAGS ${CMAKE_LIBRARY_DIR}/libgtest_main.a m gtest pthread ${CMAKE_LIBRARY_DIR}/libopenblas.a gfortran)
else()
    set(TEST_LIBS_FLAGS gtest openblas gfortran)
endif()

# BUILD OPENBLAS
add_custom_target(
    openblas
    COMMAND make && make PREFIX=${PROJECT_SOURCE_DIR} install
    WORKING_DIRECTORY ${CMAKE_LIBRARY_DIR}/OpenBLAS)

# OUTPUT LIBRARIES
add_library(${PROJECT_NAME} SHARED ${PROJECT_SRCS})
target_link_libraries(${PROJECT_NAME} ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(${PROJECT_NAME} openblas)

add_library(${PROJECT_NAME}_static STATIC ${PROJECT_SRCS})
target_link_libraries(${PROJECT_NAME}_static ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(${PROJECT_NAME}_static PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(${PROJECT_NAME}_static openblas)

# OUTPUT EXECUTABLES

## Basic tests
add_executable(test_utilities 
    ${PROJECT_SOURCE_DIR}/test/test_utilities.cpp 
    ${PROJECT_SOURCE_DIR}/src/utilities.cpp)
target_link_libraries(test_utilities ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_utilities PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(test_utilities gtest_main)
add_executable(test_kdtree 
    ${PROJECT_SOURCE_DIR}/test/test_kdtree.cpp 
    ${PROJECT_SOURCE_DIR}/src/utilities.cpp 
    ${PROJECT_SOURCE_DIR}/src/kd-tree.cpp)
target_link_libraries(test_kdtree ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_kdtree PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(test_kdtree gtest_main)
add_executable(test_kmeans 
    ${PROJECT_SOURCE_DIR}/test/test_kmeans.cpp 
    ${PROJECT_SRCS} )
target_link_libraries(test_kmeans ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_kmeans PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(test_kmeans gtest_main)

## A test for LAPACKE
add_executable(test_lapacke 
    ${PROJECT_SOURCE_DIR}/test/test_lapacke.cpp 
    ${PROJECT_SRCS} )
target_link_libraries(test_lapacke ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_lapacke PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(test_lapacke gtest_main)

## A test for okmeans
add_executable(test_okmeans 
    ${PROJECT_SOURCE_DIR}/test/test_okmeans.cpp 
    ${PROJECT_SRCS} )
target_link_libraries(test_okmeans ${TEST_LIBS_FLAGS})
if(MSVC)
    set_target_properties(test_okmeans PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
endif()
add_dependencies(test_okmeans gtest_main)

# Only build this example when found OpenCV
# if(OpenCV_FOUND AND NOT OpenCV_VERSION VERSION_LESS 2.4.0)
#    # OpenCV paths
#    set(OPENCV_INCLUDES "")
#    foreach(dir ${OpenCV_INCLUDE_DIRS})
#        include_directories(dir)
#    endforeach()
#    add_executable(test_kmeans 
#        ${PROJECT_SOURCE_DIR}/test/test_kmeans.cpp 
#        ${PROJECT_SRCS} )
#    target_link_libraries(test_kmeans ${TEST_LIBS_FLAGS} ${OpenCV_LIBRARIES})
#	if(MSVC)
#    	set_target_properties(test_kmeans PROPERTIES COMPILE_FLAGS "/MT ${OpenMP_CXX_FLAGS}")
#	endif()
#    add_dependencies(test_kmeans test_kdtree gtest_main)
#else()
#    message(WARNING "You do not have OpenCV >= 2.4.0! Your OpenCV version is ${OpenCV_VERSION}")
#endif()
