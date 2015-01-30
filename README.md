Simple Clustering
===============================

[![Build Status](https://travis-ci.org/marker68/simple-cluster.svg?branch=master)](https://travis-ci.org/marker68/simple-cluster)

## What is this?

This is a simple implementation of the state-of-the-art clustering methods such as k-means, EM algorithm, ... 
The implementations based on many papers that are mentioned below. This project is under GNU GPL v3 License. Just read the license statement.
If you need any other licenses, feel free to contact me at [t_nguyen@hal.t.u-tokyo.ac.jp](mailto:t_nguyen@hal.t.u-tokyo.ac.jp).
The purpose of this project is to provide a more flexible interface for clustering methods.

## Features

* Supported k-means algorithm.
* Supported [CMake](http://www.cmake.org/).
* Supported only L2 metric distance.
* Supported KD-tree with ANN search.
* Supported GNU C++ Compiler and clang compiler.

## Installation

### Prerequisites

* An Unix or Windows([Cygwin](https://www.cygwin.com/) not native Windows with an MSVC compiler) Operating System. Tested on Mac OS X 10.9 and Ubuntu 14.04.
* [CMake](http://www.cmake.org/) 2.8 or newer. For UNIX users, check your CMake version in terminal by `cmake -version`.
* An GNU C++ Compiler or clang compiler that support C++ 11.
* [OpenCV](http://opencv.org/downloads.html) 2.4.8 or newer. This library is used for testing only. Just for comparing the performance with our implementations.

### Build
We use CMake as the build system. On terminal,
```bash
$ git clone git@github.com:marker68/simple-cluster.git simple-cluster
$ cd ./simple-cluster
$ cmake -H. -Bbuild && cmake --build build -- -j3
```
This script will create the following files in your `bin/` directory.

```bash
$ ls  bin
CMakeFiles  libgtest_main.a      test_kdtree  test_utilities
libgtest.a  libsimplecluster.so  test_kmeans
```
* `libgtest_main.a,libgtest.a`: [Google Testing Framework](https://code.google.com/p/googletest/)'s static library files. We are using Google Test for testing.
* `libsimplecluster.{so,dylib,dll}`: dynamic library(shared library) files for SimpleCluster.
* `test_utilities, test_kdtree, test_kmeans`: test programs. See their source code in `test/` directory to know the basics of this library.

## Samples

Samples for the usage of this library could be found at: https://github.com/marker68/example-simplecluster

## Documentation

This project uses Doxygen to generate its documentation. You cand find it in `doc/` or an online version at http://simplecluster.tech-codes.com/
To generate the documentation:
```bash
$ cd ./doc
$ doxygen config.doxygen
```
The documentation files will be generated in HTML and LaTeX format.

## Change log

* **At version 1.0(2014/10/05)**:
    * Supported k-means algorithm.
    * Supported [CMake](http://www.cmake.org/).
    * Supported only L2 metric distance.
    * Supported KD-tree with ANN search.
    * Supported GNU C++ Compiler and clang compiler.
