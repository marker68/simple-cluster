Simple Clusters
===============================

## What is this?

This is a simple implementation of the state-of-the-art clustering methods such as k-means, EM algorithm, ... 
The implementations based on many papers that are mentioned below. This project is under GNU GPL v3 License. Just read the license statement.
If you need any other licenses, feel free to contact me at [t_nguyen@hal.t.u-tokyo.ac.jp](mailto:t_nguyen@hal.t.u-tokyo.ac.jp)
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
$ git clone git@github.com:marker68/simple-k-means.git simple-cluster
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

## Change log

* **At version 1.0(2014/10/05)**:
    * Supported k-means algorithm.
    * Supported [CMake](http://www.cmake.org/).
    * Supported only L2 metric distance.
    * Supported KD-tree with ANN search.
    * Supported GNU C++ Compiler and clang compiler.
    
## References

[1] S. P. Lloyd, "Least squares quantization in PCM,"  IEEE Trans. Inform. Theory,  vol. IT-28,  no. 2, pp. 129 -137, Mar. 1982.
   
[2] D. Arthur et al., "k-means++: the advantages of careful seeding,"  SODA '07 Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms, pp. 1027-1035, 2007.

[3] D.T. Lee et al., "Worst-case analysis for region and partial region searches in multidimensional binary search trees and balanced quad trees," Acta Informatica, vol. 9, issue 1, pp. 23-29, 1977.

[4] G. Hamerly., "Making k-means even faster," Proc. SDM, pp. 130-140, 2010.
## License
```
    SIMPLE CLUSTERS: A simple library for clustering works.
    Copyright (C) 2014 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
```
