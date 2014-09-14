Simple Clusters
===============================

## What is this?

This is a simple implementation of the state-of-the-art clustering methods such as k-means, EM algorithm, ... 
The implementations based on many papers that are mentioned below. This project is under GNU GPL v3 License. Just read the license statement.

The purpose of this project is to provide a more flexible interface for clustering methods. In OpenCV, the implementations of k-means algorithm is limited. Thus, I hope to this one will prodive an implementation with more choices for k-means. Sometimes, reducing the distortion is more important than improving the runtime performance. Other clustering methods will also be considered and be implemented. The input of many methods now are vectors, I hope to improve this issue soon, i.e let you input files.

## Usage

### Prerequisites

* An Unix or Windows Operating System. Tested on Mac OS X 10.9 and Ubuntu 14.04.
* [CMake](http://www.cmake.org/) 2.8 or newer. For UNIX users, check your CMake version in terminal by `cmake -version`.
* An C++ compiler that supports C++ 11.
* [OpenCV](http://opencv.org/downloads.html) 2.4.8 or newer. This library is used for testing only. Just for comparing the performance with our implementations.

### Build
* To build this project, just type:
```bash
$ cmake -H. -B<where_you_want_to_save_CMake_build_files>
### Example
$ cmake -H. -Bbuild
```
on your terminal. This will generate CMake files under `build/` and binaries under `bin/`. Then type:
```bash
$ cmake --build <the_location_of_CMake_build_files> -- -j3
```
to generate the binaries.

* The above step would generate in `bin/` folder the following files: `libsimplecluster.*`,`test_*`. 
* The `test_*` are unit test programs. You can review the source of unit tests in `test/`. They are also self-explained documentation for the common usage.
* An example to show you how to use this library: https://github.com/marker68/example-simplecluster

## Contributing
I welcome any contributions. First, read some materials like [How to contribute to an open source project?](https://guides.github.com/activities/contributing-to-open-source/) or [an excellent example from FuelPHP](https://github.com/fuelphp/fuelphp/blob/master/CONTRIBUTING.md). I will work on this term later.

## Documentation

This project uses Doxygen to generate its documentation. You cand find it in `doc/` or an online version at http://simplecluster.tech-codes.com/

## References

[1] S. P. Lloyd, "Least squares quantization in PCM",  IEEE Trans. Inform. Theory,  vol. IT-28,  no. 2, pp. 129 -137, Mar. 1982
   
[2] D. Arthur et al., "k-means++: the advantages of careful seeding",  SODA '07 Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms, pp. 1027-1035, 2007

[3] D.T. Lee et al., "Worst-case analysis for region and partial region searches in multidimensional binary search trees and balanced quad trees", Acta Informatica, vol. 9, issue 1, pp. 23-29, 1977

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
