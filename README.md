Simple Clustering
===============================

## What is this?

This is a simple implementation of the state-of-the-art clustering methods such as k-means, EM algorithm, ... 
The implementations based on many papers that are mentioned below. This project is under GNU GPL v3 License. Just see the license statement.

## References

[1] S. P. Lloyd, "Least squares quantization in PCM",  IEEE Trans. Inform. Theory,  vol. IT-28,  no. 2, pp. 129 -137, Mar. 1982
   

[2] D. Arthur et al., "k-means++: the advantages of careful seeding",  SODA '07 Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete algorithms, pp. 1027-1035, 2007

## Usage

### Prerequisites

* [CMake](http://www.cmake.org/) 2.8 or newer. For UNIX users, check your CMake version in terminal by `cmake -version`.
* An C++ compiler that supports C++ 11.

### Build
* To build this project, just type:
```bash
$ cmake -H. -B<where_you_want_to_save_CMake_build_files>
```
on your terminal. This will generate CMake files under `cmake/`. Then type:
```bash
$ cmake --build <the_location_of_CMake_build_files> -- -j3
```
to generate the binaries.

* The above step would generate in `bin/` folder the following files: `libsimple_cluster.*`,`test_kmeans`. 
* The `test_kmeans` is an unit test program.

## Documentation


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