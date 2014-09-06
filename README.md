Simple Clustering
===============================

## What is this?

This is a simple implementation of the state-of-the-art clustering methods such as k-means, EM algorithm, ... 
The implementations based on many papers that are mentioned below. You can use the source code in your project just for private uses as it will be stated in the license statement.

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
The MIT License (MIT)

Copyright (c) 2014, Nguyen Anh Tuan<t_nguyen@hal.t.u-tokyo.ac.jp>

Permission is hereby granted, free of charge, to any person obtaining a copy of 
this software and associated documentation files (the "Software"), to deal in 
the Software without restriction, including without limitation the rights to 
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of 
the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all 
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
FOR ANY CLAIM,DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT 
OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
OTHER DEALINGS IN THE SOFTWARE.
```