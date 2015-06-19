#!/bin/sh

# Check the number of arguments
if [ "$#" -ne 1 ] || ! [ -d "$1" ]; then
	echo "Usage: $0 DIRECTORY" >&2
	echo "Note: DIRECTORY should be an absolute path."
	echo "The current directory is '`pwd`'."
	exit 1
fi

HOME=$1

# Clean the top directory
rm -rf ${HOME}/bin ${HOME}/build
rm -f ${HOME}/lib/*.so
rm -f ${HOME}/lib/*.dll
rm -f ${HOME}/lib/*.dylib
rm -f ${HOME}/lib/*.lib
rm -f ${HOME}/lib/*.a

rm -f ${HOME}/include/lapacke*
rm -f ${HOME}/include/f77blas.h
rm -f ${HOME}/include/cblas.h
rm -f ${HOME}/include/openblas_config.h

# Clean the submodules' directory
cd ${HOME}/lib/OpenBLAS && make clean
