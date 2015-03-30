/*
 *  Copyright (C) 2015 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 *  mat_utilities.h
 *
 *  Created on: 2015/03/29
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef INCLUDE_MAT_UTILITIES_H_
#define INCLUDE_MAT_UTILITIES_H_

#include <iostream>
#include <vector>
#include <exception>
#include <cmath>
#include <cstring>
#include <cstdio>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

#include <cblas.h>
#include <lapacke.h>

using namespace std;

namespace SimpleCluster {
/**
 * Find the mean vector of all values in each row or column of
 * a matrix A.
 */
template<typename DataType>
void mean(
		DataType * A,
		int m,
		int n,
		int type,
		DataType *& B,
		bool verbose) {
	if(m <= 0 || n <= 0) return;
	int size = 0;
	int i, j;
	if(type == 1)
		size = n;
	else if(type == 2)
		size = m;
	else size = m * n;

	B = (DataType *)::operator new(size * sizeof(DataType));

	if(type < 1 || type > 2) {
		memcpy(B,A,size * sizeof(DataType));
		return;
	}

	if(type == 2) {
		// TODO Use OpenMP?
		for(i = 0; i < m; i++) {
			B[i] = 0;
			for(j = 0; j < n; j++) {
				B[i] += A[0];
				A++;
			}
			B[i] /= n;
		}
		return;
	}

	if(type == 1) {
		for(i = 0; i < n; i++) B[i] = 0;
		// TODO Use OpenMP?
		for(i = 0; i < m; i++) {
			for(j = 0; j < n; j++) {
				B[j] += A[0];
				A++;
			}
		}
		for(i = 0; i < n; i++) B[i] /= m;
		return;
	}
}

template<typename DataType>
void transpose(
		DataType * a,
		DataType *& b,
		int m,
		int n,
		bool verbose) {
	int i, j;
	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			b[j * m + i] = a[i * n + j];
		}
	}
}

void mabs(
		float * a,
		float *& b,
		int size,
		int nthread,
		bool verbose);
void sign(
		float * a,
		float *& b,
		int size,
		int nthread,
		bool verbose);
}



#endif /* INCLUDE_MAT_UTILITIES_H_ */
