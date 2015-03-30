/*
 *  SIMPLE CLUSTERS: A simple library for clustering works.
 *  Copyright (C) 2014 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
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
 *  mat_utilities.cpp
 *
 *  Created on: 2015/03/30
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <cstdlib>
#include <vector>
#include <exception>
#include <cmath>
#include <ctime>

#ifdef _WIN32
#include <windows.h>
#else
#include <sys/time.h>
#endif

#include "mat_utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace std;

namespace SimpleCluster {
/**
 * Get sign of each elements in a
 */
void sign(
		float * a,
		float *& b,
		int size,
		int nthread,
		bool verbose) {
	if(size <= 0) return;
	int i, i0;
	int blk = size / nthread;
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#pragma omp parallel
	{
#pragma omp for private(i, i0)
#endif
		for(i0 = 0; i0 < nthread; i0++) {
			int start = i0 * blk;
			int end = start + blk;
			if(end > size) end = size;
			float * tmp = a + start;
			float * tmp2 = b + start;
			for(i = start; i < end; i++) {
				if(*tmp > 0.0f) *(tmp2++) = 1.0f;
				else if(*tmp < 0.0f) *(tmp2++) = -1.0f;
				else *(tmp2++) = 0.0f;
			}
		}
#ifdef _OPENMP
	}
#endif
}

void abs(
		float * a,
		float *& b,
		int size,
		int nthread,
		bool verbose) {
	if(size <= 0) return;
	int i, i0;
	int blk = size / nthread;
#ifdef _OPENMP
	omp_set_num_threads(nthread);
#pragma omp parallel
	{
#pragma omp for private(i, i0)
#endif
		for(i0 = 0; i0 < nthread; i0++) {
			int start = i0 * blk;
			int end = start + blk;
			if(end > size) end = size;
			float * tmp = a + start;
			float * tmp2 = b + start;
			for(i = start; i < end; i++) {
				*(tmp2++) = fabs(*(tmp++));
			}
		}
#ifdef _OPENMP
	}
#endif
}
}


