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
 * k-means.cpp
 *
 *  Created on: 2014/09/04
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <cmath>
#include "k-means.h"

#ifdef _OPENMP
#include <omp.h>
#define SET_THREAD_NUM omp_set_num_threads(n_thread)
#else
#define SET_THREAD_NUM 0 // disable multi-thread
#endif

namespace SimpleCluster {
/**
 * Update the bounds
 * @param moved the distances that centers moved
 * @param label the labels of posize_t data
 * @param upper
 * @param lower
 * @param N
 * @param k
 * @param d
 * @param n_thread the number of threads
 */
void update_bounds(
		double * moved,
		size_t * label,
		double *& upper,
		double *& lower,
		size_t N,
		size_t k,
		int n_thread) {
	size_t r = 0;
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	double max = 0.0, max2 = 0.0, sub;
	for(i = 0; i < k; i++) {
		if(max <= moved[i]) {
			max2 = max;
			max = moved[i];
			r = i;
		} else {
			if(max2 <= moved[i]) {
				max2 = moved[i];
			}
		}
	}
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < N; i++) {
			upper[i] += moved[label[i]];
			sub = 0.0;
			if(r == label[i]) {
				sub = max2;
			} else {
				sub = max;
			}
			lower[i] = fabs(lower[i] - sub);
		}
	}
}
}
