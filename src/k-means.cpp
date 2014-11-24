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
#endif

namespace SimpleCluster {

/**
 * Update the centers
 * @param sum the sum vector of all points in the cluster
 * @param size the size of each cluster
 * @param centers the centers of clusters
 * @param moved the distances that centers moved
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param k the number of clusters
 * @param d the number of dimensions
 * @param n_thread the number of threads
 * @return nothing
 */
void update_center(
		float ** sum,
		int * size,
		float **& centers,
		float *& moved,
		DistanceType d_type,
		int k,
		int d,
		int n_thread) {
	float * c_tmp;
	init_array<float>(c_tmp,d);
	int i;
	for(i = 0; i < k; i++) {
		copy_array<float>(centers[i],c_tmp,d);
		for(int j = 0; j < d; j++) {
			centers[i][j] = static_cast<float>(sum[i][j] / size[i]);
		}
		if(d_type == DistanceType::NORM_L2)
			moved[i] = distance_l2<float>(c_tmp,centers[i],d);
		else if(d_type == DistanceType::NORM_L1)
			moved[i] = distance_l1<float>(c_tmp,centers[i],d);
	}
}

/**
 * Update the bounds
 * @param moved the distances that centers moved
 * @param label the labels of point data
 * @param upper
 * @param lower
 * @param N
 * @param k
 * @param d
 * @param n_thread the number of threads
 */
void update_bounds(
		float * moved,
		int * label,
		float *& upper,
		float *& lower,
		int N,
		int k,
		int n_thread) {
	int r = 0;
	int i;

	float max = 0.0, max2 = 0.0, sub;
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
#ifdef _OPENMP
	omp_set_num_threads(n_thread);
#pragma omp parallel
	{
#pragma omp for private(i)
#endif
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
#ifdef _OPENMP
	}
#endif
}
}
