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
 *  utilities.cpp
 *
 *  Created on: 2014/09/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <exception>
#include <math.h>
#include <sys/time.h>
#include "utilities.h"

using namespace std;

namespace SimpleCluster {
/**
 * Calculate the distance between two vectors
 * @param x
 * @param y
 * @param d
 * @return the distance between x and y in d dimensional space
 */
double distance(
		double * x,
		double * y,
		size_t d) {
	size_t i;
	double dis = 0.0, tmp = 0.0;
	try {
		for(i = 0; i < d; i++) {
			tmp = x[i] - y[i];
			dis += tmp * tmp;
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		exit(1);
	}

	return sqrt(dis);
}

/**
 * Calculate the square of distance between two vectors
 * @param x
 * @param y
 * @param d
 * @return the square of distance between x and y in d dimensional space
 */
double distance_square(
		double * x,
		double * y,
		size_t d) {
	size_t i;
	double dis = 0.0, tmp = 0.0;
	try {
		for(i = 0; i < d; i++) {
			tmp = x[i] - y[i];
			dis += tmp * tmp;
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		exit(1);
	}

	return dis;
}

/**
 * Calculate all the mean vectors
 * @param data
 * @param label
 * @param size
 * @param d
 * @param N
 * @param k
 * @param centroids
 * @return nothing
 */
void all_mean_vector(
		double ** data,
		int * label,
		size_t * size,
		double **& centroids,
		double *& moved,
		size_t d,
		size_t N,
		size_t k) {
	size_t i, j, t;
	double * tmp;
	double ** c_tmp;
	init_array_2<double>(c_tmp,k,d);
	for(i = 0; i < k; i++) {
		for(j = 0; j < d; j++) {
			c_tmp[i][j] = 0.0;
		}
	}
	for(i = 0; i < N; i++) {
		t = label[i];
		tmp = data[t];
		for(j = 0; j < d; j++) {
			c_tmp[t][j] += tmp[j];
		}
	}
	for(i = 0; i < k; i++) {
		t = size[i];
		if(t > 0) {
			for(j = 0; j < d; j++) {
				c_tmp[i][j] /= t;
			}
			moved[i] = distance_square(c_tmp[i],centroids[i],d);
			copy_array<double>(c_tmp[i],centroids[i],d);
		} else moved[i] = 0.0;
	}
}

/**
 * Calculate the mean of a cluster
 * @param data
 * @param index
 * @param d
 * @param size
 * @param centroid
 * @return the mean posize_t of a cluster
 */
double * mean_vector(
		double ** data,
		const int * index,
		double * centroid,
		size_t d,
		size_t size) {
	size_t i, j = 0, k;
	if(size <= 0) {
		return centroid;
	}

	double * tmp = (double *)calloc(d, sizeof(double));
	if(tmp == nullptr) {
		cerr << "Cannot allocate memory" << endl;
		exit(1);
	}

	try {
		for(i = 0; i < size; i++) {
			j = index[i];
			for(k = 0; k < d; k++) {
				tmp[k] += data[j][k];
			}
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		exit(1);
	}

	try {
		for(i = 0; i < d; i++) {
			tmp[i] /= static_cast<double>(size);
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		exit(1);
	}

	return tmp;
}

/**
 * Calculate the mean of a cluster
 * @param data
 * @param index a vector of integers
 * @param d
 * @return the mean posize_t of a cluster
 */
double * mean_vector(
		double ** data,
		const i_vector index,
		double * centroid,
		size_t d) {
	size_t i, j = 0, k;
	if(index.size() <= 0) {
		return centroid;
	}
	size_t size = index.size();

	double * tmp = (double *)calloc(d,sizeof(double));
	if(tmp == nullptr) {
		cerr << "Cannot allocate memory" << endl;
		exit(1);
	}

	try {
		for(i = 0; i < size; i++) {
			j = index[i];
			for(k = 0; k < d; k++) {
				tmp[k] += data[j][k];
			}
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		exit(1);
	}

	try {
		for(i = 0; i < d; i++) {
			tmp[i] /= static_cast<double>(size);
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		exit(1);
	}

	return tmp;
}

/**
 * Get system time in milliseconds
 */
unsigned long get_millisecond_time() {
	struct timeval tv;
	if(gettimeofday(&tv, nullptr) != 0) return 0;
	return (unsigned long)((tv.tv_sec * 1000ul) + (tv.tv_usec / 1000ul));
}

/**
 * Utilities for printing vector
 * @param data the input data
 * @param d the number of dimensions
 * @param N the size of the input
 */
void print_vector(
		double ** data,
		size_t d,
		size_t N) {
	size_t i, j;
	try {
		for(i = 0; i < N; i++) {
			cout << "Data point " << i << ":";
			for(j = 0; j < d; j++)
				cout << data[i][j] << " ";
			cout << endl;
		}
	} catch(exception& e) {
		cerr << "Got an exception: " << e.what() << endl;
		return;
	}
}
}


