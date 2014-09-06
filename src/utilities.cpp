/*
 * utilities.cpp
 *
 *  Created on: 2014/09/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <math.h>
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
long long double distance(d_vector x, d_vector y, int d) {
	int i;
	long long double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += tmp * tmp;
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
long long double distance_square(d_vector x, d_vector y, int d) {
	int i;
	long long double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += tmp * tmp;
	}

	return dis;
}

void allocate(vector<i_vector> ivec, int size) {
	int i;
	i_vector tmp;
	for(i = 0; i < size; i++) {
		ivec.push_back(tmp);
	}
}

/**
 * Calculate the mean of a cluster
 * @param data
 * @param index
 * @param d
 * @return the mean point of a cluster
 */
d_vector mean_vector(vector<d_vector> data, i_vector index, int d) {
	int i, j, size  = index.size();
	d_vector tmp, d_tmp;
	tmp.reserve(d);

	for(i = 0; i < d; i++)
		tmp[i] = 0.0;

	for(i = 0; i < size; i++) {
		j = index[i];
		d_tmp = data[j];
		for(j = 0; j < d; j++) {
			tmp[j] += d_tmp[j];
		}
	}

	for(i = 0; i < d; i++)
		tmp[i] /= static_cast<double>(size);

	return tmp;
}
}


