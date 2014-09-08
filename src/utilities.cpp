/*
 * utilities.cpp
 *
 *  Created on: 2014/09/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
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
double distance(d_vector x, d_vector y, size_t d) {
	if(x.size() < d || y.size() < d) {
		cerr << "distance: Your vector have not enough dimensions!" << endl;
		cerr << "x has " << x.size() << " dimensions!" << endl;
		cerr << "y has " << y.size() << " dimensions!" << endl;
		exit(1);
	}
	size_t i;
	double dis = 0.0, tmp = 0.0;
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
double distance_square(d_vector x, d_vector y, size_t d) {
	if(x.size() < d || y.size() < d) {
		cerr << "distance_square: Your vector have not enough dimensions!" << endl;
		cerr << "x has " << x.size() << " dimensions!" << endl;
		cerr << "y has " << y.size() << " dimensions!" << endl;
		exit(1);
	}
	size_t i;
	double dis = 0.0, tmp = 0.0;
	for(i = 0; i < d; i++) {
		tmp = x[i] - y[i];
		dis += tmp * tmp;
	}

	return dis;
}

/**
 * Calculate the mean of a cluster
 * @param data
 * @param index
 * @param d
 * @return the mean posize_t of a cluster
 */
d_vector mean_vector(vector<d_vector> data, i_vector index, size_t d, d_vector centroid) {
	size_t i, j, size  = index.size();
	if(size <= 0) {
		return centroid;
	}
	d_vector d_tmp;
	double tmp[d];

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
	d_tmp.clear();
	for(i = 0; i < d; i++)
		d_tmp.push_back(tmp[i]);
	return d_tmp;
}

/**
 * Get system time in milliseconds
 */
unsigned long get_millisecond_time() {
	struct timeval tv;
	if(gettimeofday(&tv, NULL) != 0) return 0;
	return (unsigned long)((tv.tv_sec * 1000ul) + (tv.tv_usec / 1000ul));
}

/**
 * Utilities for printing vector
 */
void print_vector(vector<d_vector> data, size_t d) {
	size_t i, j, size = data.size();
	for(i = 0; i < size; i++) {
		cout << "Data point " << i << ":";
		for(j = 0; j < d; j++)
			cout << data[i][j] << " ";
		cout << endl;
	}
}
}


