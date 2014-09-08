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
double distance(double * x, double * y, size_t d) {
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
double distance_square(double * x, double * y, size_t d) {
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
double * mean_vector(vector<double *> data, const int * index, size_t d, size_t size, double * centroid) {
	size_t i, j;
	if(size <= 0) {
		return centroid;
	}
	double * d_tmp = (double *)malloc(d * sizeof(double));
	double * tmp = (double *)malloc(d * sizeof(double));
	if(d_tmp == NULL || tmp == NULL) {
		cerr << "Cannot allocate memory" << endl;
		delete d_tmp;
		delete tmp;
		exit(1);
	}

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
	for(i = 0; i < d; i++)
		d_tmp[i] = tmp[i];
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
void print_vector(vector<double *> data, size_t d) {
	size_t i, j, size = data.size();
	for(i = 0; i < size; i++) {
		cout << "Data point " << i << ":";
		for(j = 0; j < d; j++)
			cout << data[i][j] << " ";
		cout << endl;
	}
}
}


