/*
 * utilities.h
 *
 *  Created on: 2014/09/05
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef UTILITIES_H_
#define UTILITIES_H_

#include <iostream>
#include <vector>

using namespace std;

namespace SimpleCluster {
typedef vector<double> d_vector;
typedef vector<int> i_vector;

/**
 * Calculate the distance between two vectors
 */
double distance(double *, double *, size_t);
double distance_square(double *, double *, size_t);
double * mean_vector(vector<double *>, const int *, size_t, size_t, double *);
/**
 * Get system time in milliseconds
 */
unsigned long get_millisecond_time();
void print_vector(vector<double *>, size_t);
}

#endif /* UTILITIES_H_ */
