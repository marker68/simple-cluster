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
double distance(d_vector, d_vector, size_t);
double distance_square(d_vector, d_vector, size_t);
void allocate(vector<i_vector>, size_t);
d_vector mean_vector(vector<d_vector>, i_vector, size_t);
}



#endif /* UTILITIES_H_ */
