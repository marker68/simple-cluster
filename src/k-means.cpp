/*
 * k-means.cpp
 *
 *  Created on: 2014/09/04
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include "k-means.h"

using namespace std;

namespace cluster {
double simple-k-means(KmeansType type, size_t N, size_t k,int iterations,
		vector<double> data, vector<double> centroids,
		vector<vector<double>> clusters, vector<double> seeds) {
	if(type == KmeansType::RANDOM_SEEDS) {

	}
}
}


