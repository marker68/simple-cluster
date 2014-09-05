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

/**
 * Random seeding method
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param seeds seeds will be stored here
 * @return this method return nothing
 */
void random_seeds(size_t N, size_t k, vector<double> data, vector<double> seeds) {

}

/**
 * K-means++ seeding method
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param seeds seeds will be stored here
 * @return this method return nothing
 */
void kmeans_pp_seeds(size_t N, size_t k, vector<double> data, vector<double> seeds) {

}

/**
 * The k-means method
 * @param N the number of data
 * @param k the number of clusters
 * @param iters the number of iterations. If provided value is non-positive, it will be set to 10000.
 * @param data the data
 * @param centers the centers after the execution finished.
 * @param clusters the clusters that labeled by the centers' indices.
 * @param seeds seeds will be stored here.
 * @return this method return the final distortion value.
 */
double simple_k_means(KmeansType type, size_t N, size_t k, int iters,
		vector<double> data, vector<double> centers,
		vector<vector<double>> clusters, vector<double> seeds) {
	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds(N,k,data,seeds);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds(N,k,data,seeds);
	}

	// Iterations
	int _iters = 10000;
	if (iters > 0) _iters = iters;

	// Clear all vectors before doing anything
	centers.clear();
	clusters.clear();

	while (_iters > 0) {
		--_iters;
	}
}
}

