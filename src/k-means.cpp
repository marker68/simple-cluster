/*
 * k-means.cpp
 *
 *  Created on: 2014/09/04
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include "k-means.h"
#include "utilities.h"

using namespace std;

namespace SimpleCluster {

/**
 * Random seeding method
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param seeds seeds will be stored here
 * @return this method return nothing
 */
void random_seeds(size_t d, size_t N, size_t k, vector<d_vector> data, vector<d_vector> seeds) {

}

/**
 * K-means++ seeding method
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param seeds seeds will be stored here
 * @return this method return nothing
 */
void kmeans_pp_seeds(size_t d, size_t N, size_t k, vector<d_vector> data, vector<d_vector> seeds) {

}

/**
 * Assign the data points to clusters
 * The execution time would be O(N*k*d)
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @param clusters
 */
void assign_to_closest_centroid(size_t d, size_t N, size_t k, vector<d_vector> data,
		vector<d_vector> centroids, vector<i_vector> clusters) {
	int i, j, tmp;
	double min = 0.0, temp = 0.0;
	d_vector d_tmp;

	for(i = 0; i < N; i++) {
		d_tmp = data[i];
		// Find the minimum distances between d_tmp and a centroid
		min = SimpleCluster::distance(d_tmp, centroids[0],d);
		tmp = 0;
		for(j = 1; j < k; j++) {
			temp = SimpleCluster::distance(d_tmp,centroids[j],d);
			if(min > temp) {
				min = temp;
				tmp = j;
			}
		}
		// Assign the data[i] into cluster j
		clusters[j].push_back(i);
	}
}

/**
 * The k-means method
 * @param N the number of data
 * @param k the number of clusters
 * @param criteria the term of accuracy and maximum of iterations.
 * @param d the number of dimensions
 * @param data the data
 * @param centers the centers after the execution finished.
 * @param clusters the clusters that labeled by the centers' indices.
 * @param seeds seeds will be stored here.
 */
void simple_k_means(KmeansType type, size_t N, size_t k, KmeansCriteria criteria,size_t d,
		vector<d_vector> data, vector<d_vector> centroids,
		vector<i_vector> clusters, vector<d_vector> seeds) {
	// Pre-check conditions
	if (N < k) {
		cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds(d,N,k,data,seeds);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds(d,N,k,data,seeds);
	}

	// Criteria's setup
	int iters = criteria.iterations, i = 0;
	double error = criteria.accuracy, e = 0.0;

	// Clear all vectors before doing anything
	centroids.clear();
	clusters.clear();
	SimpleCluster::allocate(clusters, k);
	vector<d_vector> c_tmp;
	double tmp = 0.0;

	// Initialize the centers
	centroids = seeds;

	while (i < iters && e > error) {
		// Clear the last clusters contents
		for(int j = 0; j < k; j++)
			clusters[j].clear();
		// Assign the data points to clusters
		assign_to_closest_centroid(d,N,k,data,centroids,clusters);

		// Recalculate the centroids
		c_tmp.clear();
		for(int j = 0; j < k; j++) {
			c_tmp.push_back(SimpleCluster::mean_vector(data,clusters[j],d));
		}

		// Calculate the distortion
		e = 0.0;
		for(int j = 0; j < k; j++) {
			tmp = SimpleCluster::distance_square(centroids[j],c_tmp[j],d);
			e += tmp;
		}
		e = sqrt(e);

		centroids = c_tmp;
		i++;
	}
}

/**
 * Calculate the distortion of a set of clusters
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @param clusters
 */
double distortion(size_t d, size_t N, size_t k,
		vector<d_vector> data, vector<d_vector> centroids, vector<i_vector> clusters) {
	double tmp = 0.0, e = 0.0;
	int i, j;
	d_vector d_tmp;
	i_vector i_tmp;

	for(i = 0; i < k; i++) {
		i_tmp = clusters[i];
		d_tmp = centroids[i];
		i_vector::iterator it = i_tmp.begin();
		i_vector::iterator ie = i_tmp.end();
		while(it != ie) {
			e += SimpleCluster::distance_square(d_tmp, data[(*it)], d);
			++it;
		}
		d_tmp.clear();
		i_tmp.clear();
	}

	return sqrt(e);
}
}

