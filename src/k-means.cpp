/*
 * k-means.cpp
 *
 *  Created on: 2014/09/04
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <random>
#include <stdlib.h>
#include <math.h>
#include "k-means.h"
#include "utilities.h"

using namespace std;

namespace SimpleCluster {

/**
 * Random seeding method.
 * Just pick up randomly k distinctive posize_ts from the input data.
 * We use the Reservoir sampling algorithm for this implementation.
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @return this method return nothing
 */
vector<d_vector> random_seeds(size_t d, size_t N, size_t k, vector<d_vector> data) {
	vector<d_vector> seeds;
	size_t i, j;
	int tmp[k];

	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	size_t size = static_cast<size_t>(data.size());
	for(i = 0; i < k; i++)
		tmp[i] = static_cast<int>(i);
	for(i = k; i < size; i++) {
		uniform_int_distribution<> int_dis(i,size-1);
		j = int_dis(gen);
		if(j < k)
			tmp[j] = static_cast<int>(i);
	}

	for(i = 0; i < k; i++)
		seeds.push_back(data[i]);

	return seeds;
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
vector<d_vector> kmeans_pp_seeds(size_t d, size_t N, size_t k, vector<d_vector> data) {
	vector<d_vector> seeds;
	size_t size = static_cast<size_t>(data.size());

	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	uniform_int_distribution<> int_dis(0,size-1);
	size_t tmp = int_dis(gen);

	d_vector d_tmp = data[tmp];
	seeds.push_back(d_tmp);
	d_vector distances;
	distances.reserve(size);
	size_t i, j;
	for(i = 0; i < size; i++) {
		distances[i] = SimpleCluster::distance_square(data[i],d_tmp,d);
	}

	while(seeds.size() < k) {
		double sum = 0.0;
		for(i = 0; i < size; i++)
			sum += distances[i];
		uniform_real_distribution<> real_dis(0, sum);
		double pivot = real_dis(gen);
		sum = 0.0;
		for(i = 0; i < size - 1; i++) {
			sum += distances[i];
			if(sum < pivot && pivot <= sum + distances[i+1])
				break;
		}
		seeds.push_back(data[i+1]);
		// Update the distances
		for(i = 0; i < size; i++) {
			d_tmp = data[i];
			double d_min = SimpleCluster::distance_square(d_tmp,seeds[0],d);
			double tmp2;
			for(j = 1; j < seeds.size(); j++) {
				tmp2 = SimpleCluster::distance_square(d_tmp,seeds[j],d);
				if(tmp2 < d_min) {
					d_min = tmp2;
				}
			}
			distances[i] = d_min;
			d_tmp.clear();
		}
	}

	return seeds;
}

/**
 * Assign the data posize_ts to clusters
 * The execution time would be O(N*k*d)
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @return the clusters
 */
vector<i_vector> assign_to_closest_centroid(size_t d, size_t N, size_t k,
		vector<d_vector> data, vector<d_vector> centroids) {
	vector<i_vector> clusters;
	size_t i, j, tmp;
	i_vector i_tmp;
	for(i = 0; i < k; i++) {
		clusters.push_back(i_tmp);
		i_tmp.clear();
	}

	double min = 0.0, temp = 0.0;
	d_vector d_tmp;

	for(i = 0; i < N; i++) {
		d_tmp = data[i];
		// Find the minimum distances between d_tmp and a centroid
		min = SimpleCluster::distance(d_tmp, centroids[0],d);
		tmp = 0;
		for(j = 1; j < k; j++) {
			temp = SimpleCluster::distance(d_tmp,centroids[j],d);
			if(min < temp) {
				min = temp;
				tmp = j;
			}
		}
		// Assign the data[i] into cluster tmp
		clusters[tmp].push_back(static_cast<int>(i));
	}

	return clusters;
}

/**
 * The k-means method: a description of the method can be found at
 * http://home.deib.polimi.it/matteucc/Clustering/tutorial_html/kmeans.html
 * @param N the number of data
 * @param k the number of clusters
 * @param criteria the term of accuracy and maximum of iterations.
 * @param d the number of dimensions
 * @param data the data
 * @param centers the centers after the execution finished.
 * @param clusters the clusters that labeled by the centers' indices.
 * @param seeds seeds will be stored here.
 */
void simple_k_means(KmeansType type, size_t N, size_t k, KmeansCriteria criteria, size_t d,
		vector<d_vector> data, vector<d_vector>& centroids,
		vector<i_vector>& clusters, vector<d_vector> seeds) {
	// Pre-check conditions
	if (N < k) {
		cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		seeds.clear();
		seeds = random_seeds(d,N,k,data);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		seeds.clear();
		seeds = kmeans_pp_seeds(d,N,k,data);
	}

	if(seeds.size() < k) {
		seeds.clear();
		seeds = kmeans_pp_seeds(d,N,k,data);
	}

	// Criteria's setup
	size_t iters = criteria.iterations, i = 0;
	double error = criteria.accuracy, e = error + 1.0;

	vector<d_vector> c_tmp;
	d_vector d_tmp;
	double tmp = 0.0;

	// Initialize the centers
	centroids = seeds;

	while (i < iters && e > error) {
		// Clear the last clusters contents
		clusters.clear();
		// Assign the data posize_ts to clusters
		clusters = assign_to_closest_centroid(d,N,k,data,centroids);
		// Recalculate the centroids
		c_tmp.clear();
		for(size_t j = 0; j < k; j++) {
			d_tmp = SimpleCluster::mean_vector(data,clusters[j],d,centroids[j]);
			c_tmp.push_back(d_tmp);
		}
		// Calculate the distortion
		e = 0.0;
		for(size_t j = 0; j < k; j++) {
			tmp = SimpleCluster::distance_square(centroids[j],c_tmp[j],d);
			e += tmp;
		}
		e = sqrt(e);

		centroids = c_tmp;
		i++;
	}

	cout << "Finished clustering with error is " << e << endl;
}

/**
 * Calculate the distortion of a set of clusters.
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @param clusters
 */
double distortion(size_t d, size_t N, size_t k,
		vector<d_vector> data, vector<d_vector> centroids, vector<i_vector> clusters) {
	if(centroids.size() < k || clusters.size() < k) {
		cerr << "You don't have enough clusters!" << endl;
		cerr << "You have " << centroids.size() << " centroids" << endl;
		cerr << "You have " << clusters.size() << " clusters" << endl;
		exit(1);
	}
	double e = 0.0;
	size_t i;
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

