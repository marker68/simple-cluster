/*
 * k-means.cpp
 *
 *  Created on: 2014/09/04
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <cstring>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include "k-means.h"
#include "utilities.h"
#include "kd-tree.h"

using namespace std;

namespace SimpleCluster {


/*
 * Global variables
 */
size_t * range;

void random_seeds(size_t d, size_t N, size_t k, double ** data, double ** seeds) {
	size_t i, j;
	int tmp[k];

	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	for(i = 0; i < k; i++)
		tmp[i] = static_cast<int>(i);
	for(i = k; i < N; i++) {
		uniform_int_distribution<> int_dis(i,N-1);
		j = int_dis(gen);
		if(j < k)
			tmp[j] = static_cast<int>(i);
	}

	for(i = 0; i < k; i++)
		seeds[i] = data[tmp[i]];
}

void kmeans_pp_seeds(size_t d, size_t N, size_t k, double ** data, double ** seeds) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	uniform_int_distribution<size_t> int_dis(0, N - 1);
	size_t tmp = int_dis(gen);

	double * d_tmp = data[tmp];
	seeds[0] = d_tmp;
	double * distances = (double *)::operator new(N * sizeof(double));
	double * sum_distances = (double *)::operator new(N * sizeof(double));
	size_t i;
	for(i = 0; i < N; i++) {
		::new(distances + i) double;
		::new(sum_distances + i) double;
		distances[i] = SimpleCluster::distance_square(data[i], d_tmp, d);
		sum_distances[i] = 0.0;
	}

	double sum, tmp2, sum1, sum2, pivot;
	size_t count = 1;
	while(count < k) {
		sum = 0.0;
		for(i = 0; i < N; i++) {
			sum += distances[i];
			sum_distances[i] = sum;
		}
		uniform_real_distribution<double> real_dis(0, sum);
		pivot = real_dis(gen);
		for(i = 0; i < N - 1; i++) {
			sum1 = sum_distances[i];
			sum2 = sum_distances[i + 1];
			if(sum1 < pivot && pivot <= sum2)
				break;
		}
		d_tmp = data[i + 1];
		seeds[count++] = d_tmp;
		// Update the distances
		if(count < k) {
			for(i = 0; i < N; i++) {
				tmp2 = SimpleCluster::distance_square(d_tmp,data[i],d);
				if(distances[i] > tmp2) distances[i] = tmp2;
			}
		}
	}
}

/**
 * Assign the data posize_ts to clusters
 * The execution time would be O(N*k*d)
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @return
 */
void assign_to_closest_centroid(size_t d, size_t N, size_t k,
		double ** data, double ** centroids, int ** clusters) {
	size_t i, j, tmp;
	i_vector i_tmp;
	for(i = 0; i < k; i++) {
		clusters[i] = new int();
	}

	double min = 0.0, temp = 0.0;
	for(i = 0; i < k; i++) {
		::new(range + i) int;
		range[i] = 0;
	}

	for(i = 0; i < N; i++) {
		// Find the minimum distances between d_tmp and a centroid
		min = SimpleCluster::distance(data[i],centroids[0],d);
		tmp = 0;
		for(j = 1; j < k; j++) {
			temp = SimpleCluster::distance(data[i],centroids[j],d);
			if(min < temp) {
				min = temp;
				tmp = j;
			}
		}
		// Assign the data[i] into cluster tmp
		clusters[tmp][range[tmp]++] = static_cast<int>(i);
	}
}

KDNode * convert_data_to_kd_nodes(double ** data, size_t N, size_t d) {
	KDNode * tree = (KDNode *)::operator new(N * sizeof(KDNode));

	size_t i;
	for(i = 0; i < N; i++) {
		::new(tree + i) KDNode;
		tree[i].data = data[i];
		tree[i].dim = d;
	}

	return tree;
}

/**
 * Assign the data posize_ts to clusters
 * The execution time would be O(N*k*d)
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @return
 */
void assign_to_closest_centroid_2(size_t d, size_t N, size_t k,
		double ** data, double ** centroids, int ** clusters) {
	size_t i, tmp;
	KDNode * tree = convert_data_to_kd_nodes(centroids,k,d);
	KDNode * root = make_tree(tree,k,0,d);
	if(root == NULL || tree == NULL) return;
	root->dim = d;
	KDNode node, * found = new KDNode();
	for(i = 0; i < k; i++) {
		clusters[i] = new int();
	}

	double min = DBL_MAX;

	for(i = 0; i < k; i++) {
		::new(range + i) int;
		range[i] = 0;
	}

	for(i = 0; i < N; i++) {
		// Find the minimum distances between d_tmp and a centroid
		node.data = data[i];
		node.dim = d;
		find_nearest(root,&node,*found,min,0);
		tmp = found - tree;
		// Assign the data[i] into cluster tmp
		clusters[tmp][range[tmp]++] = static_cast<int>(i);
	}
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
void simple_k_means(KmeansType type, size_t N, size_t k, KmeansCriteria criteria,size_t d,
		double ** data, double ** centroids,
		int ** clusters, double ** seeds) {
	// Pre-check conditions
	if (N < k) {
		cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	// Seeding
	// In case of defective seeds from users, we should overwrite it by kmeans++ seeds
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds(d,N,k,data,seeds);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds(d,N,k,data,seeds);
	}

	cout << "Finished seeding" << endl;

	// Criteria's setup
	size_t iters = criteria.iterations, i = 0;
	double error = criteria.accuracy, e = error, e_prev = 0.0;

	double ** c_tmp;
	range = (size_t *)::operator new(k * sizeof(size_t));
	if(range == NULL) {
		cerr << "Cannot allocate memory" << endl;
		exit(1);
	}

	double tmp = 0.0;

	// Initialize the centroids
	centroids = seeds;

	while (1) {
		cout << "Loop " << i << endl;
		// Assign the data points to clusters
		assign_to_closest_centroid(d,N,k,data,centroids,clusters);
		// Recalculate the centroids
		for(size_t j = 0; j < k; j++) {
			double * d_tmp = SimpleCluster::mean_vector(data,clusters[j],d,range[j],centroids[j]);
			c_tmp[j] = d_tmp;
		}
		cout << "Loop " << i << endl;
		// Calculate the distortion
		e_prev = e;
		e = 0.0;
		for(size_t j = 0; j < k; j++) {
			tmp = SimpleCluster::distance_square(centroids[j],c_tmp[j],d);
			e += tmp;
		}
		e = sqrt(e);

		centroids = c_tmp;
		i++;
		if(i >= iters ||
//				(e - e_prev < error && e - e_prev > -error)) break;
				(e < error && e > -error)) break;
	}

	cout << "Finished clustering with error is " << e << " after " << i << " iterations." << endl;
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
		double ** data, double ** centroids, int ** clusters) {
	double e = 0.0;
	size_t i, j = 0;
	double * d_tmp;
	int * i_tmp;

	for(i = 0; i < k; i++) {
		i_tmp = clusters[i];
		d_tmp = centroids[i];
		while(j < range[i]) {
			e += SimpleCluster::distance_square(d_tmp, data[i_tmp[j]], d);
			++j;
		}
	}

	return sqrt(e);
}
}

