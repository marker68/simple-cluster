/*
 *  SIMPLE CLUSTERS: A simple library for clustering works.
 *  Copyright (C) 2014 Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
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

/**
 * Create random seeds for k-means
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param seeds the seeds
 * @param verbose for debugging
 */
void random_seeds(size_t d, size_t N, size_t k,
		double ** data, double **& seeds, bool verbose) {
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

	for(i = 0; i < k; i++) {
		if(!copy_array<double>(data[tmp[i]],seeds[i],d)) {
			if(verbose)
				cerr << "Cannot copy array! The program exit now!" << endl;
			exit(1);
		}
	}
}

/**
 * Create seeds for k-means++
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param seeds the seeds
 * @param verbose for debugging
 */
void kmeans_pp_seeds(size_t d, size_t N, size_t k,
		double ** data, double **& seeds, bool verbose) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	uniform_int_distribution<size_t> int_dis(0, N - 1);
	size_t tmp = int_dis(gen);

	double * d_tmp;
	init_array<double>(d_tmp,d);
	copy_array<double>(data[tmp],seeds[0],d);
	double * distances;
	init_array<double>(distances,N);
	double * sum_distances;
	init_array<double>(sum_distances,N);
	size_t i;
	for(i = 0; i < N; i++) {
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
		copy_array<double>(data[i+1],d_tmp,d);
		copy_array<double>(d_tmp,seeds[count++],d);
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
 * After having a set of centroids,
 * we need to assign data into each cluster respectively.
 * This solution uses linear search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param verbose for debugging
 */
void assign_to_closest_centroid(size_t d, size_t N, size_t k,
		double ** data, double ** centroids, vector<i_vector>& clusters, bool verbose) {
	size_t i, tmp;

	double min = 0.0;

	for(i = 0; i < N; i++) {
		// Find the minimum distances between d_tmp and a centroid
		linear_search(centroids,data[i],tmp,min,k,d,verbose);
		// Assign the data[i] into cluster tmp
		clusters[tmp].push_back(static_cast<int>(i));
	}
}

/**
 * After having a set of centroids,
 * we need to assign data into each cluster respectively.
 * This solution uses kd-tree search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param verbose for debugging
 */
void assign_to_closest_centroid_2(size_t d, size_t N, size_t k,
		double ** data, double ** centroids, vector<i_vector>& clusters, bool verbose) {
	size_t i, tmp;
	KDNode<double> * root = NULL;
	make_balanced_tree(root,centroids,k,d,0,0,verbose);
	if(root == NULL) return;

	KDNode<double> query;

	for(i = 0; i < N; i++) {
		KDNode<double> * nn = NULL;
		double min = DBL_MAX;
		query.add_data(data[i],d);
		// Find the minimum distances between d_tmp and a centroid
		nn_search(root,&query,nn,min,d,0,verbose);
		tmp = nn->id;
		// Assign the data[i] into cluster tmp
		clusters[tmp].push_back(static_cast<int>(i));
		query.clear_data();
	}
}

/**
 * After having a set of centroids,
 * we need to assign data into each cluster respectively.
 * This solution uses kd-tree search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param verbose for debugging
 */
void assign_to_closest_centroid_3(size_t d, size_t N, size_t k,
		double ** data, double ** centroids, vector<i_vector>& clusters, double alpha, bool verbose) {
	size_t i, tmp;
	KDNode<double> * root = NULL;
	make_random_tree(root,centroids,k,d,0,0,verbose);
	if(root == NULL) return;

	KDNode<double> query;

	for(i = 0; i < N; i++) {
		KDNode<double> * nn = NULL;
		double min = DBL_MAX;
		query.add_data(data[i],d);
		// Find the minimum distances between d_tmp and a centroid
		ann_search(root,&query,nn,min,alpha,d,0,verbose);
		tmp = nn->id;
		// Assign the data[i] into cluster tmp
		clusters[tmp].push_back(static_cast<int>(i));
		query.clear_data();
	}
}

/**
 * The k-means method: a description of the method can be found at
 * http://home.deib.polimi.it/matteucc/Clustering/tutorial_html/kmeans.html
 * @param type the type of seeding method
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param criteria the criteria
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param seeds the initial centroids = the seeds
 * @param verbose for debugging
 */
void simple_k_means(KmeansType type, KmeansAssignType assign,
		size_t N, size_t k, KmeansCriteria criteria,size_t d,
		double ** data, double **& centroids,
		vector<i_vector>& clusters, double **& seeds, bool verbose) {
	// Pre-check conditions
	if (N < k) {
		if(verbose)
			cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	if(seeds == NULL) {
		init_array_2<double>(seeds,k,d);
	}

	// Seeding
	// In case of defective seeds from users, we should overwrite it by kmeans++ seeds
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds(d,N,k,data,seeds,verbose);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds(d,N,k,data,seeds,verbose);
	}

	if(verbose)
		cout << "Finished seeding" << endl;

	// Criteria's setup
	size_t iters = criteria.iterations, i = 0;
	double error = criteria.accuracy, e = error, e_prev = 0.0;
	double alpha = criteria.alpha;

	double ** c_tmp;
	init_array_2<double>(c_tmp,k,d);
	init_vector<i_vector>(clusters,k);

	double tmp = 0.0;

	// Initialize the centroids
	copy_array_2<double>(seeds,centroids,k,d);

	while (1) {
		// Assign the data points to clusters
		if(assign == KmeansAssignType::LINEAR)
			assign_to_closest_centroid(d,N,k,data,centroids,clusters,verbose);
		else if(assign == KmeansAssignType::NN_KD_TREE)
			assign_to_closest_centroid_2(d,N,k,data,centroids,clusters,verbose);
		else
			assign_to_closest_centroid_3(d,N,k,data,centroids,clusters,alpha,verbose);
		// Recalculate the centroids
		for(size_t j = 0; j < k; j++) {
			double * d_tmp = SimpleCluster::mean_vector(data,clusters[j],
					d,centroids[j]);
			copy_array<double>(d_tmp,c_tmp[j],d);
		}
		// Calculate the distortion
		e_prev = e;
		e = 0.0;
		for(size_t j = 0; j < k; j++) {
			tmp = SimpleCluster::distance_square(centroids[j],c_tmp[j],d);
			e += tmp;
		}
		e = sqrt(e);

		copy_array_2<double>(c_tmp,centroids,k,d);
		i++;
		if(i >= iters ||
				(e - e_prev < error && e - e_prev > -error)) break;
		//				(e < error && e > -error)) break;
	}

	if(verbose)
		cout << "Finished clustering with error is " <<
		e << " after " << i << " iterations." << endl;

	dealloc_array_2<double>(c_tmp,k);
}

/**
 * Calculate the distortion of a set of clusters.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param verbose for debugging
 */
double distortion(size_t d, size_t N, size_t k,
		double ** data, double ** centroids,
		vector<i_vector> clusters, bool verbose) {
	double e = 0.0;
	size_t i, j = 0;
	i_vector i_tmp;

	for(i = 0; i < k; i++) {
		i_tmp = clusters[i];
		while(j < i_tmp.size()) {
			e += SimpleCluster::distance_square(centroids[i], data[i_tmp[j]], d);
			++j;
		}
	}

	i_tmp.clear();

	return sqrt(e);
}
}

