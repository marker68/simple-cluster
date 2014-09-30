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
#include <exception>
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
 * @param data input data
 * @param seeds the seeds
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param verbose for debugging
 */
void random_seeds(
		double ** data,
		double **& seeds,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
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
void kmeans_pp_seeds(
		double ** data,
		double **& seeds,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	uniform_int_distribution<size_t> int_dis(0, N - 1);
	size_t tmp = int_dis(gen);

	double * d_tmp;
	if(!init_array<double>(d_tmp,d)) {
		cerr << "d_tmp: Cannot allocate the memory" << endl;
		exit(1);
	}
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
		copy_array<double>(data[(i+1)%N],d_tmp,d);
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
void linear_assign(
		double ** data,
		double ** centroids,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
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
 * @param k the number of clusters
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param verbose for debugging
 */
void kd_nn_assign(
		double ** data,
		double ** centroids,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
	size_t i, tmp;
	KDNode<double> * root = nullptr;
	make_random_tree(root,centroids,k,d,0,verbose);
	if(root == nullptr) return;

	KDNode<double> query(d);

	for(i = 0; i < N; i++) {
		KDNode<double> * nn = nullptr;
		double min = DBL_MAX;
		size_t visited = 0;
		query.add_data(data[i]);
		// Find the minimum distances between d_tmp and a centroid
		nn_search(root,&query,nn,min,d,0,visited,verbose);
		tmp = nn->id;
		// Assign the data[i] into cluster tmp
		clusters[tmp].push_back(static_cast<int>(i));
	}
}

/**
 * After having a set of centroids,
 * we need to assign data into each cluster respectively.
 * This solution uses ANN kd-tree search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param verbose for debugging
 */
void kd_ann_assign(
		double ** data,
		double ** centroids,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		double alpha,
		bool verbose) {
	size_t i, tmp;
	KDNode<double> * root = nullptr;
	make_random_tree(root,centroids,k,d,0,verbose);
	if(root == nullptr) return;

	KDNode<double> query(d);

	for(i = 0; i < N; i++) {
		KDNode<double> * nn = nullptr;
		double min = DBL_MAX;
		size_t visited = 0;
		query.add_data(data[i]);
		// Find the minimum distances between d_tmp and a centroid
		ann_search(root,&query,nn,min,alpha,d,0,visited,verbose);
		tmp = nn->id;
		// Assign the data[i] into cluster tmp
		clusters[tmp].push_back(static_cast<int>(i));
	}
}

/**
 * The k-means method: a description of the method can be found at
 * http://home.deib.polimi.it/matteucc/Clustering/tutorial_html/kmeans.html
 * @param type the type of seeding method
 * @param assign the type of assigning method
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param criteria the criteria
 * @param data input data
 * @param centroids the centroids
 * @param label the labels of data points
 * @param seeds the initial centroids = the seeds
 * @param verbose for debugging
 */
void simple_k_means(
		double ** data,
		double **& centroids,
		int *& label,
		double **& seeds,
		KmeansType type,
		KmeansAssignType assign,
		KmeansCriteria criteria,
		size_t N,
		size_t k,
		size_t d,
		bool verbose) {
	// Pre-check conditions
	if (N < k) {
		if(verbose)
			cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	if(seeds == nullptr) {
		init_array_2<double>(seeds,k,d);
	}

	// Seeding
	// In case of defective seeds from users, we should overwrite it by kmeans++ seeds
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds(data,seeds,d,N,k,verbose);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds(data,seeds,d,N,k,verbose);
	}

	if(verbose)
		cout << "Finished seeding" << endl;

	// Criteria's setup
	size_t iters = criteria.iterations, i = 0;
	double error = criteria.accuracy, e = error, e_prev = 0.0;
	double alpha = criteria.alpha;

//	double ** c_tmp;
//	init_array_2<double>(c_tmp,k,d);

	// Variables for Greg's method
	double * moved;
	double * closest;
	double * upper;
	double * lower;
	size_t * size;

	init_array<double>(moved,k);
	init_array<double>(closest,k);
	init_array<double>(upper,N);
	init_array<double>(lower,N);
	init_array<size_t>(size,k);

	// Initialize the centroids
	copy_array_2<double>(seeds,centroids,k,d);

	while (1) {
		// Initialize
		for(size_t j = 0; j < k; j++) {
			size[j] = 0;
		}
		// Assign the data points to clusters
		size_t tmp, visited;
		double min = 0.0;
		KDNode<double> * root = nullptr;
		KDNode<double> query(d);
		if(assign != KmeansAssignType::LINEAR) {
			make_random_tree(root,centroids,k,d,0,verbose);
			if(root == nullptr) return;
		}

		for(size_t j = 0; j < N; j++) {
			// Assign the data to clusters
			if(assign == KmeansAssignType::LINEAR) {
				linear_search(centroids,data[j],tmp,min,k,d,verbose);
			} else {
				KDNode<double> * nn = nullptr;
				min = DBL_MAX;
				visited = 0;
				query.add_data(data[j]);
				if(assign == KmeansAssignType::NN_KD_TREE) {
					nn_search(root,&query,nn,min,d,0,visited,verbose);
				} else {
					ann_search(root,&query,nn,min,alpha,d,0,visited,verbose);
				}
				tmp = nn->id;
			}
			// Assign the data[i] into cluster tmp
			label[j] = tmp;
			size[tmp]++;
		}

		// Recalculate the centroids
		all_mean_vector(data,label,size,centroids,moved,d,N,k);
		// Calculate the distortion
		e_prev = e;
		e = 0.0;
		for(size_t j = 0; j < k; j++) {
			e += moved[j];
		}
		e = sqrt(e);
		i++;
		if(i >= iters ||
//				(e - e_prev < error && e - e_prev > -error) ||
				(e < error && e > -error)
				) break;
	}

//	if(verbose)
		cout << "Finished clustering with error is " <<
		e << " after " << i << " iterations." << endl;
}

/**
 * Calculate the distortion of a set of clusters.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param centroids the centroids
 * @param clusters the clusters
 * @param verbose for debugging
 */
double distortion(
		double ** data,
		double ** centroids,
		int * label,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
	double e = 0.0;
	for(size_t i = 0; i < N; i++) {
		e += distance_square(data[i],centroids[label[i]],d);
	}
	return sqrt(e);
}
}

