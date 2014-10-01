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
 * After having a set of centers,
 * we need to assign data into each cluster respectively.
 * This solution uses linear search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param verbose for debugging
 */
void linear_assign(
		double ** data,
		double ** centers,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
	size_t i, tmp;

	double min = 0.0;

	for(i = 0; i < N; i++) {
		// Find the minimum distances between d_tmp and a centroid
		linear_search(centers,data[i],tmp,min,k,d,verbose);
		// Assign the data[i] into cluster tmp
		clusters[tmp].push_back(static_cast<int>(i));
	}
}

/**
 * After having a set of centers,
 * we need to assign data into each cluster respectively.
 * This solution uses kd-tree search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param verbose for debugging
 */
void kd_nn_assign(
		double ** data,
		double ** centers,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
	size_t i, tmp;
	KDNode<double> * root = nullptr;
	make_random_tree(root,centers,k,d,0,verbose);
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
 * After having a set of centers,
 * we need to assign data into each cluster respectively.
 * This solution uses ANN kd-tree search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the numbe rof clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param verbose for debugging
 */
void kd_ann_assign(
		double ** data,
		double ** centers,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		double alpha,
		bool verbose) {
	size_t i, tmp;
	KDNode<double> * root = nullptr;
	make_random_tree(root,centers,k,d,0,verbose);
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
 * Initialize for Greg's method
 * @param data the point data
 * @param centers the centers of clusters
 * @param sum vector sum of all points in each cluster
 * @param upper upper bound on the distance between point data and its assigned center
 * @param lower lower bound on the distance between point data and its second closest center
 * @param label the label of point data that specify its assigned cluster
 * @param N the size of data set
 * @param k the number of clusters
 * @param d the number of dimensions
 * @param verbose enable it to see the log
 * @return nothing
 */
void greg_initialize(
		double ** data,
		double ** centers,
		double **& sum,
		double *& upper,
		double *& lower,
		int *& label,
		size_t *& size,
		size_t N,
		size_t k,
		size_t d,
		bool verbose) {
	size_t tmp = -1;
	double min = DBL_MAX;
	double min2 = DBL_MAX;
	double d_tmp;

	// Initializing size and vector sum
	for(size_t i = 0; i < k; i++) {
		size[i] = 0;
		for(size_t j = 0; j < d; j++) {
			sum[i][j] = 0.0;
		}
	}

	for(size_t i = 0; i < N; i++) {
		for(size_t j = 0; j < k; j++) {
			d_tmp = SimpleCluster::distance(centers[j],data[i],d);
			if(min >= d_tmp) {
				min2 = min;
				min = d_tmp;
				tmp = j;
			} else {
				if(min2 > d_tmp) min2 = d_tmp;
			}
		}
		label[i] = tmp; // Update the label
		upper[i] = min; // Update the upper bound on this distance
		lower[i] = min2; // Update the lower bound on this distance

		// Update the size
		size[tmp]++;

		// Update the vector sum
		for(size_t j = 0; j < d; j++) {
			sum[tmp][j] += data[i][j];
		}
	}
}

/**
 * Update the centers
 * @param sum vector sum of all points in the cluster
 * @param size the size of each cluster
 * @param centers the centers of clusters
 * @param moved the distances that centers moved
 * @param k the number of clusters
 * @param d the number of dimensions
 * @return nothing
 */
void update_center(
		double ** sum,
		size_t * size,
		double **& centers,
		double *& moved,
		size_t k,
		size_t d) {
	double * c_tmp;
	init_array<double>(c_tmp,d);
	for(size_t i = 0; i < k; i++) {
		copy_array<double>(centers[i],c_tmp,d);
		for(size_t j = 0; j < d; j++) {
			if(size[i] > 0) centers[i][j] /= static_cast<double>(size[i]);
		}
		moved[i] = SimpleCluster::distance(c_tmp,centers[i],d);
	}
}

/**
 * Update the bounds
 * @param moved the distances that centers moved
 * @param label the labels of point data
 * @param upper
 * @param lower
 * @param N
 * @param k
 * @param d
 */
void update_bounds(
		double * moved,
		int * label,
		double *& upper,
		double *& lower,
		size_t N,
		size_t k) {
	size_t r = 0, r2 = 0;
	double max = 0.0, max2 = 0.0;
	for(size_t i = 0; i < k; i++) {
		if(max <= moved[i]) {
			max2 = max;
			max = moved[i];
			r2 = r;
			r = i;
		} else {
			if(max2 <= moved[i]) {
				max2 = moved[i];
				r2 = i;
			}
		}
	}
	for(size_t i = 0; i < N; i++) {
		upper[i] += moved[label[i]];
		double sub = 0.0;
		if(r == label[i]) {
			sub = moved[r2];
		} else {
			sub = moved[r];
		}
		lower[i] = fabs(lower[i] - sub);
//		cout << "lower " << i << ":" << lower[i] << endl;
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
 * @param centers the centers
 * @param label the labels of data points
 * @param seeds the initial centers = the seeds
 * @param verbose for debugging
 */
void simple_k_means(
		double ** data,
		double **& centers,
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

	// Variables for Greg's method
	double ** c_sum;
	double * moved;
	double * closest;
	double * upper;
	double * lower;
	size_t * size;

	init_array_2<double>(c_sum,k,d);
	init_array<double>(moved,k);
	init_array<double>(closest,k);
	init_array<double>(upper,N);
	init_array<double>(lower,N);
	init_array<size_t>(size,k);

	// Initialize the centers
	copy_array_2<double>(seeds,centers,k,d);
	greg_initialize(data,centers,c_sum,upper,lower,label,size,N,k,d,verbose);

	while (1) {
		// Assign the data points to clusters
		size_t tmp, visited;
		double min, min2, min_tmp, d_tmp, m;
		// Update the closest distances
		for(size_t j = 0; j < k; j++) {
			min2 = min = DBL_MAX;
			for(size_t t = 0; t < k, t != j; t++) {
				min_tmp = SimpleCluster::distance(centers[j],centers[t],d);
				if(min > min_tmp) min = min_tmp;
			}
			closest[j] = min;
		}

//		KDNode<double> * root = nullptr;
//		KDNode<double> query(d);
//		if(assign != KmeansAssignType::LINEAR) {
//			make_random_tree(root,centers,k,d,0,verbose);
//			if(root == nullptr) return;
//		}

		for(size_t j = 0; j < N; j++) {
			// Update m for bound test
			d_tmp = closest[label[j]]/2.0;
			m = MAX(d_tmp,lower[j]);
			// First bound test
			if(upper[j] > m) {
				// We need to tighten the upper bound
				upper[j] = SimpleCluster::distance(data[j],centers[label[j]],d);
				// Second bound test
				if(upper[j] > m) {
					size_t l = label[j];
					// Assign the data to clusters
					//					if(assign == KmeansAssignType::LINEAR) {
					for(size_t t = 0; t < k; t++) {
						d_tmp = SimpleCluster::distance(centers[t],data[j],d);
						if(min >= d_tmp) {
							min2 = min;
							min = d_tmp;
							tmp = j;
						} else {
							if(min2 > d_tmp) min2 = d_tmp;
						}
					}
					//					} else {
					//						KDNode<double> * nn = nullptr;
					//						min = DBL_MAX;
					//						visited = 0;
					//						query.add_data(data[j]);
					//						if(assign == KmeansAssignType::NN_KD_TREE) {
					//							nn_search(root,&query,nn,min,d,0,visited,verbose);
					//						} else {
					//							ann_search(root,&query,nn,min,alpha,d,0,visited,verbose);
					//						}
					//						tmp = nn->id;
					//					}

					// Assign the data[i] into cluster tmp
					label[j] = tmp; // Update the label
					upper[j] = min; // Update the upper bound on this distance
					lower[j] = min2; // Update the lower bound on this distance

					if(l != tmp) {
						size[tmp]++;
						size[l]--;
						for(size_t t = 0; t < d; t++) {
							c_sum[tmp][t] += data[j][t];
							c_sum[l][t] -= data[j][t];
						}
					}
				}
			}
		}

		// Move the centers
		update_center(c_sum,size,centers,moved,k,d);

		// Update the bounds
		update_bounds(moved,label,upper,lower,N,k);

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
				(e < error && e > -error)) break;
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
 * @param centers the centers
 * @param clusters the clusters
 * @param verbose for debugging
 */
double distortion(
		double ** data,
		double ** centers,
		int * label,
		size_t d,
		size_t N,
		size_t k,
		bool verbose) {
	double e = 0.0;
	for(size_t i = 0; i < N; i++) {
		e += distance_square(data[i],centers[label[i]],d);
	}
	return sqrt(e);
}
}

