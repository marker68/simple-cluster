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
#include <climits>
#include <stdlib.h>
#include <float.h>
#include <cmath>
#include "k-means.h"
#include "utilities.h"
#include "kd-tree.h"

#ifdef _OPENMP
#include <omp.h>
#define SET_THREAD_NUM omp_set_num_threads(n_thread)
#else
#define SET_THREAD_NUM 0 // disable multi-thread
#endif

#ifndef FLT_MAX
#define FLT_MAX 3.40282347E+38F
#endif

using namespace std;

namespace SimpleCluster {

/**
 * Create random seeds for k-means
 * @param data input data
 * @param seeds the seeds
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
void random_seeds(
		float ** data,
		float **& seeds,
		size_t d,
		size_t N,
		size_t k,
		int n_thread,
		bool verbose) {
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	size_t j;
#ifdef _WIN32
	size_t * tmp;
	SimpleCluster::init_array(tmp,k);
#else
	size_t tmp[k];
#endif

	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	for(i = 0; i < k; i++)
		tmp[i] = static_cast<size_t>(i);
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = k; i < N; i++) {
			uniform_int_distribution<> size_t_dis(i,N-1);
			j = size_t_dis(gen);
			if(j < k)
				tmp[j] = static_cast<size_t>(i);
		}
	}

	for(i = 0; i < k; i++) {
		if(!copy_array<float>(data[tmp[i]],seeds[i],d)) {
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
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
void kmeans_pp_seeds(
		float ** data,
		float **& seeds,
		size_t d,
		size_t N,
		size_t k,
		int n_thread,
		bool verbose) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	uniform_int_distribution<size_t> size_t_dis(0, N - 1);
	size_t tmp = size_t_dis(gen);

	float * d_tmp;
	if(!init_array<float>(d_tmp,d)) {
		cerr << "d_tmp: Cannot allocate the memory" << endl;
		exit(1);
	}

	copy_array<float>(data[tmp],seeds[0],d);
	float * distances;
	init_array<float>(distances,N);
	float * sum_distances;
	init_array<float>(sum_distances,N);
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < N; i++) {
			distances[i] = SimpleCluster::distance_square(data[i], d_tmp, d);
			sum_distances[i] = 0.0;
		}
	}

	float sum, tmp2, sum1, sum2, pivot;
	size_t count = 1;
	for(count = 1; count < k; count++) {
		sum = 0.0;
		for(i = 0; i < N; i++) {
			sum += distances[i];
			sum_distances[i] = sum;
		}
		uniform_real_distribution<float> real_dis(0, sum);
		pivot = real_dis(gen);

		for(i = 0; i < N - 1; i++) {
			sum1 = sum_distances[i];
			sum2 = sum_distances[i + 1];
			if(sum1 < pivot && pivot <= sum2)
				break;
		}
		copy_array<float>(data[(i+1)%N],d_tmp,d);
		copy_array<float>(d_tmp,seeds[count],d);
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
 * we need to assign data size_to each cluster respectively.
 * This solution uses linear search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
void linear_assign(
		float ** data,
		float ** centers,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		int n_thread,
		bool verbose) {
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	size_t tmp;

	float min = 0.0;
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < N; i++) {
			// Find the minimum distances between d_tmp and a centroid
			linear_search(centers,data[i],tmp,min,k,d,verbose);
			// Assign the data[i] size_to cluster tmp
			clusters[tmp].push_back(static_cast<size_t>(i));
		}
	}
}

/**
 * After having a set of centers,
 * we need to assign data size_to each cluster respectively.
 * This solution uses kd-tree search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
void kd_nn_assign(
		float ** data,
		float ** centers,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		int n_thread,
		bool verbose) {
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	size_t tmp;
	KDNode<float> * root = nullptr;
	make_random_tree(root,centers,k,d,0,verbose);
	if(root == nullptr) return;

	KDNode<float> query(d);
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < N; i++) {
			KDNode<float> * nn = nullptr;
			float min = FLT_MAX;
			size_t visited = 0;
			query.add_data(data[i]);
			// Find the minimum distances between d_tmp and a centroid
			nn_search(root,&query,nn,min,d,0,visited,verbose);
			tmp = nn->id;
			// Assign the data[i] size_to cluster tmp
			clusters[tmp].push_back(static_cast<size_t>(i));
		}
	}
}

/**
 * After having a set of centers,
 * we need to assign data size_to each cluster respectively.
 * This solution uses ANN kd-tree search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
void kd_ann_assign(
		float ** data,
		float ** centers,
		vector<i_vector>& clusters,
		size_t d,
		size_t N,
		size_t k,
		int n_thread,
		float alpha,
		bool verbose) {
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	size_t tmp;
	KDNode<float> * root = nullptr;
	make_random_tree(root,centers,k,d,0,verbose);
	if(root == nullptr) return;

	KDNode<float> query(d);

	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < N; i++) {
			KDNode<float> * nn = nullptr;
			float min = FLT_MAX;
			size_t visited = 0;
			query.add_data(data[i]);
			// Find the minimum distances between d_tmp and a centroid
			ann_search(root,&query,nn,min,alpha,d,0,visited,verbose);
			tmp = nn->id;
			// Assign the data[i] size_to cluster tmp
			clusters[tmp].push_back(static_cast<size_t>(i));
		}
	}
}


/**
 * Initialize for Greg's method
 * @param data the point data
 * @param centers the centers of clusters
 * @param sum vector sum of all posize_ts in each cluster
 * @param upper upper bound on the distance between posize_t data and its assigned center
 * @param lower lower bound on the distance between posize_t data and its second closest center
 * @param label the label of posize_t data that specify its assigned cluster
 * @param N the size of data set
 * @param k the number of clusters
 * @param d the number of dimensions
 * @param n_thread the number of threads
 * @param verbose enable it to see the log
 * @return nothing
 */
void greg_initialize(
		float ** data,
		float ** centers,
		float **& sum,
		float *& upper,
		float *& lower,
		size_t *& label,
		size_t *& size,
		size_t N,
		size_t k,
		size_t d,
		int n_thread,
		bool verbose) {
#ifdef _WIN32
	int i, j;
#else
	size_t i, j;
#endif
	// Initializing size and vector sum
	for(i = 0; i < k; i++) {
		size[i] = 0;
		for(j = 0; j < d; j++) {
			sum[i][j] = 0.0;
		}
	}

	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i,j)
		for(i = 0; i < N; i++) {
			float min = FLT_MAX;
			float min2 = FLT_MAX;
			float d_tmp;
			size_t tmp = -1;
			for(j = 0; j < k; j++) {
				d_tmp = SimpleCluster::distance(centers[j],data[i],d);
				if(min >= d_tmp) {
					min2 = min;
					min = d_tmp;
					tmp = j;
				} else {
					if(min2 >= d_tmp) min2 = d_tmp;
				}
			}

			label[i] = tmp; // Update the label
			upper[i] = min; // Update the upper bound on this distance
			lower[i] = min2; // Update the lower bound on this distance

			// Update the size
			size[tmp]++;

			// Update the vector sum
			for(j = 0; j < d; j++) {
				sum[tmp][j] += data[i][j];
			}
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
 * @param n_thread the number of threads
 * @return nothing
 */
void update_center(
		float ** sum,
		size_t * size,
		float **& centers,
		float *& moved,
		size_t k,
		size_t d,
		int n_thread) {
	float * c_tmp;
	init_array<float>(c_tmp,d);
	size_t i;
	for(i = 0; i < k; i++) {
		copy_array<float>(centers[i],c_tmp,d);
		for(size_t j = 0; j < d; j++) {
			centers[i][j] = (sum[i][j] + centers[i][j]) / static_cast<float>(size[i] + 1);
		}
		moved[i] = SimpleCluster::distance(c_tmp,centers[i],d);
	}
}

/**
 * Update the bounds
 * @param moved the distances that centers moved
 * @param label the labels of posize_t data
 * @param upper
 * @param lower
 * @param N
 * @param k
 * @param d
 * @param n_thread the number of threads
 */
void update_bounds(
		float * moved,
		size_t * label,
		float *& upper,
		float *& lower,
		size_t N,
		size_t k,
		int n_thread) {
	size_t r = 0;
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	float max = 0.0, max2 = 0.0, sub;
	for(i = 0; i < k; i++) {
		if(max <= moved[i]) {
			max2 = max;
			max = moved[i];
			r = i;
		} else {
			if(max2 <= moved[i]) {
				max2 = moved[i];
			}
		}
	}
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < N; i++) {
			upper[i] += moved[label[i]];
			sub = 0.0;
			if(r == label[i]) {
				sub = max2;
			} else {
				sub = max;
			}
			lower[i] = fabs(lower[i] - sub);
		}
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
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
void simple_k_means(
		float ** data,
		float **& centers,
		size_t *& label,
		float **& seeds,
		KmeansType type,
		KmeansAssignType assign,
		KmeansCriteria criteria,
		size_t N,
		size_t k,
		size_t d,
		int n_thread,
		bool verbose) {
	// Pre-check conditions
	if (N < k) {
		if(verbose)
			cerr << "There will be some empty clusters!" << endl;
		exit(1);
	}

	if(seeds == nullptr) {
		init_array_2<float>(seeds,k,d);
	}

	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds(data,seeds,d,N,k,n_thread,verbose);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds(data,seeds,d,N,k,n_thread,verbose);
	}

	if(verbose)
		cout << "Finished seeding" << endl;

	// Criteria's setup
	size_t iters = criteria.iterations, i = 0, count = 0;
	float error = criteria.accuracy, e = error, e_prev;

	// Variables for Greg's method
	float ** c_sum;
	float * moved;
	float * closest;
	float * upper;
	float * lower;
	size_t * size;

	init_array_2<float>(c_sum,k,d);
	init_array<float>(moved,k);
	init_array<float>(closest,k);
	init_array<float>(upper,N);
	init_array<float>(lower,N);
	init_array<size_t>(size,k);

	// Initialize the centers
	copy_array_2<float>(seeds,centers,k,d);
	greg_initialize(data,centers,c_sum,upper,lower,label,size,N,k,d,n_thread,verbose);
#ifdef _WIN32
	int j;
#else
	size_t j;
#endif

	while (1) {
		// Assign the data posize_ts to clusters
		size_t tmp = 0;
		float min, min2, min_tmp, d_tmp, m;
		// Update the closest distances
		for(j = 0; j < k; j++) {
			min2 = min = FLT_MAX;
			for(size_t t = 0; t < k; t++) {
				if(t != j) {
					min_tmp = SimpleCluster::distance(centers[j],centers[t],d);
					if(min > min_tmp) min = min_tmp;
				}
			}
			closest[j] = min;
		}

		SET_THREAD_NUM;
#pragma omp parallel
		{
#pragma omp for private(j,d_tmp,m,min,min2,tmp)
			for(j = 0; j < N; j++) {
				// Update m for bound test
				d_tmp = closest[label[j]]/2.0;
				m = std::max(d_tmp,lower[j]);
				// First bound test
				if(upper[j] > m) {
					// We need to tighten the upper bound
					upper[j] = SimpleCluster::distance(data[j],centers[label[j]],d);
					// Second bound test
					if(upper[j] > m) {
						size_t l = label[j];
						min2 = min = FLT_MAX;
						tmp = -1;
						// Assign the data to clusters
						for(size_t t = 0; t < k; t++) {
							d_tmp = SimpleCluster::distance(centers[t],data[j],d);
							if(min >= d_tmp) {
								min2 = min;
								min = d_tmp;
								tmp = t;
							} else {
								if(min2 > d_tmp) min2 = d_tmp;
							}
						}

						// Assign the data[i] size_to cluster tmp
						label[j] = tmp; // Update the label
						upper[j] = min; // Update the upper bound on this distance
						lower[j] = min2; // Update the lower bound on this distance

						if(l != tmp) {
							size[tmp]++;
							size[l]--;
							if(size[l] == 0) {
								if(verbose)
									cout << "An empty cluster was found!"
									" label = " << l << endl;
							}
							for(size_t t = 0; t < d; t++) {
								c_sum[tmp][t] += data[j][t];
								c_sum[l][t] -= data[j][t];
							}
						}
					}
				}
			}
		}

		// Move the centers
		update_center(c_sum,size,centers,moved,k,d,n_thread);

		// Update the bounds
		update_bounds(moved,label,upper,lower,N,k,n_thread);

		// Calculate the distortion
		e_prev = e;
		e = 0.0;
		for(size_t j = 0; j < k; j++) {
			e += moved[j];
		}
		e = sqrt(e);
		count += (fabs(e-e_prev) < error? 1 : 0);
		if(verbose)
			cout << "Iterator " << i
			<< "-th with error = " << e
			<< " and distortion = " << distortion(data,centers,label,d,N,k,n_thread,false)
			<< endl;
		i++;
		if(i >= iters || e < error || count >= 10) break;
	}

	if(verbose)
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
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
float distortion(
		float ** data,
		float ** centers,
		size_t * label,
		size_t d,
		size_t N,
		size_t k,
		int n_thread,
		bool verbose) {
	float e = 0.0;
#ifdef _WIN32
	int i;
#else
	size_t i;
#endif
	SET_THREAD_NUM;
#pragma omp parallel
	{
#pragma omp for private(i)
		for(i = 0; i < N; i++) {
			e += distance_square(data[i],centers[label[i]],d);
		}
	}
	return sqrt(e);
}
}
