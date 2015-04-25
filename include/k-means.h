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
 *  k-means.h
 *
 *  Created on: 2014/09/04
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#ifndef K_MEANS_H_
#define K_MEANS_H_

#include <iostream>
#include <exception>
#include <algorithm>
#include <vector>
#include <random>
#include <cstring>
#include <climits>
#include <cstdlib>
#include <cfloat>
#include <cmath>
#include "utilities.h"
#include "kd-tree.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef FLT_MAX
#define FLT_MAX 3.40282346638528859812e+38F
#endif

using namespace std;

/**
 * Cluster methods' space
 */
namespace SimpleCluster {

/**
 * Types of assigning methods
 */
enum class KmeansAssignType {
	LINEAR,
	NN_KD_TREE,
	ANN_KD_TREE
};

/**
 * Types of the k-means seeding
 */
enum class KmeansType {
	RANDOM_SEEDS, // randomly generated seeds
	KMEANS_PLUS_SEEDS, // k-means++
	USER_SEEDS // take the seeds from input
};

/**
 * Empty actions: how we treat the empty clusters
 */
enum class EmptyActs {
	SINGLETON,
	SINGLETON_2,
	NONE
};

/**
 * Criteria
 */
typedef struct {
	float alpha;
	float accuracy;
	int iterations;
} KmeansCriteria;

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
template<typename DataType>
inline void random_seeds(
		DataType * data,
		float *& seeds,
		int d,
		int N,
		int k,
		int n_thread,
		bool verbose) {
	int i;
	int j;
#ifdef _WIN32
	int * tmp;
	SimpleCluster::init_array<int>(tmp,k);
#else
	int tmp[k];
#endif

	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	for(i = 0; i < k; i++)
		tmp[i] = static_cast<int>(i);
#ifdef _OPENMP
	omp_set_num_threads(n_thread);
#pragma omp parallel
	{
#pragma omp for private(i)
#endif
		for(i = k; i < N; i++) {
			uniform_int_distribution<> int_dis(i,N-1);
			j = int_dis(gen);
			if(j < k + i)
				tmp[j - i] = static_cast<int>(i);
		}
#ifdef _OPENMP
	}
#endif

	int base = 0, base1;
	for(i = 0; i < k; i++) {
		base1 = tmp[i] * d;
		for(j = 0; j < d; j++) {
			seeds[base++] = static_cast<float>(data[base1++]);
		}
	}
}

/**
 * Create seeds for k-means++
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param seeds the seeds
 * @param n_thread the number of threads
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param verbose for debugging
 */
template<typename DataType>
inline void kmeans_pp_seeds(
		DataType * data,
		float *& seeds,
		DistanceType d_type,
		int d,
		int N,
		int k,
		int n_thread,
		bool verbose) {
	// For generating random numbers
	random_device rd;
	mt19937 gen(rd());

	uniform_int_distribution<int> int_dis(0, N - 1);
	int tmp = int_dis(gen);

	size_t base = static_cast<size_t>(tmp) * static_cast<size_t>(d);
	int i, i0, start, end, p = N / n_thread;
	for(i = 0; i < d; i++) {
		seeds[i] = static_cast<float>(data[base++]);
	}

	float * distances;
	init_array<float>(distances,N);
	float * sum_distances;
	init_array<float>(sum_distances,N);
#ifdef _OPENMP
	omp_set_num_threads(n_thread);
#pragma omp parallel
	{
#pragma omp for private(i, start, end, i0)
#endif
		for(i0 = 0; i0 < n_thread; i0++) {
			start = p * i0;
			end = start + p;
			if(end >= N || i0 == n_thread - 1) end = N;
			DataType * d_tmp2 = data + static_cast<size_t>(start) * static_cast<size_t>(d);
			float * d_tmp = seeds;
			for(i = start; i < end; i++) {
				if(d_type == DistanceType::NORM_L2)
					distances[i] = distance_l2_square<DataType,float>(d_tmp2, d_tmp, d);
				else if(d_type == DistanceType::NORM_L1)
					distances[i] = distance_l1<DataType,float>(d_tmp2, d_tmp, d);
				sum_distances[i] = 0.0;
				d_tmp2 += d;
			}
		}
#ifdef _OPENMP
	}
#endif
	float sum, tmp2 = 0.0, sum1, sum2, pivot;
	int count = 1, j, t;
	size_t base1, base2;
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
		j = (i + 1) % N;
		base1 = count * d;
		base2 = static_cast<size_t>(j) * d;
		for(t = 0; t < d; t++) {
			seeds[base1++] = static_cast<float>(data[base2++]);
		}

		// Update the distances
		if(count < k) {
#ifdef _OPENMP
			omp_set_num_threads(n_thread);
#pragma omp parallel
			{
#pragma omp for private(i, start, end, i0, tmp2)
#endif
				for(i0 = 0; i0 < n_thread; i0++) {
					start = p * i0;
					end = start + p;
					if(end >= N || i0 == n_thread - 1) end = N;
					DataType * d_tmp2 = data + static_cast<size_t>(start) * static_cast<size_t>(d);
					float * d_tmp = seeds + (count - 1) * d; // We only need to compare the old closest distances with the new one
					for(i = start; i < end; i++) {
						if(d_type == DistanceType::NORM_L2)
							tmp2 = distance_l2_square<float,DataType>(d_tmp,d_tmp2,d);
						else if(d_type == DistanceType::NORM_L1)
							tmp2 = distance_l1<float,DataType>(d_tmp,d_tmp2,d);
						if(distances[i] > tmp2) distances[i] = tmp2;
						d_tmp2 += d;
					}
				}
#ifdef _OPENMP
			}
#endif
		}
		if(verbose)
			cout << "Got " << count << " centers" << endl;
	}
}

/**
 * After having a set of centers,
 * we need to assign data into each cluster respectively.
 * This solution uses linear search to assign data.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
template<typename DataType>
inline void linear_assign(
		DataType * data,
		float * centers,
		int *& labels,
		int *& size,
		float *& sum,
		DistanceType d_type,
		int d,
		int N,
		int k,
		int n_thread,
		bool verbose) {
	if(n_thread < 1) n_thread = 1;
	int i, j, m, i0, start, end;
	int tmp;
	DataType * d_tmp;
	float *  d_tmp1;

	float min = FLT_MAX, min_tmp = 0.0;
	d_tmp = data;
	int base1, base2;
	for(i = 0; i < N; i++) {
		// Find the minimum distances between d_tmp and a centroid
		min = FLT_MAX;
		d_tmp1 = centers;
		for(j = 0; j < k; j++) {
			if(d_type == DistanceType::NORM_L2)
				min_tmp = distance_l2_square<DataType,float>(d_tmp,d_tmp1,d);
			else if(d_type == DistanceType::NORM_L1)
				min_tmp = distance_l1<DataType,float>(d_tmp,d_tmp1,d);
			if(min > min_tmp) {
				min = min_tmp;
				tmp = j;
			}
			d_tmp1 += d;
		}
		// Assign the data[i] into cluster tmp
		if(labels[i] > -1) {
			size[labels[i]]--;
			base1 = labels[i] * d;
			base2 = i * d;
			for(m = 0; m < d; m++) {
				sum[base1++] -= static_cast<float>(data[base2++]);
			}
		}
		labels[i] = tmp;
		size[tmp]++;
		base1 = tmp * d;
		base2 = i * d;
		for(m = 0; m < d; m++) {
			sum[base1++] += static_cast<float>(data[base2++]);
		}
	}
}

/**
 * Update the centers
 * @param sum the sum vector of all points in the cluster
 * @param size the size of each cluster
 * @param centers the centers of clusters
 * @param moved the distances that centers moved
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param k the number of clusters
 * @param d the number of dimensions
 * @param n_thread the number of threads
 * @return nothing
 */
inline void update_center(
		float * sum,
		int * size,
		float *& centers,
		float *& moved,
		DistanceType d_type,
		int k,
		int d,
		int n_thread) {
	float * c_tmp;
	init_array<float>(c_tmp,d);
	int i, base = 0;
	c_tmp = centers;
	for(i = 0; i < k; i++) {
		for(int j = 0; j < d; j++) {
			centers[base] = static_cast<float>(sum[base++] / size[i]);
		}
		if(d_type == DistanceType::NORM_L2)
			moved[i] = distance_l2<float>(c_tmp,centers + (base - d),d);
		else if(d_type == DistanceType::NORM_L1)
			moved[i] = distance_l1<float>(c_tmp,centers + (base - d),d);
		c_tmp += d;
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
 * @param n_thread the number of threads
 */
inline void update_bounds(
		float * moved,
		int * label,
		float *& upper,
		float *& lower,
		int N,
		int k,
		int n_thread) {
	int r = 0;
	int i;

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
#ifdef _OPENMP
	omp_set_num_threads(n_thread);
#pragma omp parallel
	{
#pragma omp for private(i)
#endif
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
#ifdef _OPENMP
	}
#endif
}

/**
 * Calculate the distortion of a set of clusters.
 * @param d the dimensions of the data
 * @param N the number of the data
 * @param k the number of clusters
 * @param data input data
 * @param centers the centers
 * @param clusters the clusters
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
template<typename DataType>
inline float distortion(
		DataType * data,
		float * centers,
		int * label,
		DistanceType d_type,
		int d,
		int N,
		int k,
		bool verbose) {
	float e = 0.0;
	int j;
	DataType * tmp = data;
	for(j = 0; j < N; j++) {
		if(d_type == DistanceType::NORM_L2)
			e += distance_l2_square<DataType,float>(tmp,centers + label[j] * d,d);
		else if(d_type == DistanceType::NORM_L1)
			e += distance_l1<DataType,float>(tmp,centers + label[j] * d,d);
		tmp += d;
	}
	return sqrt(e);
}

/**
 * Update the farthest distances
 */
template<typename DataType>
inline void find_farthest(
		DataType * data,
		float * centers,
		int * labels,
		DistanceType d_type,
		int id,
		float& dfst,
		int& fst,
		int N,
		int k,
		int d,
		bool verbose) {
	int i;
	float d_tmp;
	dfst = -1.0f;
	DataType * tmp = data;
	for(i = 0; i < N; i++) {
		if(labels[i] == id) {
			if(d_type == DistanceType::NORM_L2)
				d_tmp = distance_l2_square<DataType,float>(tmp,centers,d);
			else if(d_type == DistanceType::NORM_L1)
				d_tmp = distance_l1<DataType,float>(tmp,centers,d);
			if(dfst < d_tmp) {
				dfst = d_tmp;
				fst = i;
			}
		}
		tmp += d;
	}
	dfst = sqrt(dfst);
}

/**
 * Find a lonely observer
 */
template<typename DataType>
inline void find_lonely(
		DataType * data,
		float * centers,
		int * labels,
		DistanceType d_type,
		float& dfst,
		int& fst,
		int N,
		int k,
		int d,
		bool verbose) {
	int i;
	float d_tmp;
	dfst = -1.0f;
	DataType * tmp = data;
	for(i = 0; i < N; i++) {
		if(d_type == DistanceType::NORM_L2)
			d_tmp = distance_l2_square<DataType,float>(tmp,centers + labels[i] * d,d);
		else if(d_type == DistanceType::NORM_L1)
			d_tmp = distance_l1<DataType,float>(tmp,centers + labels[i] * d,d);
		if(dfst < d_tmp) {
			dfst = d_tmp;
			fst = i;
		}
		tmp += d;
	}
	dfst = sqrt(dfst);
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
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
template<typename DataType>
inline void greg_initialize(
		DataType * data,
		float * centers,
		float *& sum,
		float *& upper,
		float *& lower,
		int *& label,
		int *& size,
		DistanceType d_type,
		EmptyActs ea,
		int N,
		int k,
		int d,
		int n_thread,
		bool verbose) {
	size_t base = 0, p = N / n_thread;
	// Initializing size and vector sum
	for(int i = 0; i < k; i++) {
		size[i] = 0;
		for(int j = 0; j < d; j++) {
			sum[base++] = 0.0;
		}
	}

	float min, min2, d_tmp;
	int tmp;
#ifdef _OPENMP
	omp_set_num_threads(n_thread);
#pragma omp parallel
	{
#pragma omp for private(d_tmp,min,min2,tmp)
#endif
		for(int i0 = 0; i0 < n_thread; i0++) {
			size_t start = p * i0;
			size_t end = start + p;
			size_t base1 = 0, base2 = 0;
			if(end > N || i0 == n_thread - 1) end = N;
			DataType * dt = data + start * static_cast<size_t>(d);
			for(size_t i = start; i < end; i++) {
				min = FLT_MAX;
				min2 = FLT_MAX;
				d_tmp = 0.0;
				tmp = -1;
				for(size_t j = 0; j < k; j++) {
					if(d_type == DistanceType::NORM_L2) {
						d_tmp = distance_l2_square<float,DataType>(centers + j * d,dt,d);
					} else if(d_type == DistanceType::NORM_L1) {
						d_tmp = distance_l1<float,DataType>(centers + j * d,dt,d);
					}
					if(min >= d_tmp) {
						min2 = min;
						min = d_tmp;
						tmp = j;
					} else {
						if(min2 >= d_tmp) min2 = d_tmp;
					}
				}

				label[i] = tmp; // Update the label
				upper[i] = sqrt(min); // Update the upper bound on this distance
				lower[i] = sqrt(min2); // Update the lower bound on this distance

				// Update the size
				size[tmp]++;

				// Update the vector sum
				base1 = tmp * d;
				base2 = i * d;
				for(size_t j = 0; j < d; j++) {
					sum[base1++] += static_cast<float>(data[base2++]);
				}
				dt += d;
			}
		}
#ifdef _OPENMP
	}
#endif
	
    size_t s_max, l_tmp, base3, base4;
    int fst;
	float dfst;
	// Check for empty clusters
	if(ea != EmptyActs::NONE) {
		for(int i = 0; i < k; i++) {
			if(size[i] <= 0) {
				if(ea == EmptyActs::SINGLETON_2) {
					l_tmp = 0;
					for(int j = 0; j < k; j++) {
						if(l_tmp < size[j]) {
							l_tmp = size[j];
							s_max = j;
						}
					}
				}
				// Move the centers
				base = i * d;
				if(ea == EmptyActs::SINGLETON)
					find_lonely<DataType>(data,centers,label,d_type,
							dfst,fst,N,k,d,verbose);
				else if(ea == EmptyActs::SINGLETON_2)
					find_farthest<DataType>(data,centers + base,label,d_type,
							s_max,dfst,fst,N,k,d,verbose);
				base3 = static_cast<size_t>(fst) * d;
				base4 = label[fst] * d;
				for(int j = 0; j < d; j++) {
					centers[base] = static_cast<float>(data[base3++]);
					sum[base] = centers[base];
					sum[base4++] -= centers[base++];
				}
				size[i] = 1;
				size[label[fst]]--;
				label[fst] = i;
			}
		}
	}
}

template<typename DataType>
inline void greg_kmeans(
		DataType * data,
		float *& centers,
		int *& label,
		float *& seeds,
		KmeansType type,
		KmeansCriteria criteria,
		DistanceType d_type,
		EmptyActs ea,
		int N,
		int k,
		int d,
		int n_thread,
		bool verbose) {
	// Pre-check conditions
	if (N < k) {
		if(verbose)
			cerr << "There will be some empty clusters!" << endl;
		// Create a new infinite point
		float * inf = (float *)::operator new(d * sizeof(float));
		fill(inf, inf + d, FLT_MAX);
		for(int i = 0; i < k; i++) {
			label[i] = i;
			if(i < N) {
				memcpy(centers + i * d, data + i * d, d * sizeof(float));
			} else {
				memcpy(centers + i * d, inf, d * sizeof(float));
			}
		}

		return;
	}

	if(seeds == nullptr) {
		init_array<float>(seeds,k * d);
	}

	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds<DataType>(data,seeds,d,N,k,n_thread,verbose);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds<DataType>(data,seeds,d_type,d,N,k,n_thread,verbose);
	}

	if(verbose)
		cout << "Finished seeding" << endl;

	// Criteria's setup
	int iters = criteria.iterations, it = 0, count = 0;
	float error = criteria.accuracy, e = error, e_prev;

	// Variables for Greg's method
	float * c_sum;
	float * moved;
	float * closest;
	float * upper;
	float * lower;
	int * size;

	init_array<float>(c_sum,k * d);
	init_array<float>(moved,k);
	init_array<float>(closest,k);
	init_array<float>(upper,N);
	init_array<float>(lower,N);
	init_array<int>(size,k);

	int i0, i, j, s_max, l_tmp, fst, base, base0, base1, base2;
	size_t p = N / n_thread;
	float * fpt1, * fpt2;
	DataType * dpt = data;
	int tmp = 0;
	float min, min2, min_tmp = 0.0,
			d_tmp = 0.0, m, dfst;

	// Initialize the centers
	copy_array<float>(seeds,centers,k * d);
	greg_initialize<DataType>(data,centers,c_sum,upper,lower,
			label,size,d_type,ea,N,k,d,n_thread,verbose);
	if(verbose)
		cout << "Finished initialization" << endl;

	while (1) {
		// Update the closest distances
		fpt1 = centers;
		for(i = 0; i < k; i++) {
			min2 = min = FLT_MAX;
			fpt2 = centers;
			for(j = 0; j < k; j++) {
				if(j != i) {
					if(d_type == DistanceType::NORM_L2)
						min_tmp = distance_l2_square<float>(fpt1,fpt2,d);
					else if(d_type == DistanceType::NORM_L1)
						min_tmp = distance_l1<float>(fpt1,fpt2,d);
					if(min > min_tmp) min = min_tmp;
				}
				fpt2 += d;
			}
			closest[i] = sqrt(min);
			fpt1 += d;
		}

#ifdef _OPENMP
		omp_set_num_threads(n_thread);
#pragma omp parallel
		{
#pragma omp for private(i,j,d_tmp,m,min,min2,tmp)
#endif
			for(i0 = 0; i0 < n_thread; i0++) {
				size_t start = p * i0;
				size_t end = start + p;
				if(end > N || i0 == n_thread - 1) end = N;
				int l;
				for(i = start; i < end; i++) {
					// Update m for bound test
					d_tmp = closest[label[i]]/2.0;
					m = std::max(d_tmp,lower[i]);
					// First bound test
					if(upper[i] > m) {
						// We need to tighten the upper bound
						if(d_type == DistanceType::NORM_L2)
							upper[i] = distance_l2<DataType,float>(data + i * d,centers + label[i] * d,d);
						else if(d_type == DistanceType::NORM_L1)
							upper[i] = distance_l1<DataType,float>(data + i * d,centers + label[i] * d,d);
						// Second bound test
						if(upper[i] > m) {
							l = label[i];
							min2 = min = FLT_MAX;
							tmp = -1;
							// Assign the data to clusters
							fpt1 = centers;
							for(j = 0; j < k; j++) {
								if(d_type == DistanceType::NORM_L2)
									d_tmp = distance_l2_square<float,DataType>(fpt1,data + i * d,d);
								else if(d_type == DistanceType::NORM_L1)
									d_tmp = distance_l1<float,DataType>(fpt1,data + i * d,d);
								if(min >= d_tmp) {
									min2 = min;
									min = d_tmp;
									tmp = j;
								} else {
									if(min2 > d_tmp) min2 = d_tmp;
								}
								fpt1 += d;
							}

							// Assign the data[i] into cluster tmp
							label[i] = tmp; // Update the label
							upper[i] = sqrt(min); // Update the upper bound on this distance
							lower[i] = sqrt(min2); // Update the lower bound on this distance

							if(l != tmp) {
								size[tmp]++;
								size[l]--;
								if(size[l] == 0) {
									if(verbose)
										cout << "An empty cluster was found!"
										" label = " << l << endl;
								}
								base = i * d;
								base0 = tmp * d;
								base1 = l * d;
								for(j = 0; j < d; j++) {
									c_sum[base0++] += static_cast<float>(data[base]);
									c_sum[base1++] -= static_cast<float>(data[base++]);
								}
							}
						}
					}
				}
			}
#ifdef _OPENMP
		}
#endif
		// Check for empty clusters
		if(ea != EmptyActs::NONE) {
			for(i = 0; i < k; i++) {
				if(size[i] <= 0) {
					if(ea == EmptyActs::SINGLETON_2) {
						l_tmp = 0;
						for(j = 0; j < k; j++) {
							if(l_tmp < size[j]) {
								l_tmp = size[j];
								s_max = j;
							}
						}
					}
					// Move the centers
					base = i * d;
					if(ea == EmptyActs::SINGLETON)
						find_lonely<DataType>(data,centers,label,d_type,
								dfst,fst,N,k,d,verbose);
					else if(ea == EmptyActs::SINGLETON_2)
						find_farthest<DataType>(data,centers + base,label,d_type,
								s_max,dfst,fst,N,k,d,verbose);
					base1 = fst * d;
					base2 = label[fst] * d;
					for(j = 0; j < d; j++) {
						centers[base] = static_cast<float>(data[base1++]);
						c_sum[base] = centers[base];
						c_sum[base2++] -= centers[base++];
					}
					size[i] = 1;
					size[label[fst]]--;
					label[fst] = i;
				}
			}
		}
		// Move the centers
		update_center(c_sum,size,centers,moved,d_type,k,d,n_thread);
		// Update the bounds
		update_bounds(moved,label,upper,lower,N,k,n_thread);

		if(verbose) {
			cout << "Spec " << it << ":";
			for(i = 0; i < k; i++) {
				cout << size[i] << " ";
			}
			cout << endl;
		}

		// Calculate the distortion
		e_prev = e;
		e = 0.0;
		for(i = 0; i < k; i++) {
			e += moved[i];
		}
		e = sqrt(e);
		count += (fabs(e-e_prev) < error? 1 : 0);
		if(verbose)
			cout << "Iterator " << it
			<< "-th with error = " << e
			<< " and distortion = "
			<< distortion(data,centers,label,d_type,d,N,k,false)
			<< endl;
		it++;
		if(it >= iters || e < error || count >= 10) break;
	}

	if(verbose)
		cout << "Finished clustering with error is " <<
		e << " after " << it << " iterations." << endl;
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
 * @param d_type the type of distance. Available options are NORM_L1, NORM_L2, HAMMING
 * @param n_thread the number of threads
 * @param verbose for debugging
 */
template<typename DataType>
inline void simple_kmeans(
		DataType * data,
		float *& centers,
		int *& labels,
		float *& seeds,
		KmeansType type,
		KmeansAssignType assign,
		KmeansCriteria criteria,
		DistanceType d_type,
		EmptyActs ea,
		int N,
		int k,
		int d,
		int n_thread,
		bool verbose) {
	// Pre-check conditions
	if (N < k) {
		if(verbose)
			cerr << "There will be some empty clusters!" << endl;
		// Create a new infinite point
		float * inf = (float *)::operator new(d * sizeof(float));
		fill(inf, inf + d, FLT_MAX);
		for(int i = 0; i < k; i++) {
			labels[i] = i;
			if(i < N) {
				memcpy(centers + i * d, data + i * d, d * sizeof(float));
			} else {
				memcpy(centers + i * d, inf, d * sizeof(float));
			}
		}

		return;
	}

	if(seeds == nullptr) {
		init_array<float>(seeds,k*d);
	}

	// Seeding
	if (type == KmeansType::RANDOM_SEEDS) {
		random_seeds<DataType>(data,seeds,d,N,k,n_thread,verbose);
	} else if(type == KmeansType::KMEANS_PLUS_SEEDS) {
		kmeans_pp_seeds<DataType>(data,seeds,d_type,d,N,k,n_thread,verbose);
	}

	if(verbose)
		cout << "Finished seeding" << endl;

	// Criteria's setup
	int iters = criteria.iterations, it = 0, count = 0;
	float error = criteria.accuracy, e = error, e_prev,
			d_tmp, dfst;
	int i, j, f_tmp, fst,s_max, l_tmp;

	// Initialize the centers
	copy_array<float>(seeds,centers,k*d);
	init_array<int>(labels,N);
	for(i = 0; i < N; i++) labels[i] = -1;
	int * size;
	init_array<int>(size,k);
	float * sum;
	float * c_tmp;
	init_array<float>(c_tmp,d);
	init_array<float>(sum,k*d);
	int base = 0, base1, base2;
	for(i = 0; i < k; i++) {
		for(j = 0; j < d; j++)
			sum[base++] = 0.0;
		size[i] = 0;
	}
	if(verbose)
		cout << "Finished initialization" << endl;

	while (1) {
		// Assigning
		linear_assign<DataType>(data,centers,labels,size,sum,
				d_type,d,N,k,n_thread,verbose);
		// Check for empty clusters
		if(ea != EmptyActs::NONE) {
			for(i = 0; i < k; i++) {
				if(size[i] <= 0) {
					if(ea == EmptyActs::SINGLETON_2) {
						l_tmp = 0;
						for(j = 0; j < k; j++) {
							if(l_tmp < size[j]) {
								l_tmp = size[j];
								s_max = j;
							}
						}
					}
					// Move the centers
					base = i * d;
					if(ea == EmptyActs::SINGLETON)
						find_lonely<DataType>(data,centers,labels,d_type,
								dfst,fst,N,k,d,verbose);
					else if(ea == EmptyActs::SINGLETON_2)
						find_farthest<DataType>(data,centers + base,labels,d_type,
								s_max,dfst,fst,N,k,d,verbose);
					base1 = fst * d;
					base2 = labels[fst] * d;
					for(j = 0; j < d; j++) {
						centers[base] = static_cast<float>(data[base1++]);
						sum[base] = centers[base];
						sum[base2++] -= centers[base++];
					}
					size[i] = 1;
					size[labels[fst]]--;
					labels[fst] = i;
				}
			}
		}

		if(verbose) {
			cout << "Spec " << it << ":";
			for(i = 0; i < k; i++) {
				cout << size[i] << " ";
			}
			cout << endl;
		}

		// Update centers
		e_prev = e;
		e = 0.0;
		base = 0;
		for(i = 0; i < k; i++) {
			for(j = 0; j < d; j++) {
				c_tmp[j] = centers[base];
				centers[base] = sum[base++] / size[i];
			}
			e += distance_l2_square<float>(c_tmp,centers + (base - d),d);
		}
		e = sqrt(e);
		count += (fabs(e-e_prev) < error? 1 : 0);

		if(verbose)
			cout << "Iterator " << it
			<< "-th with error = " << e
			<< " and distortion = " <<
			distortion(data,centers,labels,d_type,
					d,N,k,false)
					<< endl;
		it++;

		if(
				it >= iters ||
				fabs(e-e_prev) < error ||
				count >= 10
		) break;
	}

	if(verbose)
		cout << "Finished clustering with error is " <<
		e << " after " << it << " iterations." << endl;
}
}

#endif /* K_MEANS_H_ */
