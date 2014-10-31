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
#include <vector>

using namespace std;

/**
 * Cluster methods' space
 */
namespace SimpleCluster {

typedef vector<float> d_vector;
typedef vector<size_t> i_vector;

/**
 * Criteria
 */
typedef struct {
	float alpha;
	float accuracy;
	int iterations;
} KmeansCriteria;

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

void random_seeds(
		float **,
		float **&,
		size_t,
		size_t,
		size_t,
		int,
		bool);
void kmeans_pp_seeds(
		float **,
		float **&,
		size_t,
		size_t,
		size_t,
		int,
		bool);
void linear_assign(
		float **,
		float **,
		vector<i_vector>&,
		size_t,
		size_t,
		size_t,
		int,
		bool);
void kd_nn_assign(
		float **,
		float **,
		vector<i_vector>&,
		size_t,
		size_t,
		size_t,
		int,
		bool);
void kd_ann_assign(
		float **,
		float **,
		vector<i_vector>&,
		size_t,
		size_t,
		size_t,
		int,
		float,
		bool);
void greg_initialize(
		float **,
		float **,
		float **&,
		float *&,
		float *&,
		size_t *&,
		size_t *&,
		size_t,
		size_t,
		size_t,
		int,
		bool
		);
void update_center(
		float **,
		size_t *,
		float **&,
		float *&,
		size_t,
		size_t,
		int);
void update_bounds(
		float *,
		size_t *,
		float *&,
		float *&,
		size_t,
		size_t,
		int);
void simple_k_means(
		float **,
		float **&,
		size_t *&,
		float **&,
		KmeansType,
		KmeansAssignType,
		KmeansCriteria,
		size_t,
		size_t,
		size_t,
		int,
		bool);
float distortion(
		float **,
		float **,
		size_t *,
		size_t,
		size_t,
		size_t,
		int,
		bool);
}

#endif /* K_MEANS_H_ */
