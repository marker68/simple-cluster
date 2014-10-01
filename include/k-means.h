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

typedef vector<double> d_vector;
typedef vector<int> i_vector;

/**
 * Criteria
 */
typedef struct {
	double alpha;
	double accuracy;
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
		double **,
		double **&,
		size_t,
		size_t,
		size_t,
		bool);
void kmeans_pp_seeds(
		double **,
		double **&,
		size_t,
		size_t,
		size_t,
		bool);
void linear_assign(
		double **,
		double **,
		vector<i_vector>&,
		size_t,
		size_t,
		size_t,
		bool);
void kd_nn_assign(
		double **,
		double **,
		vector<i_vector>&,
		size_t,
		size_t,
		size_t,
		bool);
void kd_ann_assign(
		double **,
		double **,
		vector<i_vector>&,
		size_t,
		size_t,
		size_t,
		double,
		bool);
void greg_initialize(
		double **,
		double **,
		double **&,
		double *&,
		double *&,
		int *&,
		size_t *&,
		size_t,
		size_t,
		size_t,
		bool
		);
void simple_k_means(
		double **,
		double **&,
		int *&,
		double **&,
		KmeansType,
		KmeansAssignType,
		KmeansCriteria,
		size_t,
		size_t,
		size_t,
		bool);
double distortion(
		double **,
		double **,
		int *,
		size_t,
		size_t,
		size_t,
		bool);
}

#endif /* K_MEANS_H_ */
