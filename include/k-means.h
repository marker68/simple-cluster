/*
 * k-means.h
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

/**
 * Criteria
 */
typedef struct {
	double accuracy;
	int iterations;
} KmeansCriteria;

/**
 * Types of the k-means seeding
 */
enum class KmeansType {
	RANDOM_SEEDS, // randomly generated seeds
	KMEANS_PLUS_SEEDS, // k-means++
	USER_SEEDS // take the seeds from input
};

void random_seeds(size_t, size_t, size_t, double **, double **&, bool);
void kmeans_pp_seeds(size_t, size_t, size_t, double **, double **&, bool);
void assign_to_closest_centroid(size_t, size_t, size_t, double **,
		double **, int **&, bool);
void assign_to_closest_centroid_2(size_t, size_t, size_t, double **,
		double **, int **&, bool);
void simple_k_means(KmeansType, size_t, size_t, KmeansCriteria,size_t,
		double **, double **,
		int **, double **, bool);
double distortion(size_t d, size_t N, size_t k,
		double **, double **, int **);
}

#endif /* K_MEANS_H_ */
