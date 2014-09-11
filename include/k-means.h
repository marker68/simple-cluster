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

void random_seeds(size_t d, size_t N, size_t k, double ** data, double ** seeds);
void kmeans_pp_seeds(size_t d, size_t N, size_t k, double ** data, double ** seeds);
void assign_to_closest_centroid(size_t d, size_t N, size_t k, double ** data,
		double ** centroids, int ** clusters);
void assign_to_closest_centroid_2(size_t d, size_t N, size_t k, double ** data,
		double ** centroids, int ** clusters);
void simple_k_means(KmeansType type, size_t N, size_t k, KmeansCriteria criteria,size_t d,
		double ** data, double ** centroids,
		int ** clusters, double ** seeds);
double distortion(size_t d, size_t N, size_t k,
		double ** data, double ** centroids, int ** clusters);
}

#endif /* K_MEANS_H_ */
