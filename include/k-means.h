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

typedef vector<double> d_vector;
typedef vector<int> i_vector;

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

/**
 * Random seeding method
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param seeds seeds will be stored here
 * @return
 */
void random_seeds(size_t d, size_t N, size_t k, vector<d_vector> data, vector<d_vector>& seeds);

/**
 * K-means++'s seeding method
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param seeds seeds will be stored here
 * @return
 */
void kmeans_pp_seeds(size_t d, size_t N, size_t k, vector<d_vector> data, vector<d_vector>& seeds);

/**
 * Assign the data points to clusters
 * The execution time would be O(N*k*d)
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @param clusters
 */
void assign_to_closest_centroid(size_t d, size_t N, size_t k, vector<d_vector> data,
		vector<d_vector> centroids, vector<i_vector>& clusters);

/**
 * The k-means method
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
		vector<d_vector> data, vector<d_vector>& centroids,
		vector<i_vector>& clusters, vector<d_vector> seeds);

/**
 * Calculate the distortion of a set of clusters
 * @param d the number of dimensions
 * @param N the number of data
 * @param k the number of clusters
 * @param data the data
 * @param centroids
 * @param clusters
 */
double distortion(size_t d, size_t N, size_t k,
		vector<d_vector> data, vector<d_vector> centroids, vector<i_vector> clusters);
}

#endif /* K_MEANS_H_ */
