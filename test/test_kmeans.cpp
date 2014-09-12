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
 *  test_kmeans.cpp
 *
 *  Created on: 2014/09/12
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */


#include <iostream>
#include <vector>
#include <random>
#include <float.h>
#include <gtest/gtest.h>
#include <math.h>
#include <stdlib.h>
#include "k-means.h"
#include "utilities.h"

using namespace std;
using namespace SimpleCluster;

/**
 * Customized test case for testing
 */
class KmeansTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		N = 10000;
		d = 128;
		k = 256;
		int i, j;

		// For generating random numbers
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<double> real_dis(0.0, static_cast<double>(N));

		if(!init_array_2<double>(data,N,d)) {
			cerr << "Cannot allocate memory for test data!" << endl;
			exit(1);
		}
		for(i = 0; i < N; i++) {
			for(j = 0; j < d; j++) {
				data[i][j] = real_dis(gen);
			}
		}

		if(!init_array_2<double>(seeds,k,d)) {
			cerr << "Cannot allocate memory for seeds data!" << endl;
			exit(1);
		}
		if(!init_array_2<double>(centroids,k,d)) {
			cerr << "Cannot allocate memory for centroids data!" << endl;
			exit(1);
		}
		if(!init_array<int *>(clusters,k)) {
			cerr << "Cannot allocate memory for clusters data!" << endl;
			exit(1);
		}
		//		print_vector(data,d,N);
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		if(!dealloc_array_2<double>(data,N)) {
			cerr << "data: Deallocating failed!" << endl;
			exit(1);
		}
		if(!dealloc_array_2<double>(seeds,k)) {
			cerr << "seeds: Deallocating failed!" << endl;
			exit(1);
		}
		if(!dealloc_array_2<double>(centroids,k)) {
			cerr << "centroids: Deallocating failed!" << endl;
			exit(1);
		}
		if(!dealloc_array_2<int>(clusters,k)) {
			cerr << "clusters: Deallocating failed!" << endl;
			exit(1);
		}
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static double ** data;
	static double ** seeds;
	static double ** centroids;
	static int ** clusters;
	static int N, d, k;
};

double ** KmeansTest::data;
double ** KmeansTest::seeds;
double ** KmeansTest::centroids;
int ** KmeansTest::clusters;
int KmeansTest::N;
int KmeansTest::d;
int KmeansTest::k;

TEST_F(KmeansTest, test1) {
	random_seeds(d,N,k,data,seeds,true);
}

TEST_F(KmeansTest, test2) {
	kmeans_pp_seeds(d,N,k,data,seeds,true);
}

TEST_F(KmeansTest, test3) {
	random_seeds(d,N,k,data,seeds,true);
	assign_to_closest_centroid_2(d,N,k,data,seeds,clusters,true);
}

TEST_F(KmeansTest, DISABLED_test4) {
	KmeansCriteria criteria = {1.0,10};
	simple_k_means(KmeansType::RANDOM_SEEDS,N,k,criteria,d,
			data,centroids,clusters,seeds,true);
}
