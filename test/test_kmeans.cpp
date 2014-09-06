/*
 * test_kmeans.cpp
 *
 *  Created on: 2014/09/06
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
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
		k = 16;
		int i, j, scale;
		d_vector tmp;
		double x;
		for(i = 0; i < N; i++) {
			for(j = 0; j < d; j++) {
				scale = rand() % N;
				x = sqrt(scale);
				tmp.push_back(x);
			}
			data.push_back(tmp);
			tmp.clear();
		}
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		data.clear();
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static vector<d_vector> data;
	static int N;
	static int d;
	static int k;
};

vector<d_vector> KmeansTest::data;
int KmeansTest::N;
int KmeansTest::d;
int KmeansTest::k;

TEST_F(KmeansTest, test1) {
	vector<d_vector> centroids, seeds;
	vector<i_vector> clusters;
	KmeansCriteria criteria = {11.0,10};
	unsigned long int t1, t2;
	t1 = get_millisecond_time();
	simple_k_means(KmeansType::RANDOM_SEEDS,N,k,criteria,d,data,centroids,clusters,seeds);
	t2 = get_millisecond_time();
	cout << "Finished in " << t2-t1 << "[ms]" << endl;
	cout << "The distortion is " << distortion(d,N,k,data,centroids,clusters) << endl;
}

/*int main(int argc, char * argv[]) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}*/
