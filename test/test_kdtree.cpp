/*
 * test_kdtree.cpp
 *
 *  Created on: 2014/09/09
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

/*
 * test_utilities.cpp
 *
 *  Created on: 2014/09/09
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <random>
#include <float.h>
#include <gtest/gtest.h>
#include <math.h>
#include <stdlib.h>
#include "kd-tree.h"
#include "utilities.h"

using namespace std;
using namespace SimpleCluster;

/**
 * Customized test case for testing
 */
class KDTreeTest : public ::testing::Test {
protected:
	// Per-test-case set-up.
	// Called before the first test in this test case.
	// Can be omitted if not needed.
	static void SetUpTestCase() {
		N = 10000;
		d = 128;
		int i, j;

		// For generating random numbers
		random_device rd;
		mt19937 gen(rd());
		uniform_real_distribution<double> real_dis(0.0, static_cast<double>(N));

		data = (double **)::operator new(N * sizeof(double *));
		for(i = 0; i < N; i++) {
			::new(data + i) double *;
			data[i] = (double *)::operator new(d * sizeof(double));
			for(j = 0; j < d; j++) {
				::new(data[i] + j) double;
				data[i][j] = real_dis(gen);
			}
		}
//		print_vector(data,d,N);
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		size_t i;
		for(i = 0; i < N; i++)
			::delete[] data[i];
		::delete[] data;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static double ** data;
	static int N, d;
};

double ** KDTreeTest::data;
int KDTreeTest::N;
int KDTreeTest::d;

TEST_F(KDTreeTest, test1) {
	KDNode<double> * a, * b;
	a = ::new KDNode<double>;
	b = ::new KDNode<double>;
	for(int i = 0; i < d; i++) {
		a->add_data(data[0][i]);
		b->add_data(data[1][i]);
	}
	EXPECT_LT(0.0, kd_distance(a,b,d,false));
}

TEST_F(KDTreeTest, test2) {
	KDNode<double> * root = NULL;
	double d[][2] = {{2.0,3.0},{4.0,3.0},{7.0,9.0}};
	for(int i = 0; i < 3; i++) {
		kd_insert(root,d[i],2,0,i,false);
	}
	kd_travel(root,2,0);
}

TEST_F(KDTreeTest, test3) {
	double ** _d;
	_d = new double*[10000];
	for(int i = 0; i < 10000; i++) {
		_d[i] = new double[2];
		_d[i][0] = _d[i][1] = static_cast<double>(i);
	}
	EXPECT_EQ(5000,find_median(_d,10000,2,1,true));
//	cout << find_median(data,N,d,64) << endl;
}

TEST_F(KDTreeTest, test4) {
	KDNode<double> * root = NULL;
	make_balanced_tree(root,data,N,d,0,0,false);
//	kd_travel(root,d,0);
}

TEST_F(KDTreeTest, test5) {
	KDNode<double> * root = NULL;
	make_random_tree(root,data,N,d,0,0,false);
//	kd_travel(root,d,0);
}

TEST_F(KDTreeTest, test6) {
	KDNode<double> * root = NULL;
	make_random_tree(root,data,N,d,0,0,false);
//	kd_travel<double>(root,d,0);
	KDNode<double> * result = NULL;
	double best_dist = DBL_MAX;
	nn_search(root,data[33],result,best_dist,d,0,false);
	cout << "best:" << best_dist << " at " << result->id << endl;
	EXPECT_EQ(0.0,best_dist);
}

TEST_F(KDTreeTest, test7) {
	KDNode<double> * root = NULL;
	make_balanced_tree(root,data,N,d,0,0,false);
//	kd_travel<double>(root,d,0);
	KDNode<double> * result = NULL;
	double best_dist = DBL_MAX;
	nn_search(root,data[33],result,best_dist,d,0,false);
	cout << "best:" << best_dist << " at " << result->id << endl;
	EXPECT_EQ(0.0,best_dist);
}
