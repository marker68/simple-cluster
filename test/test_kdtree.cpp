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
 *  test_kdtree.cpp
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
		N = 1000;
		d = 2;
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
			::delete data[i];
		::delete data;
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

TEST_F(KDTreeTest, DISABLED_test1) {
	KDNode<double> * a, * b;
	a = ::new KDNode<double>(d);
	b = ::new KDNode<double>(d);
	a->add_data(data[0]);
	b->add_data(data[1]);
	EXPECT_LT(0.0, kd_distance(a,b,false));
}

TEST_F(KDTreeTest, DISABLED_test2) {
	KDNode<double> * root = nullptr;
	double d[][2] = {{2.0,3.0},{4.0,3.0},{7.0,9.0}};
	for(int i = 0; i < 3; i++) {
		kd_insert(root,d[i],2,0,i,false);
	}
	kd_travel(root,2,0);
}

TEST_F(KDTreeTest, DISABLED_test3) {
	double ** _d;
	_d = new double*[10000];
	for(int i = 0; i < 10000; i++) {
		_d[i] = new double[2];
		_d[i][0] = _d[i][1] = static_cast<double>(i);
	}
	EXPECT_EQ(5000,find_median(_d,10000,2,1,true));
	//	cout << find_median(data,N,d,64) << endl;
}

TEST_F(KDTreeTest, DISABLED_test4) {
	KDNode<double> * root = nullptr;
	make_balanced_tree(root,data,N,d,0,0,false);
	kd_travel(root,d,0);
}

TEST_F(KDTreeTest, DISABLED_test5) {
	KDNode<double> * root = nullptr;
	make_random_tree(root,data,N,d,0,0,false);
	kd_travel(root,d,0);
}

TEST_F(KDTreeTest, test6) {
	size_t loop = 1;
	KDNode<double> * root = nullptr;
	make_random_tree(root,data,N,d,0,0,false);
	//	kd_travel<double>(root,d,0);
	KDNode<double> * query = ::new KDNode<double>(d);
	size_t pos = 10;
	query->add_data(data[pos]);
	for(int i = 0; i < loop; i++) {
		KDNode<double> * result = nullptr;
		double best_dist = DBL_MAX;
		nn_search(root,query,result,best_dist,d,0,false);
		EXPECT_EQ(0.0,best_dist);
	}
}

TEST_F(KDTreeTest, DISABLED_test7) {
	size_t loop = 1;
	KDNode<double> * root = nullptr;
	make_balanced_tree(root,data,N,d,0,0,false);
	//	kd_travel<double>(root,d,0);
	KDNode<double> * query = ::new KDNode<double>(d);
	size_t pos = 10;
	query->add_data(data[pos]);
	for(int i = 0; i < loop; i++) {
		KDNode<double> * result = nullptr;
		double best_dist = DBL_MAX;
		nn_search(root,query,result,best_dist,d,0,false);
		EXPECT_EQ(0.0,best_dist);
	}
}

TEST_F(KDTreeTest, test8) {
	size_t loop = 1;
	KDNode<double> * root = nullptr;
	make_balanced_tree(root,data,N,d,0,0,false);
	size_t pos = 10;
	for(int i = 0; i < loop; i++) {
		size_t best = 0;
		double best_dist = DBL_MAX;
		linear_search(data,data[pos],best,best_dist,N,d,true);
		EXPECT_EQ(0.0,best_dist);
	}
}

TEST_F(KDTreeTest, DISABLED_test9) {
	KDNode<double> * a, * b;
	a = ::new KDNode<double>(d);
	b = ::new KDNode<double>(d);
	a->add_data(data[0]);
	b->add_data(data[1]);
	unsigned long int t1, t2, t3;
	t1 = get_millisecond_time();
	double d1 = kd_distance(a,b,true);
	t2 = get_millisecond_time();
	double d2 = distance(data[0],data[1],d);
	t3 = get_millisecond_time();
	cout << "kd_distance's time:" << t2-t1 << "[ms]" << endl;
	cout << "distance's time:" << t3-t2 << "[ms]" << endl;
	EXPECT_EQ(d1,d2);
}

TEST_F(KDTreeTest, test10) {
	size_t loop = 1;
	KDNode<double> * root = nullptr;
	make_random_tree(root,data,N,d,0,0,false);
	//	kd_travel<double>(root,d,0);
	KDNode<double> * query = ::new KDNode<double>(d);
	size_t pos = 10;
	query->add_data(data[pos]);
	for(int i = 0; i < loop; i++) {
		KDNode<double> * result = nullptr;
		double best_dist = DBL_MAX;
		ann_search(root,query,result,best_dist,1.5,d,0,false);
		cout << "best distances " << best_dist << endl;
	}
}
