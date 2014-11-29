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
#include <cfloat>
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
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
		uniform_real_distribution<float> real_dis(0.0, static_cast<float>(N));

		data = (float **)::operator new(N * sizeof(float *));
		for(i = 0; i < N; i++) {
			::new(data + i) float *;
			data[i] = (float *)::operator new(d * sizeof(float));
			for(j = 0; j < d; j++) {
				::new(data[i] + j) float;
				data[i][j] = real_dis(gen);
			}
		}
	}

	// Per-test-case tear-down.
	// Called after the last test in this test case.
	// Can be omitted if not needed.
	static void TearDownTestCase() {
		int i;
		for(i = 0; i < N; i++)
			::delete data[i];
		::delete data;
	}

	// You can define per-test set-up and tear-down logic as usual.
	virtual void SetUp() { }
	virtual void TearDown() {}

public:
	// Some expensive resource shared by all tests.
	static float ** data;
	static int N, d;
};

float ** KDTreeTest::data;
int KDTreeTest::N;
int KDTreeTest::d;

TEST_F(KDTreeTest, DISABLED_test1) {
	KDNode<float> * a, * b;
	a = ::new KDNode<float>(d);
	b = ::new KDNode<float>(d);
	a->add_data(data[0]);
	b->add_data(data[1]);
	EXPECT_LT(0.0, kd_distance<float>(a,b,DistanceType::NORM_L2,false));
}

TEST_F(KDTreeTest, DISABLED_test2) {
	KDNode<float> * root = nullptr;
	float d[][2] = {{2.0,3.0},{4.0,3.0},{7.0,9.0}};
	for(int i = 0; i < 3; i++) {
		kd_insert<float>(root,d[i],2,0,i,false);
	}
	kd_travel<float>(root,2,0);
}

TEST_F(KDTreeTest, DISABLED_test3) {
	float ** _d;
	_d = new float*[10000];
	for(int i = 0; i < 10000; i++) {
		_d[i] = new float[2];
		_d[i][0] = _d[i][1] = static_cast<float>(i);
	}
	EXPECT_EQ(5000,find_median<float>(_d,10000,2,1,true));
	//	cout << find_median(data,N,d,64) << endl;
}

TEST_F(KDTreeTest, DISABLED_test4) {
	KDNode<float> * root = nullptr;
	make_balanced_tree<float>(root,data,N,d,0,0,false);
	kd_travel<float>(root,d,0);
}

TEST_F(KDTreeTest, DISABLED_test5) {
	KDNode<float> * root = nullptr;
	make_random_tree<float>(root,data,N,d,0,false);
	kd_travel<float>(root,d,0);
}

TEST_F(KDTreeTest, test6) {
	unsigned long int t1, t2, t3;
	t1 = get_millisecond_time();
	KDNode<float> * root = nullptr;
	make_random_tree<float>(root,data,N,d,0,false);
	t2 = get_millisecond_time();
	//	kd_travel<float>(root,d,0);
	KDNode<float> * query = ::new KDNode<float>(d);
	int pos = 10;
	query->add_data(data[pos]);
	KDNode<float> * result = nullptr;
	double best_dist = DBL_MAX;
	int visited = 0;
	nn_search<float>(root,query,result,DistanceType::NORM_L2,best_dist,d,0,visited,false);
	t3 = get_millisecond_time();
	EXPECT_EQ(0.0,best_dist);
	cout << "Tree build time is " << t2-t1 << "[ms]" << endl;
	cout << "Search time is " << t3-t2 << "[ms]" << endl;
	cout << "Visited " << visited << " nodes" << endl;
}

TEST_F(KDTreeTest, test7) {
	unsigned long int t1, t2, t3;
	t1 = get_millisecond_time();
	KDNode<float> * root = nullptr;
	make_balanced_tree<float>(root,data,N,d,0,0,false);
	t2 = get_millisecond_time();
	cout << "Tree build time is " << t2-t1 << "[ms]" << endl;
	KDNode<float> * query = ::new KDNode<float>(d);
	int pos = 10;
	query->add_data(data[pos]);
	KDNode<float> * result = nullptr;
	double best_dist = DBL_MAX;
	int visited = 0;
	nn_search<float>(root,query,result,DistanceType::NORM_L2,best_dist,d,0,visited,false);
	t3 = get_millisecond_time();
	EXPECT_EQ(0.0,best_dist);
	cout << "Search time is " << t3-t2 << "[ms]" << endl;
	cout << "Visited " << visited << " nodes" << endl;
}

TEST_F(KDTreeTest, DISABLED_test8) {
	KDNode<float> * root = nullptr;
	make_balanced_tree<float>(root,data,N,d,0,0,false);
	int pos = 10;
	int best = 0;
	double best_dist = DBL_MAX;
	linear_search<float>(data,data[pos],DistanceType::NORM_L2,best,best_dist,N,d,true);
	EXPECT_EQ(0.0,best_dist);
}

TEST_F(KDTreeTest, DISABLED_test9) {
	KDNode<float> * a, * b;
	a = ::new KDNode<float>(d);
	b = ::new KDNode<float>(d);
	a->add_data(data[0]);
	b->add_data(data[1]);
	unsigned long int t1, t2, t3;
	t1 = get_millisecond_time();
	double d1 = kd_distance<float>(a,b,DistanceType::NORM_L2,true);
	t2 = get_millisecond_time();
	double d2 = distance_l2<float>(data[0],data[1],d);
	t3 = get_millisecond_time();
	cout << "kd_distance's time:" << t2-t1 << "[ms]" << endl;
	cout << "distance's time:" << t3-t2 << "[ms]" << endl;
	EXPECT_EQ(d1,d2);
}

TEST_F(KDTreeTest, test10) {
	KDNode<float> * root = nullptr;
	make_random_tree(root,data,N,d,0,false);
	//	kd_travel<float>(root,d,0);
	KDNode<float> * query = ::new KDNode<float>(d);
	int pos = 10;
	query->add_data(data[pos]);
	KDNode<float> * result = nullptr;
	double best_dist = DBL_MAX;
	int visited = 0;
	ann_search<float>(root,query,result,DistanceType::NORM_L2,best_dist,1.5,d,0,visited,false);
	cout << "best distances " << best_dist << endl;
	cout << "Visited " << visited << " nodes" << endl;
}

int main(int argc, char * argv[])
{
	/*The method is initializes the Google framework and must be called before RUN_ALL_TESTS */
	::testing::InitGoogleTest(&argc, argv);

	/*RUN_ALL_TESTS automatically detects and runs all the tests defined using the TEST macro.
	It's must be called only once in the code because multiple calls lead to conflicts and,
	therefore, are not supported.
	*/
	return RUN_ALL_TESTS();
}
