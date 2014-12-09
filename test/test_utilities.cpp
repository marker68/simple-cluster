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
 *  test_utilities.cpp
 *
 *  Created on: 2014/09/09
 *      Author: Nguyen Anh Tuan <t_nguyen@hal.t.u-tokyo.ac.jp>
 */

#include <iostream>
#include <vector>
#include <random>
#include <gtest/gtest.h>
#include <cmath>
#include <cstdlib>
#include "utilities.h"

#ifdef _OPENMP
#include <omp.h>
#define SET_THREAD_NUM omp_set_num_threads(n_thread)
#else
#define SET_THREAD_NUM 0 // disable multi-thread
#endif

using namespace std;
using namespace SimpleCluster;

/**
 * A comparator
 */
int compare_float(const float * _a, const float * _b) {
	if(*_a > *_b) return 1;
	else if(*_a < *_b) return -1;
	else return 0;
}

/**
 * Customized test case for testing
 */
class UtilTest : public ::testing::Test {
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
			::new(data +i) float *;
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
	static int N;
	static int d;
};

float ** UtilTest::data;
int UtilTest::N;
int UtilTest::d;

TEST_F(UtilTest, test1) {
	float arr[] = {1.0, 3.0, 5.0, 7.0, 9.0};
	float m = quick_select_k(arr,5,3,compare_float);
	EXPECT_EQ(7.0,m);
}

TEST_F(UtilTest, test2) {
	float arr[10000];
	for(int i = 0; i < 10000; i++)
		arr[i] = i;
	float m = quick_select_k(arr,10000,5000,compare_float);
	EXPECT_EQ(5000.0,m);
}

TEST_F(UtilTest, test3) {
	float * t;
	EXPECT_TRUE(init_array<float>(t,1000) && t!=nullptr);
}

TEST_F(UtilTest, test4) {
	float * t;
	EXPECT_TRUE(init_array<float>(t,d) && t != nullptr);
	EXPECT_TRUE(copy_array<float>(data[0],t,d));
	for(int i = 0; i < d; i++) {
		EXPECT_TRUE(data[0][i] == t[i]);
	}
}

TEST_F(UtilTest, test5) {
	float ** t;
	EXPECT_TRUE(init_array_2<float>(t,1000,200) && t!=nullptr);
}

TEST_F(UtilTest, test7) {
	unsigned char x[3] = {1,2,3};
	unsigned char y[3] = {128,225,123};
	EXPECT_LT(0.0,distance_l2<unsigned char>(x,y,3));
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
